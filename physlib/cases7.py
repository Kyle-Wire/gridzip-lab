"""PROTOCOL-7 part D: the real cases again with the v2 engine (grammar and guards on). Exploratory, not a test.

python -m physlib.cases7 [w1 w2 storm]       writes results/physlib/v7_D_<case>.json
  w1     AIJ Case C, LBM D/10 and RANS (results/physics/aij_c_arrays.npz)
  w2     Layer 1 storm gusts, thunderstorm and ordinary (results/physics/layer1.json)
  storm  WindSim protocol 6's 903 thunderstorm gusts with cold-pool, pressure-jump and dryness inputs (spent data, post
         hoc): STORM_TERMS=<storm_terms.json> (gridzip-lab branch claude/project-thread-fnbdyl, 75980a8)
"""
import json, math, os, sys, time, warnings
from pathlib import Path

import numpy as np

from . import engine as E, terms as T, guards as Q
from .validate7 import _safe

ROOT = Path(__file__).resolve().parents[1]
PHYS = ROOT / "results" / "physics"; OUT = ROOT / "results" / "physlib"
sys.path.insert(0, str(ROOT / "physics"))


def summary(res, before=None):
    gd = res["guards"]
    out = {"winner": res["winner"]["terms"], "coef": res["winner"]["coef"], "alpha": res["winner"]["alpha"], "code": res["winner"]["code"],
           "bits_saved_per_reading": res["bits_saved_per_reading"], "sd_before": res["sd_before"], "sd_after": res["sd_after"],
           "transfers": res["transfers"], "holdout": res["holdout"], "library_size": res["library_size"], "grammar_size": res["grammar_size"],
           "verdict": Q.verdict(res), "guards": {k: v for k, v in gd.items() if k != "null_bits_saved"}, "n": res["n"]}
    if "null_bits_saved" in gd: out["null_bits_saved"] = gd["null_bits_saved"]
    if before is not None: out["findings4_winner"] = before
    return out


def w1():
    import aij_c as A
    d = np.load(PHYS / "aij_c_arrays.npz", allow_pickle=True); out = {}; arr = {}
    old = json.loads((OUT / "w1.json").read_text()) if (OUT / "w1.json").exists() else {}
    g = E.Grid(0.01, -4.0, 4.0)
    for eng in ("lbm", "rans"):
        x, y, m, p, hc, lay = (d[f"{eng}_{k}"] for k in ("x", "y", "m", "p", "hc", "layout"))
        r = np.log(m / p)
        sx = np.array([A.street(v) for v in x]); sy = np.array([A.street(v) for v in y])
        dx = np.maximum(np.abs(x) - 0.1, 0); dy = np.maximum(np.abs(y) - 0.1, 0)
        ctx = {"lp": np.log(p / A.U_APP), "fetch": (x + 0.5) / A.D, "along": (sy & ~sx).astype(float), "cross": (sx & ~sy).astype(float),
               "hc": hc.astype(float), "corner": ((np.hypot(dx, dy) <= 0.1 + 1e-9) & (hc > 0)).astype(float),
               "lateral": np.abs(y) / A.D, "x_pos": x / A.D}
        folds = [(f"leave out {L}", np.where(lay != L)[0], np.where(lay == L)[0]) for L in ("0D", "1D", "2D")]
        t = time.time()
        res = Q.analyse(T.WIND, ctx, r, g, folds=folds, groups=lay, bootstrap_B=30, null_N=19, name=f"AIJ Case C {eng.upper()} D/10",
                        exclude=(), boot_by_group=False)
        out[eng] = summary(res, old.get(eng, {}).get("winner", {}).get("terms"))
        arr[eng] = {"y": res["_y"], "mu": res["_mu"], "mu0": res["_mu0"], "x": x, "yy": y}
        print("W1", eng, out[eng]["winner"], round(out[eng]["bits_saved_per_reading"], 3), f"{time.time() - t:.0f}s", flush=True)
    return out, arr


def w2():
    L1 = json.loads((PHYS / "layer1.json").read_text()); out = {}; arr = {}
    old = json.loads((OUT / "w2.json").read_text()) if (OUT / "w2.json").exists() else {}
    g = E.Grid(0.005, -3.0, 3.0)
    for q, lab in (("Q1_thunderstorm", "thunderstorm"), ("Q2_ordinary", "ordinary")):
        D = L1[q]["data"]; y = np.array(D["y"], float); n = len(y)
        st = np.array(D["station"]); pl = np.array(D["place"]); key = np.char.add(np.char.add(pl, "|"), st)
        loo = np.zeros(n)
        for k in set(key.tolist()):
            ii = np.where(key == k)[0]
            if len(ii) > 1: loo[ii] = (y[ii].sum() - y[ii]) / (len(ii) - 1)
        ctx = {"stn_loo": loo, "dz": np.array(D["dz_elev"], float), "dirn": np.array(D["dirn"], float),
               "hour": np.array(D["hl"], float), "month": np.array(D["month"], float)}
        folds = [(f"leave out {p}", np.where(pl != p)[0], np.where(pl == p)[0]) for p in sorted(set(pl.tolist()))]
        t = time.time()
        res = Q.analyse(T.WIND, ctx, y, g, folds=folds, groups=key, bootstrap_B=30, null_N=19, name=f"Layer 1 {lab}", exclude=("stn_loo",))
        out[lab] = summary(res, old.get(lab, {}).get("winner", {}).get("terms"))
        print("W2", lab, out[lab]["winner"], round(out[lab]["bits_saved_per_reading"], 3), f"{time.time() - t:.0f}s", flush=True)
    return out, arr


def storm():
    path = Path(os.environ.get("STORM_TERMS", PHYS / "storm_terms.json"))
    st = json.loads(path.read_text()); L1 = json.loads((PHYS / "layer1.json").read_text())["Q1_thunderstorm"]["data"]
    key = {(p, s, t): i for i, (p, s, t) in enumerate(zip(L1["place"], L1["station"], L1["storm"]))}
    rows = []
    for k, meta in enumerate(st["meta"]):
        r = st["pairs"].get(str(k)); i = key.get((meta["place"], meta["station"], meta["storm"]))
        if r is None or i is None: continue
        rows.append((meta, r, i))
    f = lambda nm: np.array([np.nan if r.get(nm) is None else r[nm] for _, r, _ in rows], float)
    dT, dp, dtd, dry, T0 = f("dT"), f("dp"), f("dTd"), f("dryness"), f("T0")
    keep = np.flatnonzero(~(np.isnan(dT) & np.isnan(dp) & np.isnan(dtd) & np.isnan(dry)))
    fill = lambda v: np.where(np.isnan(v), np.nanmean(v), v)[keep]
    idx = np.array([i for _, _, i in rows])[keep]
    y = np.array(L1["y"], float)[idx]; n = len(y)
    pred = np.array([m["pred"] for m, _, _ in rows], float)[keep]
    station = np.array([f'{m["place"]}|{m["station"]}' for m, _, _ in rows])[keep]
    place = np.array([m["place"] for m, _, _ in rows])[keep]; group = np.array([m["group"] for m, _, _ in rows])[keep]
    loo = np.zeros(n)
    for s in set(station.tolist()):
        ii = np.where(station == s)[0]
        if len(ii) > 1: loo[ii] = (y[ii].sum() - y[ii]) / (len(ii) - 1)
    ctx = {"dT": fill(dT), "Tk": fill(T0) + 273.15, "dp": 100 * fill(dp), "u_model": pred, "dTd": fill(dtd), "dryness": fill(dry),
           "stn_loo": loo, "dz": np.array(L1["dz_elev"], float)[idx], "dirn": np.array(L1["dirn"], float)[idx],
           "hour": np.array(L1["hl"], float)[idx], "month": np.array(L1["month"], float)[idx]}
    folds = [("calibration -> blind", np.where(group == "cal")[0], np.where(group == "blind")[0])]
    folds += [(f"leave out {p}", np.where(place != p)[0], np.where(place == p)[0]) for p in sorted(set(place.tolist()))]
    g = E.Grid(0.005, -3.0, 3.0); t = time.time()
    res = Q.analyse(T.WIND_FULL, ctx, y, g, folds=folds, groups=station, bootstrap_B=30, null_N=19, name="storm table (post hoc)",
                    exclude=("stn_loo",))
    out = {"storm": summary(res, "FINDINGS-6: station offsets only (cold_pool +0.52 per sqrt(dT/T), below MDE)")}
    out["storm"]["n_pairs"] = n; out["storm"]["source"] = str(path.name)
    print("storm", out["storm"]["winner"], round(out["storm"]["bits_saved_per_reading"], 3), f"{time.time() - t:.0f}s", flush=True)
    return out, {}


if __name__ == "__main__":
    warnings.filterwarnings("ignore"); OUT.mkdir(parents=True, exist_ok=True)
    for c in (sys.argv[1:] or ["w1", "w2", "storm"]):
        o, a = {"w1": w1, "w2": w2, "storm": storm}[c]()
        (OUT / f"v7_D_{c}.json").write_text(json.dumps(_safe(o), indent=1))
        if a: np.savez_compressed(OUT / f"v7_D_{c}_arrays.npz", **{f"{k}_{kk}": vv for k, v in a.items() for kk, vv in v.items()})
