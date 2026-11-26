"""Protocol 13 part H (A5): is the rate error left after the bridged fix a fading-signal (sawtooth) bias?

    python -m physics.insar_m2h [synth|H1|H2|all]

Reads data/insar_m2/fix_rerun/v2_fix_points.npz, writes results/physics/insar_m2_partH.json.
Written fresh for SPARC / gridzip-lab.
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

from physics import insar_m2 as M
from physlib import guards as Q
from physlib.insar import INSAR_RATE, SAWTOOTH, sawtooth_amplitude

FIX = M.ROOT / "data/insar_m2/fix_rerun/v2_fix_points.npz"
OUT = M.ROOT / "results/physics/insar_m2_partH.json"
TERMS = INSAR_RATE + [SAWTOOTH]


def load():
    M.NPZ = FIX
    D = M.load()
    D["valid"] = D["valid"] & ~D["bridged"]
    ok = np.all(np.isfinite(D["series_mm"]), axis=1)
    Y = np.where(ok[:, None], D["series_mm"], 0.0)
    amp, P = sawtooth_amplitude(Y, D["t_yr"], D["ministack"].astype(int), ref_rows=np.where(D["valid"])[0])
    amp[~ok] = np.nan
    D["sawtooth"] = amp; D["P"] = P
    return D


def ctx(D, rows):
    c = M.ctx_for(D, rows); c["sawtooth"] = np.asarray(D["sawtooth"][rows], float); return c


def analyse(D, rows, y, seed, B=10, N=9):
    c = ctx(D, rows); tiles = D["tile"][rows]; t = time.time()
    res = Q.analyse(TERMS, c, y, M.GRID, grammar=True, folds=M.folds_for(c), groups=tiles, bootstrap_B=B, null_N=N,
                    blocks=M.blocks_for(tiles), name=f"V2-fix part H seed {seed}", seed=seed, boot_by_group=True)
    out = M.summary(res); out["seconds"] = time.time() - t; out["n"] = len(rows); return out


def synth(D, seeds=(0, 1, 2, 3, 4)):
    rows = M.sample(D, 20, n=3000); c = ctx(D, rows); base = np.asarray(D["rate_mm_yr"][rows], float)
    o = {"null": [], "plant": []}
    for s in seeds:
        rng = np.random.default_rng(200 + s); y0 = rng.permutation(base)
        o["null"].append(analyse(D, rows, y0, s, B=0, N=0)["winner"])
        a = M._c(c["sawtooth"]); k = np.sqrt(0.10 / 0.90) * np.std(y0) / np.std(a)
        o["plant"].append(analyse(D, rows, y0 + k * a, s, B=0, N=0)["winner"])
    o["pass_null"] = sum(1 for w in o["null"] if not w) >= 4
    o["pass_plant"] = sum(1 for w in o["plant"] if any("sawtooth" in t for t in w)) >= 4
    return o


def stations(D):
    t = D["t_yr"]; out = {}
    for i, s in enumerate(D["stations"]):
        r, g = D["station_harness_radar_mm"][i], D["station_gnss_los_mm"][i]; m = np.isfinite(r) & np.isfinite(g)
        if m.sum() < 10: continue
        near = (np.hypot((D["lat"] - D["station_lat"][i]) * 110.57, (D["lon"] - D["station_lon"][i]) * 111.32 * np.cos(np.radians(D["station_lat"][i]))) < 0.2) & D["valid"]
        res_rate = np.linalg.lstsq(M.design(t[m]), (r - g)[m], rcond=None)[0][1]
        out[str(s)] = {"n_dates": int(m.sum()), "resid_rate": float(res_rate), "sawtooth_200m": float(np.nanmean(D["sawtooth"][near])), "n_pts": int(near.sum())}
    mu = np.mean([v["resid_rate"] for v in out.values()])
    for v in out.values(): v["resid_rate_minus_mean"] = v["resid_rate"] - mu
    xs = [v["sawtooth_200m"] for v in out.values()]; ys = [v["resid_rate_minus_mean"] for v in out.values()]
    return {"stations": out, "corr": float(np.corrcoef(xs, ys)[0, 1])}


def main(part="all"):
    D = load(); res = json.loads(OUT.read_text()) if OUT.exists() else {}
    v = D["valid"] & np.isfinite(D["sawtooth"])
    res["pattern"] = {"P_mm": D["P"].tolist(), "sawtooth_median": float(np.median(D["sawtooth"][v])),
                      "sawtooth_rsd": M.rsd(D["sawtooth"][v]), "corr_rate": float(np.corrcoef(D["sawtooth"][v], D["rate_mm_yr"][v])[0, 1])}
    print("pattern", {k: x for k, x in res["pattern"].items() if k != "P_mm"}, flush=True)
    if part in ("synth", "all"):
        res["H0"] = synth(D); print("H0", res["H0"], flush=True)
    if part in ("H2", "all"):
        res["H2"] = stations(D); print("H2", res["H2"], flush=True)
    if part in ("H1", "all"):
        r1 = M.sample(D, 30, n=3000); r2 = M.sample(D, 31, n=3000, exclude=r1)
        s1 = analyse(D, r1, np.asarray(D["rate_mm_yr"][r1], float), 30); print("H1a", s1["winner"], s1["sd_before"], s1["sd_after"], flush=True)
        s2 = analyse(D, r2, np.asarray(D["rate_mm_yr"][r2], float), 31); print("H1b", s2["winner"], s2["sd_before"], s2["sd_after"], flush=True)
        res["H1"] = {"a": s1, "b": s2, "named": M.named_A(s1, s2)}; print("named", res["H1"]["named"], flush=True)
    OUT.write_text(json.dumps(M._clean(res), indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "all")
