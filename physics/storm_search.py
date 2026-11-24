"""Protocol 6, step 1a: do storm terms explain WindSim's thunderstorm gust misses beyond the station offsets?

python physics/storm_search.py [--synthetic-only]
Reads results/physics/storm_terms.json (physics/storm_terms.py) and layer1.json (protocol 3's features).
Writes results/physics/storm_search.json and physics/figures/storm_*.png.
"""
import argparse, json, math, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "physics"))
from physlib import engine as E, terms as T  # noqa: E402
from physlib.engine import Term  # noqa: E402

G = E.Grid(0.005, -3.0, 3.0)
BASE = ("station_bias", "elevation", "direction", "stability", "season")
STORM = [
    Term("cold_pool", "wind", "b · sqrt(ΔT⁺ / T₀)", "the cold pool behind the gust front pushes harder than the chain knows (density current)",
         ("cold_pool",), lambda c: c["cold_pool"]),
    Term("pressure_jump", "wind", "b · Δp⁺", "dense outflow air arrives with a pressure jump", ("dp",), lambda c: c["dp"]),
    Term("dew_drop", "wind", "b · ΔTd⁺", "rain-cooled downdraft air reaches the ground", ("dtd",), lambda c: c["dtd"]),
    Term("dryness", "wind", "b · (T − Td) before", "dry air under the cloud: evaporation drives stronger downdrafts", ("dry",), lambda c: c["dry"]),
    Term("gust_front", "wind", "δ · [ΔT ≥ 3 K and Δp ≥ 1 hPa]", "a gust front passed the station", ("gf",), lambda c: c["gf"]),
    Term("heat", "wind", "b · T₀", "hot afternoons: stronger convection", ("heat",), lambda c: c["heat"]),
]
STORM_NAMES = [t.name for t in STORM]


def library():
    return [t for t in T.WIND if t.name in BASE] + STORM


def load():
    st = json.loads((ROOT / "results/physics/storm_terms.json").read_text())
    L1 = json.loads((ROOT / "results/physics/layer1.json").read_text())["Q1_thunderstorm"]["data"]
    key = {(p, s, t): i for i, (p, s, t) in enumerate(zip(L1["place"], L1["station"], L1["storm"]))}
    rows = []
    for k, meta in enumerate(st["meta"]):
        r = st["pairs"].get(str(k)); i = key.get((meta["place"], meta["station"], meta["storm"]))
        if r is None or i is None: continue
        y = math.log(meta["obs"] / meta["pred"])
        assert abs(y - L1["y"][i]) < 1e-3, (meta, y, L1["y"][i])
        rows.append((meta, r, i))
    n = len(rows); f = lambda name: np.array([np.nan if r.get(name) is None else r[name] for _, r, _ in rows], float)
    dT, dp, dtd, dry, T0 = f("dT"), f("dp"), f("dTd"), f("dryness"), f("T0")
    have = ~(np.isnan(dT) & np.isnan(dp) & np.isnan(dtd) & np.isnan(dry))
    fill = lambda v: np.where(np.isnan(v), np.nanmean(v), v)
    ctx = {"cold_pool": fill(np.sqrt(np.maximum(dT, 0) / (T0 + 273.15))), "dp": fill(np.maximum(dp, 0)), "dtd": fill(np.maximum(dtd, 0)),
           "dry": fill(dry), "gf": ((dT >= 3) & (dp >= 1)).astype(float), "heat": fill(T0)}
    idx = [i for _, _, i in rows]
    for k_src, k_dst in (("stn_loo", "stn_loo"), ("dz_elev", "dz"), ("dirn", "dirn"), ("hl", "hour"), ("month", "month")):
        ctx[k_dst] = np.array(L1[k_src], float)[idx]
    y = np.array(L1["y"], float)[idx]
    place = np.array([m["place"] for m, _, _ in rows]); group = np.array([m["group"] for m, _, _ in rows])
    station = np.array([m["station"] for m, _, _ in rows])
    info = {"pairs_in_file": len(st["meta"]), "pairs_used": int(have.sum()), "matched_rate": float(np.mean([r.get("matched", False) for _, r, _ in rows])),
            "missing": {k: int(np.isnan(v).sum()) for k, v in (("dT", dT), ("dp", dp), ("dTd", dtd), ("dryness", dry))},
            "gust_front_share": float(ctx["gf"].mean())}
    keep = np.where(have)[0]
    ctx = {k: v[keep] for k, v in ctx.items()}
    return y[keep], ctx, place[keep], group[keep], station[keep], info, {"dT": dT[keep], "dp": dp[keep]}


def folds(place, group):
    out = [("calibration -> blind", np.where(group == "cal")[0], np.where(group == "blind")[0])]
    out += [(f"leave out {p}", np.where(place != p)[0], np.where(place == p)[0]) for p in sorted(set(place.tolist()))]
    return out


def synthetic(ctx, station, y, seeds=10):
    rng0 = np.random.default_rng(1)
    stns = sorted(set(station.tolist())); sd_stn = float(np.std([y[station == s].mean() for s in stns]))
    sd_noise = float(np.sqrt(max(np.var(y) - sd_stn ** 2, 1e-4)))
    out = {"sd_station": sd_stn, "sd_noise": sd_noise}

    def sim(seed, extra):
        rng = np.random.default_rng(seed); off = {s: rng.normal(0, sd_stn) for s in stns}
        yy = np.array([off[s] for s in station]) + rng.normal(0, sd_noise, len(y)) + extra
        c = dict(ctx); loo = np.zeros(len(y))
        for s in stns:
            ii = np.where(station == s)[0]
            if len(ii) > 1: loo[ii] = (yy[ii].sum() - yy[ii]) / (len(ii) - 1)
        c["stn_loo"] = loo
        return E.run_case(library(), c, yy, G, (), max_terms=4)["winner"]["terms"]
    picks = [sim(100 + s, 0.0) for s in range(seeds)]
    out["a offsets only"] = {"no_storm_term": sum(not set(p) & set(STORM_NAMES) for p in picks), "of": seeds, "winners": picks}
    cp = (ctx["cold_pool"] - ctx["cold_pool"].mean()) / ctx["cold_pool"].std()
    gf = ctx["gf"] - ctx["gf"].mean()
    out["power"] = {}
    for nm, v in (("cold_pool", cp), ("gust_front", gf / gf.std())):
        row = {}
        for sd in (0.02, 0.04, 0.06, 0.08, 0.10):
            row[str(sd)] = sum(nm in sim(200 + s, sd * v) for s in range(seeds))
        out["power"][nm] = {"found_of_10_by_sd": row, "mde_sd": next((float(k) for k, h in row.items() if h >= 8), None)}
    out["passed"] = out["a offsets only"]["no_storm_term"] >= 9
    return out


def coverage(y, mask_cal, mask_blind, split):
    """10-90 % range of observed/predicted from calibration, applied to blind: single, and per half of `split`."""
    r = {}
    lo, hi = np.percentile(y[mask_cal], [10, 90]); r["single"] = float(np.mean((y[mask_blind] >= lo) & (y[mask_blind] <= hi)))
    if split is not None:
        hit = np.zeros(mask_blind.sum(), bool); yb = y[mask_blind]; sb = split[mask_blind]
        for v in (False, True):
            c = mask_cal & (split == v)
            if c.sum() < 20: return r | {"split": None, "note": "fewer than 20 calibration pairs in a half"}
            lo, hi = np.percentile(y[c], [10, 90]); m = sb == v; hit[m] = (yb[m] >= lo) & (yb[m] <= hi)
        r["split"] = float(hit.mean())
    return r


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--synthetic-only", action="store_true"); a = ap.parse_args()
    y, ctx, place, group, station, info, raw = load()
    out = {"protocol": "physics/PROTOCOL-6.md", "data": info}
    print("data", info, flush=True)
    syn = synthetic(ctx, station, y); out["synthetic"] = syn
    print("synthetic", json.dumps({k: v for k, v in syn.items() if k != "a offsets only"}), syn["a offsets only"]["no_storm_term"], flush=True)
    res_path = ROOT / "results/physics/storm_search.json"
    if a.synthetic_only or not syn["passed"]:
        res_path.write_text(json.dumps(out, indent=1, default=float)); return
    fo = folds(place, group)
    res = E.run_case(library(), ctx, y, G, fo, max_terms=4, name="Layer 1 thunderstorm gusts + storm terms")
    out["search"] = E.public(res)
    built = E.assemble(library(), ctx, len(y)); b0 = E.exact_bits([], y, G, len(built))["total_bits"]
    base = [b for b in built if b.name == "station_bias"]
    bb = E.exact_bits(base, y, G, len(built))["total_bits"]
    out["single_term_savings_bits"] = {b.name: float(b0 - E.exact_bits([b], y, G, len(built))["total_bits"]) for b in built}
    out["savings_beyond_station_bits"] = {b.name: float(bb - E.exact_bits(base + [b], y, G, len(built))["total_bits"])
                                          for b in built if b.name != "station_bias"}
    out["coef_beyond_station"] = {}
    for b in built:
        if b.name in STORM_NAMES:
            th, al, _, _ = E.fit(base + [b], y); out["coef_beyond_station"][b.name] = float(th[-1])
    win = res["winner"]["terms"]; ho = res["holdout"]
    cb = next((h for h in ho if h["fold"] == "calibration -> blind"), None)
    n_place = sum(h["transfers"] for h in ho if h["fold"].startswith("leave out"))
    out["contributes"] = {t: bool(t in win and cb and cb["transfers"] and n_place >= 10) for t in STORM_NAMES}
    out["place_folds_transferring"] = n_place
    split = None
    if out["contributes"].get("gust_front"): split = ctx["gf"] > 0.5
    else:
        c = [t for t in STORM_NAMES if out["contributes"][t] and t != "gust_front"]
        if c:
            key = {"cold_pool": "cold_pool", "pressure_jump": "dp", "dew_drop": "dtd", "dryness": "dry", "heat": "heat"}[c[0]]
            split = ctx[key] > np.median(ctx[key][group == "cal"])
    out["coverage_blind"] = coverage(y, group == "cal", group == "blind", split)
    out["gust_front_vs_not"] = {"n_front": int(ctx["gf"].sum()), "median_ratio_front": float(np.exp(np.median(y[ctx["gf"] > 0.5]))),
                                "median_ratio_other": float(np.exp(np.median(y[ctx["gf"] < 0.5]))),
                                "p10_p90_front": np.exp(np.percentile(y[ctx["gf"] > 0.5], [10, 90])).tolist(),
                                "p10_p90_other": np.exp(np.percentile(y[ctx["gf"] < 0.5], [10, 90])).tolist()}
    print("winner", win, res["winner"]["coef"], "cal->blind", cb, "places", n_place, flush=True)
    print("beyond station", out["savings_beyond_station_bits"], out["coef_beyond_station"], flush=True)
    print("coverage", out["coverage_blind"], out["gust_front_vs_not"], flush=True)
    np.savez(ROOT / "results/physics/storm_arrays.npz", y=y, mu=res["_mu"], gf=ctx["gf"], cold=ctx["cold_pool"], dT=raw["dT"], dp=raw["dp"],
             place=place, group=group, stn_loo=ctx["stn_loo"])
    res_path.write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
