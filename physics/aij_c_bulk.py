"""Protocol 5: rank the suspects for WindSim's ~1.3x street excess on AIJ Case C 1D (physics/PROTOCOL-5.md).

python physics/aij_c_bulk.py [--synthetic-only] [--windsim ../WindSim]
Reads WindSim's committed per-point and profile files (bf9e352); writes results/physics/aij_c_bulk.json and
physics/figures/bulk_*.png. CPU only, a few minutes.
"""
import argparse, json, math, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "physics"))
from physlib import engine as E, terms as T  # noqa: E402
from physlib.engine import Term  # noqa: E402
from physlib.cases import neighbour, spatial_ctx  # noqa: E402
import aij_c as A  # noqa: E402

D, U_APP = 0.2, 2.43
G = E.Grid(0.01, -4.0, 4.0)
RUNS = {  # folder: settings (cells per D, tunnel z0, walls, inflow length 0.5 D, Cs 0.17)
    "d16":          dict(n=16, z0=0, walls=0, L05=0, cs17=0),
    "d16_z0tunnel": dict(n=16, z0=1, walls=0, L05=0, cs17=0),
    "d20":          dict(n=20, z0=0, walls=0, L05=0, cs17=0),
    "d24":          dict(n=24, z0=0, walls=0, L05=0, cs17=0),
    "d32_z0tunnel": dict(n=32, z0=1, walls=0, L05=0, cs17=0),
    "d32_walls":    dict(n=32, z0=1, walls=1, L05=0, cs17=0),
    "d32_L0.5D":    dict(n=32, z0=1, walls=0, L05=1, cs17=0),
    "d32_cs0.17":   dict(n=32, z0=1, walls=0, L05=0, cs17=1),
}
D32 = [k for k in RUNS if k.startswith("d32")]
BLOCKS = [(bx, by) for bx in (-0.4, 0.0, 0.4) for by in (-0.4, 0.0, 0.4)]
SD_FLOOR = 0.063          # protocol 1's mirror-pair sd (ln)

S5_NAMES = ("contrast", "fetch_decay", "fetch_linear", "canyon_along", "canyon_cross", "corner", "contrast_along", "fetch_decay_along")
LATERAL = Term("lateral", "wind", "b · |y| / D", "outer streets (nearer the model's side boundaries) are faster: a domain or blockage effect",
               ("lat",), lambda c: c["lat"])
NEAR_WALL = Term("near_wall", "wind", "A · exp(−d_wall / L)", "the miss is concentrated next to building faces (wall drag)",
                 ("dwall",), lambda c: (lambda L: np.exp(-c["dwall"] / L)), family=(0.05, 2.0, True), unit="L (D)")
POOL = [Term("point_pattern", "wind", "δ · (the same point's mean miss in the other runs)",
             "the pattern is fixed by something every run shares (measurements, inflow, domain, scoring)", ("pp",), lambda c: c["pp"]),
        Term("res", "wind", "b · ln(cells per D)", "the level moves with grid resolution", ("res",), lambda c: c["res"]),
        Term("z0_tunnel", "wind", "δ · [tunnel ground z0]", "the ground roughness setting moves the level", ("z0",), lambda c: c["z0"]),
        Term("walls", "wind", "δ · [building wall model]", "the wall model moves the level", ("walls",), lambda c: c["walls"]),
        Term("L05", "wind", "δ · [inflow length 0.5 D]", "the inflow eddy size moves the level", ("L05",), lambda c: c["L05"]),
        Term("cs17", "wind", "δ · [Cs 0.17]", "the subgrid constant moves the level", ("cs17",), lambda c: c["cs17"])]


def s5_terms():
    return [t for t in T.WIND if t.name in S5_NAMES] + [LATERAL, NEAR_WALL]


# ------------------------------------------------------------------------------------------------- data
def load(vdir):
    base = vdir / "lbm_case_c_1D_resolution"; runs = {}
    for k in RUNS:
        d = np.genfromtxt(base / k / "points_1D.csv", delimiter=",", names=True); o = np.argsort(d["id"])
        runs[k] = {c: d[c][o] for c in ("id", "x", "y", "vs_obs", "vs_pred")}
    ids = runs["d32_z0tunnel"]["id"]
    for k, v in runs.items():
        assert np.array_equal(v["id"], ids) and np.allclose(v["vs_obs"], runs["d32_z0tunnel"]["vs_obs"]), k
    prof = {k: json.loads((base / k / "profiles_1D.json").read_text()) for k in D32 if (base / k / "profiles_1D.json").exists()}
    rd = np.genfromtxt(vdir / "rans_case_c_d20" / "results_points.csv", delimiter=",", names=True); o = np.argsort(rd["id"])
    rans = {c: rd[c][o] for c in ("id", "x", "y", "vs_obs", "pred_vector")}
    return runs, prof, rans


def geometry(x, y):
    sx = np.array([A.street(v) for v in x]); sy = np.array([A.street(v) for v in y])
    dwall = np.full(len(x), np.inf)
    for bx, by in BLOCKS:
        dx = np.maximum(np.abs(x - bx) - D / 2, 0); dy = np.maximum(np.abs(y - by) - D / 2, 0)
        dwall = np.minimum(dwall, np.hypot(dx, dy))
    ddx = np.maximum(np.abs(x) - 0.1, 0); ddy = np.maximum(np.abs(y) - 0.1, 0)
    return {"fetch": (x + 0.5) / D, "along": (sy & ~sx).astype(float), "cross": (sx & ~sy).astype(float),
            "corner": (np.hypot(ddx, ddy) <= 0.1 + 1e-9).astype(float), "lat": np.abs(y) / D, "dwall": dwall / D}


# ------------------------------------------------------------------------------------------------- tests
def s5(x, y, r, pred, name):
    ctx = geometry(x, y) | {"lp": np.log(pred / U_APP)}
    folds = [("fit y<=0, code y>0", np.where(y <= 0)[0], np.where(y > 0)[0]), ("fit y>=0, code y<0", np.where(y >= 0)[0], np.where(y < 0)[0])]
    return E.run_case(s5_terms(), ctx, r, G, folds, max_terms=4, name=name)


def single_term_savings(terms, ctx, r):
    built = E.assemble(terms, ctx, len(r)); b0 = E.exact_bits([], r, G, len(built))["total_bits"]
    return {b.name: float(b0 - E.exact_bits([b], r, G, len(built))["total_bits"]) for b in built}


def pool(runs, r_by):
    keys = list(RUNS); n = len(r_by[keys[0]]); R = np.array([r_by[k] for k in keys])
    y, ctx = [], {k: [] for k in ("pp", "res", "z0", "walls", "L05", "cs17")}; run_of = []
    for j, k in enumerate(keys):
        others = np.delete(R, j, axis=0).mean(0)
        y.append(R[j]); ctx["pp"].append(others); ctx["res"].append(np.full(n, math.log(RUNS[k]["n"])))
        for s in ("z0", "walls", "L05", "cs17"): ctx[s].append(np.full(n, float(RUNS[k][s])))
        run_of.append(np.full(n, j))
    y = np.concatenate(y); ctx = {k: np.concatenate(v) for k, v in ctx.items()}; run_of = np.concatenate(run_of)
    folds = [(f"leave out {k}", np.where(run_of != j)[0], np.where(run_of == j)[0]) for j, k in enumerate(keys)]
    return y, ctx, folds, run_of


def ledgers(prof, level):
    out = {"level_ln": level, "level_factor": math.exp(-level)}
    up02, upD, slopes, over = [], [], {}, []
    for k, P in prof.items():
        z = np.array(P["z_m"]); cols = P["columns"]; up = np.array(cols["upstream_x-0.8"])
        up02.append(float(np.interp(0.02, z, up))); upD.append(float(np.interp(D, z, up)))
        i = int(np.searchsorted(z, 0.02)); lo, hi = i - 1, i
        for c in ("street_behind_first_row_x-0.2", "side_street_x0_y0.2"):
            u = np.array(cols[c]); slopes.setdefault(c, []).append(float(np.log(u[hi] / u[lo]) / np.log(z[hi] / z[lo])))
        am = np.array(P["array_mean"]); slopes.setdefault("array_mean", []).append(float(np.log(am[hi] / am[lo]) / np.log(z[hi] / z[lo])))
        sel = (z >= 1.25 * D - 1e-9) & (z <= 2 * D + 1e-9)
        over.append(float(np.mean(np.log(np.array(cols["side_street_x0_y0.2"])[sel] / up[sel]))))
    tot = -level
    s1 = math.log(np.mean(up02) / U_APP)
    half = 0.5 * D / 32; smax = max(abs(np.mean(v)) for v in slopes.values())
    s2 = smax * math.log((0.02 + half) / 0.02)
    s3_geo = -math.log(1 - 3 * D * D / (1.8 * 1.2))
    s3_over = float(np.mean(over))
    out.update({
        "S1_inflow": {"U_model_up_0.02": up02, "U_measured_0.02": U_APP, "ln": s1, "share": s1 / tot, "U_model_up_D": upD},
        "S2_height": {"slope_dlnU_dlnz_at_0.02": {k: float(np.mean(v)) for k, v in slopes.items()}, "height_error_m": half,
                      "ln_bound": s2, "share_bound": s2 / tot},
        "S3_domain": {"frontal_blockage": 3 * D * D / (1.8 * 1.2), "ln_blockage_bound": s3_geo, "share_blockage_bound": s3_geo / tot,
                      "over_array_ln_vs_upstream_1.25_to_2D": s3_over, "share_over_array": s3_over / tot}})
    out["made_inside_array_share"] = 1 - out["S1_inflow"]["share"] - out["S3_domain"]["share_blockage_bound"]
    return out


def verdict_share(s):
    return "ruled out" if s < 0.10 else ("a part" if s <= 0.50 else "a main cause")


# ------------------------------------------------------------------------------------------------- synthetic check
def synthetic(x, y, pred, seeds=10):
    rng0 = np.random.default_rng(5); out = {}
    ctx = geometry(x, y) | {"lp": np.log(pred / U_APP)}
    r0 = lambda s: -0.26 + np.random.default_rng(100 + s).normal(0, SD_FLOOR, len(x))
    picks = [sorted(E.run_case(s5_terms(), ctx, r0(s), G, (), max_terms=4)["winner"]["terms"]) for s in range(seeds)]
    out["a uniform only"] = {"hits": sum(w == [] for w in picks), "of": seeds, "winners": picks}
    shapes = {"lateral": ctx["lat"], "near_wall": np.exp(-ctx["dwall"] / 0.25), "fetch_decay": np.exp(-ctx["fetch"] / 1.0),
              "canyon_along": ctx["along"]}
    out["power"] = {}
    for nm, v in shapes.items():
        v = (v - v.mean()) / v.std(); row = {}
        for sd in (0.01, 0.015, 0.02, 0.03, 0.04, 0.06):
            h = 0
            for s in range(seeds):
                res = E.run_case(s5_terms(), ctx, r0(s) + sd * v, G, (), max_terms=4)
                h += nm in res["winner"]["terms"]
            row[str(sd)] = h
        mde = next((float(k) for k, h in row.items() if h >= 8), None)
        out["power"][nm] = {"found_of_10_by_sd": row, "mde_sd_ln": mde}
    # d: the 8-run pool, a shared point pattern plus a resolution shift
    hits = 0; picks = []
    for s in range(seeds):
        rng = np.random.default_rng(200 + s); pat = rng.normal(0, 0.10, len(x)) + 0.05 * rng0.normal(size=len(x))
        r_by = {k: -0.26 + pat + 0.06 / math.log(2) * (math.log(v["n"]) - math.log(32)) + rng.normal(0, 0.04, len(x))
                for k, v in RUNS.items()}
        yy, cc, _, _ = pool(None, r_by)
        res = E.run_case(POOL, cc, yy, G, (), max_terms=3); w = sorted(res["winner"]["terms"]); picks.append(w)
        hits += set(["point_pattern", "res"]) <= set(w) and len(w) == 2
    out["d pool: pattern + resolution"] = {"hits": int(hits), "of": seeds, "winners": picks}
    out["passed"] = out["d pool: pattern + resolution"]["hits"] >= 8 and out["a uniform only"]["hits"] >= 9
    return out


# ------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--windsim", default=str(ROOT.parent / "WindSim"))
    ap.add_argument("--synthetic-only", action="store_true"); a = ap.parse_args()
    runs, prof, rans = load(Path(a.windsim) / "docs" / "validation")
    x, y, m = runs["d32_z0tunnel"]["x"], runs["d32_z0tunnel"]["y"], runs["d32_z0tunnel"]["vs_obs"]
    p32 = np.exp(np.mean([np.log(runs[k]["vs_pred"]) for k in D32], axis=0))
    out = {"protocol": "physics/PROTOCOL-5.md", "windsim_files": "bf9e352"}
    syn = synthetic(x, y, p32); out["synthetic"] = syn
    print("synthetic", json.dumps({k: (v.get("hits") if "hits" in v else v) if isinstance(v, dict) else v for k, v in syn.items()}, default=str), flush=True)
    if a.synthetic_only or not syn["passed"]:
        (ROOT / "results/physics/aij_c_bulk.json").write_text(json.dumps(out, indent=1, default=float)); return
    r_by = {k: np.log(m / runs[k]["vs_pred"]) for k in RUNS}
    out["run_levels"] = {k: {"median_ratio_model_over_measured": float(np.exp(-np.median(v))), "mean_ln": float(v.mean()), **RUNS[k]} for k, v in r_by.items()}
    r32 = np.mean([r_by[k] for k in D32], axis=0); level = float(np.mean(r32))
    # ledgers
    out["ledgers"] = led = ledgers(prof, level)
    # S4: the 8-run pool
    yy, cc, folds, run_of = pool(runs, r_by)
    rp = E.run_case(POOL, cc, yy, G, folds, max_terms=4, name="AIJ C 1D, 8 LBM runs pooled")
    out["S4_pool"] = E.public(rp); out["S4_pool"]["single_term_savings_bits"] = single_term_savings(POOL, cc, yy)
    print("S4", rp["winner"]["terms"], rp["winner"]["coef"], rp["transfers"], flush=True)
    # S5 on the D/32 mean
    r5 = s5(x, y, r32, p32, "AIJ C 1D D/32, mean of 4 runs")
    out["S5_d32"] = E.public(r5)
    out["S5_d32"]["single_term_savings_bits"] = single_term_savings(s5_terms(), geometry(x, y) | {"lp": np.log(p32 / U_APP)}, r32)
    print("S5", r5["winner"]["terms"], r5["winner"]["coef"], r5["transfers"], flush=True)
    # S5 on each D/32 run alone (is the answer the same run by run?)
    out["S5_each"] = {k: {"winner": (rr := s5(x, y, r_by[k], runs[k]["vs_pred"], k))["winner"]["terms"], "transfers": rr["transfers"],
                          "bits_saved_per_reading": rr["bits_saved_per_reading"]} for k in D32}
    # S7: neighbours on S5's remainder
    order, cx = spatial_ctx(np.c_[x, y], 0.075)
    rem = r5["_y"] - r5["_mu"]
    out["S7_neighbours"] = {"before": neighbour(r32 - r32.mean(), order, cx, G, [np.arange(len(x))], n_perm=199),
                            "after_S5": neighbour(rem, order, cx, G, [np.arange(len(x))], n_perm=199)}
    # S6: LBM vs RANS point patterns
    assert np.array_equal(rans["id"], runs["d32_z0tunnel"]["id"])
    rr_ = np.log(rans["vs_obs"] / rans["pred_vector"]); pl, pr = r32 - r32.mean(), rr_ - rr_.mean()
    c = float(np.corrcoef(pl, pr)[0, 1]); shared_noise = SD_FLOOR ** 2 / (pl.std() * pr.std())
    # the LBM pattern's own repeatability: mean correlation between D/32 runs' patterns
    pats = [r_by[k] - r_by[k].mean() for k in D32]
    rep = float(np.mean([np.corrcoef(pats[i], pats[j])[0, 1] for i in range(4) for j in range(i + 1, 4)]))
    out["S6_rans"] = {"corr_lbm_rans": c, "corr_from_shared_measurement_noise_alone": float(shared_noise),
                      "corr_between_d32_runs": rep, "rans_level_factor": float(np.exp(-rr_.mean())), "lbm_level_factor": float(np.exp(-level)),
                      "sd_pattern_lbm": float(pl.std()), "sd_pattern_rans": float(pr.std())}
    print("S6", out["S6_rans"], flush=True)
    out["ranking"] = rank(out)
    np.savez(ROOT / "results/physics/aij_c_bulk_arrays.npz", x=x, y=y, r32=r32, mu=r5["_mu"], rans=rr_,
             **{f"r_{k}": v for k, v in r_by.items()})
    (ROOT / "results/physics/aij_c_bulk.json").write_text(json.dumps(out, indent=1, default=float))
    for row in out["ranking"]: print(row)


def rank(o):
    L = o["ledgers"]; rows = []
    rows.append({"suspect": "S1 approach flow near the ground", "evidence": f"model upstream {np.mean(L['S1_inflow']['U_model_up_0.02']):.2f} vs measured 2.43 m/s",
                 "share": L["S1_inflow"]["share"], "verdict": verdict_share(L["S1_inflow"]["share"])})
    rows.append({"suspect": "S2 evaluation height 0.02 m", "evidence": f"slope ≤ {max(abs(v) for v in L['S2_height']['slope_dlnU_dlnz_at_0.02'].values()):.2f}, half-cell error",
                 "share": L["S2_height"]["share_bound"], "verdict": verdict_share(L["S2_height"]["share_bound"])})
    s5w = o["S5_d32"]["winner"]["terms"]; lat_ok = "lateral" in s5w and o["S5_d32"]["transfers"]
    rows.append({"suspect": "S3 domain blockage / fixed top", "evidence": f"frontal blockage bound; lateral term {'pays' if lat_ok else 'does not pay'}",
                 "share": L["S3_domain"]["share_blockage_bound"], "verdict": verdict_share(L["S3_domain"]["share_blockage_bound"]) + " (upper bound)"})
    pw = o["S4_pool"]["winner"]["terms"]; tr = o["S4_pool"]["transfers"]
    for t in ("res", "z0_tunnel", "walls", "L05", "cs17"):
        rows.append({"suspect": f"S4 setting: {t}", "evidence": f"{o['S4_pool']['single_term_savings_bits'].get(t, 0):.1f} bits alone",
                     "share": None, "verdict": "contributes" if (t in pw and tr) else "does not pay"})
    rows.append({"suspect": "S5 street-scale mechanism", "evidence": f"winner {s5w or 'one factor only'}",
                 "share": None, "verdict": "contributes" if (s5w and o["S5_d32"]["transfers"]) else "does not pay"})
    s6 = o["S6_rans"]
    rows.append({"suspect": "S6 shared with RANS (measurements, probe, geometry)",
                 "evidence": f"pattern corr {s6['corr_lbm_rans']:.2f} (noise alone {s6['corr_from_shared_measurement_noise_alone']:.2f}); levels LBM ×{s6['lbm_level_factor']:.2f}, RANS ×{s6['rans_level_factor']:.2f}",
                 "share": None, "verdict": "pattern shared" if s6["corr_lbm_rans"] > 0.5 * s6["corr_between_d32_runs"] else "not shared"})
    s7 = o["S7_neighbours"]["after_S5"]
    rows.append({"suspect": "S7 unnamed structure", "evidence": f"{s7['structure_left_bits']:.2f} bit/point left, p {s7['p']:.3f}",
                 "share": None, "verdict": "present" if s7["p"] < 0.01 else "none detected"})
    rows.append({"suspect": "made inside the array (what is left)", "evidence": "1 − S1 − S3 bound", "share": L["made_inside_array_share"], "verdict": "the main cause"})
    return rows


if __name__ == "__main__":
    main()
