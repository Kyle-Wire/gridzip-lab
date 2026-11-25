"""Protocol 10: does the flow over the roofs name the LBM's Case C pattern or its level? (WindSim F3 step 6)

python physics/roof_fields.py [--synthetic-only] [--seeds 10]
Reads results/physics/f3s6_R{0,1}_fields.npz and f3s6_R{0,1}_points.csv (Case C 1D, D/16, L 0.5 D; R0 standard
section, R1 3.6 x 2.0 m). Writes results/physics/roof_fields.json and roof_fields_arrays.npz.
"""
import argparse, json, math, sys, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "physics"))
from physlib import engine as E, terms as T, guards as Q  # noqa: E402
from physlib.engine import Term  # noqa: E402
import lateral_fields as LF  # noqa: E402  (protocol 9's helpers and field terms, unchanged)

PHYS = ROOT / "results" / "physics"
D, U_APP = LF.D, LF.U_APP
X_APP = -0.8                      # upstream column of the step 6 approach check (x = -0.8, y = 0)
FOLD = np.array([1, -1, 1, 1, 1, 1, 1, -1, 1, -1])   # u v w uu vv ww |u| uv uw vw under y -> -y

ROOF = [
    Term("roof_speed", "wind", "b · u(1.25 D) / U_app(1.25 D)", "the flow just over the roofs runs fast or slow above this street",
         ("roof_speed",), lambda c: c["roof_speed"]),
    Term("roof_stress", "wind", "b · −u'w'(D) / U_ref²", "momentum carried down into the street at roof height",
         ("roof_stress",), lambda c: c["roof_stress"]),
    Term("roof_k", "wind", "b · k(D) / U_ref²", "turbulence energy at roof height", ("roof_k",), lambda c: c["roof_k"]),
    Term("roof_shear", "wind", "b · ∂u/∂z(D) · D / U_ref", "the strength of the shear layer over the street", ("roof_shear",),
         lambda c: c["roof_shear"]),
    Term("coupling", "wind", "b · u(0.02 m) / u(1.25 D)", "how much of the roof flow reaches the probe", ("coupling",),
         lambda c: c["coupling"]),
    Term("ti_full", "wind", "b · sqrt(2k/3) / |U_h| at 0.02 m", "local turbulence from the full stress tensor", ("ti_full",),
         lambda c: c["ti_full"]),
]
ROOF_NAMES = [t.name for t in ROOF]
DIFF_KEYS = ROOF_NAMES + ["street_flux", "side_flux", "downwash"]
DIFF = [Term(f"d_{k}", "wind", f"b · Δ{k} (R0 − R1)", f"the wider section changes {k} at this point", (f"d_{k}",),
             (lambda kk: (lambda c: c[f"d_{kk}"]))(k)) for k in DIFF_KEYS]


def load(run):
    d = np.load(PHYS / f"f3s6_{run}_fields.npz")
    st = d["mean_sums"].astype(np.float64)
    st = 0.5 * (st + st[:, ::-1] * FOLD)                      # mirror fold, as WindSim scores it
    solid = d["solid"] | d["solid"][:, ::-1]
    p = np.genfromtxt(PHYS / f"f3s6_{run}_points.csv", delimiter=",", names=True)
    p = p[np.argsort(p["id"])]
    return st, solid, float(d["dx_m"]), d["origin"].astype(float), p


def at(st, solid, o, dx, pts, cols):
    """Trilinear sample of the given columns; a sample whose nearest cell is solid moves to the nearest fluid cell
    at the same height (protocol 10). Returns values and the indices of moved points."""
    pts = np.asarray(pts, float).copy(); moved = []
    xc, yc, zc = LF.centres(o, dx, solid.shape)
    for k, (x, y, z) in enumerate(pts):
        i, j, kz = (int(np.argmin(np.abs(c - v))) for c, v in ((xc, x), (yc, y), (zc, z)))
        if solid[i, j, kz]:
            fi, fj = np.where(~solid[:, :, kz]); m = np.argmin((xc[fi] - x) ** 2 + (yc[fj] - y) ** 2)
            pts[k, :2] = xc[fi[m]], yc[fj[m]]; moved.append(k)
    return LF.interp(st[..., cols], o, dx, pts), moved


def ctx_build(run):
    st, solid, dx, o, p = load(run)
    xc, yc, zc = LF.centres(o, dx, solid.shape)
    x, y = p["x"], p["y"]; n = len(x)
    col = lambda z: np.c_[x, y, np.full(n, z)]
    # protocol 9 field terms on this run, same definitions
    f, mv0 = at(st, solid, o, dx, np.c_[x, y, p["z"]], list(range(7)))
    mean, var = f[:, :3], np.maximum(f[:, 3:6] - f[:, :3] ** 2, 0); uh = np.hypot(mean[:, 0], mean[:, 1])
    w05 = at(st, solid, o, dx, col(0.05), [2])[0][:, 0]
    ctx = {"turn": np.abs(mean[:, 1]) / np.maximum(uh, 1e-9), "reverse": (mean[:, 0] < 0).astype(float),
           "ti_local": np.sqrt(var.sum(1) / 3) / np.maximum(uh, 1e-9), "downwash": w05 / U_APP,
           "side_flux": mean[:, 1] * np.sign(y) / U_APP}
    sf = np.zeros(n)
    for k in range(n):
        ay = abs(y[k])
        if 0.1 - 1e-9 <= ay <= 0.3 + 1e-9:
            lo, hi = (0.1, 0.3) if y[k] >= 0 else (-0.3, -0.1)
            sf[k] = LF.flux_u(st, solid, xc, yc, zc, x[k], lo, hi) / LF.flux_u(st, solid, xc, yc, zc, LF.X_UP, lo, hi)
        else:
            lo, hi = (0.1, 0.3) if x[k] >= 0 else (-0.3, -0.1)
            ref = LF.flux_u(st, solid, xc, yc, zc, LF.X_UP, -0.1, 0.1) * (hi - lo) / 0.2
            sf[k] = LF.flux_v(st, solid, xc, yc, zc, y[k], lo, hi) / ref if y[k] != 0 else 0.0
    ctx["street_flux"] = sf
    # roof terms
    app = lambda z, c=0: float(LF.interp(st[..., [c]], o, dx, np.array([[X_APP, 0.0, z]]))[0, 0])
    u_ref = app(D)
    u125, mv1 = at(st, solid, o, dx, col(1.25 * D), [0]); u075, mv2 = at(st, solid, o, dx, col(0.75 * D), [0])
    r1, mv3 = at(st, solid, o, dx, col(D), [0, 1, 2, 3, 4, 5, 8])
    uw = r1[:, 6] - r1[:, 0] * r1[:, 2]; k1 = 0.5 * np.maximum(r1[:, 3:6] - r1[:, :3] ** 2, 0).sum(1)
    ctx.update({"roof_speed": u125[:, 0] / app(1.25 * D), "roof_stress": -uw / u_ref ** 2, "roof_k": k1 / u_ref ** 2,
                "roof_shear": (u125[:, 0] - u075[:, 0]) / (0.5 * D) * D / u_ref,
                "coupling": mean[:, 0] / u125[:, 0],
                "ti_full": np.sqrt(2 * 0.5 * var.sum(1) / 3) / np.maximum(uh, 1e-9)})
    street = lambda v: 0.1 < round(abs(v), 3) < 0.3
    sx = np.array([street(v) for v in x]); sy = np.array([street(v) for v in y])
    dxb = np.maximum(np.abs(x) - 0.1, 0); dyb = np.maximum(np.abs(y) - 0.1, 0)
    ctx.update({"lp": np.log(p["vs_pred"] / U_APP), "fetch": (x + 0.5) / D, "along": (sy & ~sx).astype(float),
                "cross": (sx & ~sy).astype(float), "hc": np.ones(n), "corner": (np.hypot(dxb, dyb) <= 0.1 + 1e-9).astype(float),
                "lateral": np.abs(y) / D, "x_pos": x / D})
    r = np.log(p["vs_pred"] / p["vs_obs"])
    # Q3 ledger
    Q0 = LF.flux_u(st, solid, xc, yc, zc, LF.X_UP, yc.min(), yc.max())
    arr = (xc >= -0.5) & (xc <= 0.5); arry = (yc >= -0.5) & (yc <= 0.5); kD = int(np.argmin(np.abs(zc - D)))
    wD = st[:, :, kD, 2][np.ix_(arr, arry)]; fD = ~solid[:, :, kD][np.ix_(arr, arry)]
    ledger = {"Q_over_Q_upstream": {nm: LF.flux_u(st, solid, xc, yc, zc, xs, yc.min(), yc.max()) / Q0
                                    for nm, xs in (("front", -0.4), ("middle", 0.0), ("back", 0.4))},
              "array_mean_roof_speed": float(ctx["roof_speed"].mean()), "array_mean_roof_stress": float(ctx["roof_stress"].mean()),
              "plan_mean_w_at_D_over_Uapp": float((wD * fD).sum() / fD.sum() / U_APP), "u_ref_upstream_at_D": u_ref,
              "approach_minus_uw_over_uref2": {f"{z / D:.2f}D": -(app(z, 8) - app(z, 0) * app(z, 2)) / u_ref ** 2
                                               for z in (0.5 * D, D, 1.25 * D, 1.5 * D, 2 * D)},
              "moved_samples": {"z0.02": len(mv0), "1.25D": len(mv1), "0.75D": len(mv2), "D": len(mv3)}}
    return ctx, r, x, y, p["id"], ledger


def named_rule(res, boot, keys):
    inc = res["guards"].get("inclusion_by_input", {}); q = LF._quad_transfers(res)
    return {k: {"inclusion": inc.get(k, 0), "bootstrap_share_any_term": boot.get(k, 0), "quadrant_folds": q,
                "pays": bool(inc.get(k, 0) >= 0.8 and boot.get(k, 0) >= 0.6 and q >= 3)} for k in keys}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--synthetic-only", action="store_true"); ap.add_argument("--seeds", type=int, default=10)
    a = ap.parse_args()
    c0, y0, x, y, id0, led0 = ctx_build("R0"); c1, y1, _, _, id1, led1 = ctx_build("R1")
    assert (id0 == id1).all(); n = len(y0); fo = LF.folds(x, y)
    out = {"protocol": "physics/PROTOCOL-10.md", "n": n, "level": {"R0": math.exp(y0.mean()), "R1": math.exp(y1.mean())},
           "ledger": {"R0": led0, "R1": led1},
           "roof_term_summary": {k: {"mean": float(c0[k].mean()), "sd": float(c0[k].std()), "min": float(c0[k].min()),
                                     "max": float(c0[k].max())} for k in ROOF_NAMES}}
    print("ledger", json.dumps(out["ledger"], default=float), flush=True)
    FIELD = LF.FIELD
    t = time.time(); A = LF.search(T.WIND + FIELD, {k: c0[k] for k in LF.BASE_KEYS + tuple(LF.FIELD_NAMES)}, y0, fo, full=True,
                                   name="Case C 1D D/16 R0: FINDINGS-7 inputs + protocol 9 terms")
    print("A", A["winner"]["terms"], round(A["sd_after"], 4), f"{time.time() - t:.0f}s", flush=True)
    sdA, muA = float(A["sd_after"]), A["_mu"]
    syn = {"noise_sd": sdA}
    def named(yy):
        w = LF.search(T.WIND + FIELD + ROOF, c0, yy)["winner"]["terms"]
        return w, [t for t in w if any(f in t for f in ROOF_NAMES)]
    picks = []
    for s in range(a.seeds):
        rng = np.random.default_rng(1000 + s); w, f = named(muA + rng.normal(0, sdA, n)); picks.append({"winner": w, "roof": f})
    syn["a_none"] = {"no_roof_term": sum(not p["roof"] for p in picks), "of": a.seeds, "runs": picks}
    print("syn a", syn["a_none"]["no_roof_term"], flush=True)
    syn["power"] = {}
    for nm in ("roof_stress", "roof_speed"):
        z = (c0[nm] - c0[nm].mean()) / c0[nm].std(); row = {}
        for sd in (0.02, 0.04, 0.06, 0.08, 0.10):
            hits = 0
            for s in range(a.seeds):
                rng = np.random.default_rng(1100 + s); w, f = named(muA + rng.normal(0, sdA, n) + sd * z); hits += any(nm in t for t in f)
            row[str(sd)] = hits
        syn["power"][nm] = {"found_by_sd": row, "mde_sd": next((float(k) for k, h in row.items() if h >= 8), None)}
        print("syn power", nm, row, flush=True)
    syn["passed"] = syn["a_none"]["no_roof_term"] >= 9
    out["synthetic"] = syn; out["A"] = LF._sum(A)
    if a.synthetic_only:
        (PHYS / "roof_fields.json").write_text(json.dumps(out, indent=1, default=float)); return
    t = time.time(); B = LF.search(T.WIND + FIELD + ROOF, c0, y0, fo, full=True, name="Case C 1D D/16 R0: + roof terms")
    print("B", B["winner"]["terms"], round(B["sd_after"], 4), f"{time.time() - t:.0f}s", flush=True)
    out["B"] = LF._sum(B)
    built = E.assemble(T.WIND + FIELD + ROOF, c0, n); L = E.library_info(built); b0 = E.exact_bits([], y0, E.Grid(0.01, -4, 4), L)["total_bits"]
    out["single_term_savings_bits"] = {b.name: float(b0 - E.exact_bits([b], y0, LF.G, L)["total_bits"]) for b in built if b.name in ROOF_NAMES}
    out["single_term_coef"] = {b.name: float(E.fit([b], y0)[0][-1]) for b in built if b.name in ROOF_NAMES}
    bA, bB = LF.boot_inputs(A, y0), LF.boot_inputs(B, y0); out["bootstrap_by_input"] = {"A": bA, "B": bB}
    incA = A["guards"].get("inclusion_by_input", {}); incB = B["guards"].get("inclusion_by_input", {})
    roof = named_rule(B, bB, ROOF_NAMES)
    for k in roof: roof[k]["names_lateral"] = bool(roof[k]["pays"] and incB.get("lateral", 0) < 0.5)
    out["Q1"] = {"lateral_recurs_at_D16": bool(incA.get("lateral", 0) >= 0.8 and LF._quad_transfers(A) >= 3),
                 "lateral_inclusion_A": incA.get("lateral", 0), "lateral_inclusion_B": incB.get("lateral", 0),
                 "quadrant_folds_A": LF._quad_transfers(A), "quadrant_folds_B": LF._quad_transfers(B), "roof_terms": roof,
                 "named": [k for k, v in roof.items() if v["names_lateral"]], "adds_to": [k for k, v in roof.items() if v["pays"] and not v["names_lateral"]]}
    # Q2: what the wider section changes
    dy = y0 - y1; cd = dict(c0); cd.update({f"d_{k}": c0[k] - c1[k] for k in DIFF_KEYS})
    t = time.time(); C = LF.search(T.WIND + FIELD + ROOF + DIFF, cd, dy, fo, full=True, name="dy = y(R0) - y(R1)")
    print("C", C["winner"]["terms"], round(C["sd_after"], 4), f"{time.time() - t:.0f}s", flush=True)
    bC = LF.boot_inputs(C, dy); keysC = ROOF_NAMES + LF.FIELD_NAMES + [f"d_{k}" for k in DIFF_KEYS] + list(LF.BASE_KEYS)
    nrC = named_rule(C, bC, keysC); p = C["guards"].get("null_bits_saved", {}).get("p", 1.0)
    out["Q2"] = {"dy_mean": float(dy.mean()), "dy_sd": float(dy.std()), "null_p": p, "search": LF._sum(C), "bootstrap_by_input": bC,
                 "named": [k for k, v in nrC.items() if v["pays"]],
                 "reading": "uniform shift" if (p > 0.05 or not [k for k, v in nrC.items() if v["pays"]]) else "structured: blockage footprint"}
    out["Q3"] = {"level_R0": out["level"]["R0"], "level_R1": out["level"]["R1"],
                 "roof_speed_R0": led0["array_mean_roof_speed"], "roof_speed_R1": led1["array_mean_roof_speed"],
                 "roof_stress_R0": led0["array_mean_roof_stress"], "roof_stress_R1": led1["array_mean_roof_stress"]}
    (PHYS / "roof_fields.json").write_text(json.dumps(out, indent=1, default=float))
    np.savez(PHYS / "roof_fields_arrays.npz", y0=y0, y1=y1, muA=muA, muB=B["_mu"], muC=C["_mu"], x=x, y=y,
             **{f"R0_{k}": c0[k] for k in ROOF_NAMES + LF.FIELD_NAMES}, **{f"R1_{k}": c1[k] for k in ROOF_NAMES + LF.FIELD_NAMES})
    print("done", json.dumps({"Q1": out["Q1"]["named"], "adds": out["Q1"]["adds_to"], "Q2": out["Q2"]["reading"]}), flush=True)


if __name__ == "__main__":
    main()
