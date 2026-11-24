"""Protocol 9: do terms from the LBM's own near-ground flow name its lateral Case C pattern?

python physics/lateral_fields.py [--synthetic-only] [--seeds 10]
Reads results/physics/f3s5_fields.npz and f3s5_points.csv (WindSim F3 step 5, Case C 1D, D/32, L 0.5 D).
Writes results/physics/lateral_fields.json and results/physics/lateral_fields_arrays.npz.
"""
import argparse, json, math, sys, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from physlib import engine as E, terms as T, guards as Q  # noqa: E402
from physlib.engine import Term  # noqa: E402

PHYS = ROOT / "results" / "physics"
D, U_APP, Z_EVAL, Z_TOP = 0.2, 2.43, 0.02, 0.1
G = E.Grid(0.01, -4.0, 4.0)
X_UP = -0.7          # reference section: 1 D upstream of the array front (x = -0.5), inside x < -2 D

FIELD = [
    Term("turn", "wind", "b · |sin θ|", "the model's flow at the point turns away from the approach direction", ("turn",), lambda c: c["turn"]),
    Term("reverse", "wind", "δ · [mean u < 0]", "the point sits in the model's recirculation", ("reverse",), lambda c: c["reverse"]),
    Term("ti_local", "wind", "b · σ_u / |U_h|", "the model's local turbulence at the point", ("ti_local",), lambda c: c["ti_local"]),
    Term("downwash", "wind", "b · w(0.05 m) / U_app", "air brought down into the street (negative w) or lifted", ("downwash",), lambda c: c["downwash"]),
    Term("side_flux", "wind", "b · v·sign(y) / U_app", "flow pushed outward, away from the centreline (side channelling)", ("side_flux",), lambda c: c["side_flux"]),
    Term("street_flux", "wind", "b · Q_street / Q_upstream", "the street carries more (or less) air than the same width upstream: blockage / continuity",
         ("street_flux",), lambda c: c["street_flux"]),
]
FIELD_NAMES = [t.name for t in FIELD]
BASE_KEYS = ("lp", "fetch", "x_pos", "lateral", "hc", "corner", "along", "cross")


def interp(field, origin, dx, pts):
    """Trilinear interpolation of a cell-centred (nx, ny, nz, k) field at (m, 3) points."""
    g = (np.asarray(pts, float) - np.asarray(origin, float)) / dx - 0.5
    i0 = np.floor(g).astype(int); t = g - i0; out = 0.0
    for c in np.ndindex(2, 2, 2):
        w = np.prod(np.where(np.array(c), t, 1 - t), axis=1); ii = i0 + np.array(c)
        out = out + w[:, None] * field[ii[:, 0], ii[:, 1], ii[:, 2]]
    return out


def load():
    d = np.load(PHYS / "f3s5_fields.npz")
    st = d["mean_sums"].astype(np.float64)
    st = 0.5 * (st + st[:, ::-1] * np.array([1, -1, 1, 1, 1, 1, 1]))   # mirror fold, as WindSim scores it
    solid = d["solid"] | d["solid"][:, ::-1]
    dx, o = float(d["dx_m"]), d["origin"].astype(float)
    p = np.genfromtxt(PHYS / "f3s5_points.csv", delimiter=",", names=True)
    return st, solid, dx, o, p


def centres(o, dx, shape):
    return [o[k] + (np.arange(shape[k]) + 0.5) * dx for k in range(3)]


def flux_u(st, solid, xc, yc, zc, x, ylo, yhi):
    """Streamwise flux below Z_TOP across y in [ylo, yhi] at x (m^3/s per unit... cell sums times dx^2)."""
    i = int(np.argmin(np.abs(xc - x))); jy = (yc >= ylo) & (yc <= yhi); kz = zc < Z_TOP
    u = st[i][:, :, 0][np.ix_(jy, kz)]; f = ~solid[i][np.ix_(jy, kz)]
    return float((u * f).sum())


def flux_v(st, solid, xc, yc, zc, y, xlo, xhi):
    j = int(np.argmin(np.abs(yc - y))); ix = (xc >= xlo) & (xc <= xhi); kz = zc < Z_TOP
    v = st[:, j][:, :, 1][np.ix_(ix, kz)] * np.sign(y if y != 0 else 1); f = ~solid[:, j][np.ix_(ix, kz)]
    return float((v * f).sum())


def ctx_build():
    st, solid, dx, o, p = load()
    xc, yc, zc = centres(o, dx, solid.shape)
    x, y = p["x"], p["y"]; xyz = np.c_[x, y, p["z"]]; n = len(x)
    f = interp(st, o, dx, xyz); mean, var = f[:, :3], np.maximum(f[:, 3:6] - f[:, :3] ** 2, 0)
    uh = np.hypot(mean[:, 0], mean[:, 1])
    w05 = interp(st[..., 2:3], o, dx, np.c_[x, y, np.full(n, 0.05)])[:, 0]
    ctx = {"turn": np.abs(mean[:, 1]) / np.maximum(uh, 1e-9), "reverse": (mean[:, 0] < 0).astype(float),
           "ti_local": np.sqrt(var.sum(1) / 3) / np.maximum(uh, 1e-9), "downwash": w05 / U_APP,
           "side_flux": mean[:, 1] * np.sign(y) / U_APP}
    sf = np.zeros(n)
    for k in range(n):
        ay, ax = abs(y[k]), abs(x[k])
        if 0.1 - 1e-9 <= ay <= 0.3 + 1e-9:             # along-wind street band (or its intersections / face lines)
            lo, hi = (0.1, 0.3) if y[k] >= 0 else (-0.3, -0.1)
            sf[k] = flux_u(st, solid, xc, yc, zc, x[k], lo, hi) / flux_u(st, solid, xc, yc, zc, X_UP, lo, hi)
        else:                                           # cross street in front of / behind the central block
            lo, hi = (0.1, 0.3) if x[k] >= 0 else (-0.3, -0.1)
            ref = flux_u(st, solid, xc, yc, zc, X_UP, -0.1, 0.1) * (hi - lo) / 0.2
            sf[k] = flux_v(st, solid, xc, yc, zc, y[k], lo, hi) / ref if y[k] != 0 else 0.0
    ctx["street_flux"] = sf
    street = lambda v: 0.1 < round(abs(v), 3) < 0.3
    sx = np.array([street(v) for v in x]); sy = np.array([street(v) for v in y])
    dxb = np.maximum(np.abs(x) - 0.1, 0); dyb = np.maximum(np.abs(y) - 0.1, 0)
    hc = np.ones(n)                                     # 1D layout: every point has the 1 D-taller centre block
    ctx.update({"lp": np.log(p["vs_pred"] / U_APP), "fetch": (x + 0.5) / D, "along": (sy & ~sx).astype(float),
                "cross": (sx & ~sy).astype(float), "hc": hc, "corner": (np.hypot(dxb, dyb) <= 0.1 + 1e-9).astype(float),
                "lateral": np.abs(y) / D, "x_pos": x / D})
    r = np.log(p["vs_pred"] / p["vs_obs"])
    # level ledger
    Q0 = flux_u(st, solid, xc, yc, zc, X_UP, yc.min(), yc.max())
    led = {f"{xs:+.1f}": flux_u(st, solid, xc, yc, zc, xs, yc.min(), yc.max()) / Q0 for xs in (-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6)}
    arr = (xc >= -0.5) & (xc <= 0.5); kt = len(zc) - 1
    wtop = st[:, :, kt, 2][np.ix_(arr, (yc >= -0.5) & (yc <= 0.5))]; ftop = ~solid[:, :, kt][np.ix_(arr, (yc >= -0.5) & (yc <= 0.5))]
    ledger = {"Q_over_Q_upstream_by_x": led, "z_top_layer_m": float(zc[kt]),
              "plan_mean_w_top_over_Uapp": float((wtop * ftop).sum() / ftop.sum() / U_APP),
              "note": "Q counts fluid cells below 0.1 m across the whole 1.8 m section"}
    return ctx, r, x, y, ledger


def folds(x, y):
    q = (np.where(x >= 0, 1, 0) * 2 + np.where(y >= 0, 1, 0)).astype(int); lab = {0: "x<0,y<0", 1: "x<0,y>=0", 2: "x>=0,y<0", 3: "x>=0,y>=0"}
    out = [(f"leave out {lab[k]}", np.where(q != k)[0], np.where(q == k)[0]) for k in range(4)]
    out.append(("upstream -> downstream", np.where(x < 0)[0], np.where(x > 0)[0]))
    return out


def ctx_A(ctx): return {k: ctx[k] for k in BASE_KEYS}


def search(lib, ctx, r, fo=(), full=False, name=""):
    if full:
        return Q.analyse(lib, ctx, r, G, folds=fo, bootstrap_B=30, null_N=19, name=name, boot_by_group=False)
    pool = list(lib) + Q.G.candidates(ctx, len(r))
    return E.run_case(pool, ctx, r, G, folds=fo, name=name)


def by_input_share(res):
    """Bootstrap selection share per input: a resample counts for x if any selected term reads x."""
    built = {b.name: b for b in res["_built"]}; gd = res["guards"]
    return gd.get("bootstrap", {}).get("selected", {}), {k: v for k, v in gd.get("inclusion_by_input", {}).items()}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--synthetic-only", action="store_true"); ap.add_argument("--seeds", type=int, default=10)
    a = ap.parse_args()
    ctx, r, x, y, ledger = ctx_build(); fo = folds(x, y); n = len(r)
    out = {"protocol": "physics/PROTOCOL-9.md", "n": n, "level_ln_mean": float(r.mean()), "level_factor": float(math.exp(r.mean())),
           "ledger": ledger, "field_term_summary": {k: {"mean": float(ctx[k].mean()), "sd": float(ctx[k].std()),
                                                      "min": float(ctx[k].min()), "max": float(ctx[k].max())} for k in FIELD_NAMES}}
    print("ledger", json.dumps(ledger), flush=True)
    # run A (needed for the synthetic noise level)
    t = time.time(); A = search(T.WIND, ctx_A(ctx), r, fo, full=True, name="Case C 1D D/32 L0.5D: FINDINGS-7 inputs")
    print("A", A["winner"]["terms"], round(A["sd_after"], 3), f"{time.time() - t:.0f}s", flush=True)
    sdA = float(A["sd_after"]); muA = A["_mu"]
    # synthetic check
    rng0 = np.random.default_rng(8); syn = {"noise_sd": sdA}
    def named(yy):
        w = search(T.WIND + FIELD, ctx, yy)["winner"]["terms"]
        return w, [t for t in w if t in FIELD_NAMES or any(f in t for f in FIELD_NAMES)]
    picks = []
    for s in range(a.seeds):
        rng = np.random.default_rng(800 + s); w, f = named(muA + rng.normal(0, sdA, n)); picks.append({"winner": w, "field": f})
    syn["a_none"] = {"no_field_term": sum(not p["field"] for p in picks), "of": a.seeds, "runs": picks}
    print("syn a", syn["a_none"]["no_field_term"], flush=True)
    syn["power"] = {}
    for nm in ("side_flux", "street_flux"):
        z = (ctx[nm] - ctx[nm].mean()) / ctx[nm].std(); row = {}
        for sd in (0.02, 0.04, 0.06, 0.08, 0.10):
            hits = 0
            for s in range(a.seeds):
                rng = np.random.default_rng(900 + s); w, f = named(muA + rng.normal(0, sdA, n) + sd * z)
                hits += any(nm in t for t in f)
            row[str(sd)] = hits
        syn["power"][nm] = {"found_by_sd": row, "mde_sd": next((float(k) for k, h in row.items() if h >= 8), None)}
        print("syn power", nm, row, flush=True)
    syn["passed"] = syn["a_none"]["no_field_term"] >= 9
    out["synthetic"] = syn
    out["A"] = _sum(A)
    if a.synthetic_only or not syn["passed"]:
        (PHYS / "lateral_fields.json").write_text(json.dumps(out, indent=1, default=float)); return
    t = time.time(); B = search(T.WIND + FIELD, ctx, r, fo, full=True, name="Case C 1D D/32 L0.5D: + field terms")
    print("B", B["winner"]["terms"], round(B["sd_after"], 3), f"{time.time() - t:.0f}s", flush=True)
    out["B"] = _sum(B)
    built = E.assemble(T.WIND + FIELD, ctx, n); L = E.library_info(built); b0 = E.exact_bits([], r, G, L)["total_bits"]
    out["single_term_savings_bits"] = {b.name: float(b0 - E.exact_bits([b], r, G, L)["total_bits"]) for b in built if b.name in FIELD_NAMES}
    out["single_term_coef"] = {b.name: float(E.fit([b], r)[0][-1]) for b in built if b.name in FIELD_NAMES}
    out["bootstrap_by_input"] = {"A": boot_inputs(A, r), "B": boot_inputs(B, r)}
    out["readings"] = readings(A, B, out["bootstrap_by_input"]["B"])
    (PHYS / "lateral_fields.json").write_text(json.dumps(out, indent=1, default=float))
    np.savez(PHYS / "lateral_fields_arrays.npz", r=r, muA=muA, muB=B["_mu"], x=x, y=y, **{k: ctx[k] for k in FIELD_NAMES})


def _sum(res):
    gd = res["guards"]
    return {"winner": res["winner"]["terms"], "coef": res["winner"]["coef"], "sd_before": res["sd_before"], "sd_after": res["sd_after"],
            "bits_saved_per_reading": res["bits_saved_per_reading"], "holdout": res["holdout"],
            "inclusion": gd.get("inclusion"), "inclusion_by_input": gd.get("inclusion_by_input"),
            "bootstrap": gd.get("bootstrap"), "null_p": gd.get("null_bits_saved", {}).get("p"), "neighbour": gd.get("neighbour"),
            "verdict": Q.verdict(res)}


def _quad_transfers(res):
    return sum(h["transfers"] for h in res["holdout"] if h["fold"].startswith("leave out"))


def _boot_input(res, key):
    built = {b.name: b for b in res["_built"]}; sel = res["guards"].get("bootstrap", {}).get("selected", {})
    # selection shares per term are not additive over terms; the largest share among terms reading `key` is a lower bound
    return max([v for t, v in sel.items() if t in built and key in built[t].term.needs] + [0.0])


def boot_inputs(res, r, B=30, seed=0):
    """The guards' bootstrap replayed (same seed, so the same resamples) to count, per input, the resamples whose
    winner reads that input in any term (a library term or a grammar shape of it)."""
    built = res["_built"]; L = E.library_info(built); rng = np.random.default_rng(seed); by = {b.name: b for b in built}; cnt = {}
    for _ in range(B):
        rows = Q._resample(rng, len(r), None)
        w, _, _ = E.search(built, r, G, library_size=L, rows=rows, robust=True)
        for k in {x for t in w["terms"] for x in by[t].term.needs}: cnt[k] = cnt.get(k, 0) + 1
    return {k: v / B for k, v in sorted(cnt.items(), key=lambda kv: -kv[1])}


def readings(A, B, bootB):
    incA = A["guards"].get("inclusion_by_input", {}); incB = B["guards"].get("inclusion_by_input", {})
    q1 = incA.get("lateral", 0) >= 0.8 and _quad_transfers(A) >= 3
    named = {}
    for f in FIELD_NAMES:
        ok = incB.get(f, 0) >= 0.8 and bootB.get(f, 0) >= 0.6 and _quad_transfers(B) >= 3
        named[f] = {"inclusion": incB.get(f, 0), "bootstrap_share_any_term": bootB.get(f, 0),
                    "bootstrap_share_max_single_term": _boot_input(B, f), "pays": bool(ok),
                    "names_lateral": bool(ok and incB.get("lateral", 0) < 0.5)}
    return {"Q1_lateral_recurs": bool(q1), "lateral_inclusion_A": incA.get("lateral", 0), "lateral_inclusion_B": incB.get("lateral", 0),
            "quadrant_folds_transferring_A": _quad_transfers(A), "quadrant_folds_transferring_B": _quad_transfers(B),
            "Q2_field_terms": named, "Q2_named": [f for f, v in named.items() if v["names_lateral"]],
            "Q2_adds_to": [f for f, v in named.items() if v["pays"] and not v["names_lateral"]]}


if __name__ == "__main__":
    main()
