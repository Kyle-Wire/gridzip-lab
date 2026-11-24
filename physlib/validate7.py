"""PROTOCOL-7 parts A to C: the v2 engine on planted laws, on noise, and against its guards' failure modes.

python -m physlib.validate7 [suite ...] [--quick]     suites: A Areg B C1 C2 C3 C4 C5 (default: all)
Writes results/physlib/v7_<suite>.json (every run and the summary against the protocol's rules).
"""
import json, math, sys, time, warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from . import engine as E, synth as S, terms as T, grammar as G, guards as Q
from .engine import Term

OUT = Path(__file__).resolve().parents[1] / "results" / "physlib"
WG = (0.01, -4, 4)
QUICK = "--quick" in sys.argv or "--dev" in sys.argv
DEV = "--dev" in sys.argv          # development smoke runs on seeds + 90000, never the protocol's seeds
SEED_AT = {"A": 1, "Areg": 2, "B": 1, "C1": 2, "C2": 1, "C3": 0, "C4": 1, "C5": 1}


def _col(term, ctx, alpha=None):
    f = term.build(ctx); v = np.asarray(f(alpha) if term.family else f, float)
    if v.ndim == 2: v = v[:, 0]
    return v - v.mean()


def plant_snr(cols, weights, sd, snr=1.0):
    """Unit-sd columns times weights, the sum scaled so its sd is snr x noise sd. Returns (signal, scale per column)."""
    U = [c / c.std() for c in cols]; s = sum(w * u for w, u in zip(weights, U)); k = snr * sd / s.std()
    return k * s, [k * w / c.std() for w, c in zip(weights, cols)]


def _named(found, planted, lib):
    grp = {t.name: (t.group or t.name) for t in lib}
    gg = lambda t: grp.get(t, t)
    return sorted(gg(t) for t in found) == sorted(gg(t) for t in planted)


def _r2(clean, mu):
    c = clean - clean.mean(); f = mu - mu.mean(); return 1 - float(((c - f) ** 2).sum() / (c ** 2).sum())


def _false(found, planted, lib):
    grp = {t.name: (t.group or t.name) for t in lib}; pg = {grp.get(p, p) for p in planted}
    return [t for t in found if grp.get(t, t) not in pg]


# ================================================================================================== A
A_LAWS = {
    "A1 cold pool": ("conv", [("cold_pool", None, 1.0)]),
    "A2 cold pool + pressure jump": ("conv", [("cold_pool", None, 1.0), ("pressure_jump", None, 1.0)]),
    "A3 MO stability": ("conv", [("mo_stability", None, 1.0)]),
    "A4 downdraft + station offset": ("conv", [("downdraft", None, 1.0), ("station_bias", None, 1.0)]),
    "A5 decay after landfall": ("trop", [("land_decay", 0.095, 1.0)]),
    "A6 motion asymmetry": ("trop", [("motion_asymmetry", None, 1.0)]),
}


def run_A(args):
    law, seed = args; kind, spec = A_LAWS[law]; rng = np.random.default_rng(seed)
    ctx = S.conv_ctx(rng) if kind == "conv" else S.trop_ctx(rng); n = len(ctx["hour"])
    lib = T.WIND_FULL; by = {t.name: t for t in lib}
    cols = [_col(by[nm], ctx, a) for nm, a, _ in spec]; sig, k = plant_snr(cols, [w for *_, w in spec], 0.15)
    y = 0.1 + sig + rng.normal(0, 0.15, n)
    r = E.run_case(lib, ctx, y, E.Grid(*WG), robust=True)
    found = r["winner"]["terms"]; planted = [nm for nm, _, _ in spec]
    return {"suite": "A", "law": law, "seed": seed, "found": found, "named": _named(found, planted, lib), "r2": _r2(sig, r["_mu"]),
            "false": _false(found, planted, lib), "alpha": r["winner"]["alpha"], "library_size": r["library_size"], "code": r["winner"]["code"]}


# ================================================================================================== A-reg
def _extend(dom, ctx, rng, n):
    if dom == "wind":   # convective and tropical inputs with no effect, so the full wind library is in play
        c2 = S.conv_ctx(rng, n); c3 = S.trop_ctx(rng, n)
        for k, v in {**c3, **c2}.items():
            if k not in ctx: ctx[k] = v
    if dom == "heat":
        ctx.update(zi=np.exp(rng.normal(np.log(1200), 0.3, n)), adv=rng.normal(0, 0.5, n), swc=rng.uniform(0.1, 0.4, n))
    return ctx


def run_Areg(args):
    dom, law, seed = args
    terms, ctx, y, g, planted = S.case(dom, law, seed)
    rng = np.random.default_rng(10_000 + seed); ctx = _extend(dom, ctx, rng, len(y))
    lib = T.WIND_FULL if dom == "wind" else (T.HEAT + T.HEAT_V2 if dom == "heat" else terms)
    r = E.run_case(lib, ctx, y, E.Grid(*g), prev=S.prev_of(dom, len(y)), robust=True)
    found = r["winner"]["terms"]
    return {"suite": "Areg", "law": f"{dom}: {law}", "seed": seed, "found": found, "named": _named(found, planted, lib),
            "false": _false(found, planted, lib), "library_size": r["library_size"]}


# ================================================================================================== B
B_LAWS = ("B1 hc·lp", "B2 corner·decay(fetch)", "B3 Td·cloud (heat)", "B4 hinge(z)", "B5 off-grid decay", "B6 control")


def run_B(args):
    law, seed = args; rng = np.random.default_rng(seed)
    if law.startswith("B3"):
        ctx = S.heat_ctx(rng); n = len(ctx["hour"]); lib = list(T.HEAT); sd = 0.8; g = (0.1, -20, 20); prev = S.prev_of("heat", n)
        col = ctx["Td"] * ctx["cloud"]
    elif law.startswith("B6"):
        terms, ctx, y, g, planted = S.case("wind", "contrast + entrance decay", seed)
        r = E.run_case(list(terms) + G.candidates(ctx, len(y)), ctx, y, E.Grid(*g), robust=True)
        found = r["winner"]["terms"]
        return {"suite": "B", "law": law, "seed": seed, "found": found, "named": _named(found, planted, terms), "false": _false(found, planted, terms)}
    else:
        ctx = S.wind_ctx(rng); n = len(ctx["lp"]); lib = list(T.WIND); sd = 0.15; g = WG; prev = None
        f = ctx["fetch"]; lmed = float(np.quantile(f[f > 0], 0.5))
        col = {"B1": ctx["hc"] * ctx["lp"], "B2": ctx["corner"] * np.exp(-f / lmed), "B4": np.maximum(0, ctx["z"] - np.quantile(ctx["z"], 0.5)),
               "B5": np.exp(-f / 3.1)}[law[:2]]
        if law.startswith("B5"): lib = [t for t in lib if t.name != "fetch_decay"]
    sig, _ = plant_snr([col - col.mean()], [1.0], sd); y = 0.1 + sig + rng.normal(0, sd, n)
    r = E.run_case(lib + G.candidates(ctx, n), ctx, y, E.Grid(*g), robust=True, prev=prev)
    built = {b.name: b for b in r["_built"]}; found = r["winner"]["terms"]
    u = Q._unitv(col); best = 0.0
    for t in found:
        b = built[t]
        v = Q.term_column(b, r["winner"]) if b.fn is not None else b.cols[:, 0]
        best = max(best, abs(float(Q._unitv(v) @ u)))
    return {"suite": "B", "law": law, "seed": seed, "found": found, "best_r": best, "r2": _r2(sig, r["_mu"]),
            "grammar_size": r["grammar_size"]}


# ================================================================================================== C1
def run_C1(args):
    dom, noise, seed, n_small = args; rng = np.random.default_rng(seed)
    if n_small:
        ctx = S.wind_ctx(rng, n=n_small); y = 0.1 + rng.normal(0, 0.15, n_small); terms = T.WIND; g = WG; prev = None
    else:
        terms, ctx, y, g, _ = S.case(dom, None, seed, noise="ar1" if noise == "ar1" else "iid")
        prev = S.prev_of(dom, len(y))
        if noise == "t3":
            sd = S.SCENARIOS[dom]["sd"]; y = 0.1 + sd * rng.standard_t(3, len(y)) / math.sqrt(3.0)
    r = E.run_case(list(terms) + G.candidates(ctx, len(y)), ctx, y, E.Grid(*g), robust=True, prev=prev)
    return {"suite": "C1", "set": f"{dom} {noise}" + (f" n={n_small}" if n_small else ""), "seed": seed, "found": r["winner"]["terms"],
            "code": r["winner"].get("code")}


# ================================================================================================== C2
def run_C2(args):
    rho, seed = args; rng = np.random.default_rng(seed)
    ctx = S.wind_ctx(rng); n = len(ctx["lp"]); a = rng.normal(0, 1, n); b = rho * a + math.sqrt(1 - rho ** 2) * rng.normal(0, 1, n)
    ctx.update(ta=a, tb=b)
    lib = list(T.WIND) + [Term("twin_a", "wind", "b · A", "planted", ("ta",), lambda c: c["ta"]),
                          Term("twin_b", "wind", "b · B", "look-alike of A", ("tb",), lambda c: c["tb"])]
    sig, _ = plant_snr([a - a.mean()], [1.0], 0.15); y = 0.1 + sig + rng.normal(0, 0.15, n)
    r = Q.analyse(lib, ctx, y, E.Grid(*WG), neighbour=False)
    found = r["winner"]["terms"]; gd = r["guards"]; la = gd["lookalikes"]; inc = gd["inclusion"]
    ii = gd["inclusion_by_input"]; bt = {b.name: b for b in r["_built"]}
    reads = lambda t, k: k in bt[t].term.needs
    chosen = [t for t in found if reads(t, "ta") or reads(t, "tb")]
    cant_tell = bool(chosen) and any(any(reads(z["term"], "tb" if reads(c, "ta") else "ta") for z in la.get(c, [])) or
                                     ii.get("ta" if reads(c, "ta") else "tb", 0) < 0.9 for c in chosen)
    return {"suite": "C2", "rho": rho, "seed": seed, "found": found, "lookalikes": la, "inclusion_a": ii.get("ta", 0.0),
            "inclusion_b": ii.get("tb", 0.0), "cant_tell": cant_tell, "a_sure": any(reads(t, "ta") for t in found) and ii.get("ta", 0) >= 0.9}


# ================================================================================================== C3
def run_C3(args):
    (seed,) = args; rng = np.random.default_rng(seed)
    ctx = S.wind_ctx(rng); n = len(ctx["lp"]); stn = np.repeat(np.arange(6), n // 6)
    off = rng.normal(0, 0.15, 6); attr = rng.normal(0, 1, 6)
    ctx.update(station=stn, site_attr=attr[stn])
    lib = list(T.WIND) + [Term("station_offsets", "wind", "δ_s per station", "siting offset per station", ("station",), lambda c: T._onehot(c["station"])),
                          Term("site_attr", "wind", "b · a_s", "a station-level attribute (no effect planted)", ("site_attr",), lambda c: c["site_attr"])]
    y = 0.1 + off[stn] + rng.normal(0, 0.15, n)
    r = Q.analyse(lib, ctx, y, E.Grid(*WG), groups=stn, neighbour=False, exclude=("station",))
    found = r["winner"]["terms"]; gc = r["guards"].get("group_confound", {}); bt = {b.name: b for b in r["_built"]}
    false = [t for t in found if t != "station_offsets"]
    unflagged = [t for t in false if not gc.get(t, {}).get("flag", False)]
    attr = [t for t in found if "site_attr" in bt[t].term.needs]       # A1: the library term or any grammar shape of it
    return {"suite": "C3", "seed": seed, "found": found, "group_confound": gc, "attr_named": bool(attr),
            "attr_flagged": bool(attr) and all(gc.get(t, {}).get("flag", False) for t in attr), "unflagged_false": unflagged}


# ================================================================================================== C4
def run_C4(args):
    n, seed = args; rng = np.random.default_rng(seed); ctx = S.conv_ctx(rng, n); lib = T.WIND_FULL; by = {t.name: t for t in lib}
    col = _col(by["cold_pool"], ctx); sig, k = plant_snr([col], [1.0], 0.15); y = 0.1 + sig + rng.normal(0, 0.15, n)
    B = 30 if n == 120 else 0
    r = Q.analyse(lib, ctx, y, E.Grid(*WG), bootstrap_B=B, neighbour=False, seed=seed)
    found = r["winner"]["terms"]; out = {"suite": "C4", "n": n, "seed": seed, "found": found, "named": _named(found, ["cold_pool"], lib),
                                         "false": _false(found, ["cold_pool"], lib), "true_coef": k[0]}
    if B and "cold_pool" in found:
        iv = r["guards"]["bootstrap"]["coef_90"]["cold_pool"]; out.update(coef_90=iv, covered=bool(iv[0] <= k[0] <= iv[1]),
                                                                        selected=r["guards"]["bootstrap"]["selected"])
    return out


# ================================================================================================== C5
def run_C5(args):
    kind, seed = args; rng = np.random.default_rng(seed)
    if kind == "a":
        terms, ctx, y, g, planted = S.case("wind", "contrast + entrance decay", seed)
    else:
        ctx = S.wind_ctx(rng); n = len(ctx["lp"]); terms = T.WIND; g = WG
        if kind == "b":
            col = np.sin(3 * ctx["lp"]) * np.cos(ctx["fetch"]); sig, _ = plant_snr([col - col.mean()], [1.0], 0.15)
        else: sig = np.zeros(n)
        y = 0.1 + sig + rng.normal(0, 0.15, n)
    r = Q.analyse(terms, ctx, y, E.Grid(*g), seed=seed)
    nb = r["guards"]["neighbour"]
    return {"suite": "C5", "kind": kind, "seed": seed, "found": r["winner"]["terms"], "p": nb["p"], "where": nb["where"],
            "structure_left": nb["structure_left"], "sd_after": r["sd_after"], "sd_ratio": r["sd_after"] / 0.15}


# ================================================================================================== jobs, summary
def jobs(suite):
    R = 4 if QUICK else 20
    if suite == "A": return run_A, [(l, s) for l in A_LAWS for s in range(R)]
    if suite == "Areg": return run_Areg, [(d, l, s) for d, sc in S.SCENARIOS.items() for l in sc["laws"] for s in range(R)]
    if suite == "B": return run_B, [(l, s) for l in B_LAWS for s in range(R)]
    if suite == "C1":
        k = (lambda m: 3 if QUICK else m)
        J = [(d, "iid", 1000 + s, 0) for d in ("wind", "heat", "flood", "ground") for s in range(k(50))]
        J += [(d, "ar1", 2000 + s, 0) for d in ("heat", "ground") for s in range(k(30))]
        J += [(d, "t3", 4000 + s, 0) for d in ("wind", "flood") for s in range(k(50))]
        J += [("wind", "iid", 5000 + s, 40) for s in range(k(50))]
        return run_C1, J
    if suite == "C2": return run_C2, [(rho, s) for rho in (0.998, 0.97, 0.6) for s in range(R)]
    if suite == "C3": return run_C3, [(s,) for s in range(R)]
    if suite == "C4": return run_C4, [(n, s) for n in (30, 60, 120) for s in range(R)]
    if suite == "C5": return run_C5, [(k, (6000 if k == "a" else 7000 if k == "b" else 8000) + s) for k in "abc" for s in range(R)]
    raise KeyError(suite)


def frac(xs): xs = list(xs); return float(np.mean(xs)) if xs else float("nan")


def summarise(suite, rows):
    P4 = json.loads((OUT / "synthetic.json").read_text())["summary"]["planted"] if (OUT / "synthetic.json").exists() else {}
    s = {}
    if suite == "A":
        for law in A_LAWS:
            rr = [r for r in rows if r["law"] == law]; nm, ex = frac(r["named"] for r in rr), frac(r["r2"] >= 0.9 for r in rr)
            s[law] = {"named": nm, "explained": ex, "pass": nm >= 0.8 and ex >= 0.8, "found": _count(rr)}
    elif suite == "Areg":
        for law in sorted({r["law"] for r in rows}):
            rr = [r for r in rows if r["law"] == law]; nm = frac(r["named"] for r in rr); was = P4.get(law, {}).get("named_exactly")
            s[law] = {"named": nm, "protocol4": was, "pass": was is None or nm >= was - 0.10 - 1e-9, "found": _count(rr), "library_size": rr[0]["library_size"]}
    elif suite == "B":
        for law in B_LAWS:
            rr = [r for r in rows if r["law"] == law]
            if law.startswith("B6"):
                nm = frac(r["named"] for r in rr); s[law] = {"named": nm, "pass": nm >= 0.8, "found": _count(rr)}
            elif law.startswith("B5"):
                ex = frac(r["r2"] >= 0.8 for r in rr); s[law] = {"explained_0.8": ex, "pass": ex >= 0.7, "found": _count(rr)}
            else:
                hit, ex = frac(r["best_r"] >= 0.95 for r in rr), frac(r["r2"] >= 0.9 for r in rr)
                s[law] = {"shape_found": hit, "explained": ex, "pass": hit >= 0.7 and ex >= 0.8, "found": _count(rr)}
    elif suite == "C1":
        for st in sorted({r["set"] for r in rows}):
            rr = [r for r in rows if r["set"] == st]; a = frac(len(r["found"]) > 0 for r in rr)
            s[st] = {"runs": len(rr), "any_term": a, "pass": a <= 0.05, "found": _count(rr), "codes": _count([(r["code"],) for r in rr], raw=True)}
    elif suite == "C2":
        rr = [r for r in rows if r["rho"] == 0.998]; ct = frac(r["cant_tell"] for r in rr)
        s["twins r=0.998"] = {"says_cant_tell": ct, "pass": ct >= 0.8, "found": _count(rr)}
        rr = [r for r in rows if r["rho"] == 0.97]
        s["pair r=0.97 (reported)"] = {"says_cant_tell": frac(r["cant_tell"] for r in rr), "a_named_sure": frac(r["a_sure"] for r in rr), "found": _count(rr)}
        rr = [r for r in rows if r["rho"] == 0.6]; su = frac(r["a_sure"] for r in rr)
        s["pair r=0.6"] = {"a_named_sure": su, "pass": su >= 0.8, "found": _count(rr)}
    elif suite == "C3":
        named = [r for r in rows if r["attr_named"]]; fl = frac(r["attr_flagged"] for r in named) if named else 1.0
        uf = frac(len(r["unflagged_false"]) > 0 for r in rows)
        s["group confounding"] = {"attr_named_runs": len(named), "flagged_when_named": fl, "unflagged_false_share": uf,
                                  "pass": fl >= 0.9 and uf <= 0.05, "found": _count(rows)}
    elif suite == "C4":
        for n in (30, 60, 120):
            rr = [r for r in rows if r["n"] == n]; fpr = frac(len(r["false"]) for r in rr)
            d = {"false_per_run": fpr, "named": frac(r["named"] for r in rr), "pass_false": fpr <= 0.2, "found": _count(rr)}
            if n == 120:
                cv = [r["covered"] for r in rr if "covered" in r]; d.update(coverage=frac(cv), coverage_runs=len(cv), pass_coverage=(frac(cv) >= 0.8) if cv else None)
            s[f"n={n}"] = d
        s["pass"] = all(v["pass_false"] for k, v in s.items() if k.startswith("n=")) and s["n=120"].get("pass_coverage") is not False
    elif suite == "C5":
        a = [r for r in rows if r["kind"] == "a"]; b = [r for r in rows if r["kind"] == "b"]; c = [r for r in rows if r["kind"] == "c"]
        s["a law named, no structure left"] = {"flagged": frac(r["structure_left"] for r in a), "pass": frac(r["structure_left"] for r in a) <= 0.10}
        s["b form outside grammar"] = {"flagged": frac(r["structure_left"] for r in b), "pass": frac(r["structure_left"] for r in b) >= 0.70,
                                       "where": _count([(r["where"],) for r in b], raw=True), "found": _count(b)}
        s["c noise: sd within 10 %"] = {"share": frac(abs(r["sd_ratio"] - 1) <= 0.10 for r in c), "pass": frac(abs(r["sd_ratio"] - 1) <= 0.10 for r in c) >= 0.9}
    return s


def _count(rr, raw=False):
    c = {}
    for r in rr:
        k = " + ".join(r) if raw else (" + ".join(r["found"]) or "(nothing)"); c[k] = c.get(k, 0) + 1
    return dict(sorted(c.items(), key=lambda kv: -kv[1])[:8])


def _safe(o):
    if isinstance(o, dict): return {str(k): _safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [_safe(v) for v in o]
    if isinstance(o, (np.floating, np.integer)): return o.item()
    if isinstance(o, np.bool_): return bool(o)
    if isinstance(o, np.ndarray): return o.tolist()
    return o


def main(suites):
    warnings.filterwarnings("ignore"); OUT.mkdir(parents=True, exist_ok=True)
    for su in suites:
        fn, J = jobs(su); t = time.time()
        if DEV: J = [tuple(v + 90000 if i == SEED_AT[su] else v for i, v in enumerate(j)) for j in J]
        with Pool(4) as p: rows = p.map(fn, J, chunksize=1)
        summ = summarise(su, rows)
        (OUT / f"v7_{su}{'_dev' if DEV else '_quick' if QUICK else ''}.json").write_text(json.dumps(_safe({"summary": summ, "runs": rows}), indent=1))
        print(f"== {su}: {len(J)} runs, {time.time() - t:.0f} s"); print(json.dumps(_safe(summ), indent=1)[:4000], flush=True)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    main(args or ["A", "Areg", "B", "C1", "C2", "C3", "C4", "C5"])
