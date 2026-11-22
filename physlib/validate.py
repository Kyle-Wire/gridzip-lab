"""Validate the physics compressor on planted laws and on pure noise (PROTOCOL-4, part A).

python -m physlib.validate [--quick]      writes results/physlib/synthetic.json
Scores, per planted law over seeds: named exactly (the winner's terms equal the planted set, a declared alias group
counting as a match), law explained (the winner's fitted correction reproduces the planted signal, R^2 >= 0.9),
false terms per run. On noise: the share of runs that name any term at all (the false-discovery rate), for
independent noise and for AR(1) noise (rho 0.8), which breaks the code's independence assumption on purpose.
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from . import engine as E, synth as S

OUT = Path(__file__).resolve().parents[1] / "results" / "physlib"


def _groups(terms):
    return {t.name: (t.group or t.name) for t in terms}


def one(args):
    domain, law, seed, scale, noise = args
    terms, ctx, y, g, planted = S.case(domain, law, seed, scale=scale, noise=noise)
    res = E.run_case(terms, ctx, y, E.Grid(*g), name=f"{domain}:{law}:{seed}")
    found = res["winner"]["terms"]; grp = _groups(terms)
    named = sorted(grp[t] for t in found) == sorted(grp[t] for t in planted)
    rng = np.random.default_rng(seed); sig = (y - 0.1)  # noise included; recompute the clean planted signal
    clean = S.plant(domain, ctx, S.SCENARIOS[domain]["laws"][law]) * scale if law else np.zeros(len(y))
    fitted = res["_mu"] - res["_mu"].mean()
    c = clean - clean.mean()
    r2 = 1 - float(np.sum((c - fitted) ** 2) / np.sum(c ** 2)) if np.sum(c ** 2) > 0 else None
    false_terms = [t for t in found if grp[t] not in {grp[p] for p in planted}]
    return {"domain": domain, "law": law, "seed": seed, "scale": scale, "noise": noise, "planted": planted, "found": found,
            "alpha": res["winner"]["alpha"], "named": bool(named), "r2": r2, "false_terms": false_terms,
            "bits_saved": res["bits_saved_per_reading"], "n": res["n"], "library_size": res["library_size"]}


def jobs(quick=False):
    R, RN, RA, RS = (4, 10, 6, 3) if quick else (20, 100, 50, 10)
    J = []
    for d, sc in S.SCENARIOS.items():
        for law in sc["laws"]: J += [(d, law, s, 1.0, "iid") for s in range(R)]
        J += [(d, None, 1000 + s, 1.0, "iid") for s in range(RN)]
        J += [(d, None, 2000 + s, 1.0, "ar1") for s in range(RA)]
        first = next(iter(sc["laws"]))
        J += [(d, first, 3000 + s, k, "iid") for k in (0.1, 0.2, 0.35, 0.5, 0.75, 1.0) for s in range(RS)]
    return J


def summarise(rows):
    out = {"planted": {}, "noise": {}, "sweep": {}}
    for d in S.SCENARIOS:
        for law in S.SCENARIOS[d]["laws"]:
            rr = [r for r in rows if r["domain"] == d and r["law"] == law and r["scale"] == 1.0 and r["seed"] < 1000]
            if not rr: continue
            out["planted"][f"{d}: {law}"] = {"runs": len(rr), "named_exactly": float(np.mean([r["named"] for r in rr])),
                                             "law_explained": float(np.mean([(r["r2"] or 0) >= 0.9 for r in rr])),
                                             "median_r2": float(np.median([r["r2"] or 0 for r in rr])),
                                             "false_terms_per_run": float(np.mean([len(r["false_terms"]) for r in rr])),
                                             "found_sets": _count([tuple(r["found"]) for r in rr])}
        for nz in ("iid", "ar1"):
            rr = [r for r in rows if r["domain"] == d and r["law"] is None and r["noise"] == nz]
            if rr: out["noise"][f"{d}: {nz}"] = {"runs": len(rr), "any_term": float(np.mean([len(r["found"]) > 0 for r in rr])),
                                                 "terms_named": _count([t for r in rr for t in r["found"]])}
        rr = [r for r in rows if r["domain"] == d and r["seed"] >= 3000]
        if rr:
            out["sweep"][d] = {"law": rr[0]["law"], "points": [
                {"scale": k, "named": float(np.mean([r["named"] for r in rr if r["scale"] == k])),
                 "explained": float(np.mean([(r["r2"] or 0) >= 0.9 for r in rr if r["scale"] == k]))}
                for k in sorted({r["scale"] for r in rr})]}
    return out


def _count(xs):
    c = {}
    for x in xs: k = " + ".join(x) if isinstance(x, tuple) else x; k = k or "(nothing)"; c[k] = c.get(k, 0) + 1
    return dict(sorted(c.items(), key=lambda kv: -kv[1]))


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    quick = "--quick" in sys.argv; t = time.time(); J = jobs(quick)
    with Pool() as p: rows = p.map(one, J, chunksize=4)
    OUT.mkdir(parents=True, exist_ok=True)
    summ = summarise(rows)
    (OUT / ("synthetic_quick.json" if quick else "synthetic.json")).write_text(json.dumps({"summary": summ, "runs": rows}, indent=1, default=float))
    print(json.dumps(summ, indent=1, default=float)[:6000]); print(f"{len(J)} runs, {time.time() - t:.0f} s")
