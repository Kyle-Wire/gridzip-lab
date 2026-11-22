"""Shared steps for a physics-gap case: search, coder check, hold-out folds, neighbour test (PROTOCOL-2.md)."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from codelen import lawsearch as L


def strip(d): return {k: v for k, v in d.items() if not k.startswith("_") and k not in ("rounded", "alphas")}


def case(lib, y, g, folds, exclusive=(), max_terms=3, keep=10, n_hold=5):
    top, n = L.search(lib, y, g, max_terms=max_terms, exclusive=exclusive, keep=keep)
    none = L.score(lib, [], y, g)
    hold = []
    for res in top[:n_hold]:
        rows = []
        for name, tr, te in folds:
            bf = L.holdout(lib, res["terms"], y, g, tr, te) / len(te); b0 = L.holdout(lib, [], y, g, tr, te) / len(te)
            rows.append({"fold": name, "formula_bpr": bf, "one_correction_bpr": b0})
        hold.append({"terms": res["terms"], "folds": rows, "transfers": all(r["formula_bpr"] < r["one_correction_bpr"] for r in rows)})
    return {"n": len(y), "candidates": n, "one_correction": strip(none), "top": [strip(t) for t in top],
            "coder_check": L.verify(top[0], y, g), "holdout": hold, "_mu": top[0]["_mu"], "_mu0": none["_mu"]}


def neighbour_test(y, order, ctx, g, perm_groups, n_perm=199, seed=0):
    """ctx[i]: indices (into y) of context for the i-th coded item; shuffles permute y within each perm group."""
    rng = np.random.default_rng(seed)
    obs = L.context_bits(y, order, ctx, g).sum() / len(y); null = np.zeros(n_perm)
    for t in range(n_perm):
        yp = y.copy()
        for grp in perm_groups: yp[grp] = rng.permutation(y[grp])
        null[t] = L.context_bits(yp, order, ctx, g).sum() / len(y)
    return {"bits_per_reading": float(obs), "null_mean": float(null.mean()), "null_sd": float(null.std()),
            "p": float((1 + np.sum(null <= obs)) / (n_perm + 1))}, null
