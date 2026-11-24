"""PROTOCOL-7 guards: how sure is the compressor of what it names, and what could it be confusing it with?

Written fresh for SPARC / gridzip-lab.
  lookalikes      swap each winning term for every candidate that resembles it (|r| > 0.8); any swap within 3.3 bits
                  (odds under 10:1) is listed as "can't rule out"
  posterior       approximate odds over every formula the search scored (weight 2^-bits): inclusion probability per term
  group_confound  a winning term nearly constant within groups is compared with plain group offsets
  bootstrap       how often each term is chosen on resampled data, and a 90 % interval per coefficient
  null_search     the bits a search of this size saves on shuffled readings (the whole search re-run)
  neighbour_test  is structure left that no term named? k-nearest neighbours in input space (the full space and every
                  pair of inputs, the best taken, in the shuffles too) predict the remainder; 99 shuffles
  analyse         all of the above for one case
"""
from __future__ import annotations

import math

import numpy as np

from . import engine as E, grammar as G
from .engine import Built, Term

ODDS_BITS = math.log2(10.0)      # 3.32 bits = 10:1


# ------------------------------------------------------------------------------------------------ helpers
def term_column(b: Built, res) -> np.ndarray:
    """The fitted contribution of one winning term, as one column."""
    if b.fn is not None:
        v = b.fn(res["alpha"][b.name])
        if b.term.fixed: return v.sum(1)
        c = np.atleast_1d(res["coef"].get(b.name, 1.0)); return v @ c
    c = np.atleast_1d(res["coef"].get(b.name, 1.0)); return b.cols @ c if b.cols.shape[1] == len(c) else b.cols[:, 0]


def _unitv(v): v = np.asarray(v, float) - np.mean(v); return v / (np.linalg.norm(v) + 1e-300)


def _cand_corr(b: Built, u):
    """Largest |r| between a candidate and the unit column u (over a family's grid, over a multi-column term's columns)."""
    if b.fn is not None:
        best = 0.0
        for a in b.grid[::2]:
            v = b.fn(a)
            if not np.all(np.isfinite(v)): continue
            for j in range(v.shape[1]): best = max(best, abs(float(_unitv(v[:, j]) @ u)))
        return best
    if b.cols.shape[1] == 1: return abs(float(_unitv(b.cols[:, 0]) @ u))
    Q, _ = np.linalg.qr(b.cols - b.cols.mean(0)); return float(np.linalg.norm(Q.T @ u))


# ------------------------------------------------------------------------------------------------ guards
def lookalikes(built, res, y, g, L, prev=None, robust=False, r_min=0.8, top=8, rows=None):
    by = {b.name: b for b in built}; win = res["terms"]; out = {}
    for t in win:
        u = _unitv(term_column(by[t], res)); others = [n for n in win if n != t]
        sims = []
        for b in built:
            if b.name in win: continue
            st = [by[n] for n in others] + [b]
            grp = [x.term.group for x in st if x.term.group]
            if len(grp) != len(set(grp)): continue
            r = _cand_corr(b, u)
            if r > r_min: sims.append((r, b))
        sims.sort(key=lambda z: -z[0]); found = []
        lib_sims = [z for z in sims if z[1].term.cost_bits is None][:top]; gr_sims = [z for z in sims if z[1].term.cost_bits is not None][:top]
        for r, b in lib_sims + gr_sims:
            ex = E.exact_bits([by[n] for n in others] + [b], y, g, L, rows=rows, prev=prev, robust=robust)
            d = ex["total_bits"] - res["total_bits"]
            if d < ODDS_BITS: found.append({"term": b.name, "r": round(r, 3), "extra_bits": round(d, 2), "grammar": b.term.cost_bits is not None})
        out[t] = sorted(found, key=lambda z: z["extra_bits"])
    return out


def posterior(memo, built, groups_as_one=True):
    st = [s for s in memo if isinstance(s, tuple) and np.isfinite(memo[s])]
    if not st: return {}, {}
    b = np.array([memo[s] for s in st]); w = np.exp2(-(b - b.min())); w /= w.sum()
    name = (lambda i: built[i].term.group or built[i].name) if groups_as_one else (lambda i: built[i].name)
    inc = {}
    for s, wi in zip(st, w):
        for k in {name(i) for i in s}: inc[k] = inc.get(k, 0.0) + float(wi)
    best = st[int(np.argmin(b))]
    inp = {}   # A1: inclusion per input (a grammar shape of x counts for x, as a library term reading x does)
    for s, wi in zip(st, w):
        for k in {x for i in s for x in built[i].term.needs}: inp[k] = inp.get(k, 0.0) + float(wi)
    return dict(sorted(inc.items(), key=lambda kv: -kv[1])), {"map_formula": [built[i].name for i in best], "map_weight": float(w.max()),
                                                              "inclusion_by_input": {k: v for k, v in sorted(inp.items(), key=lambda kv: -kv[1]) if v >= 0.01}}


def _onehot(labels):
    labels = np.asarray(labels); u = sorted(set(labels.tolist()))
    return np.column_stack([(labels == v).astype(float) for v in u[1:]])


def group_confound(built, res, y, g, L, groups, prev=None, robust=False, share=0.05):
    by = {b.name: b for b in built}; groups = np.asarray(groups); out = {}
    if len(set(groups.tolist())) < 2: return out
    off = E.Built(Term("group_offsets", "guard", "δ per group", "plain offsets per group (siting)", (), None), cols=_onehot(groups) - _onehot(groups).mean(0))
    for t in res["terms"]:
        v = term_column(by[t], res); v = v - v.mean(); tot = float(v @ v)
        if tot <= 0: continue
        gm = np.zeros_like(v)
        for k in set(groups.tolist()): m = groups == k; gm[m] = v[m].mean()
        within = float(((v - gm) ** 2).sum()) / tot
        if within >= share: continue
        others = [by[n] for n in res["terms"] if n != t and n != "station_offsets"]
        ex = E.exact_bits(others + [off], y, g, L, prev=prev, robust=robust)
        d = ex["total_bits"] - res["total_bits"]
        # A1: a group-level term's evidence is one value per group. Code the k group means of the remainder without it,
        # with and without the term (approx MDL: k/2 log2 of the SSE ratio against naming + 1 parameter); flag if it
        # does not pay at that level.
        without = E.exact_bits(others, y, g, L, prev=prev, robust=robust) if others else E.exact_bits([], y, g, L, prev=prev, robust=robust)
        e = y - without["_mu"]; ks = sorted(set(groups.tolist())); k = len(ks)
        rg = np.array([e[groups == q].mean() for q in ks]); vg = np.array([v[groups == q].mean() for q in ks])
        rg0 = rg - rg.mean(); X = np.c_[np.ones(k), vg]; e1 = rg - X @ np.linalg.lstsq(X, rg, rcond=None)[0]
        saved = 0.5 * k * math.log2(max(float(rg0 @ rg0), 1e-300) / max(float(e1 @ e1), 1e-300))
        cost = math.log2(E.vocab_size(L[0] if isinstance(L, tuple) else L)) + 1 + 0.5 * math.log2(max(k, 2)) + 2
        pays = saved > cost
        out[t] = {"within_group_share": round(within, 4), "offsets_extra_bits": round(d, 2), "group_level_saved_bits": round(saved, 2),
                  "group_level_cost_bits": round(cost, 2), "flag": bool(d < ODDS_BITS or not pays), "n_groups": k}
    return out


def _resample(rng, n, groups):
    if groups is None: return rng.integers(0, n, n)
    groups = np.asarray(groups); u = sorted(set(groups.tolist())); pick = rng.choice(len(u), len(u))
    return np.concatenate([np.flatnonzero(groups == u[k]) for k in pick])


def bootstrap(built, res, y, g, L, B=30, groups=None, prev=None, robust=False, seed=0, search_kw=None):
    """Selection share of each term (by group) over B re-searches, and the 5-95 % range of each winning coefficient
    with the winner's formula refitted on each resample."""
    rng = np.random.default_rng(seed); by = {b.name: b for b in built}; sel = {}; coefs = {t: [] for t in res["terms"]}
    kw = dict(search_kw or {})
    for _ in range(B):
        rows = _resample(rng, len(y), groups)
        w, _, _ = E.search(built, y, g, library_size=L, rows=rows, prev=None if groups is not None else prev, robust=robust, **kw)
        for k in {by[t].term.group or t for t in w["terms"]}: sel[k] = sel.get(k, 0) + 1
        if res["terms"]:
            th, al, _, _ = E.fit([by[t] for t in res["terms"]], y, rows)
            k = 1
            for t in res["terms"]:
                b = by[t]
                if b.fn is not None and b.term.fixed: coefs[t].append(float(al[res["terms"].index(t)])); continue
                kk = b.k; coefs[t].append(float(th[k]) if kk == 1 else float(np.linalg.norm(th[k:k + kk]))); k += kk
    iv = {t: [float(np.quantile(v, 0.05)), float(np.quantile(v, 0.95))] for t, v in coefs.items() if v}
    return {"B": B, "selected": {k: v / B for k, v in sorted(sel.items(), key=lambda kv: -kv[1])}, "coef_90": iv}


def null_search(built, y, g, L, N=19, blocks=None, prev=None, robust=False, seed=0, search_kw=None):
    rng = np.random.default_rng(seed); out = []
    for _ in range(N):
        yp = y.copy()
        if blocks is None: yp = rng.permutation(y)
        else:
            for bl in blocks: yp[bl] = rng.permutation(y[bl])
        w, ranked, _ = E.search(built, yp, g, library_size=L, prev=prev, robust=robust, **(search_kw or {}))
        none = next(r for r in ranked if not r["terms"]) if any(not r["terms"] for r in ranked) else E.exact_bits([], yp, g, L, prev=prev, robust=robust)
        out.append((none["total_bits"] - w["total_bits"]) / len(y))
    return out


def neighbour_sets(prims, k=8, pairs=True):
    """Neighbour index arrays (n, k) over rank-scaled inputs: the full space, then every pair."""
    keys = list(prims); X = np.column_stack([np.argsort(np.argsort(prims[c])) / (len(prims[c]) - 1) for c in keys])
    sets = [("all inputs", X)]
    if pairs:
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)): sets.append((f"{keys[i]} × {keys[j]}", X[:, [i, j]]))
    out = []
    for nm, Z in sets:
        d = ((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1); np.fill_diagonal(d, np.inf)
        out.append((nm, np.argpartition(d, k, axis=1)[:, :k]))
    return out


def _gain(e, nb):
    e = e - e.mean(); v = float(e @ e)
    return max(((nm, 1.0 - float(((e - e[idx].mean(1)) ** 2).sum()) / v) for nm, idx in nb), key=lambda z: z[1])


def neighbour_test(prims, e, n_perm=99, k=8, blocks=None, seed=0, nb=None):
    rng = np.random.default_rng(seed); nb = nb or neighbour_sets(prims, k)
    where, obs = _gain(e, nb); null = []
    for _ in range(n_perm):
        ep = e.copy()
        if blocks is None: ep = rng.permutation(e)
        else:
            for bl in blocks: ep[bl] = rng.permutation(e[bl])
        null.append(_gain(ep, nb)[1])
    null = np.array(null)
    return {"gain": obs, "where": where, "null_95": float(np.quantile(null, 0.95)), "p": float((1 + np.sum(null >= obs)) / (n_perm + 1)),
            "structure_left": bool((1 + np.sum(null >= obs)) / (n_perm + 1) <= 0.05)}


# ------------------------------------------------------------------------------------------------ one case
def analyse(terms, ctx, y, g, *, grammar=True, robust=True, prev=None, groups=None, folds=(), bootstrap_B=0, null_N=0,
            neighbour=True, replicate_sd=None, name="", exclude=(), blocks=None, seed=0, boot_by_group=True):
    y = np.asarray(y, float); n = len(y)
    pool = list(terms) + (G.candidates(ctx, n, exclude=exclude) if grammar else [])
    memo = {}
    res = E.run_case(pool, ctx, y, g, folds=folds, prev=prev, robust=robust, name=name, memo_out=memo)
    built = res["_built"]; L = E.library_info(built)
    win_full = E.exact_bits([next(b for b in built if b.name == t) for t in res["winner"]["terms"]], y, g, L, prev=prev, robust=robust)
    out = {"lookalikes": lookalikes(built, win_full, y, g, L, prev=prev, robust=robust)}
    inc, mp = posterior(memo, built); out["inclusion"] = {k: v for k, v in inc.items() if v >= 0.01}; out.update(mp)
    if groups is not None: out["group_confound"] = group_confound(built, win_full, y, g, L, groups, prev=prev, robust=robust)
    if bootstrap_B: out["bootstrap"] = bootstrap(built, win_full, y, g, L, B=bootstrap_B, groups=groups if boot_by_group else None, prev=prev, robust=robust, seed=seed)
    if null_N:
        nb = null_search(built, y, g, L, N=null_N, blocks=blocks, prev=prev, robust=robust, seed=seed)
        out["null_bits_saved"] = {"max": max(nb), "mean": float(np.mean(nb)), "all": nb,
                                  "p": float((1 + sum(v >= res["bits_saved_per_reading"] for v in nb)) / (len(nb) + 1))}
    if neighbour:
        prims = G.primitives(ctx, n, exclude)
        if prims: out["neighbour"] = neighbour_test(prims, y - res["_mu"], blocks=blocks, seed=seed)
    out["sd_after"] = res["sd_after"]
    if replicate_sd: out["sd_over_replicate"] = res["sd_after"] / replicate_sd
    res["guards"] = out
    return res


def verdict(res, planted_groups=None):
    """Plain-words reading per winning term: named and sure, named but can't rule out X, or group-level."""
    gd = res["guards"]; out = []
    for t in res["winner"]["terms"]:
        la = gd["lookalikes"].get(t, []); gc = gd.get("group_confound", {}).get(t)
        if gc and gc["flag"]: out.append((t, f"group-level: can't be told from plain offsets over {gc['n_groups']} groups"))
        elif la: out.append((t, "named, but can't rule out " + ", ".join(f"{z['term']} (+{z['extra_bits']:.1f} bits)" for z in la[:3])))
        else: out.append((t, "named; no look-alike within 10:1"))
    return out
