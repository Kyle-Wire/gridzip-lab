"""The physics compressor's engine: sparse formula search over a curated term library, scored by total code length.

Written fresh for SPARC / gridzip-lab (no lpaq1, cmix or other GPL code). It generalises `codelen.lawsearch` in three
ways: a term may carry several columns (station offsets, a harmonic pair), the library is named physics with its
formula and units rather than a per-case list, and the search is a beam over sparse combinations instead of every
combination of three, so a library of 20 to 40 terms stays affordable.

Code length of a formula (two-part MDL, the same parts as codelen.mdl):
  structure   Elias-delta(number of tokens + 1) + log2(|library| + 2) per token. A bigger library makes every
              term dearer to name: that is the built-in multiple-testing penalty.
  parameters  every coefficient and family parameter rounded to a step 2**-p, p chosen per parameter to minimise
              the total (zigzag Elias-gamma for p and for the integer), plus the residual sd.
  residuals   the readings binned on the grid, under N(formula, sd), with codelen's escape mixture.
`approx_bits` (search) replaces the parameter part with 0.5 log2(n) + 2 bits per parameter and uses unrounded
values; finalists are always re-scored exactly with `exact_bits`, and the winner is checked with the arithmetic
coder (`verify`).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np

from codelen import mdl as M
from codelen import scoring as sc
from codelen.lawsearch import Grid, bits_gaussian  # noqa: F401  (re-exported)


# --------------------------------------------------------------------------------------------------- the library
@dataclass
class Term:
    """One named piece of physics a model's errors might follow.

    build(ctx) returns an (n,) or (n, k) array, or for a one-parameter family a function a -> array.
    """
    name: str
    domain: str
    formula: str                  # what the term is, in symbols
    meaning: str                  # what it says about the model when it wins, in plain words
    needs: tuple                  # ctx keys it reads
    build: Callable
    family: Optional[tuple] = None   # (lo, hi, log) for a one-parameter family
    fixed: bool = False           # family enters with coefficient 1 (a law with one unknown), not a fitted slope
    group: Optional[str] = None   # terms in one group say the same thing; at most one per formula
    physics: bool = True          # False for shape-only terms (harmonics with no physics name)
    unit: str = ""                # unit of the family parameter


@dataclass
class Built:
    term: Term
    cols: Optional[np.ndarray] = None     # (n, k) for a plain term
    fn: Optional[Callable] = None         # a -> (n, k) for a family
    grid: Optional[np.ndarray] = None

    @property
    def name(self): return self.term.name

    @property
    def k(self):
        return (0 if self.term.fixed else self.fn(self.grid[0]).shape[1]) if self.fn else self.cols.shape[1]


def _as2d(v):
    v = np.asarray(v, np.float64); return v[:, None] if v.ndim == 1 else v


def assemble(terms: Sequence[Term], ctx: dict, n: int) -> list[Built]:
    """Build every term whose inputs are in ctx; drop terms that are constant or duplicate an earlier one."""
    out, seen = [], []
    for t in terms:
        if any(k not in ctx or ctx[k] is None for k in t.needs): continue
        try:
            if t.family:
                lo, hi, lg = t.family
                grid = np.geomspace(lo, hi, 33) if lg else np.linspace(lo, hi, 33)
                f = t.build(ctx)
                fn = (lambda f: (lambda a: _as2d(f(a))))(f)
                probe = fn(grid[len(grid) // 2])
                if probe.shape[0] != n or not np.all(np.isfinite(probe)): continue
                out.append(Built(t, fn=fn, grid=grid))
            else:
                cols = _as2d(t.build(ctx))
                if cols.shape[0] != n or not np.all(np.isfinite(cols)): continue
                cols = cols - cols.mean(0)                # the intercept carries the level
                keep = cols.std(0) > 1e-12 * (1 + np.abs(cols).max())
                if not keep.any(): continue
                cols = cols[:, keep]
                sig = np.round(cols[:, 0] / (np.linalg.norm(cols[:, 0]) + 1e-300), 9)
                if any(np.allclose(sig, s) or np.allclose(sig, -s) for s in seen): continue
                seen.append(sig); out.append(Built(t, cols=cols))
        except (ValueError, FloatingPointError, ZeroDivisionError, IndexError, TypeError):
            continue
    return out


# --------------------------------------------------------------------------------------------------- fitting
def _design(built: list[Built], alphas, rows):
    cols, off = [np.ones((len(rows), 1))], np.zeros(len(rows))
    for b, a in zip(built, alphas):
        if b.fn is not None:
            v = b.fn(a)[rows]
            if b.term.fixed: off = off + v.sum(1)
            else: cols.append(v)
        else:
            cols.append(b.cols[rows])
    return np.hstack(cols), off


def _lsq(X, y):
    th, *_ = np.linalg.lstsq(X, y, rcond=None); return th


def fit(built: list[Built], y, rows=None, refine=True):
    """Least squares; each family's parameter by grid search, then golden-section refinement (coordinate passes)."""
    rows = np.arange(len(y)) if rows is None else np.asarray(rows)
    yr = y[rows]
    fam = [j for j, b in enumerate(built) if b.fn is not None]
    al = [None] * len(built)
    for j in fam: al[j] = float(built[j].grid[len(built[j].grid) // 2])

    def sse(a):
        X, off = _design(built, a, rows); z = yr - off; th = _lsq(X, z); e = z - X @ th; return float(e @ e)
    if len(fam) > 1:   # A1: start the families one at a time, each at its best given those already started
        done = [j for j, b in enumerate(built) if b.fn is None]
        for j in fam:
            sub = done + [j]; gj = built[j].grid; best = None
            for v in gj:
                a = [al[i] for i in sub]; a[-1] = float(v)
                X1, o1 = _design([built[i] for i in sub], a, rows); z = yr - o1; e = z - X1 @ _lsq(X1, z); q = float(e @ e)
                if best is None or q < best[0]: best = (q, float(v))
            al[j] = best[1]; done.append(j)
    for _ in range(2 if len(fam) > 1 else 1):
        for j in fam:
            g = built[j].grid; s = []
            for v in g:
                a = al[:]; a[j] = float(v); s.append(sse(a))
            k = int(np.argmin(s)); al[j] = float(g[k])
            if refine:
                lg = built[j].term.family[2]; tr, it = (math.log, math.exp) if lg else ((lambda v: v), (lambda v: v))
                u, w = tr(g[max(k - 1, 0)]), tr(g[min(k + 1, len(g) - 1)])
                for _ in range(24):
                    m1, m2 = u + 0.382 * (w - u), u + 0.618 * (w - u)
                    a1, a2 = al[:], al[:]; a1[j], a2[j] = it(m1), it(m2)
                    if sse(a1) < sse(a2): w = m2
                    else: u = m1
                al[j] = it((u + w) / 2)
    X, off = _design(built, al, rows); th = _lsq(X, yr - off)
    return th, al, X, off


def vocab_size(library_size): return library_size + 2       # terms + intercept + '+'


def n_tokens(built): return 1 + 2 * len(built)               # c, then '+ term' per term


def structure_bits(built, library_size):
    return M.delta_bits(n_tokens(built) + 1) + n_tokens(built) * math.log2(vocab_size(library_size))


def n_params(built): return 1 + sum(b.k for b in built) + sum(1 for b in built if b.fn is not None)


def _prev_in(prev, rows):
    """Map each row's previous-reading index into positions within `rows` (-1 where it is not in rows)."""
    if prev is None: return None
    pos = np.full(len(prev), -1); pos[rows] = np.arange(len(rows)); p = np.asarray(prev)[rows]
    return np.where(p >= 0, pos[np.maximum(p, 0)], -1)


def serial_bits(yr, mu, pr, g, rho=None, sd=None):
    """A1: residuals coded as AR(1) along `pr` (position of each reading's predecessor, -1 for the first of a block).

    Each reading is coded given the one before it: N(mu + rho (y_prev - mu_prev), sd), the first of each block under the
    stationary sd / sqrt(1 - rho^2). The decoder has the earlier readings, so this is a valid sequential code. Returns
    (bits, rho, sd)."""
    e = yr - mu; has = pr >= 0; ep = np.where(has, e[np.maximum(pr, 0)], 0.0)
    if rho is None:
        rho = float(np.clip((e[has] @ ep[has]) / max(ep[has] @ ep[has], 1e-300), -0.99, 0.99)) if has.sum() > 2 else 0.0
    pred = mu + rho * ep
    if sd is None: sd = max(float(np.sqrt(np.mean((yr - pred)[has] ** 2))) if has.any() else float(np.std(e)), g.res / 2)
    sdv = np.where(has, sd, sd / math.sqrt(max(1 - rho ** 2, 1e-4)))
    return bits_gaussian(yr, pred, sdv, g), rho, sd


def approx_bits(built, y, g, library_size, rows=None, prev=None):
    """Search score: exact structure and residual bits at the least-squares fit, ~0.5 log2(n) + 2 bits per parameter."""
    rows = np.arange(len(y)) if rows is None else np.asarray(rows)
    th, al, X, off = fit(built, y, rows, refine=False)
    mu = X @ th + off; yr = y[rows]
    sd = max(float(np.sqrt(np.mean((yr - mu) ** 2))), g.res / 2)
    n = len(rows); pb = (n_params(built) + 1 + (prev is not None)) * (0.5 * math.log2(max(n, 2)) + 2)
    if prev is not None: return structure_bits(built, library_size) + pb + float(serial_bits(yr, mu, _prev_in(prev, rows), g)[0].sum())
    return structure_bits(built, library_size) + pb + float(bits_gaussian(yr, mu, sd, g).sum())


def exact_bits(built, y, g, library_size, rows=None, p_range=range(-3, 17), prev=None):
    """Total code length of the formula, every parameter rounded at its best step. Returns a result dict."""
    rows = np.arange(len(y)) if rows is None else np.asarray(rows)
    th, al, X, off = fit(built, y, rows)
    yr = y[rows]
    fam = [j for j, b in enumerate(built) if b.fn is not None]
    pr = _prev_in(prev, rows)
    # A1: a family parameter is sent relative to its range in the library (the decoder knows the range): a - lo, or
    # ln(a / lo) for a log-scale family, so a date like 2022.0 costs what its position in the record costs
    lg = {j: built[j].term.family[2] for j in fam}; lo0 = {j: float(min(built[j].grid[0], built[j].grid[-1])) for j in fam}
    enc = lambda j, a: math.log(a / lo0[j]) if lg[j] else a - lo0[j]
    dec = lambda j, u: lo0[j] * math.exp(u) if lg[j] else u + lo0[j]
    vals = list(th) + [enc(j, al[j]) for j in fam]
    nth = len(th)
    if pr is not None: vals.append(serial_bits(yr, X @ th + off, pr, g)[1])

    def total(steps):
        pb, rv = zip(*[M.param_bits(v, p) for v, p in zip(vals, steps)])
        a = al[:]
        for m, j in enumerate(fam):
            lo, hi = built[j].grid[0], built[j].grid[-1]; a[j] = min(max(dec(j, rv[nth + m]), min(lo, hi)), max(lo, hi))
        if fam: Xr, offr = _design(built, a, rows)
        else: Xr, offr = X, off
        mu = Xr @ np.array(rv[:nth]) + offr
        if not np.all(np.isfinite(mu)): return math.inf, None
        if pr is not None:
            rho = float(np.clip(rv[-1], -0.99, 0.99)); _, _, sd0 = serial_bits(yr, mu, pr, g, rho=rho)
        else:
            rho = None; sd0 = float(np.sqrt(np.mean((yr - mu) ** 2)))
        p_sd = max(0, -math.floor(math.log2(max(sd0, g.res)))) + 3
        sdb, sdr = M.param_bits(sd0, p_sd); sdr = max(sdr, g.res / 2)
        rb = float((serial_bits(yr, mu, pr, g, rho=rho, sd=sdr)[0] if pr is not None else bits_gaussian(yr, mu, sdr, g)).sum())
        return sum(pb) + sdb + rb, dict(param_bits=float(sum(pb) + sdb), residual_bits=rb, sd=sdr, mu=mu, rounded=list(rv), alphas=a, rho=rho)
    steps = [8] * len(vals); best, info = total(steps)
    for _ in range(3):
        changed = False
        for j in range(len(vals)):
            for p in p_range:
                st = steps[:j] + [p] + steps[j + 1:]; t, inf = total(st)
                if t < best - 1e-9: best, info, steps, changed = t, inf, st, True
        if not changed: break
    sb = structure_bits(built, library_size); n = len(rows)
    coef, k = {"level": float(th[0])}, 1
    for b in built:
        if b.fn is not None and b.term.fixed: continue
        kk = b.k; coef[b.name] = [float(v) for v in th[k:k + kk]] if kk > 1 else float(th[k]); k += kk
    return {"terms": [b.name for b in built], "structure_bits": sb, "param_bits": info["param_bits"],
            "residual_bits": info["residual_bits"], "total_bits": sb + best, "bits_per_reading": (sb + best) / n,
            "n": n, "sd": info["sd"], "rho": info.get("rho"), "coef": coef, "alpha": {built[j].name: float(al[j]) for j in fam},
            "_mu": info["mu"], "_steps": steps}


def verify(res, y, g):
    """Re-code the winner's residual part with codelen's arithmetic coder and decode it back."""
    r = sc.score_gaussian(y, res["_mu"], np.full(len(y), res["sd"]), g.res, g.lo, g.hi, verify=True, eps=g.eps)
    return {"ideal_residual_bits": res["residual_bits"], "coded_bits": r["coded_bits"], "roundtrip_ok": bool(r["roundtrip_ok"])}


# --------------------------------------------------------------------------------------------------- search
def _ok(state, built):
    grp = [built[i].term.group for i in state if built[i].term.group]
    return len(grp) == len(set(grp))


def search(built: list[Built], y, g, max_terms=5, beam=6, finalists=10, library_size=None, rows=None, prev=None):
    """Beam search over sparse formulas by approx bits; finalists (plus 'no term') re-scored exactly.

    Then backward elimination on the exact winner (drop a term while that shortens the code).
    Returns (winner, ranked exact results, number of formulas scored).
    """
    L = library_size or len(built)
    memo = {}

    def A(state):
        if state not in memo: memo[state] = approx_bits([built[i] for i in state], y, g, L, rows, prev)
        return memo[state]
    frontier = [()]; A(())
    for _ in range(max_terms):
        cand = set()
        for s in frontier:
            for i in range(len(built)):
                if i in s: continue
                t = tuple(sorted(s + (i,)))
                if _ok(t, built): cand.add(t)
        if not cand: break
        scored = sorted(cand, key=A)
        frontier = scored[:beam]
    ranked = sorted(memo, key=memo.get)[:finalists]
    if () not in ranked: ranked.append(())
    ex = {s: exact_bits([built[i] for i in s], y, g, L, rows, prev=prev) for s in ranked}
    best = min(ex, key=lambda s: ex[s]["total_bits"])
    improved = True
    while improved and best:
        improved = False
        for i in best:
            s = tuple(j for j in best if j != i)
            if s not in ex: ex[s] = exact_bits([built[j] for j in s], y, g, L, rows, prev=prev)
            if ex[s]["total_bits"] < ex[best]["total_bits"] - 1e-9: best, improved = s, True; break
    out = sorted(ex.values(), key=lambda r: r["total_bits"])
    return ex[best], out, len(memo)


def ladder(built, winner_terms, y, g, library_size, prev=None):
    """Add the winner's terms one at a time, each step the one that saves most: the bits waterfall."""
    by = {b.name: b for b in built}; left = list(winner_terms); cur = []; steps = [exact_bits([], y, g, library_size, prev=prev)]
    steps[0]["added"] = "one correction"
    while left:
        trial = [(exact_bits([by[t] for t in cur + [n]], y, g, library_size, prev=prev), n) for n in left]
        r, nme = min(trial, key=lambda z: z[0]["total_bits"]); cur.append(nme); left.remove(nme); r["added"] = nme; steps.append(r)
    return steps


def holdout(built, terms, y, g, train, test, library_size, prev=None):
    """Fit on train rows, then code test rows with that formula; bits per test reading, formula and none."""
    by = {b.name: b for b in built}; out = []
    for names in (terms, []):
        bb = [by[t] for t in names]
        th, al, _, _ = fit(bb, y, train)
        X, off = _design(bb, al, np.arange(len(y))); mu = X @ th + off
        if prev is not None:
            _, rho, sd = serial_bits(y[train], mu[train], _prev_in(prev, train), g)
            out.append(float(serial_bits(y[test], mu[test], _prev_in(prev, test), g, rho=rho, sd=sd)[0].sum()) / len(test))
        else:
            sd = max(float(np.sqrt(np.mean((y[train] - mu[train]) ** 2))), g.res / 2)
            out.append(float(bits_gaussian(y[test], mu[test], sd, g).sum()) / len(test))
    return {"formula_bpr": out[0], "one_correction_bpr": out[1], "transfers": out[0] < out[1]}


def run_case(terms: Sequence[Term], ctx: dict, y, g, folds=(), max_terms=5, beam=6, name="", prev=None):
    """prev (A1): index of each reading's predecessor in its own sequence (-1 for the first), to code the residuals as
    AR(1); None codes them as independent."""
    """The whole case: assemble, search, ladder, hold-out folds, coder check. Arrays prefixed _ are not serialised."""
    y = np.asarray(y, np.float64); n = len(y)
    built = assemble(terms, ctx, n)
    win, ranked, nscored = search(built, y, g, max_terms=max_terms, beam=beam, prev=prev)
    lad = ladder(built, win["terms"], y, g, len(built), prev=prev)
    ho = [dict(fold=f, **holdout(built, win["terms"], y, g, tr, te, len(built), prev=prev)) for f, tr, te in folds] if win["terms"] else []
    meta = {b.name: {"formula": b.term.formula, "meaning": b.term.meaning, "physics": b.term.physics, "unit": b.term.unit} for b in built}
    return {"name": name, "n": n, "library": [b.name for b in built], "library_size": len(built), "formulas_scored": nscored,
            "winner": strip(win), "ranked": [strip(r) for r in ranked[:8]], "ladder": [strip(r) for r in lad],
            "holdout": ho, "transfers": all(h["transfers"] for h in ho) if ho else None,
            "coder_check": verify(win, y, g) if prev is None else {"note": "AR(1) code: ideal bits only"}, "terms_meta": meta,
            "noise_model": "AR(1)" if prev is not None else "independent",
            "bits_saved_per_reading": lad[0]["bits_per_reading"] - win["bits_per_reading"],
            "sd_before": lad[0]["sd"], "sd_after": win["sd"],
            "_y": y, "_mu": win["_mu"], "_mu0": lad[0]["_mu"], "_built": built}


def strip(d): return {k: v for k, v in d.items() if not k.startswith("_")}


def public(res): return {k: v for k, v in res.items() if not k.startswith("_")}
