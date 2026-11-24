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
    cost_bits: Optional[float] = None  # v2: a grammar term's own description length (None = a library term)


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


def _sigkey(sig):
    v = np.round(sig, 6); nz = np.flatnonzero(np.abs(v) > 1e-6)
    if len(nz) and v[nz[0]] < 0: v = -v
    return (v + 0.0).tobytes()


def assemble(terms: Sequence[Term], ctx: dict, n: int) -> list[Built]:
    """Build every term whose inputs are in ctx; drop terms that are constant or duplicate an earlier one."""
    out, seen, keys = [], [], set()
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
                key = _sigkey(sig)
                if t.cost_bits is not None:        # v2 grammar term: hashed duplicate check (thousands of them)
                    if key in keys: continue
                elif any(np.allclose(sig, s) or np.allclose(sig, -s) for s in seen): continue
                seen.append(sig) if t.cost_bits is None else None; keys.add(key); out.append(Built(t, cols=cols))
        except (ValueError, FloatingPointError, ZeroDivisionError, IndexError, TypeError):
            continue
    return _drop_library_shapes(out)


def _unit(M_):
    M_ = M_ - M_.mean(0); return M_ / (np.linalg.norm(M_, axis=0) + 1e-300)


def _drop_library_shapes(out, r_max=0.995):
    """v2: a grammar term whose column the library already has (|r| > r_max against a plain library term, or against a
    family member anywhere on its grid) is dropped: the library names that shape, and the family's free parameter
    already covers the grammar's fixed constant. Grammar columns are scaled to unit sd (the decoder can redo it)."""
    gi = [i for i, b in enumerate(out) if b.term.cost_bits is not None]
    if not gi: return out
    lib = []
    for b in out:
        if b.term.cost_bits is not None: continue
        if b.fn is not None:
            for a in b.grid:
                v = b.fn(a)
                if np.all(np.isfinite(v)): lib.append(v)
        else: lib.append(b.cols)
    for i in gi:
        c = out[i].cols; out[i].cols = c / (c.std(0) + 1e-300)
    if not lib: return out
    Lm = _unit(np.hstack(lib)); Gm = _unit(np.hstack([out[i].cols for i in gi]))
    r = np.abs(Gm.T @ Lm).max(1)
    drop = {gi[k] for k in np.flatnonzero(r > r_max)}
    return [b for i, b in enumerate(out) if i not in drop]


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


def _lib(library_size):
    """library_size is an int (library only) or (library terms, grammar on) -- v2."""
    return (library_size, False) if isinstance(library_size, (int, np.integer)) else (int(library_size[0]), bool(library_size[1]))


def structure_bits(built, library_size):
    """c, then '+ term' per term. A library term costs log2(vocab); with the grammar on, one more bit says library or
    grammar, and a grammar term costs that bit plus its own expression length. With the grammar off this is the
    PROTOCOL-4 formula exactly."""
    L, flag = _lib(library_size); lv = math.log2(vocab_size(L))
    names = sum(lv + (1.0 if flag else 0.0) if b.term.cost_bits is None else 1.0 + b.term.cost_bits for b in built)
    return M.delta_bits(n_tokens(built) + 1) + (1 + len(built)) * lv + names


def library_info(built):
    """(number of library terms, grammar on) for a built pool."""
    return (sum(1 for b in built if b.term.cost_bits is None), any(b.term.cost_bits is not None for b in built))


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


MIX_W, MIX_K = (0.05, 0.10, 0.20), 4.0


def bits_mix(y, mu, sd, w, g, k=MIX_K):
    """v2 heavy-tail code: (1 - w) N(mu, sd) + w N(mu, k sd), on codelen's bins and escape."""
    from codelen.lawsearch import _ndtr
    y = np.asarray(y, np.float64); mu = np.broadcast_to(np.asarray(mu, np.float64), y.shape)
    r = np.rint(y / g.res) * g.res; out = 0.0
    for wt, s in ((1 - w, sd), (w, k * sd)):
        a, b = (r - g.res / 2 - mu) / s, (r + g.res / 2 - mu) / s
        out = out + wt * np.where(a > 0, _ndtr(-a) - _ndtr(-b), _ndtr(b) - _ndtr(a))
    return -np.log2((1 - g.eps) * out + g.eps * g.res / (g.hi - g.lo))


def robust_sd(e): return 1.4826 * float(np.median(np.abs(e - np.median(e))))


def resid_bits(yr, mu, sd_rms, g, robust, sd_rob=None):
    """Residual bits under the shortest of: Gaussian, or (robust) the two-scale mixtures; 2 bits name the choice.
    Returns (bits, choice, sd used)."""
    best = (float(bits_gaussian(yr, mu, sd_rms, g).sum()), "gauss", sd_rms)
    if not robust: return best
    sdr = max(sd_rob if sd_rob is not None else robust_sd(yr - mu), g.res / 2)
    for w in MIX_W:
        b = float(bits_mix(yr, mu, sdr, w, g).sum())
        if b < best[0]: best = (b, f"mix{w:g}", sdr)
    return best[0] + 2.0, best[1], best[2]


def approx_bits(built, y, g, library_size, rows=None, prev=None, robust=False):
    """Search score: exact structure and residual bits at the least-squares fit, ~0.5 log2(n) + 2 bits per parameter."""
    rows = np.arange(len(y)) if rows is None else np.asarray(rows)
    th, al, X, off = fit(built, y, rows, refine=False)
    mu = X @ th + off; yr = y[rows]
    sd = max(float(np.sqrt(np.mean((yr - mu) ** 2))), g.res / 2)
    n = len(rows); pb = (n_params(built) + 1 + (prev is not None)) * (0.5 * math.log2(max(n, 2)) + 2)
    if prev is not None: v = structure_bits(built, library_size) + pb + float(serial_bits(yr, mu, _prev_in(prev, rows), g)[0].sum())
    else: v = structure_bits(built, library_size) + pb + resid_bits(yr, mu, sd, g, robust)[0]
    return v if math.isfinite(v) else math.inf     # v2: a non-finite score must not scramble the beam's sort


def exact_bits(built, y, g, library_size, rows=None, p_range=range(-3, 17), prev=None, robust=False):
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
        sdb, sdr = M.param_bits(sd0, p_sd); sdr = max(sdr, g.res / 2); code = "gauss"
        if pr is not None: rb = float(serial_bits(yr, mu, pr, g, rho=rho, sd=sdr)[0].sum())
        else:
            rb = float(bits_gaussian(yr, mu, sdr, g).sum())
            if robust:   # v2: the mixture's core sd is sent at the same step as the Gaussian sd
                _, sdm = M.param_bits(max(robust_sd(yr - mu), g.res / 2), p_sd); sdm = max(sdm, g.res / 2)
                rb2, code2, _ = resid_bits(yr, mu, sdr, g, True, sd_rob=sdm)
                rg = rb + 2.0
                if rb2 < rg: rb, code, sdr = rb2, code2, sdm
                else: rb = rg
        return sum(pb) + sdb + rb, dict(param_bits=float(sum(pb) + sdb), residual_bits=rb, sd=sdr, mu=mu, rounded=list(rv), alphas=a, rho=rho, code=code)
    if not all(np.isfinite(v) for v in vals):
        return {"terms": [b.name for b in built], "total_bits": math.inf, "bits_per_reading": math.inf, "structure_bits": math.inf,
                "param_bits": math.inf, "residual_bits": math.inf, "n": len(rows), "sd": math.inf, "coef": {}, "alpha": {}, "code": "gauss",
                "_mu": np.zeros(len(rows)), "_steps": []}
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
            "n": n, "sd": info["sd"], "rho": info.get("rho"), "code": info.get("code", "gauss"), "coef": coef, "alpha": {built[j].name: float(al[j]) for j in fam},
            "_mu": info["mu"], "_steps": steps}


def verify(res, y, g):
    """Re-code the winner's residual part with codelen's arithmetic coder and decode it back."""
    r = sc.score_gaussian(y, res["_mu"], np.full(len(y), res["sd"]), g.res, g.lo, g.hi, verify=True, eps=g.eps)
    return {"ideal_residual_bits": res["residual_bits"], "coded_bits": r["coded_bits"], "roundtrip_ok": bool(r["roundtrip_ok"])}


# --------------------------------------------------------------------------------------------------- search
def _ok(state, built):
    grp = [built[i].term.group for i in state if built[i].term.group]
    return len(grp) == len(set(grp))


def search(built: list[Built], y, g, max_terms=5, beam=6, finalists=10, library_size=None, rows=None, prev=None,
           robust=False, prescreen=30, memo_out=None):
    """Beam search over sparse formulas by approx bits; finalists (plus 'no term') re-scored exactly.

    Then backward elimination on the exact winner (drop a term while that shortens the code).
    v2: grammar terms (cost_bits set) are prescreened per beam state: only the `prescreen` whose columns correlate
    best with that state's remainder are scored (library terms are always scored). memo_out, if a dict, receives
    every scored formula's approx bits (for the posterior odds in guards).
    Returns (winner, ranked exact results, number of formulas scored).
    """
    L = library_size if library_size is not None else library_info(built)
    rows = np.arange(len(y)) if rows is None else np.asarray(rows)
    memo = {}
    gi = [i for i, b in enumerate(built) if b.term.cost_bits is not None and b.fn is None]
    li = [i for i in range(len(built)) if i not in set(gi)]
    if gi:
        Z = np.column_stack([built[i].cols[rows, 0] for i in gi]); Z = Z - Z.mean(0); Z /= np.linalg.norm(Z, axis=0) + 1e-300

    def A(state):
        if state not in memo: memo[state] = approx_bits([built[i] for i in state], y, g, L, rows, prev, robust=robust)
        return memo[state]

    def pool(state):
        if not gi: return range(len(built))
        bb = [built[i] for i in state]; th, al, X, off = fit(bb, y, rows, refine=False); e = y[rows] - X @ th - off
        sc_ = np.abs(Z.T @ (e - e.mean())); top = np.argsort(-sc_)[:prescreen]
        return li + [gi[k] for k in top]
    frontier = [()]; A(())
    for _ in range(max_terms):
        cand = set()
        for s in frontier:
            for i in pool(s):
                if i in s: continue
                t = tuple(sorted(s + (i,)))
                if _ok(t, built): cand.add(t)
        if not cand: break
        scored = sorted(cand, key=A)
        frontier = scored[:beam]
    ranked = sorted(memo, key=memo.get)[:finalists]
    if gi:   # v2: approx bits undercount parameters, so bigger formulas crowd the finalists; keep the best 3 of each size
        for k in range(1, max_terms + 1):
            for s_ in sorted((s_ for s_ in memo if len(s_) == k), key=memo.get)[:3]:
                if s_ not in ranked: ranked.append(s_)
    if () not in ranked: ranked.append(())
    EX = lambda s: exact_bits([built[i] for i in s], y, g, L, rows, prev=prev, robust=robust)
    ex = {s: EX(s) for s in ranked}
    best = min(ex, key=lambda s: ex[s]["total_bits"])
    improved = True
    while improved and best:
        improved = False
        for i in best:
            s = tuple(j for j in best if j != i)
            if s not in ex: ex[s] = EX(s)
            if ex[s]["total_bits"] < ex[best]["total_bits"] - 1e-9: best, improved = s, True; break
    if memo_out is not None: memo_out.update(memo); memo_out["_exact"] = {s: r["total_bits"] for s, r in ex.items()}; memo_out["_best"] = best
    out = sorted(ex.values(), key=lambda r: r["total_bits"])
    return ex[best], out, len(memo)


def ladder(built, winner_terms, y, g, library_size, prev=None, robust=False):
    """Add the winner's terms one at a time, each step the one that saves most: the bits waterfall."""
    by = {b.name: b for b in built}; left = list(winner_terms); cur = []
    steps = [exact_bits([], y, g, library_size, prev=prev, robust=robust)]
    steps[0]["added"] = "one correction"
    while left:
        trial = [(exact_bits([by[t] for t in cur + [n]], y, g, library_size, prev=prev, robust=robust), n) for n in left]
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


def run_case(terms: Sequence[Term], ctx: dict, y, g, folds=(), max_terms=5, beam=6, name="", prev=None, robust=False,
             memo_out=None):
    """prev (A1): index of each reading's predecessor in its own sequence (-1 for the first), to code the residuals as
    AR(1); None codes them as independent."""
    """The whole case: assemble, search, ladder, hold-out folds, coder check. Arrays prefixed _ are not serialised."""
    y = np.asarray(y, np.float64); n = len(y)
    built = assemble(terms, ctx, n); L = library_info(built)
    win, ranked, nscored = search(built, y, g, max_terms=max_terms, beam=beam, prev=prev, robust=robust, memo_out=memo_out)
    lad = ladder(built, win["terms"], y, g, L, prev=prev, robust=robust)
    ho = [dict(fold=f, **holdout(built, win["terms"], y, g, tr, te, L, prev=prev)) for f, tr, te in folds] if win["terms"] else []
    meta = {b.name: {"formula": b.term.formula, "meaning": b.term.meaning, "physics": b.term.physics, "unit": b.term.unit} for b in built}
    return {"name": name, "n": n, "library": [b.name for b in built if b.term.cost_bits is None], "library_size": L[0],
            "grammar_size": sum(1 for b in built if b.term.cost_bits is not None), "formulas_scored": nscored,
            "winner": strip(win), "ranked": [strip(r) for r in ranked[:8]], "ladder": [strip(r) for r in lad],
            "holdout": ho, "transfers": all(h["transfers"] for h in ho) if ho else None,
            "coder_check": (verify(win, y, g) if win.get("code", "gauss") == "gauss" else {"note": "mixture code: ideal bits only"})
            if prev is None else {"note": "AR(1) code: ideal bits only"}, "terms_meta": {b.name: meta[b.name] for b in built if b.term.cost_bits is None or b.name in win["terms"]},
            "noise_model": "AR(1)" if prev is not None else "independent",
            "bits_saved_per_reading": lad[0]["bits_per_reading"] - win["bits_per_reading"],
            "sd_before": lad[0]["sd"], "sd_after": win["sd"],
            "_y": y, "_mu": win["_mu"], "_mu0": lad[0]["_mu"], "_built": built}


def strip(d): return {k: v for k, v in d.items() if not k.startswith("_")}


def public(res): return {k: v for k, v in res.items() if not k.startswith("_")}
