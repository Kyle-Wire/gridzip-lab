"""Search short formulas by total code length (MDL): which formula, added to a model, stores its errors in the fewest bits.

Written fresh for SPARC / gridzip-lab; no lpaq1, cmix or other GPL code.

A candidate is  y ~ N(c + sum_j theta_j * f_j(x; alpha_j), sd)  with f_j taken from a named feature library.
A feature is either a fixed column or a one-parameter family f(alpha) (e.g. exp(-xi / L)), whose alpha is fitted and
paid for like any other parameter. Total bits follow codelen.mdl exactly:
  structure   Elias-delta length + log2 |vocab| per token (one token per term, one for the intercept, one '+' per join)
  parameters  each value rounded to a step 2**-p, p chosen per parameter to minimise the total (zigzag Elias-gamma)
  residuals   the readings binned on the grid, under the Gaussian, with codelen's escape mixture
`bits_gaussian` is the vectorised ideal code length; it equals codelen.score_gaussian to ~1e-6 bits per reading (the
tests check), and the winners are re-scored with the real arithmetic coder (`verify`).
"""
import itertools
import math
import numpy as np
from . import mdl as M
from . import scoring as sc
from .dists import _ndtr, EPS


class Grid:
    def __init__(self, res, lo, hi, eps=EPS): self.res, self.lo, self.hi, self.eps = float(res), float(lo), float(hi), float(eps)


def bits_gaussian(y, mu, sd, g):
    """Per-reading ideal bits of y on grid g under N(mu, sd): codelen's bins, escape and unrenormalised mass."""
    y = np.asarray(y, np.float64); mu = np.broadcast_to(np.asarray(mu, np.float64), y.shape)
    sd = np.maximum(np.broadcast_to(np.asarray(sd, np.float64), y.shape), 1e-9)
    r = np.rint(y / g.res) * g.res
    a, b = (r - g.res / 2 - mu) / sd, (r + g.res / 2 - mu) / sd
    mass = np.where(a > 0, _ndtr(-a) - _ndtr(-b), _ndtr(b) - _ndtr(a))
    return -np.log2((1 - g.eps) * mass + g.eps * g.res / (g.hi - g.lo))


class Family:
    """A one-parameter feature f(alpha) with alpha searched on [lo, hi] (log scale when log=True)."""
    def __init__(self, fn, lo, hi, log=False, n=41, fixed=False):
        """fixed=True: the term enters with coefficient 1 (a physical law with one unknown), not a fitted slope."""
        self.fn, self.lo, self.hi, self.log, self.n, self.fixed = fn, lo, hi, log, n, fixed
    def grid(self):
        return np.geomspace(self.lo, self.hi, self.n) if self.log else np.linspace(self.lo, self.hi, self.n)


def _fixed(lib, t): return isinstance(lib[t], Family) and lib[t].fixed


def _design(lib, terms, alphas, rows=None):
    """Design matrix (intercept + fitted terms) and the offset of fixed-coefficient families."""
    cols, off = [], None
    for t, a in zip(terms, alphas):
        f = lib[t]; v = np.asarray(f.fn(a) if isinstance(f, Family) else f, np.float64)
        if _fixed(lib, t): off = v if off is None else off + v
        else: cols.append(v)
    n = _n(lib)
    X = np.column_stack([np.ones(n)] + cols); off = np.zeros(n) if off is None else off
    return (X, off) if rows is None else (X[rows], off[rows])


def _n(lib):
    for v in lib.values():
        return len(v.fn(v.lo) if isinstance(v, Family) else v)


def _lstsq(X, y):
    th, *_ = np.linalg.lstsq(X, y, rcond=None); return th


def _fit_alphas(lib, terms, y, rows):
    """Least-squares fit with each family's alpha chosen on its grid, then refined by golden section (coordinate)."""
    fams = [j for j, t in enumerate(terms) if isinstance(lib[t], Family)]
    alphas = [None] * len(terms)
    for j in fams: alphas[j] = float(np.median(lib[terms[j]].grid()))
    def sse(al):
        X, off = _design(lib, terms, al, rows); yy = y[rows] - off; th = _lstsq(X, yy); e = yy - X @ th; return float(e @ e)
    for _ in range(2):
        for j in fams:
            fam = lib[terms[j]]; best = None
            for a in fam.grid():
                al = alphas[:]; al[j] = float(a); s = sse(al)
                if best is None or s < best[0]: best = (s, float(a))
            alphas[j] = best[1]
            # golden-section refinement in the family's own scale around the best grid point
            gr = fam.grid(); k = int(np.argmin(np.abs(gr - best[1])))
            lo, hi = gr[max(k - 1, 0)], gr[min(k + 1, len(gr) - 1)]
            tr = (math.log, math.exp) if fam.log else (lambda v: v, lambda v: v)
            u, v = tr[0](lo), tr[0](hi)
            for _ in range(30):
                m1, m2 = u + 0.382 * (v - u), u + 0.618 * (v - u)
                a1, a2 = alphas[:], alphas[:]; a1[j], a2[j] = tr[1](m1), tr[1](m2)
                if sse(a1) < sse(a2): v = m2
                else: u = m1
            alphas[j] = tr[1]((u + v) / 2)
    return alphas


def _vocab(lib): return ["c", "+"] + list(lib)


def tokens(terms): return ["c"] + [x for t in terms for x in ("+", t)]


def score(lib, terms, y, g, rows=None, p_range=range(-3, 17)):
    """Total MDL bits of formula `terms` on y[rows]. Returns a dict with the fitted, rounded parameters."""
    rows = np.arange(len(y)) if rows is None else np.asarray(rows)
    alphas = _fit_alphas(lib, terms, y, rows)
    X, off = _design(lib, terms, alphas, rows); th = _lstsq(X, y[rows] - off); yr = y[rows]
    fam_idx = [j for j, t in enumerate(terms) if isinstance(lib[t], Family)]
    vals = list(th) + [alphas[j] for j in fam_idx]          # linear coefficients, then family parameters
    steps = [8] * len(vals)

    def total(st):
        pb, rv = zip(*[M.param_bits(v, p) for v, p in zip(vals, st)])
        al = alphas[:]
        for n, j in enumerate(fam_idx): al[j] = rv[len(th) + n]
        Xr, offr = _design(lib, terms, al, rows); mu = Xr @ np.array(rv[:len(th)]) + offr
        sd0 = float(np.sqrt(np.mean((yr - mu) ** 2)))
        p_sd = max(0, -math.floor(math.log2(max(sd0, g.res)))) + 3
        sdb, sdr = M.param_bits(sd0, p_sd); sdr = max(sdr, g.res / 2)
        rb = float(bits_gaussian(yr, mu, sdr, g).sum())
        return sum(pb) + sdb + rb, dict(param_bits=float(sum(pb)), sd_bits=sdb, residual_bits=rb, sd=sdr,
                                        rounded=list(rv), mu=mu, alphas=al)
    best, info = total(steps)
    for _ in range(3):
        changed = False
        for j in range(len(vals)):
            for p in p_range:
                st = steps[:j] + [p] + steps[j + 1:]; t, inf = total(st)
                if t < best - 1e-9: best, info, steps, changed = t, inf, st, True
        if not changed: break
    sb = M.structure_bits(tokens(terms), _vocab(lib))
    n = len(rows)
    return {"terms": list(terms), "structure_bits": sb, **{k: v for k, v in info.items() if k != "mu"},
            "total_bits": sb + best, "bits_per_reading": (sb + best) / n, "n": n,
            "coef": dict(zip(["c"] + [t for t in terms if not _fixed(lib, t)], info["rounded"][:len(th)])),
            "alpha": {terms[j]: info["rounded"][len(th) + n_] for n_, j in enumerate(fam_idx)},
            "fit_coef": dict(zip(["c"] + [t for t in terms if not _fixed(lib, t)], map(float, th))),
            "fit_alpha": {terms[j]: float(alphas[j]) for j in fam_idx},
            "steps": steps, "_mu": info["mu"]}


def verify(res, y, g, rows=None):
    """Re-score a result's residual part with codelen's real arithmetic coder; returns coded residual bits."""
    rows = np.arange(len(y)) if rows is None else np.asarray(rows)
    r = sc.score_gaussian(y[rows], res["_mu"], np.full(len(rows), res["sd"]), g.res, g.lo, g.hi, verify=True, eps=g.eps)
    return {"ideal_residual_bits": res["residual_bits"], "coder_ideal_bits": r["total_bits"],
            "coded_bits": r["coded_bits"], "roundtrip_ok": r["roundtrip_ok"]}


def allowed(terms, exclusive):
    """No two terms from the same exclusive group (e.g. A*xi and A*exp_xi say the same thing twice)."""
    for grp in exclusive:
        if sum(t in grp for t in terms) > 1: return False
    return True


def _quick_bits(lib, terms, y, rows):
    """Fast stand-in for the total bits of a formula without families: the Gaussian code length of the least-squares
    residuals (n/2 log2 of the residual variance) plus 0.5 log2 n per fitted value. Used only to choose which
    formulas the exact scorer sees; it never decides a ranking on its own."""
    X, off = _design(lib, terms, [None] * len(terms), rows); yr = y[rows] - off
    th = _lstsq(X, yr); r = yr - X @ th; n = len(rows)
    return 0.5 * n * math.log2(max(float(r @ r) / n, 1e-12)) + 0.5 * (len(th) + 1) * math.log2(n)


def search(lib, y, g, max_terms=3, rows=None, exclusive=(), keep=20, candidates=None, prescreen=None):
    """Every formula of up to max_terms library terms (plus intercept), ranked by total bits.

    ``prescreen`` (an int): when the library has no families, every formula is first ranked by `_quick_bits`, and only
    the best ``prescreen`` of them (and every formula of fewer than 2 terms) get the exact MDL score. Large libraries
    (tens of features and their products) are then searched in minutes; the count returned is still every formula
    considered. Checked against the exhaustive search in codelen/test_lawsearch.py (test_prescreen_matches_exhaustive)."""
    names = list(lib) if candidates is None else list(candidates)
    rows_ = np.arange(len(y)) if rows is None else np.asarray(rows)
    combos = [t for k in range(max_terms + 1) for t in itertools.combinations(names, k) if allowed(t, exclusive)]
    if prescreen and not any(isinstance(lib[t], Family) for t in names):
        q = sorted(combos, key=lambda t: _quick_bits(lib, t, y, rows_))
        chosen = set(q[:prescreen]) | {t for t in combos if len(t) < 2}
        todo = [t for t in combos if t in chosen]
    else:
        todo = combos
    out = [score(lib, terms, y, g, rows) for terms in todo]
    out.sort(key=lambda r: r["total_bits"])
    return out[:keep], len(combos)


def holdout(lib, terms, y, g, train, test):
    """Fit on train rows (alphas, coefficients and sd), then the bits to code the test rows with that formula."""
    al = _fit_alphas(lib, terms, y, train)
    X, off = _design(lib, terms, al); th = _lstsq(X[train], y[train] - off[train])
    mu = X @ th + off
    sd = max(float(np.sqrt(np.mean((y[train] - mu[train]) ** 2))), g.res / 2)
    return float(bits_gaussian(y[test], mu[test], sd, g).sum())


def context_bits(y, order, ctx, g, s0=0.3, w0=2.0):
    """Prequential bits of y coded in `order`, each predicted by the mean of its context (earlier) readings,
    sd fitted online; the same rule as codelen.ContextMean (checked in the tests)."""
    seen = {}; ss, n_err, tot, cnt = 0.0, 0, 0.0, 0
    bits = np.zeros(len(order))
    for i, j in enumerate(order):
        c = [seen[q] for q in ctx[i] if q >= 0 and q in seen]
        mu = float(np.mean(c)) if c else (tot / cnt if cnt else (g.lo + g.hi) / 2)
        sd = max(math.sqrt((w0 * s0 ** 2 + ss) / (w0 + n_err)), g.res / 2)
        yq = round(y[j] / g.res) * g.res
        bits[i] = bits_gaussian(np.array([y[j]]), mu, sd, g)[0]
        ss += (yq - mu) ** 2; n_err += 1; seen[j] = yq; tot += yq; cnt += 1
    return bits
