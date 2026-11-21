"""Predictive distributions -> integer frequency tables over a fixed grid of bins.

Conventions follow the pre-registered SPARC field-model scores (sparc-transformers docs/design/field-model.md,
"Scores for real tests", 21 Nov 2026), so these bits agree with its analytic reference src/sparc_tf/bits.py:
  bins      centred on multiples of `res`: a reading y is coded as the bin [r - res/2, r + res/2], r = round(y / res) * res,
            for r from lo to hi
  escape    each bin's probability is (1 - eps) * model mass + eps * res / (hi - lo)  (a uniform over [lo, hi]);
            the model mass is the exact bin mass and is not renormalised. A reading costs at most -log2(eps * res / (hi - lo)).
  Gaussian  bin mass as a difference of CDFs, or of survival functions in the upper tail (keeps precision)
  quantiles linear CDF between the predicted quantiles, Gaussian tails matched to the outer two quantiles each side
The coder needs integer counts: count = 1 + floor(p * room) per bin, so the coded probability is p up to ~2^-38;
model mass outside [lo, hi] is kept in a never-coded slot, not spread over the bins.
The model part lives in a window of bins around its mass; outside it every bin has the same escape count, so a
table costs O(window), not O(bins).
"""
import math
import numpy as np

try:
    from scipy.special import ndtr as _ndtr, ndtri as _ndtri
except ImportError:                                            # no scipy: the same functions from the standard library
    from statistics import NormalDist
    _erf = np.frompyfunc(math.erf, 1, 1)
    def _ndtr(z): return 0.5 * (1.0 + _erf(np.asarray(z, np.float64) / math.sqrt(2.0)).astype(np.float64))
    def _ndtri(q): return np.array([NormalDist().inv_cdf(float(v)) for v in np.atleast_1d(q)])

SCALE = 1 << 38          # nominal total of a table (escape bins get ~2e5 counts: rounding < 1e-5 bits)
EPS = 1e-3               # weight of the uniform escape (pre-registered value)
WIDTH = 12.0             # window half-width in sds (or kernel widths)


class Grid:
    def __init__(self, res, lo, hi, eps=EPS):
        assert res > 0 and hi > lo and 0 < eps < 1
        self.res, self.lo, self.hi, self.eps = float(res), float(lo), float(hi), float(eps)
        self.k0 = int(round(lo / res)); self.K = int(round(hi / res)) - self.k0 + 1
        self.u = self.res / (self.hi - self.lo)                # escape mass per bin
        assert self.K <= 1 << 24, "too many bins: widen res or narrow [lo, hi]"

    def index(self, y):
        k = np.rint(np.asarray(y, np.float64) / self.res).astype(np.int64) - self.k0
        bad = (k < 0) | (k >= self.K)
        if np.any(bad): raise ValueError(f"{int(bad.sum())} readings outside [{self.lo}, {self.hi}]")
        return k

    def value(self, k): return (np.asarray(k) + self.k0) * self.res

    def edges(self, a, b):
        """Lower edges of bins a..b-1 and the upper edge of b-1 (b - a + 1 values)."""
        return (np.arange(a, b + 1) + self.k0) * self.res - self.res / 2

    def bin_of(self, x): return int(min(max(math.floor(x / self.res + 0.5) - self.k0, 0), self.K - 1))

    def max_bits(self): return -math.log2(self.eps * self.u)


class Table:
    """Integer counts over K bins: escape count c outside [a, b), window counts W inside."""
    __slots__ = ("K", "a", "b", "c", "W", "total")

    def __init__(self, grid, a, b, mass):
        g = np.clip(np.nan_to_num(np.asarray(mass, np.float64), nan=0.0), 0.0, 1.0)
        room = SCALE - grid.K; eps = grid.eps
        self.K, self.a, self.b = grid.K, a, b
        self.c = 1 + int(eps * grid.u * room)
        w = 1 + np.floor(((1.0 - eps) * g + eps * grid.u) * room).astype(np.int64)
        self.W = np.concatenate([[0], np.cumsum(w)])
        # model mass that falls outside the window (past [lo, hi], or in far tails) is not handed to the bins: it
        # stays in a slot after the last bin that is never coded, so the score is the plain (unrenormalised) mixture
        lost = int(np.floor(max(0.0, 1.0 - float(g.sum())) * (1.0 - eps) * room))
        self.total = int(self.c * (grid.K - (b - a)) + self.W[-1]) + lost

    def cum(self, k):
        a, b, c = self.a, self.b, self.c
        if k <= a: return c * k
        if k <= b: return c * a + int(self.W[k - a])
        return c * a + int(self.W[-1]) + c * (k - b)

    def interval(self, k): return self.cum(k), self.cum(k + 1)

    def find(self, t):
        a, c = self.a, self.c
        base = c * a; ws = int(self.W[-1])
        if t < base: return t // c
        if t < base + ws: return a + int(np.searchsorted(self.W, t - base, side="right")) - 1
        return self.b + (t - base - ws) // c

    def bits(self, k):
        lo, hi = self.interval(k); return math.log2(self.total / (hi - lo))


def _gauss_mass(E, m, s):
    """Exact bin masses between consecutive edges E for N(m, s^2); survival functions in the upper tail."""
    z = (E - m) / s; lo, hi = z[:-1], z[1:]
    return np.where(lo > 0, _ndtr(-lo) - _ndtr(-hi), _ndtr(hi) - _ndtr(lo))


def mixture(grid, means, sds, weights=None):
    """Table for a Gaussian mixture (one Gaussian: pass length-1 arrays)."""
    means = np.atleast_1d(np.asarray(means, np.float64)); sds = np.maximum(np.atleast_1d(np.asarray(sds, np.float64)), 1e-9)
    wts = np.full(len(means), 1.0 / len(means)) if weights is None else np.asarray(weights, np.float64) / np.sum(weights)
    a = grid.bin_of(float(np.min(means - WIDTH * sds))); b = grid.bin_of(float(np.max(means + WIDTH * sds))) + 1
    E = grid.edges(a, b)
    mass = sum(w * _gauss_mass(E, m, s) for m, s, w in zip(means, sds, wts))
    return Table(grid, a, b, mass)


def gaussian(grid, mu, sd):
    return mixture(grid, [mu], [sd])


def samples(grid, xs, bandwidth=None):
    """Table from predictive samples: a Gaussian kernel on each sample (Silverman's bandwidth, at least one bin)."""
    xs = np.asarray(xs, np.float64); xs = xs[np.isfinite(xs)]
    if len(xs) == 0: raise ValueError("no finite samples")
    if bandwidth is None:
        sd = xs.std() if len(xs) > 1 else 0.0
        iqr = np.subtract(*np.percentile(xs, [75, 25])) / 1.349 if len(xs) > 1 else 0.0
        spread = min(sd, iqr) if min(sd, iqr) > 0 else max(sd, iqr)
        bandwidth = 0.9 * spread * len(xs) ** -0.2
    return mixture(grid, xs, np.full(len(xs), max(bandwidth, grid.res)))


class QuantileShape:
    """Fixed quantile levels qs (k >= 2, strictly inside (0, 1)); z-scores computed once."""
    def __init__(self, qs):
        self.qs = np.asarray(qs, np.float64); assert np.all(np.diff(self.qs) > 0) and 0 < self.qs[0] and self.qs[-1] < 1
        self.z = _ndtri(self.qs)


def quantiles(grid, shape, values):
    """Table from predicted quantiles `values` at shape.qs: linear CDF between them, Gaussian tails matched to the
    outer two quantiles on each side (sparc_tf.bits.quantile_cdf's convention)."""
    qs, z = shape.qs, shape.z; V = np.sort(np.asarray(values, np.float64))
    sl = max((V[1] - V[0]) / (z[1] - z[0]), 1e-6); sh = max((V[-1] - V[-2]) / (z[-1] - z[-2]), 1e-6)
    a = grid.bin_of(V[0] - WIDTH * sl); b = grid.bin_of(V[-1] + WIDTH * sh) + 1
    E = grid.edges(a, b)
    F = np.interp(E, V, qs)
    low = E <= V[0]; F[low] = _ndtr(z[0] + (E[low] - V[0]) / sl)
    up = E >= V[-1]; F[up] = _ndtr(z[-1] + (E[up] - V[-1]) / sh)
    S = 1.0 - F; S[up] = _ndtr(-(z[-1] + (E[up] - V[-1]) / sh))       # survival, exact in the upper tail
    both_up = up[:-1] & up[1:]
    mass = np.where(both_up, S[:-1] - S[1:], F[1:] - F[:-1])
    return Table(grid, a, b, mass)


HN_MEDIAN = 0.6744897501960817   # HalfNormal(1).icdf(0.5) = ndtri(0.75)


def bars(grid, borders, probs):
    """Table for a TabPFN-style bar distribution (FullSupportBarDistribution, the density it is trained with):
    inner buckets uniform (P_k / width_k); each outer bucket's mass P runs outward from its inner border as a
    half-normal with scale width / HN_MEDIAN (half its mass within one bucket width), past the outer border.
    borders (B+1,), probs (B,) summing to ~1."""
    E0 = np.asarray(borders, np.float64); P = np.clip(np.asarray(probs, np.float64), 0, None)
    assert len(E0) == len(P) + 1 and len(P) >= 3 and np.all(np.diff(E0) > 0)
    bl, br = E0[1], E0[-2]; sl = (E0[1] - E0[0]) / HN_MEDIAN; sr = (E0[-1] - E0[-2]) / HN_MEDIAN
    a = grid.bin_of(bl - WIDTH * sl); b = grid.bin_of(br + WIDTH * sr) + 1
    E = grid.edges(a, b)
    cumP = np.concatenate([[0.0], np.cumsum(P)])                        # cumP[k] = mass below borders[k]
    F = np.interp(E, E0[1:-1], cumP[1:-1])                              # inner buckets: linear
    low = E < bl; F[low] = 2 * P[0] * _ndtr(-(bl - E[low]) / sl)
    up = E > br; S = cumP[-1] - F; S[up] = 2 * P[-1] * _ndtr(-(E[up] - br) / sr)   # survival, exact in the tail
    F[up] = cumP[-1] - S[up]
    both_up = up[:-1] & up[1:]
    mass = np.where(both_up, S[:-1] - S[1:], F[1:] - F[:-1])
    return Table(grid, a, b, mass)
