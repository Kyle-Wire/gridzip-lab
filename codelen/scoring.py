"""Code length of readings under a forecaster: bits to store each quantised reading, exactly, with a real coder.

A predictor gives, for reading i, a distribution over the grid, using only what a decoder would already have:
readings 0..i-1 and side information that is not being scored. `score` sums -log2 p (the ideal length of the integer
tables the coder uses) and, with verify=True, also runs the arithmetic encoder and decoder and asserts the readings
come back. Reading order is the coding order: put context sensors' readings earlier than the readings that use them.

Predictors (objects with reset(), table(i, y_so_far) and update(i, k, y)):
  Given(grid, mu, sd)          fixed Gaussian forecasts, one per reading (a model's output)
  GivenSamples(grid, S)        fixed sample forecasts, S[i] = samples for reading i (generative model)
  GivenQuantiles(grid, qs, V)  fixed quantile forecasts, V[i] = values at levels qs (linear CDF, Gaussian tails)
  GivenBars(grid, borders, P)  fixed bar-distribution forecasts (TabPFN), half-normal outer tails
  ContextMean(grid, ctx, ...)  baseline: mean of earlier readings named in ctx[i] (indices < i, -1 = none), sd fitted
                               online from its own past errors; falls back to the running mean of all earlier readings
  Persistence(grid, group)     baseline: the last reading of the same sensor (group id), sd fitted online per sensor
"""
import math
import numpy as np
from . import coder, dists


class Given:
    def __init__(self, grid, mu, sd):
        self.grid, self.mu, self.sd = grid, np.asarray(mu, np.float64), np.asarray(sd, np.float64)
    def reset(self): pass
    def table(self, i, y): return dists.gaussian(self.grid, self.mu[i], self.sd[i])
    def update(self, i, k, y): pass


class GivenSamples(Given):
    def __init__(self, grid, S, bandwidth=None):
        self.grid, self.S, self.bw = grid, S, bandwidth
    def table(self, i, y): return dists.samples(self.grid, self.S[i], self.bw)


class GivenQuantiles(Given):
    """Fixed quantile forecasts: V[i] = values at levels qs for reading i (e.g. TabPFN at 0.01..0.99)."""
    def __init__(self, grid, qs, V):
        self.grid, self.shape, self.V = grid, dists.QuantileShape(qs), np.asarray(V, np.float64)
    def table(self, i, y): return dists.quantiles(self.grid, self.shape, self.V[i])


class GivenBars(Given):
    """Fixed bar-distribution forecasts (TabPFN predict(output_type="full")): borders (B+1,) shared or (n, B+1),
    probs (n, B) = softmax(logits). Outer buckets get TabPFN's half-normal tails."""
    def __init__(self, grid, borders, probs):
        self.grid, self.B, self.P = grid, np.asarray(borders, np.float64), np.asarray(probs, np.float64)
    def table(self, i, y): return dists.bars(self.grid, self.B if self.B.ndim == 1 else self.B[i], self.P[i])


class _OnlineSD:
    """sd from past squared errors with a prior of w0 pseudo-readings at s0; floored at half a bin."""
    def __init__(self, s0, w0, floor): self.s0, self.w0, self.floor = s0, w0, floor; self.reset()
    def reset(self): self.ss, self.n = 0.0, 0
    def sd(self): return max(math.sqrt((self.w0 * self.s0 ** 2 + self.ss) / (self.w0 + self.n)), self.floor)
    def add(self, e): self.ss += e * e; self.n += 1


class ContextMean:
    def __init__(self, grid, ctx, s0=1.0, w0=2.0, m0=None):
        self.grid, self.ctx = grid, np.asarray(ctx, np.int64).reshape(len(ctx), -1)
        assert np.all(self.ctx < np.arange(len(self.ctx))[:, None]), "context must name earlier readings only"
        self.s0, self.w0, self.m0 = s0, w0, m0; self.reset()
    def reset(self):
        self.err = _OnlineSD(self.s0, self.w0, self.grid.res / 2); self.sum, self.n = 0.0, 0
    def _mu(self, i, y):
        c = [y[j] for j in self.ctx[i] if j >= 0]
        if c: return float(np.mean(c))
        if self.n: return self.sum / self.n
        return self.m0 if self.m0 is not None else (self.grid.lo + self.grid.hi) / 2
    def table(self, i, y):
        self._last = self._mu(i, y); return dists.gaussian(self.grid, self._last, self.err.sd())
    def update(self, i, k, y):
        self.err.add(y[i] - self._last); self.sum += y[i]; self.n += 1


class Persistence:
    def __init__(self, grid, group, s0=1.0, w0=2.0, m0=None):
        self.grid, self.group = grid, np.asarray(group); self.s0, self.w0, self.m0 = s0, w0, m0; self.reset()
    def reset(self): self.last, self.err, self.sum, self.n = {}, {}, 0.0, 0
    def table(self, i, y):
        g = self.group[i]
        if g in self.last: self._mu = self.last[g]
        elif self.n: self._mu = self.sum / self.n
        else: self._mu = self.m0 if self.m0 is not None else (self.grid.lo + self.grid.hi) / 2
        e = self.err.setdefault(g, _OnlineSD(self.s0, self.w0, self.grid.res / 2))
        return dists.gaussian(self.grid, self._mu, e.sd())
    def update(self, i, k, y):
        g = self.group[i]; self.err[g].add(y[i] - self._mu); self.last[g] = y[i]; self.sum += y[i]; self.n += 1


def score(y, predictor, grid, verify=False, per_reading=False):
    """Bits to store readings y (quantised to the grid) under predictor. Returns a dict."""
    k = grid.index(y); yq = grid.value(k); n = len(k)
    predictor.reset(); enc = coder.Encoder() if verify else None
    bits = np.zeros(n); seen = np.full(n, np.nan)   # the predictor sees readings only once they are coded (no peeking)
    for i in range(n):
        t = predictor.table(i, seen)
        lo, hi = t.interval(int(k[i])); bits[i] = math.log2(t.total / (hi - lo))
        if enc: enc.encode(lo, hi, t.total)
        seen[i] = yq[i]; predictor.update(i, int(k[i]), seen)
    out = {"n": n, "resolution": grid.res, "lo": grid.lo, "hi": grid.hi, "eps": grid.eps, "bins": grid.K,
           "total_bits": float(bits.sum()), "bits_per_reading": float(bits.mean()) if n else 0.0}
    if verify:
        stream = coder.pack(enc.finish())
        back = decode(stream, n, predictor, grid)
        out.update({"coded_bytes": len(stream), "coded_bits": 8 * len(stream), "roundtrip_ok": bool(np.array_equal(back, k))})
        assert out["roundtrip_ok"], "arithmetic decode differs from the readings"
    if per_reading: out["bits"] = bits
    return out


def decode(stream, n, predictor, grid):
    """Bin indices back from a stream; the predictor sees only what has been decoded so far."""
    predictor.reset(); dec = coder.Decoder(coder.unpack(stream))
    k = np.zeros(n, np.int64); yq = np.full(n, np.nan)
    for i in range(n):
        t = predictor.table(i, yq)
        s = t.find(dec.target(t.total)); lo, hi = t.interval(s); dec.consume(lo, hi, t.total)
        k[i] = s; yq[i] = grid.value(s); predictor.update(i, s, yq)
    return k


def score_gaussian(y, mu, sd, res, lo, hi, verify=False, per_reading=False, eps=dists.EPS):
    g = dists.Grid(res, lo, hi, eps); return score(y, Given(g, mu, sd), g, verify, per_reading)


def score_samples(y, S, res, lo, hi, verify=False, per_reading=False, bandwidth=None, eps=dists.EPS):
    g = dists.Grid(res, lo, hi, eps); return score(y, GivenSamples(g, S, bandwidth), g, verify, per_reading)


def score_quantiles(y, qs, V, res, lo, hi, verify=False, per_reading=False, eps=dists.EPS):
    g = dists.Grid(res, lo, hi, eps); return score(y, GivenQuantiles(g, qs, V), g, verify, per_reading)
