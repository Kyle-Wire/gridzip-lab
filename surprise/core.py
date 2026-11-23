"""Surprise score: how many extra bits it takes to describe one unit (a house) given what its neighbours do.

For each unit i and reading s (a storm, a summer week, an InSAR epoch) there is a value y[i, s] on a common coding
scale. The neighbour model says: unit i behaves like one of its neighbours j, picked with weight w_ij, plus kernel
noise. So reading s of unit i is coded with a mixture of Gaussians at the neighbours' values y[j, s], widths
    h_ij,s^2 = h^2 + noise[i, s]^2 + noise[j, s]^2        (h fitted once over the block; noise = model noise, e.g. run
                                                        variants)
and the readings of one unit are coded in order with the chain rule: after each reading the weights move toward the
neighbours that matched it (w_ij <- w_ij * N(y_is; y_js, h_ij,s)). That is the exact code length of the joint kernel
density, so a house that is odd in the same way in every storm is not charged six times, and a dry house among dry
neighbours costs almost nothing even when half the block is wet (a single Gaussian would get both wrong).

Spatial weights: Gaussian in distance (sigma = radius / 2) over the units within `radius` (at least min_n), self
excluded. Optional similarity weights on a descriptor (e.g. ground height): w_ij also times N(d_i - d_j; 0, tau), so
"surprise given height" compares a house only with neighbours that sit at about the same height.

bits[i]     code length of unit i's readings (codelen integer tables; with verify=True the block is arithmetic-coded
            and decoded and the round trip asserted)
typical[i]  the expected code length of each of those tables (what a reading drawn from the model would cost)
surprise[i] bits - typical, in bits. Around 0 = as its neighbours suggest; each extra bit halves how often a house like
            its neighbours would look this odd. It says "unlike its neighbours", not "worse": sign[i] says whether the
            house reads wetter (+1) or drier (-1) than its weighted neighbours.
Pure numpy plus codelen; CPU only. Written fresh for gridzip-lab (no GPL code).
"""
import math
import numpy as np
import codelen
from codelen import dists
from codelen.scoring import score as _score


def neighbours(xy, radius, min_n=4):
    """For each unit: indices of the others within radius (at least min_n nearest) and Gaussian distance weights."""
    xy = np.asarray(xy, float)
    d = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1]); np.fill_diagonal(d, np.inf)
    out = []
    for i in range(len(xy)):
        idx = np.where(d[i] <= radius)[0]
        if len(idx) < min_n: idx = np.argsort(d[i])[:min_n]
        w = np.exp(-0.5 * (d[i, idx] / (radius / 2)) ** 2); out.append((idx, w / w.sum()))
    return out


class _Tables:
    """A codelen predictor that hands out prebuilt tables (the chain-rule weights depend only on the unit's own earlier
    readings, which the decoder has, so the tables can be built before coding without peeking)."""
    def __init__(self, tables): self.t = tables
    def reset(self): pass
    def table(self, i, y): return self.t[i]
    def update(self, i, k, y): pass


def _table_entropy(t):
    T = t.total; w = np.diff(t.W).astype(float) / T; w = w[w > 0]; esc = t.c / T
    return float(-(w * np.log2(w)).sum() - (t.K - (t.b - t.a)) * esc * math.log2(esc))


def _code(y, nb, grid, h, noise, extra_w=None, want_tables=False, p_self=0.3, p_wide=0.02, share=0.05):
    """Chain-rule mixture tables for every unit and reading; returns (bits, typical, tables).
    Experts for reading s of unit i: each neighbour j (kernel at y[j, s]); from the second reading on, the unit itself
    (its previous reading plus the neighbours' mean change, width h); and a wide expert at the neighbour mean with the
    block's spread for that reading, so a house unlike every neighbour pays a bounded price instead of the escape.
    Prior weights: neighbours share 1 - p_self - p_wide (p_self = 0 on the first reading); all are updated by Bayes,
    with a fixed share of 5 % returned to the priors after each reading."""
    U, S = y.shape; yq = grid.value(grid.index(y)); tables = []; H = np.zeros((U, S)); B = np.zeros((U, S))
    wide = np.maximum(y.std(0), 4 * h)
    for i, (idx, w0) in enumerate(nb):
        wn = w0 * (extra_w[i] if extra_w is not None else 1.0); wn = wn / wn.sum()
        n = len(idx); post = None; mu_prev = None
        for s in range(S):
            sd_n = np.sqrt(h ** 2 + noise[i, s] ** 2 + noise[idx, s] ** 2); mu_n = wn @ y[idx, s]
            means = [y[idx, s]]; sds = [sd_n]
            if s > 0:
                means.append([yq[i, s - 1] + mu_n - mu_prev]); sds.append([math.sqrt(2 * h ** 2 + noise[i, s] ** 2 + noise[i, s - 1] ** 2)])
            means.append([mu_n]); sds.append([wide[s]])
            means = np.concatenate([np.atleast_1d(m) for m in means]); sds = np.concatenate([np.atleast_1d(v) for v in sds])
            if post is None or s == 1:   # priors (the self expert joins at s = 1, taking its share from the others)
                ps = p_self if s > 0 else 0.0
                prior = np.concatenate([wn * (1 - ps - p_wide), [ps] if s > 0 else [], [p_wide]])
                if post is not None:      # carry what reading 0 taught about the neighbours and the wide expert
                    c = post / post.sum(); prior = np.concatenate([c[:n] * (1 - ps), [ps], c[n:] * (1 - ps)])
                post = prior
            wts = post / post.sum(); keep = wts > 1e-9 * wts.max()
            t = dists.mixture(grid, means[keep], sds[keep], wts[keep])
            lo, hi = t.interval(int(grid.index(yq[i, s]))); B[i, s] = math.log2(t.total / (hi - lo)); H[i, s] = _table_entropy(t)
            if want_tables: tables.append(t)
            lik = np.exp(-0.5 * ((yq[i, s] - means) / sds) ** 2) / sds
            post = wts * lik
            post = post / post.sum() if post.sum() > 0 else wts
            # fixed share (Herbster and Warmuth 1998): a little weight goes back to the priors after every reading, so a
            # house can switch from "like its neighbours" to "unlike them" in one storm without paying the escape
            base = np.concatenate([wn * (1 - p_self - p_wide), [p_self], [p_wide]]) if s + 1 < S else None
            if base is not None and len(base) == len(post): post = (1 - share) * post + share * base
            elif base is not None: post = post  # reading 0: the self expert is added at reading 1
            mu_prev = mu_n
    return B, H, tables


def fit_h(y, nb, grid, noise, hs, extra_w=None):
    """The kernel width that gives the block the shortest code (one parameter, log2(len(hs)) bits to send)."""
    best = None
    for h in hs:
        B, _, _ = _code(y, nb, grid, h, noise, extra_w)
        if best is None or B.sum() < best[1]: best = (h, B.sum())
    return best[0]


def score(y, xy, radius, res, lo, hi, noise=None, hs=None, descriptor=None, taus=None, verify=True):
    """Surprise per unit. y: (units, readings) on the coding scale; noise: (units, readings) model-noise sd or None;
    descriptor: (units,) e.g. ground height, for "surprise given the descriptor". Returns a dict of arrays."""
    y = np.asarray(y, float); U, S = y.shape; grid = dists.Grid(res, lo, hi, dists.EPS)
    nb = neighbours(xy, radius); noise = np.zeros_like(y) if noise is None else np.asarray(noise, float)
    hs = hs if hs is not None else res * np.array([1, 1.5, 2, 3, 4, 6, 8, 12, 16])
    h = fit_h(y, nb, grid, noise, hs)
    B, H, T = _code(y, nb, grid, h, noise, want_tables=verify)
    out = {"h": float(h), "bits": B.sum(1), "typical": H.sum(1), "surprise": B.sum(1) - H.sum(1), "per_reading": B - H}
    mu = np.array([[w @ y[idx, s] for s in range(S)] for idx, w in nb]); out["mu"] = mu
    out["sign"] = np.sign(np.round((y - mu).mean(1), 6))
    if verify:
        r = _score(y.ravel(), _Tables(T), grid, verify=True)   # coding order: unit, then reading
        out["roundtrip_ok"] = r["roundtrip_ok"]; out["total_bits"] = r["total_bits"]
        assert abs(r["total_bits"] - B.sum()) < 1e-6 * max(1.0, B.sum())
    if descriptor is not None:
        d = np.asarray(descriptor, float); taus = taus if taus is not None else [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, np.inf]
        best = None
        for tau in taus:
            ew = [np.ones(len(idx)) if not np.isfinite(tau) else np.exp(-0.5 * ((d[i] - d[idx]) / tau) ** 2) + 1e-12 for i, (idx, _) in enumerate(nb)]
            B2, H2, _ = _code(y, nb, grid, h, noise, ew)
            if best is None or B2.sum() < best[1]: best = (tau, B2.sum(), B2, H2)
        tau, _, B2, H2 = best
        out.update({"tau": float(tau), "surprise_given": B2.sum(1) - H2.sum(1), "bits_given": B2.sum(1)})
        out["explained"] = out["surprise"] - out["surprise_given"]
    return out
