"""GridZip v3: coordinates coded with the last K cells as context, blended
PPM-style from the longest context down to order 0.

p_k(s) = (count_k(s) + w_k * p_{k-1}(s)) / (n_k + w_k),  w_k = a*distinct_k + b
p_-1   = uniform over the grid

All counts are "seen so far" (adaptive: the decoder can rebuild them), and
the code length is computed exactly, vectorised. ORDER=0..4.
"""
import math, sys, time
import numpy as np


def occ_before(keys):
    order = np.argsort(keys, kind="stable")
    k = keys[order]
    start = np.r_[0, np.nonzero(k[1:] != k[:-1])[0] + 1]
    grp = np.repeat(start, np.diff(np.r_[start, len(k)]))
    out = np.empty(len(keys), dtype=np.int64)
    out[order] = np.arange(len(k)) - grp
    return out


def distinct_before(ctx, first):
    order = np.argsort(ctx, kind="stable")
    f = first[order]
    cs = np.cumsum(f) - f
    cv = ctx[order]
    start = np.r_[0, np.nonzero(cv[1:] != cv[:-1])[0] + 1]
    base = np.repeat(cs[start], np.diff(np.r_[start, len(cv)]))
    out = np.empty(len(ctx), dtype=np.int64)
    out[order] = cs - base
    return out


def context_keys(t, k, m):
    """unique int key per (previous k tokens); positions < k get shorter
    contexts, marked by the sentinel value m."""
    key = np.zeros(len(t), dtype=np.int64)
    for j in range(1, k + 1):
        prev = np.r_[np.full(j, m), t[:-j]]
        key = key * (m + 1) + prev
    return key


def code_bytes(t, m, order, a=3.0, b=0.5, per_symbol=False):
    t = np.asarray(t, dtype=np.int64)
    p = np.full(len(t), 1.0 / m)
    for k in range(0, order + 1):
        ctx = context_keys(t, k, m) if k else np.zeros(len(t), dtype=np.int64)
        key = ctx * (m + 1) + t
        if k >= 3:  # avoid int64 overflow on big grids
            key = np.unique(ctx, return_inverse=True)[1].astype(np.int64) * (m + 1) + t
            ctx = key // (m + 1)
        c = occ_before(key)
        n = occ_before(ctx)
        d = distinct_before(ctx, (c == 0).astype(np.int64))
        w = a * d + b
        p = (c + w * p) / (n + w)
    bits = -np.log2(p)
    return bits if per_symbol else float(bits.sum() / 8)
