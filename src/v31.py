"""GridZip v3.1: learned mixing of context lengths.

v3 blends contexts with a fixed rule. v3.1 keeps v3's chain (q_k = the
blended estimate using contexts up to length k) and adds a mixer that learns
online how much to trust each q_k:

    P(s) = sum_k softmax(z[set])_k * q_k(s)

Each q_k is a normalised distribution, so P is too. The weight set is chosen
by how often the two longest contexts have been seen, which tells the mixer
how mature each one is. Weights are updated after every symbol by gradient
descent on the code length. Everything is causal, so a decoder can replay it.
"""
import numpy as np
from numba import njit
from v3 import occ_before, distinct_before


def chain(t, m, order, a=3.0, b=0.5):
    t = np.asarray(t, dtype=np.int64)
    Q = np.empty((order + 1, len(t)))
    Nc = np.empty((order + 1, len(t)), dtype=np.int64)
    p = np.full(len(t), 1.0 / m)
    ctx = np.zeros(len(t), dtype=np.int64)
    for k in range(order + 1):
        if k:
            prev = np.r_[np.full(k, m), t[:-k]]
            ctx = np.unique(ctx * (m + 1) + prev, return_inverse=True)[1].astype(np.int64)
        key = ctx * (m + 1) + t
        c = occ_before(key)
        n = occ_before(ctx)
        d = distinct_before(ctx, (c == 0).astype(np.int64))
        w = a * d + b
        p = (c + w * p) / (n + w)
        Q[k] = p
        Nc[k] = n
    return Q, Nc


@njit(cache=True)
def mix(Q, Nc, lr):
    K, N = Q.shape
    nb = 12
    Z = np.zeros((nb * nb, K))
    bits = 0.0
    lam = np.empty(K)
    for i in range(N):
        b1 = min(nb - 1, int(np.log2(Nc[K - 1, i] + 1)))
        b2 = min(nb - 1, int(np.log2(Nc[K - 2, i] + 1))) if K > 1 else 0
        s = b1 * nb + b2
        z = Z[s]
        mx = z.max()
        tot = 0.0
        for k in range(K):
            lam[k] = np.exp(z[k] - mx); tot += lam[k]
        P = 0.0
        for k in range(K):
            lam[k] /= tot; P += lam[k] * Q[k, i]
        bits -= np.log2(P)
        for k in range(K):
            z[k] += lr * lam[k] * (Q[k, i] - P) / P
    return bits / 8


def code_bytes(t, m, order, lr=0.02, a=3.0, b=0.5):
    Q, Nc = chain(t, m, order, a, b)
    return mix(Q, Nc, lr)
