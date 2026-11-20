"""A map shaped for CODE LENGTH, not for distance.

For each grid cell g (the longest cell ending at the previous byte, exactly
what the compressor sees), count what byte comes next: N[g, b]. Group cells
into R regions so that coding the next byte from the REGION costs as few bits
as possible:

    minimise  sum_g sum_b N[g,b] * -log2 q_r(g)(b),   q_r = pooled next-byte
                                                      distribution of region r

This is k-means with cross-entropy (KL divergence) as the distance, the
information-theoretic clustering of Dhillon et al. (2003). Every step lowers
the code length, the quantity we actually care about.
"""
import numpy as np
from numba import njit
from v4 import build_automaton


@njit(cache=True)
def _counts(buf, G, out, m):
    N = np.zeros((m, 256), dtype=np.float64)
    st = 0
    g = 0
    for i in range(len(buf)):
        N[g, buf[i]] += 1.0
        st = G[st, buf[i]]
        g = out[st]
    return N


def code_bits(N, lab, R, alpha=0.5):
    Q = np.zeros((R, 256))
    np.add.at(Q, lab, N)
    Q = (Q + alpha) / (Q + alpha).sum(1, keepdims=True)
    return float(-(N * np.log2(Q[lab])).sum())


def codelen_map(data, cells, R, iters=40, seed=0, alpha=0.5):
    G, out, lens = build_automaton(cells)
    buf = np.frombuffer(data, dtype=np.uint8).astype(np.int64)
    m = len(cells)
    N = _counts(buf, G, out, m)
    rng = np.random.default_rng(seed)
    # init: the R busiest cells seed the regions, everything else joins its nearest
    tot = N.sum(1)
    seeds = np.argsort(-tot)[:R]
    Q = (N[seeds] + alpha) / (N[seeds] + alpha).sum(1, keepdims=True)
    lab = (N @ np.log2(Q).T).argmax(1)
    hist = []
    for _ in range(iters):
        Q = np.zeros((R, 256)); np.add.at(Q, lab, N)
        empty = Q.sum(1) == 0
        if empty.any():  # revive empty regions with the costliest cells
            worst = np.argsort(-(tot))[: empty.sum()]
            Q[empty] = N[worst]
        Q = (Q + alpha) / (Q + alpha).sum(1, keepdims=True)
        new = (N @ np.log2(Q).T).argmax(1)  # min cross-entropy = max log-likelihood
        hist.append(code_bits(N, new, R, alpha))
        if np.array_equal(new, lab):
            break
        lab = new
    exact = code_bits(N, np.arange(m), m, alpha)
    rnd = code_bits(N, rng.integers(0, R, m), R, alpha)
    one = code_bits(N, np.zeros(m, dtype=np.int64), 1, alpha)
    info = dict(R=R, iters=len(hist), bits_per_byte_region=hist[-1] / len(buf),
                bits_per_byte_exact_cell=exact / len(buf), bits_per_byte_random=rnd / len(buf),
                bits_per_byte_no_context=one / len(buf))
    return lab, info
