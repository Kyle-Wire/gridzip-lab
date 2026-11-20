"""Reorganise, then compress: the Burrows-Wheeler transform on BYTES vs on
GRID COORDINATES.

The BWT sorts every rotation of the sequence and keeps the last column. It
is a permutation of the input that groups symbols by what follows them, and
it is undone from ONE number (the row of the original sequence). After it:
move-to-front ranks, coded with an adaptive context model (v3's coder).

Real inverse transform; the round trip is checked.
"""
import math, sys, time, bz2
import numpy as np
from numba import njit
from gridzip import learn_grid, apply_grid, map_bytes
import v3


def suffix_array(t):
    """prefix doubling; t gets a unique smallest sentinel appended."""
    s = np.r_[t + 1, 0].astype(np.int64)
    n = len(s)
    rank = s.copy()
    k = 1
    sa = np.argsort(rank, kind="stable")
    while True:
        r2 = np.where(np.arange(n) + k < n, rank[np.minimum(np.arange(n) + k, n - 1)], -1)
        sa = np.lexsort((r2, rank))
        key_r, key_2 = rank[sa], r2[sa]
        diff = np.r_[1, (key_r[1:] != key_r[:-1]) | (key_2[1:] != key_2[:-1])]
        new = np.empty(n, dtype=np.int64)
        new[sa] = np.cumsum(diff) - 1
        rank = new
        if rank.max() == n - 1:
            return sa, s
        k *= 2


def bwt(t):
    sa, s = suffix_array(t)
    last = s[(sa - 1) % len(s)]
    primary = int(np.nonzero(sa == 0)[0][0])
    return last, primary  # last contains the sentinel 0 once


@njit(cache=True)
def _walk(last, lf, primary, n):
    out = np.empty(n, dtype=np.int64)
    j = primary
    for i in range(n - 1, -1, -1):
        out[i] = last[j]
        j = lf[j]
    return out


def ibwt(last, primary):
    order = np.argsort(last, kind="stable")      # first column positions
    lf = np.empty(len(last), dtype=np.int64)
    lf[order] = np.arange(len(last))
    full = _walk(last, lf, primary, len(last))   # rotation ending in sentinel
    # full is s rotated so that it ends just before row `primary`; s ends with 0
    return full[:-1] - 1 if full[-1] == 0 else None


@njit(cache=True)
def mtf(x, m):
    table = np.arange(m)
    out = np.empty(len(x), dtype=np.int64)
    for i in range(len(x)):
        v = x[i]
        j = 0
        while table[j] != v:
            j += 1
        out[i] = j
        while j > 0:
            table[j] = table[j - 1]; j -= 1
        table[0] = v
    return out


def size_bwt(t, m, order):
    last, primary = bwt(t)
    ranks = mtf(last, m + 1)
    return v3.code_bytes(ranks, m + 1, order) + 8, last, primary  # +8: the primary index


if __name__ == "__main__":
    path, n = sys.argv[1], int(sys.argv[2])
    data = open(path, "rb").read()[:n]
    print(f"n={n:,} bzip2={len(bz2.compress(data, 9)):,}", flush=True)
    raw = np.frombuffer(data, np.uint8).astype(np.int64)
    for cells in [0, 256, 1024, 4096]:
        if cells:
            M, _ = learn_grid(data[:2_000_000], cells)
            t = apply_grid(data, M); mp = map_bytes(M)
        else:
            t, mp = raw, 0.0
        t0 = time.time()
        for order in [1, 2]:
            sz, last, primary = size_bwt(t, 256 + cells, order)
            print(f"cells={cells:5d} ctx={order}: BWT+MTF {sz + mp:12,.0f}  ({time.time()-t0:.0f}s)", flush=True)
        back = ibwt(last, primary)
        assert back is not None and np.array_equal(back, t), "inverse BWT failed"
    print("inverse BWT verified for every grid size", flush=True)
