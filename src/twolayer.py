"""The "screenshot" idea as a lossless two-part code.

Layer A (the screenshot): the text with every letter lowercased. Lossy on its
own, and smaller to code.
Layer B (the correction): one case bit per letter.

The point: when B is coded, ALL of A is already known, including what comes
AFTER each letter. So the case of a letter can be predicted from the whole
word, the next word, and what follows ("Mr. Smith", "New York", sentence
starts), which a one-pass compressor only sees in the past.

Layer B is coded with a small context-mixing model; its size is the exact
adaptive code length (sum of -log2 p as the model learns online).
"""
import os, subprocess, sys
import numpy as np
from numba import njit


@njit(cache=True)
def hsh(a, b):
    h = a * 0x5BD1E995 + (b ^ (b >> 15)) * 0x27D4EB2F + 0x165667B1
    return (h ^ (h >> 29)) & 0x7FFFFFFFFFFF


@njit(cache=True)
def case_bits(A, upper, use_future):
    n = len(A)
    TB = 22; size = 1 << TB; mask = size - 1
    NM = 6
    P = np.full((NM, size), 0.5, dtype=np.float32)
    Cn = np.zeros((NM, size), dtype=np.uint16)
    W = np.full((3 * 2 * 2, NM + 1), 0.3)
    # word hash for every position (whole word: past AND future letters)
    word = np.zeros(n, dtype=np.int64); wpos = np.zeros(n, dtype=np.int64)
    prevw = np.zeros(n, dtype=np.int64); nextw = np.zeros(n, dtype=np.int64)
    i = 0; last = 0
    while i < n:
        if 97 <= A[i] <= 122:
            j = i; h = 0
            while j < n and 97 <= A[j] <= 122:
                h = hsh(h, A[j]); j += 1
            for k in range(i, j):
                word[k] = h; wpos[k] = k - i; prevw[k] = last
            last = h; i = j
        else:
            i += 1
    nxt = 0
    for k in range(n - 1, -1, -1):
        if 97 <= A[k] <= 122:
            if k + 1 < n and not (97 <= A[k + 1] <= 122):
                pass
            nextw[k] = nxt
            if k == 0 or not (97 <= A[k - 1] <= 122):
                nxt = word[k]
    bits = 0.0
    x = np.zeros(NM + 1)
    idx = np.zeros(NM, dtype=np.int64)
    prev_case = 0; first_case = 0
    for i in range(n):
        c = A[i]
        if not (97 <= c <= 122):
            continue
        b = upper[i]
        wp = wpos[i]
        if wp == 0:
            first_case = 0
        p1 = A[i - 1] if i >= 1 else 0
        p2 = A[i - 2] if i >= 2 else 0
        p3 = A[i - 3] if i >= 3 else 0
        n1 = A[i + 1] if (use_future and i + 1 < n) else 0
        n2 = A[i + 2] if (use_future and i + 2 < n) else 0
        wd = word[i] if use_future else hsh(prevw[i], wp)  # past-only: no whole word
        nw = nextw[i] if use_future else 0
        ctx = np.empty(NM, dtype=np.int64)
        ctx[0] = hsh(1, hsh(wd, min(wp, 3)))
        ctx[1] = hsh(2, p1 | (p2 << 8) | (p3 << 16) | (n1 << 24) | (n2 << 32))
        ctx[2] = hsh(3, hsh(prevw[i], wd) * 4 + min(wp, 3))
        ctx[3] = hsh(4, min(wp, 3) * 4 + prev_case * 2 + first_case)
        ctx[4] = hsh(5, hsh(wd, nw) * 4 + min(wp, 3))
        ctx[5] = hsh(6, hsh(p1 | (p2 << 8), wd) * 4 + first_case * 2 + (1 if wp == 0 else 0))
        for k in range(NM):
            idx[k] = ctx[k] & mask
            q = P[k, idx[k]]
            if q < 0.001: q = 0.001
            if q > 0.999: q = 0.999
            x[k] = np.log(q / (1 - q))
        x[NM] = 0.3
        s = min(wp, 2) * 4 + prev_case * 2 + first_case
        dot = 0.0
        for k in range(NM + 1): dot += W[s, k] * x[k]
        if dot > 30: dot = 30.0
        if dot < -30: dot = -30.0
        p = 1 / (1 + np.exp(-dot))
        if p < 1e-6: p = 1e-6
        if p > 1 - 1e-6: p = 1 - 1e-6
        bits -= np.log2(p if b else 1 - p)
        err = b - p
        for k in range(NM + 1): W[s, k] += 0.02 * err * x[k]
        for k in range(NM):
            cc = Cn[k, idx[k]]
            P[k, idx[k]] += (b - P[k, idx[k]]) / (cc + 1.5)
            if cc < 1023: Cn[k, idx[k]] = cc + 1
        prev_case = b
        if wp == 0: first_case = b
    return bits / 8


if __name__ == "__main__":
    src, out = sys.argv[1], sys.argv[2]
    raw = np.frombuffer(open(src, "rb").read(), dtype=np.uint8)
    upper = ((raw >= 65) & (raw <= 90)).astype(np.int64)
    A = raw.copy(); A[upper == 1] += 32
    open(out, "wb").write(A.tobytes())
    A64 = A.astype(np.int64)
    future = case_bits(A64, upper, True)
    past = case_bits(A64, upper, False)
    print(f"letters={int(((A >= 97) & (A <= 122)).sum()):,} uppercase={int(upper.sum()):,}")
    print(f"case bits, past-only context : {past:,.0f} B")
    print(f"case bits, with FUTURE context: {future:,.0f} B")
