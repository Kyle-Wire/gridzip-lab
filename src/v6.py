"""GridZip v6 (= v5 with the map at two zoom levels: coarse regA, fine regB): bitwise context-mixing compressor with GRID-LOCATION contexts.

Each bit is predicted by several models; a learned mixer combines them; a
binary arithmetic coder writes the bits. Real encoder AND decoder: run()
in decode mode must rebuild the input exactly.

Byte models   : order 1,2,3,4,6 contexts, a word model, a long-match model.
Grid models   : "where am I on the grid": the longest grid cell that ends at
                the last byte (found causally with an Aho-Corasick automaton
                built from the grid), alone and paired with the cell before it.

Flags let any group be switched off, so the value of the grid signal is
measured on top of a strong byte model, not against a weak one.
"""
import sys, time
from collections import deque
import numpy as np
from numba import njit

TB = 22  # log2 table size per model


def build_automaton(cells):
    """cells: list of bytes. Returns goto[S,256], longest-cell-ending-here[S], cell lengths."""
    goto = [dict()]
    term = [-1]
    for cid, s in enumerate(cells):
        st = 0
        for ch in s:
            nx = goto[st].get(ch)
            if nx is None:
                nx = len(goto); goto[st][ch] = nx; goto.append(dict()); term.append(-1)
            st = nx
        term[st] = cid
    S = len(goto)
    G = np.zeros((S, 256), dtype=np.int32)
    fail = np.zeros(S, dtype=np.int32)
    out = np.full(S, -1, dtype=np.int32)
    q = deque()
    for ch in range(256):
        nx = goto[0].get(ch, 0)
        G[0, ch] = nx
        if nx:
            fail[nx] = 0; q.append(nx)
    while q:
        st = q.popleft()
        out[st] = term[st] if term[st] >= 0 else out[fail[st]]
        for ch in range(256):
            nx = goto[st].get(ch)
            if nx is None:
                G[st, ch] = G[fail[st], ch]
            else:
                fail[nx] = G[fail[st], ch]; G[st, ch] = nx; q.append(nx)
    lens = np.array([len(c) for c in cells], dtype=np.int32)
    return G, out, lens


@njit(cache=True)
def squash(x):
    if x > 20.0: x = 20.0
    if x < -20.0: x = -20.0
    return 1.0 / (1.0 + np.exp(-x))


@njit(cache=True)
def stretch(p):
    return np.log(p / (1.0 - p))


@njit(cache=True)
def hsh(a, b):
    h = a * 0x5BD1E995 + (b ^ (b >> 15)) * 0x27D4EB2F + 0x165667B1
    return (h ^ (h >> 29)) & 0x7FFFFFFFFFFF


@njit(cache=True)
def run(buf, n, decode, code, G, out, lens, use_bytes, use_grid, lr0, limit, use_space, reg, regB):
    NM = 16
    mask = (1 << 22) - 1
    P = np.full((NM, 1 << 22), 0.5, dtype=np.float32)
    C = np.zeros((NM, 1 << 22), dtype=np.uint16)
    W = np.zeros((256 * 3, NM + 2))
    for s in range(W.shape[0]):
        for k in range(NM + 2):
            W[s, k] = 0.3
    apm = np.empty((256 * 33,), dtype=np.float32)
    for c in range(256):
        for j in range(33):
            apm[c * 33 + j] = squash((j - 16) / 2.0)
    MM = 1 << 22
    mtab = np.zeros(MM, dtype=np.int64)
    mprob = np.full(64 * 2, 0.5)
    mcnt = np.zeros(64 * 2)
    gh = np.zeros(n + 1, dtype=np.int32)   # longest cell ending at i-1
    st = 0
    base = np.zeros(NM, dtype=np.int64)
    active = np.zeros(NM, dtype=np.bool_)
    idx = np.zeros(NM, dtype=np.int64)
    x = np.zeros(NM + 2)
    x1 = np.uint64(0); x2 = np.uint64(0xFFFFFFFF); xv = np.uint64(0)
    outp = 0; inp = 0
    if decode:
        for _ in range(4):
            xv = (xv << np.uint64(8)) | np.uint64(code[inp]); inp += 1
    word = 0; word1 = 0
    mptr = 0; mlen = 0
    bits = 0.0
    for i in range(n):
        # ---- byte-boundary contexts
        c1 = buf[i - 1] if i >= 1 else 0
        c2 = buf[i - 2] if i >= 2 else 0
        c3 = buf[i - 3] if i >= 3 else 0
        c4 = buf[i - 4] if i >= 4 else 0
        c5 = buf[i - 5] if i >= 5 else 0
        c6 = buf[i - 6] if i >= 6 else 0
        for k in range(NM):
            active[k] = False
        if use_bytes:
            base[0] = hsh(1, c1); active[0] = True
            base[1] = hsh(2, c1 | (c2 << 8)); active[1] = True
            base[2] = hsh(3, c1 | (c2 << 8) | (c3 << 16)); active[2] = True
            base[3] = hsh(4, c1 | (c2 << 8) | (c3 << 16) | (c4 << 24)); active[3] = True
            base[4] = hsh(5, hsh(c1 | (c2 << 8) | (c3 << 16), c4 | (c5 << 8) | (c6 << 16))); active[4] = True
            base[5] = hsh(6, hsh(word, word1)); active[5] = True
            base[6] = hsh(7, word); active[6] = True
        base[7] = hsh(8, 0); active[7] = True  # order 0
        if use_grid:
            g1 = gh[i]
            g2 = gh[i - lens[g1]] if i - lens[g1] >= 0 else 0
            base[8] = hsh(9, g1); active[8] = True
            base[9] = hsh(10, g1 * 8192 + g2); active[9] = True
            base[10] = hsh(11, hsh(g1, c1 | (c2 << 8) | (c3 << 16))); active[10] = True
        if use_space:
            g1 = gh[i]
            j2 = i - lens[g1]
            g2 = gh[j2] if j2 >= 0 else 0
            j3 = j2 - lens[g2]
            g3 = gh[j3] if j3 >= 0 else 0
            r1 = reg[g1]; r2 = reg[g2]; r3 = reg[g3]
            base[11] = hsh(13, r1 * 65536 + r2 * 256 + r3); active[11] = True
            base[12] = hsh(14, g1 * 65536 + r2 * 256 + r3); active[12] = True
            base[13] = hsh(15, hsh(r1 * 256 + r2, c1)); active[13] = True
            f1 = regB[g1]; f2 = regB[g2]; f3 = regB[g3]
            base[14] = hsh(16, f1 * 1048576 + f2 * 1024 + f3); active[14] = True
            base[15] = hsh(17, hsh(f1 * 1024 + f2, r3 * 256 + c1)); active[15] = True
        c0 = 1
        expected = buf[mptr] if (use_bytes and mlen > 0) else -1
        byte = 0
        for bitpos in range(8):
            ni = 0
            for k in range(NM):
                if active[k]:
                    idx[k] = (base[k] + c0 * 0x2545F491) & mask
                    p = P[k, idx[k]]
                    if p < 0.0005: p = 0.0005
                    if p > 0.9995: p = 0.9995
                    x[k] = stretch(p)
                else:
                    x[k] = 0.0
            # match input
            mi = -1
            if expected >= 0 and ((expected + 256) >> (8 - bitpos)) == c0:
                eb = (expected >> (7 - bitpos)) & 1
                mi = min(mlen, 31) * 2 + eb
                pm = mprob[mi]
                if pm < 0.001: pm = 0.001
                if pm > 0.999: pm = 0.999
                x[NM] = stretch(pm)
            else:
                x[NM] = 0.0
            x[NM + 1] = 0.25
            ws = c0 + 256 * (0 if mi < 0 else (1 if mlen < 16 else 2))
            dot = 0.0
            for k in range(NM + 2):
                dot += W[ws, k] * x[k]
            pmix = squash(dot)
            # APM on order-1
            s = stretch(min(max(pmix, 0.0001), 0.9999))
            pos = (s + 8.0) * 2.0
            if pos < 0: pos = 0.0
            if pos > 31.999: pos = 31.999
            j = int(pos); fr = pos - j
            ai = (c1 * 0 + c0) * 33 + j  # context: partial byte
            pa = apm[ai] * (1 - fr) + apm[ai + 1] * fr
            pf = 0.25 * pmix + 0.75 * pa
            p12 = int(pf * 4096)
            if p12 < 1: p12 = 1
            if p12 > 4095: p12 = 4095
            # ---- code the bit
            xmid = x1 + ((x2 - x1) >> np.uint64(12)) * np.uint64(p12)
            if decode:
                bit = 1 if xv <= xmid else 0
            else:
                bit = (buf[i] >> (7 - bitpos)) & 1
            if bit:
                x2 = xmid
            else:
                x1 = xmid + np.uint64(1)
            while ((x1 ^ x2) & np.uint64(0xFF000000)) == 0:
                if not decode:
                    code[outp] = np.uint8(x2 >> np.uint64(24)); outp += 1
                x1 = (x1 << np.uint64(8)) & np.uint64(0xFFFFFFFF)
                x2 = ((x2 << np.uint64(8)) & np.uint64(0xFFFFFFFF)) | np.uint64(255)
                if decode:
                    xv = ((xv << np.uint64(8)) & np.uint64(0xFFFFFFFF)) | np.uint64(code[inp] if inp < len(code) else 0); inp += 1
            bits -= np.log2(p12 / 4096.0 if bit else 1 - p12 / 4096.0)
            # ---- update
            err = bit - pmix
            lr = lr0
            for k in range(NM + 2):
                W[ws, k] += lr * err * x[k]
            for k in range(NM):
                if active[k]:
                    cc = C[k, idx[k]]
                    P[k, idx[k]] += (bit - P[k, idx[k]]) / (cc + 1.5)
                    if cc < limit: C[k, idx[k]] = cc + 1
            if mi >= 0:
                mcnt[mi] = min(mcnt[mi] + 1, 1000.0)
                mprob[mi] += (bit - mprob[mi]) / (mcnt[mi] + 1.5)
            apm[ai] += (bit - apm[ai]) * (1 - fr) * 0.02
            apm[ai + 1] += (bit - apm[ai + 1]) * fr * 0.02
            c0 = (c0 << 1) | bit
            byte = (byte << 1) | bit
        if decode:
            buf[i] = byte
        # ---- after the byte: advance or (re)find a long match, THEN store
        if mlen > 0:
            if expected == byte:
                mlen += 1; mptr += 1
            else:
                mlen = 0
        if i + 1 >= 6:
            h = hsh(12, byte | (c1 << 8) | (c2 << 16) | (c3 << 24) | (c4 << 32) | (c5 << 40)) & (MM - 1)
            if mlen == 0 and mtab[h] > 0:
                mptr = mtab[h]; mlen = 1
            mtab[h] = i + 1
        if (byte >= 65 and byte <= 90) or (byte >= 97 and byte <= 122):
            word = hsh(word, byte | 32)
        elif word != 0:
            word1 = word; word = 0
        st = G[st, byte]
        gh[i + 1] = out[st]
    if not decode:
        for _ in range(4):
            code[outp] = np.uint8(x1 >> np.uint64(24)); outp += 1
            x1 = (x1 << np.uint64(8)) & np.uint64(0xFFFFFFFF)
    return outp, bits / 8


def compress(data, cells, use_bytes=True, use_grid=True, lr=0.004, limit=255, verify=False,
             use_space=False, reg=None, regB=None):
    G, out, lens = build_automaton(cells)
    reg = np.zeros(len(cells), dtype=np.int64) if reg is None else np.asarray(reg, dtype=np.int64)
    regB = reg.copy() if regB is None else np.asarray(regB, dtype=np.int64)
    buf = np.frombuffer(data, dtype=np.uint8).astype(np.int64).copy()
    code = np.zeros(len(data) + 1024, dtype=np.uint8)
    nbytes, ideal = run(buf, len(data), False, code, G, out, lens, use_bytes, use_grid, lr, limit, use_space, reg, regB)
    if verify:
        back = np.zeros(len(data), dtype=np.int64)
        run(back, len(data), True, code[:nbytes], G, out, lens, use_bytes, use_grid, lr, limit, use_space, reg, regB)
        assert bytes(back.astype(np.uint8)) == data, "round trip failed"
    return nbytes, ideal


def spatial_map(tokens, n_cells, regions, dims=16, seed=0):
    """Place every cell on a map: cells that are followed and preceded by
    similar cells land near each other. Returns region id per cell."""
    t = np.asarray(tokens, dtype=np.int64)
    m = n_cells
    A = np.zeros((m, m), dtype=np.float32)
    np.add.at(A, (t[:-1], t[1:]), 1.0)
    F = np.hstack([A, A.T])  # what follows and what precedes each cell
    F = np.log1p(F)
    F /= np.linalg.norm(F, axis=1, keepdims=True) + 1e-9
    U, S, _ = np.linalg.svd(F - F.mean(0), full_matrices=False)
    X = U[:, :dims] * S[:dims]
    X /= np.linalg.norm(X, axis=1, keepdims=True) + 1e-9
    rng = np.random.default_rng(seed)
    Cn = X[rng.choice(m, regions, replace=False)]
    for _ in range(30):  # spherical k-means
        lab = (X @ Cn.T).argmax(1)
        for r in range(regions):
            pts = X[lab == r]
            if len(pts):
                c = pts.sum(0); Cn[r] = c / (np.linalg.norm(c) + 1e-9)
    return lab
