"""GridZip v8: a stronger context-mixing engine carrying all the grid ideas.

Changes from v7 (each was measured, see README):
  * tables 2^TB slots per model (default 2^24, v4-v7 used 2^22)
  * more byte contexts: order 5, order 8, sparse (skip-1, skip-2), word
    unigram/bigram, and a "bracket" context (last structural char of wiki
    markup: [ { < | = ' plus order 1)
  * two-layer mixing: three mixers selected by different things (partial
    byte; match state x bit position; how many high orders have been seen),
    combined by a small final mixer
  * two APM/SSE stages (order 0 and order 1), averaged
Grid, spatial map (coarse + fine), symmetry-orbit contexts as in v7; an
optional extra map regC (e.g. the code-length-shaped map) adds 2 contexts.
Real arithmetic coder; decode() rebuilds the input.
"""
import numpy as np
from numba import njit
from v4 import build_automaton

NMAX = 32


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
def orbit(c):
    if c >= 65 and c <= 90:
        return c + 32
    if c >= 48 and c <= 57:
        return 48
    return c


@njit(cache=True)
def run(buf, n, decode, code, G, out, lens, reg, regB, regC, flags, TB, lr0, limit):
    use_grid, use_space, use_orbit, use_c, use_extra = flags[0], flags[1], flags[2], flags[3], flags[4]
    size = 1 << TB
    mask = size - 1
    NM = NMAX
    P = np.full((NM, size), 0.5, dtype=np.float32)
    C = np.zeros((NM, size), dtype=np.uint8)
    NI = NM + 2
    W1 = np.full((256, NI), 0.2)
    W2 = np.full((3 * 8 * 2, NI), 0.2)
    W3 = np.full((11 * 8, NI), 0.2)
    WF = np.full((256, 4), 0.33)
    apm0 = np.empty(256 * 33, dtype=np.float32)
    apm1 = np.empty(65536 * 33, dtype=np.float32)
    for j in range(33):
        v = squash((j - 16) / 2.0)
        for c in range(256): apm0[c * 33 + j] = v
        for c in range(65536): apm1[c * 33 + j] = v
    MM = 1 << 22
    mtab = np.zeros(MM, dtype=np.int64)
    mprob = np.full(64 * 2, 0.5); mcnt = np.zeros(64 * 2)
    gh = np.zeros(n + 1, dtype=np.int32)
    st = 0
    base = np.zeros(NM, dtype=np.int64)
    active = np.zeros(NM, dtype=np.bool_)
    idx = np.zeros(NM, dtype=np.int64)
    x = np.zeros(NI)
    x1 = np.uint64(0); x2 = np.uint64(0xFFFFFFFF); xv = np.uint64(0)
    outp = 0; inp = 0
    if decode:
        for _ in range(4):
            xv = (xv << np.uint64(8)) | np.uint64(code[inp]); inp += 1
    word = 0; word1 = 0; brk = 0
    mptr = 0; mlen = 0
    for i in range(n):
        c1 = buf[i - 1] if i >= 1 else 0
        c2 = buf[i - 2] if i >= 2 else 0
        c3 = buf[i - 3] if i >= 3 else 0
        c4 = buf[i - 4] if i >= 4 else 0
        c5 = buf[i - 5] if i >= 5 else 0
        c6 = buf[i - 6] if i >= 6 else 0
        c7 = buf[i - 7] if i >= 7 else 0
        c8 = buf[i - 8] if i >= 8 else 0
        for k in range(NM): active[k] = False
        o3 = c1 | (c2 << 8) | (c3 << 16)
        o4 = o3 | (c4 << 24)
        base[0] = hsh(1, c1)
        base[1] = hsh(2, c1 | (c2 << 8))
        base[2] = hsh(3, o3)
        base[3] = hsh(4, o4)
        base[4] = hsh(5, o4 | (c5 << 32))
        base[5] = hsh(6, hsh(o3, c4 | (c5 << 8) | (c6 << 16)))
        base[6] = hsh(7, hsh(o4, c5 | (c6 << 8) | (c7 << 16) | (c8 << 24)))
        base[7] = hsh(8, hsh(word, word1))
        base[8] = hsh(9, word)
        base[9] = hsh(10, 0)
        for k in range(10): active[k] = True
        if use_extra:
            base[10] = hsh(11, c2 | (c3 << 8)); active[10] = True            # skip-1
            base[11] = hsh(12, c3 | (c4 << 8)); active[11] = True            # skip-2
            base[12] = hsh(13, brk * 256 + c1); active[12] = True            # markup
            base[13] = hsh(14, hsh(word, c1)); active[13] = True
        g1 = gh[i]
        j2 = i - lens[g1]
        g2 = gh[j2] if j2 >= 0 else 0
        j3 = j2 - lens[g2]
        g3 = gh[j3] if j3 >= 0 else 0
        if use_grid:
            base[14] = hsh(15, g1); active[14] = True
            base[15] = hsh(16, g1 * 8192 + g2); active[15] = True
            base[16] = hsh(17, hsh(g1, o3)); active[16] = True
        if use_space:
            r1 = reg[g1]; r2 = reg[g2]; r3 = reg[g3]
            f1 = regB[g1]; f2 = regB[g2]; f3 = regB[g3]
            base[17] = hsh(18, r1 * 65536 + r2 * 256 + r3); active[17] = True
            base[18] = hsh(19, g1 * 65536 + r2 * 256 + r3); active[18] = True
            base[19] = hsh(20, hsh(r1 * 256 + r2, c1)); active[19] = True
            base[20] = hsh(21, f1 * 1048576 + f2 * 1024 + f3); active[20] = True
            base[21] = hsh(22, hsh(f1 * 1024 + f2, r3 * 256 + c1)); active[21] = True
        if use_orbit:
            q1 = orbit(c1); q2 = orbit(c2); q3 = orbit(c3)
            q4 = orbit(c4); q5 = orbit(c5); q6 = orbit(c6)
            base[22] = hsh(23, q1 | (q2 << 8) | (q3 << 16)); active[22] = True
            base[23] = hsh(24, hsh(q1 | (q2 << 8) | (q3 << 16), q4 | (q5 << 8) | (q6 << 16))); active[23] = True
            base[24] = hsh(25, q1 | (q2 << 8)); active[24] = True
        if use_c:
            k1 = regC[g1]; k2 = regC[g2]
            base[25] = hsh(26, k1 * 4096 + k2); active[25] = True
            base[26] = hsh(27, hsh(k1, c1 | (c2 << 8))); active[26] = True
        # how many of the long contexts have been seen before (mixer selector)
        c0 = 1
        expected = buf[mptr] if mlen > 0 else -1
        byte = 0
        seen = 0
        for bitpos in range(8):
            for k in range(NM):
                if active[k]:
                    idx[k] = (base[k] + c0 * 0x2545F491) & mask
                    p = P[k, idx[k]]
                    if p < 0.0005: p = 0.0005
                    if p > 0.9995: p = 0.9995
                    x[k] = stretch(p)
                else:
                    x[k] = 0.0
            if bitpos == 0:
                seen = 0
                for k in range(2, 7):
                    if C[k, idx[k]] > 0: seen += 1
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
            s1 = c0
            ms = 0 if mi < 0 else (1 if mlen < 16 else 2)
            s2 = (ms * 8 + bitpos) * 2 + (0 if mi < 0 else (mi & 1))
            s3 = (seen + 5 * (1 if mi >= 0 else 0)) * 8 + bitpos
            d1 = 0.0; d2 = 0.0; d3 = 0.0
            for k in range(NI):
                d1 += W1[s1, k] * x[k]; d2 += W2[s2, k] * x[k]; d3 += W3[s3, k] * x[k]
            p1 = squash(d1); p2 = squash(d2); p3 = squash(d3)
            y1 = stretch(min(max(p1, 1e-4), 1 - 1e-4))
            y2 = stretch(min(max(p2, 1e-4), 1 - 1e-4))
            y3 = stretch(min(max(p3, 1e-4), 1 - 1e-4))
            df = WF[c0, 0] * y1 + WF[c0, 1] * y2 + WF[c0, 2] * y3 + WF[c0, 3] * 0.25
            pmix = squash(df)
            s = stretch(min(max(pmix, 1e-4), 1 - 1e-4))
            pos = (s + 8.0) * 2.0
            if pos < 0: pos = 0.0
            if pos > 31.999: pos = 31.999
            jj = int(pos); fr = pos - jj
            a0 = c0 * 33 + jj
            a1 = (c1 * 256 + c0) * 33 + jj
            pa0 = apm0[a0] * (1 - fr) + apm0[a0 + 1] * fr
            pa1 = apm1[a1] * (1 - fr) + apm1[a1 + 1] * fr
            pf = 0.4 * pmix + 0.3 * pa0 + 0.3 * pa1
            p12 = int(pf * 4096)
            if p12 < 1: p12 = 1
            if p12 > 4095: p12 = 4095
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
            # updates
            e1 = bit - p1; e2 = bit - p2; e3 = bit - p3
            for k in range(NI):
                W1[s1, k] += lr0 * e1 * x[k]
                W2[s2, k] += lr0 * e2 * x[k]
                W3[s3, k] += lr0 * e3 * x[k]
            ef = bit - pmix
            WF[c0, 0] += 0.002 * ef * y1; WF[c0, 1] += 0.002 * ef * y2
            WF[c0, 2] += 0.002 * ef * y3; WF[c0, 3] += 0.002 * ef * 0.25
            for k in range(NM):
                if active[k]:
                    cc = C[k, idx[k]]
                    P[k, idx[k]] += (bit - P[k, idx[k]]) / (cc + 1.5)
                    if cc < limit: C[k, idx[k]] = cc + 1
            if mi >= 0:
                mcnt[mi] = min(mcnt[mi] + 1, 1000.0)
                mprob[mi] += (bit - mprob[mi]) / (mcnt[mi] + 1.5)
            apm0[a0] += (bit - apm0[a0]) * (1 - fr) * 0.02
            apm0[a0 + 1] += (bit - apm0[a0 + 1]) * fr * 0.02
            apm1[a1] += (bit - apm1[a1]) * (1 - fr) * 0.02
            apm1[a1 + 1] += (bit - apm1[a1 + 1]) * fr * 0.02
            c0 = (c0 << 1) | bit
            byte = (byte << 1) | bit
        if decode:
            buf[i] = byte
        if mlen > 0:
            if expected == byte:
                mlen += 1; mptr += 1
            else:
                mlen = 0
        if i + 1 >= 6:
            h = hsh(30, byte | (c1 << 8) | (c2 << 16) | (c3 << 24) | (c4 << 32) | (c5 << 40)) & (MM - 1)
            if mlen == 0 and mtab[h] > 0:
                mptr = mtab[h]; mlen = 1
            mtab[h] = i + 1
        if (byte >= 65 and byte <= 90) or (byte >= 97 and byte <= 122):
            word = hsh(word, byte | 32)
        elif word != 0:
            word1 = word; word = 0
        if byte == 91 or byte == 123 or byte == 60 or byte == 124 or byte == 61 or byte == 39 or byte == 10:
            brk = byte
        st = G[st, byte]
        gh[i + 1] = out[st]
    if not decode:
        for _ in range(4):
            code[outp] = np.uint8(x1 >> np.uint64(24)); outp += 1
            x1 = (x1 << np.uint64(8)) & np.uint64(0xFFFFFFFF)
    return outp


def compress(data, cells, reg=None, regB=None, regC=None, grid=True, space=True, orbit_=True,
             extra=True, TB=24, lr=0.004, limit=255, verify=False):
    G, out, lens = build_automaton(cells)
    m = len(cells)
    z = np.zeros(m, dtype=np.int64)
    reg = z if reg is None else np.asarray(reg, dtype=np.int64)
    regB = reg if regB is None else np.asarray(regB, dtype=np.int64)
    use_c = regC is not None
    regC = z if regC is None else np.asarray(regC, dtype=np.int64)
    flags = np.array([grid, space and reg is not z, orbit_, use_c, extra], dtype=np.bool_)
    buf = np.frombuffer(data, dtype=np.uint8).astype(np.int64).copy()
    code = np.zeros(len(data) + 1024, dtype=np.uint8)
    nb = run(buf, len(data), False, code, G, out, lens, reg, regB, regC, flags, TB, lr, limit)
    if verify:
        back = np.zeros(len(data), dtype=np.int64)
        run(back, len(data), True, code[:nb], G, out, lens, reg, regB, regC, flags, TB, lr, limit)
        assert bytes(back.astype(np.uint8)) == data, "round trip failed"
    return nb
