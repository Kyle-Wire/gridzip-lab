"""GridZip: Kyle's "coordinate grid of functions" compressor.

Grid row 0 holds the 256 raw bytes. Every new cell is a function of two
earlier cells: cell(k) = concat(cell(a), cell(b)), learned by repeatedly
taking the most frequent adjacent pair (this is byte-pair encoding).
The file is then stored as a list of grid coordinates.

Size is counted honestly:
  map    = 2 coordinates per learned cell (fixed width)  [0 if pre-shared]
  stream = the coordinates, either fixed width or arithmetic coded with an
           adaptive frequency model (real encoder + decoder, round-trip checked)
"""
import math
import numpy as np


# ---------------------------------------------------------------- grid ---
def learn_grid(data: bytes, n_cells: int):
    t = np.frombuffer(data, dtype=np.uint8).astype(np.int64)
    merges = []
    nxt = 256
    for _ in range(n_cells):
        if len(t) < 2:
            break
        codes = t[:-1] * 1_000_000 + t[1:]
        vals, counts = np.unique(codes, return_counts=True)
        i = counts.argmax()
        if counts[i] < 2:
            break
        a, b = divmod(int(vals[i]), 1_000_000)
        t = _apply(t, a, b, nxt)
        merges.append((a, b))
        nxt += 1
    return merges, t


def _apply(t, a, b, new):
    pos = np.nonzero((t[:-1] == a) & (t[1:] == b))[0]
    if len(pos) == 0:
        return t
    if a == b:  # drop overlaps (aaa -> [aa]a), greedy left to right
        keep, last = [], -2
        for p in pos.tolist():
            if p != last + 1:
                keep.append(p)
                last = p
        pos = np.array(keep, dtype=np.int64)
    out = t.copy()
    out[pos] = new
    mask = np.ones(len(t), dtype=bool)
    mask[pos + 1] = False
    return out[mask]


def apply_grid(data: bytes, merges):
    t = np.frombuffer(data, dtype=np.uint8).astype(np.int64)
    for k, (a, b) in enumerate(merges):
        t = _apply(t, a, b, 256 + k)
    return t


def expand(tokens, merges):
    table = [bytes([i]) for i in range(256)]
    for a, b in merges:
        table.append(table[a] + table[b])
    return b"".join(table[x] for x in tokens.tolist())


# ------------------------------------------------- arithmetic coding ---
class Fenwick:
    def __init__(self, n, init=1):
        self.n = n
        self.tree = [0] * (n + 1)
        for i in range(n):
            self.add(i, init)
        self.total = n * init

    def add(self, i, d):
        if i < self.n and i >= 0:
            pass
        j = i + 1
        while j <= self.n:
            self.tree[j] += d
            j += j & -j

    def prefix(self, i):  # sum of [0, i)
        s = 0
        while i > 0:
            s += self.tree[i]
            i -= i & -i
        return s

    def find(self, target):  # largest i with prefix(i) <= target
        pos, rem = 0, target
        step = 1 << self.n.bit_length()
        while step:
            nx = pos + step
            if nx <= self.n and self.tree[nx] <= rem:
                pos = nx
                rem -= self.tree[nx]
            step >>= 1
        return pos


BITS = 32
FULL = (1 << BITS) - 1
HALF = 1 << (BITS - 1)
QTR = 1 << (BITS - 2)
INC = 32


def ac_encode(tokens, m):
    f = Fenwick(m)
    low, high, pending, out = 0, FULL, 0, []

    def emit(bit):
        nonlocal pending
        out.append(bit)
        out.extend([1 - bit] * pending)
        pending = 0

    for s in tokens:
        rng = high - low + 1
        cl, ch = f.prefix(s), f.prefix(s + 1)
        high = low + rng * ch // f.total - 1
        low = low + rng * cl // f.total
        while True:
            if high < HALF:
                emit(0)
            elif low >= HALF:
                emit(1); low -= HALF; high -= HALF
            elif low >= QTR and high < HALF + QTR:
                pending += 1; low -= QTR; high -= QTR
            else:
                break
            low, high = low * 2, high * 2 + 1
        f.add(s, INC); f.total += INC
    pending += 1
    emit(0 if low < QTR else 1)
    return out


def ac_decode(bits, n, m):
    f = Fenwick(m)
    low, high, val, idx = 0, FULL, 0, 0
    for _ in range(BITS):
        val = (val << 1) | (bits[idx] if idx < len(bits) else 0); idx += 1
    res = []
    for _ in range(n):
        rng = high - low + 1
        target = ((val - low + 1) * f.total - 1) // rng
        s = f.find(target)
        cl, ch = f.prefix(s), f.prefix(s + 1)
        high = low + rng * ch // f.total - 1
        low = low + rng * cl // f.total
        while True:
            if high < HALF:
                pass
            elif low >= HALF:
                low -= HALF; high -= HALF; val -= HALF
            elif low >= QTR and high < HALF + QTR:
                low -= QTR; high -= QTR; val -= QTR
            else:
                break
            low, high = low * 2, high * 2 + 1
            val = (val << 1) | (bits[idx] if idx < len(bits) else 0); idx += 1
        res.append(s)
        f.add(s, INC); f.total += INC
    return res


# ------------------------------------------------------------ sizing ---
def map_bytes(merges):
    bits = sum(2 * math.ceil(math.log2(256 + k)) for k in range(len(merges)))
    return bits / 8


def compress_report(data, n_cells, shared=None, verify=True):
    if shared is None:
        merges, toks = learn_grid(data, n_cells)
        mapcost = map_bytes(merges)
    else:
        merges = shared
        toks = apply_grid(data, merges)
        mapcost = 0.0
    m = 256 + len(merges)
    fixed = len(toks) * math.ceil(math.log2(m)) / 8
    tl = toks.tolist()
    bits = ac_encode(tl, m)
    if verify:
        back = ac_decode(bits, len(tl), m)
        assert back == tl, "arithmetic round trip failed"
        assert expand(np.array(back), merges) == data, "grid round trip failed"
    coded = len(bits) / 8
    return dict(cells=len(merges), tokens=len(toks), map=mapcost,
                fixed_total=mapcost + fixed, coded_total=mapcost + coded)
