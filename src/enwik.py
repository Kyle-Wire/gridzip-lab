"""GridZip v2 on enwik8 (Hutter Prize data), vectorised so 100 MB is feasible.

Grid learned on the first SAMPLE bytes (map cost counted), applied to all of
enwik8, coordinates coded with the same adaptive order-1 + order-0 escape
model as exp1b. Code length is computed exactly (vectorised), not estimated.
"""
import gzip, math, os, sys, time
import numpy as np
from gridzip import learn_grid, apply_grid, map_bytes

def occ_before(keys):
    """for each position, how many times keys[i] occurred before i"""
    order = np.argsort(keys, kind="stable")
    k = keys[order]
    start = np.r_[0, np.nonzero(k[1:] != k[:-1])[0] + 1]
    grp = np.repeat(start, np.diff(np.r_[start, len(k)]))
    out = np.empty(len(keys), dtype=np.int64)
    out[order] = np.arange(len(k)) - grp
    return out

def ctx_code_bytes(t, m):
    t = t.astype(np.int64)
    prev = np.r_[-1, t[:-1]] + 1                    # context id (0 = start)
    c0 = occ_before(t)
    n0 = np.arange(len(t))
    p0 = (c0 + 0.5) / (n0 + 0.5 * m)
    pair = prev * (m + 1) + t
    c1 = occ_before(pair)
    n1 = occ_before(prev)
    first = (c1 == 0).astype(np.int64)              # new follower of this context
    # distinct followers seen before i in context prev[i]
    order = np.argsort(prev, kind="stable")
    f = first[order]
    cs = np.cumsum(f) - f
    pv = prev[order]
    start = np.r_[0, np.nonzero(pv[1:] != pv[:-1])[0] + 1]
    base = np.repeat(cs[start], np.diff(np.r_[start, len(pv)]))
    d = np.empty(len(t), dtype=np.int64); d[order] = cs - base
    w = d + 1
    p = (c1 + w * p0) / (n1 + w)
    return float(-np.log2(p).sum() / 8)

if __name__ == "__main__":
    path, cells, sample = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    data = open(path, "rb").read()
    if len(sys.argv) > 4:
        data = data[: int(sys.argv[4])]
    t0 = time.time()
    merges, _ = learn_grid(data[:sample], cells)
    print(f"learned {len(merges)} cells on {sample:,} B in {time.time()-t0:.0f}s", flush=True)
    toks = apply_grid(data, merges)
    print(f"applied: {len(data):,} B -> {len(toks):,} coordinates in {time.time()-t0:.0f}s", flush=True)
    stream = ctx_code_bytes(toks, 256 + len(merges))
    mp = map_bytes(merges)
    print(f"RESULT file={len(data):,} cells={len(merges)} map={mp:,.0f} stream={stream:,.0f} total={mp+stream:,.0f} "
          f"ratio={(mp+stream)/len(data):.4f} time={time.time()-t0:.0f}s", flush=True)
