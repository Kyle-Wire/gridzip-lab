"""Round 6 B: an adaptive region tree for the map (the "adaptive mesh refinement" idea).

The champion map clusters its cells twice, independently: 16 coarse and 128 fine k-means regions over the
cells' 16-dim embedding (src/write_map.py). Here instead one tree, grown by always splitting the leaf with the
largest frequency-weighted spread (the heaviest, most varied part of the text gets the finest resolution).
The cuts at 16, 64 and 256 leaves are nested, so lpaq1g can back off from fine to coarse (GFLAGS 262144).

Same cells, tokens and embedding as write_map.py (first 5 MB of the chain, 1,024 merges).
Writes ~/gz_exp/map_tree.bin = work/map.bin + int32 3 + three int32[m] label arrays.
Cost a decoder must be sent on top of the map: the 256-leaf label (8 bits per cell) and the split order
(one leaf index, 8 bits, per split): about 1,536 B for 1,280 cells.
Run in WSL from the repo root: ~/venvs/gridzip/bin/python work/exp/make_tree_map.py   (KMEANS=1: the flat k-means control, map_km.bin;
BITS=~/gz_exp/cost_train5_host.bin: weighted by bits spent, map_treebits.bin)"""
import json, os, struct, sys, time
import numpy as np

R = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(R, "src"))
from gridzip import learn_grid
from geometry import features, pca_coords

OUT = os.path.expanduser("~/gz_exp/map_tree.bin")
t0 = time.time()
data = open(os.path.join(R, "work", "enwik8.chain"), "rb").read()[:5_000_000]
M, toks = learn_grid(data, 1024)
table = [bytes([i]) for i in range(256)]
for a, b in M: table.append(table[a] + table[b])
m = len(table)
print(f"learn_grid: {len(M)} merges, {m} cells, {time.time() - t0:.0f} s", flush=True)

# the existing map must have exactly these cells (then the labels line up with lpaq1g's cell ids)
mp = open(os.path.join(R, "work", "map.bin"), "rb").read()
pos = 4; assert struct.unpack_from("<i", mp, 0)[0] == m
for c in table:
    l = mp[pos]; assert mp[pos + 1:pos + 1 + l] == c, "map.bin cells differ"; pos += 1 + l
pos += 8 * m                                    # coarse and fine regions
nt = struct.unpack_from("<i", mp, pos)[0]; pos += 4 + 4 * m
assert pos == len(mp), (pos, len(mp))          # map.bin = cells, regions, topics; nothing after

X = pca_coords(features(toks, m))
X = X / (np.linalg.norm(X, axis=1, keepdims=True) + 1e-9)
w = np.bincount(toks, minlength=m).astype(np.float64) + 1.0   # how often each cell occurs
if os.environ.get("BITS"):   # round 7b (the user's LISA hotspot idea): weight cells by the bits lpaq1g spends on them,
    # not by how often they occur, so the tree refines where prediction is poor. BITS = a DUMPCOST file of these 5 MB
    # (per byte: uint16 cost in 1/256 bit, uint16 cell); each parse token gets the cost of its bytes.
    cost = np.fromfile(os.environ["BITS"], "<u2").reshape(-1, 2)[:, 0].astype(np.float64) / 256
    assert len(cost) == len(data)
    lens = np.array([len(c) for c in table])[toks]; ends = np.cumsum(lens)
    cb = np.concatenate([[0.0], np.cumsum(cost)])
    tokbits = cb[ends] - cb[ends - lens]
    w = np.bincount(toks, weights=tokbits, minlength=m) + 1.0
    OUT = OUT.replace("map_tree", "map_treebits")


def spread(idx):
    c = (w[idx, None] * X[idx]).sum(0) / w[idx].sum()
    return float((w[idx] * ((X[idx] - c) ** 2).sum(1)).sum())


def split(idx):
    """Weighted 2-means, started from a split at the weighted median along the main direction."""
    c = (w[idx, None] * X[idx]).sum(0) / w[idx].sum()
    Y = (X[idx] - c) * np.sqrt(w[idx, None])
    v = np.linalg.svd(Y, full_matrices=False)[2][0]
    p = (X[idx] - c) @ v
    o = np.argsort(p, kind="stable"); cw = np.cumsum(w[idx][o])
    k = int(np.searchsorted(cw, cw[-1] / 2)); k = min(max(k, 0), len(idx) - 2)
    lab = np.zeros(len(idx), bool); lab[o[k + 1:]] = True
    for _ in range(10):
        c0 = (w[idx][~lab, None] * X[idx][~lab]).sum(0) / w[idx][~lab].sum()
        c1 = (w[idx][lab, None] * X[idx][lab]).sum(0) / w[idx][lab].sum()
        new = ((X[idx] - c1) ** 2).sum(1) < ((X[idx] - c0) ** 2).sum(1)
        if new.all() or not new.any() or (new == lab).all(): break
        lab = new
    return idx[~lab], idx[lab]


leaves = [np.arange(m)]
cuts, order = {}, []
while len(leaves) < 256:
    j = max((i for i in range(len(leaves)) if len(leaves[i]) > 1), key=lambda i: spread(leaves[i]))
    a, b = split(leaves[j]); order.append(j)
    leaves[j:j + 1] = [a, b]
    if len(leaves) in (16, 64, 256):
        lab = np.empty(m, np.int32)
        for i, L in enumerate(leaves): lab[L] = i
        cuts[len(leaves)] = lab
if os.environ.get("KMEANS"):   # control: flat, unnested k-means at the same three sizes (geometry.kmeans, as write_map.py)
    from geometry import kmeans
    P = pca_coords(features(toks, m))
    cuts = {16: np.asarray(kmeans(P, 16), np.int32), 64: np.asarray(kmeans(P, 64, seed=2), np.int32), 256: np.asarray(kmeans(P, 256, seed=3), np.int32)}
    OUT = OUT.replace("map_tree", "map_km")
with open(OUT, "wb") as f:
    f.write(mp); f.write(struct.pack("<i", 3))
    for k in (16, 64, 256): f.write(cuts[k].astype("<i4").tobytes())
sizes = {k: sorted(np.bincount(cuts[k]).tolist(), reverse=True) for k in cuts}
info = {"cells": m, "tree_cost_bytes": m + len(order), "leaf_sizes_16": sizes[16], "largest_256": sizes[256][:10],
        "singletons_256": sum(1 for s in sizes[256] if s == 1), "seconds": round(time.time() - t0)}
json.dump(info, open(os.path.join(R, "work", "exp", "km_map.json" if os.environ.get("KMEANS") else "treebits_map.json" if os.environ.get("BITS") else "tree_map.json"), "w"), indent=1)
print(json.dumps(info), flush=True)
print(f"wrote {OUT}: {os.path.getsize(OUT):,} B")
