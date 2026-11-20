"""Does the SPATIAL arrangement of the grid carry information?
Compare: no space contexts / space contexts with a learned map / with a RANDOM
map (same number of regions, same extra models, no meaningful layout).
Map cost = grid cells (2 coords each) + region id per cell."""
import json, math, sys, time
import numpy as np
from gridzip import learn_grid, apply_grid, map_bytes
import v5

path, n, cells, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
regions = [int(x) for x in sys.argv[5].split(",")]
data = open(path, "rb").read()[:n]
sample = data[:min(n, 2_000_000)]
M, toks = learn_grid(sample, cells)
table = [bytes([i]) for i in range(256)]
for a, b in M: table.append(table[a] + table[b])
gmap = map_bytes(M)
res = []
def log(r):
    res.append(r); print(json.dumps(r), flush=True); json.dump(res, open(out, "w"), indent=1)
nb, _ = v5.compress(data, table, True, True)
log(dict(n=n, cells=cells, space="none", size=nb + gmap))
for R in regions:
    for kind in ["learned", "random"]:
        reg = v5.spatial_map(toks, 256 + len(M), R) if kind == "learned" else \
              np.random.default_rng(1).integers(0, R, 256 + len(M))
        rcost = (256 + len(M)) * math.log2(R) / 8
        t = time.time()
        nb, _ = v5.compress(data, table, True, True, use_space=True, reg=reg)
        log(dict(n=n, cells=cells, space=kind, regions=R, size=nb + gmap + rcost, secs=round(time.time() - t)))
