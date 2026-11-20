"""v7: symmetry (orbit) contexts on top of v6 (16 coarse / 128 fine regions)."""
import json, math, os, sys
import numpy as np
from gridzip import learn_grid, map_bytes
import v5, v7
path, n, cells, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
modes = [int(x) for x in sys.argv[5].split(",")]
data = open(path, "rb").read()[:n]
M, toks = learn_grid(data[:min(n, int(os.environ.get("SAMPLE", 2_000_000)))], cells)
table = [bytes([i]) for i in range(256)]
for a, b in M: table.append(table[a] + table[b])
m = 256 + len(M); res = []
A, B = 16, 128
for om in modes:
    ra, rb = v5.spatial_map(toks, m, A), v5.spatial_map(toks, m, B, seed=1)
    cost = map_bytes(M) + m * (math.log2(A) + math.log2(B)) / 8
    nb, _ = v7.compress(data, table, True, True, use_space=True, reg=ra, regB=rb, verify=n <= 200_000, omode=om)
    res.append(dict(n=n, cells=cells, orbit_mode=om, size=nb + cost))
    print(json.dumps(res[-1]), flush=True); json.dump(res, open(out, "w"), indent=1)
