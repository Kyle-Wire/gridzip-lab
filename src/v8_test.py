"""v8 tests. Map cost = grid + coarse(4 bits) + fine(7 bits) [+ regC bits] per cell."""
import json, math, os, sys, time
import numpy as np
from gridzip import learn_grid, map_bytes
from geometry import features, pca_coords, kmeans
import v8, v9, v10

path, n, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
variants = sys.argv[4].split(",")
data = open(path, "rb").read()[:n]
M, toks = learn_grid(data[:min(n, int(os.environ.get("SAMPLE", 2_000_000)))], 1024)
table = [bytes([i]) for i in range(256)]
for a, b in M: table.append(table[a] + table[b])
m = 256 + len(M)
P = pca_coords(features(toks, m))
ra, rb = kmeans(P, 16), kmeans(P, 128, seed=1)
base_cost = map_bytes(M) + m * 11 / 8
regC = None
if any("C" in v for v in variants):
    from codelen_map import codelen_map
    regC, info = codelen_map(data[:min(n, int(os.environ.get("SAMPLE", 2_000_000)))], table, int(os.environ.get("RC", 64)))
    print("code-length map:", info, flush=True)
res = []
for v in variants:
    kw = dict(TB=int(os.environ.get("TB", 24)), verify=n <= 200_000)
    if v.startswith("noextra"): kw["extra"] = False
    if "C" in v: kw["regC"] = regC
    cost = base_cost + (m * math.log2(int(os.environ.get("RC", 64))) / 8 if "C" in v else 0)
    t = time.time()
    eng = {"v8": v8, "v9": v9, "v10": v10}[os.environ.get("ENGINE", "v8")]
    if eng is v10 and int(os.environ.get("NT", 1)) > 1:
        kw["tlab"] = v10.topic_map(toks, m, int(os.environ["NT"]))
        kw["decay"] = float(os.environ.get("DECAY", 0.95))
        cost += m * math.log2(int(os.environ["NT"])) / 8
    nb = eng.compress(data, table, ra, rb, **kw)
    res.append(dict(n=n, variant=v, TB=kw["TB"], size=nb + cost, secs=round(time.time() - t)))
    print(json.dumps(res[-1]), flush=True); json.dump(res, open(out, "w"), indent=1)
