"""GridZip v2: coordinates coded using the PREVIOUS coordinate as context.

"Where you just were on the grid" predicts where you go next. Code length is
the adaptive ideal (sum of -log2 p as the model learns online), which an
arithmetic coder reaches within a few bytes; the v1 run already proved the
real coder round-trips. Map cost is counted as in v1.
"""
import json, math, os, sys
from collections import defaultdict
import numpy as np
from gridzip import learn_grid, map_bytes

D = os.path.join(os.path.dirname(__file__), "..", "data")
SIZE = 500_000

def ctx_code_bytes(toks, m, lam=None):
    c0 = defaultdict(int); n0 = 0
    c1 = defaultdict(lambda: defaultdict(int)); n1 = defaultdict(int)
    bits, prev = 0.0, -1
    for s in toks:
        p0 = (c0[s] + 0.5) / (n0 + 0.5 * m)
        d = len(c1[prev])  # distinct followers seen: PPM-style escape weight
        w = d + 1
        p = (c1[prev][s] + w * p0) / (n1[prev] + w)
        bits -= math.log2(p)
        c0[s] += 1; n0 += 1; c1[prev][s] += 1; n1[prev] += 1
        prev = s
    return bits / 8

rng = np.random.default_rng(0)
sets = {
    "english (Wikipedia)": open(f"{D}/test.txt", "rb").read()[:SIZE],
    "flood/heat docs": open(f"{D}/docs_test.txt", "rb").read()[:SIZE],
    "sparse bits (1% ones)": np.packbits(rng.random(SIZE * 8) < 0.01).tobytes(),
}
out = {}
for name, data in sets.items():
    out[name] = {}
    for k in [0, 256, 1024, 2048, 4096]:
        merges, toks = learn_grid(data, k)
        size = map_bytes(merges) + ctx_code_bytes(toks.tolist(), 256 + len(merges))
        out[name][k] = size
        print(f"{name:24s} cells={len(merges):5d} context-coded total={size:9.0f}", flush=True)
json.dump(out, open(os.path.join(os.path.dirname(__file__), "exp1b_results.json"), "w"), indent=1)
