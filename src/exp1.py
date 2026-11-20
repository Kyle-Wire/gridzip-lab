import bz2, gzip, json, lzma, math, os, sys, time
import numpy as np
from collections import Counter
from gridzip import compress_report, learn_grid

D = os.path.join(os.path.dirname(__file__), "..", "data")
CELLS = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [0, 256, 1024, 4096]
SIZE = 500_000  # bytes per dataset, so every row is comparable

rng = np.random.default_rng(0)
sets = {
    "english (Wikipedia)": open(f"{D}/test.txt", "rb").read()[:SIZE],
    "flood/heat docs": open(f"{D}/docs_test.txt", "rb").read()[:SIZE],
    "random bytes": rng.integers(0, 256, SIZE, dtype=np.uint8).tobytes(),
    "sparse bits (1% ones)": np.packbits(rng.random(SIZE * 8) < 0.01).tobytes(),
}
source_entropy = {"random bytes": SIZE,
                  "sparse bits (1% ones)": SIZE * 8 * (-(0.01 * math.log2(0.01) + 0.99 * math.log2(0.99))) / 8}

def order0(b):
    c = Counter(b); n = len(b)
    return sum(-v * math.log2(v / n) for v in c.values()) / 8

t0 = time.time()
shared, _ = learn_grid(open(f"{D}/train.txt", "rb").read()[:2_000_000], max(CELLS))
print(f"shared grid ({len(shared)} cells) learned on 2 MB of other Wikipedia text in {time.time()-t0:.0f}s", flush=True)

rows = []
for name, data in sets.items():
    r = dict(dataset=name, original=len(data),
             gzip=len(gzip.compress(data, 9)), bz2=len(bz2.compress(data, 9)),
             xz=len(lzma.compress(data, preset=9 | lzma.PRESET_EXTREME)),
             order0_entropy=order0(data), source_entropy=source_entropy.get(name))
    for k in CELLS:
        g = compress_report(data, k)
        r[f"grid{k}_fixed"] = g["fixed_total"]; r[f"grid{k}_coded"] = g["coded_total"]
        r[f"grid{k}_map"] = g["map"]; r[f"grid{k}_cells"] = g["cells"]
        print(f"  {name:24s} cells={g['cells']:5d} map={g['map']:9.0f} fixed={g['fixed_total']:9.0f} coded={g['coded_total']:9.0f}", flush=True)
    for k in [c for c in CELLS if c]:
        g = compress_report(data, None, shared=shared[:k])
        r[f"shared{k}_coded"] = g["coded_total"]
    rows.append(r)
    print(json.dumps({k: (round(v) if isinstance(v, float) else v) for k, v in r.items()}), flush=True)

json.dump(rows, open(os.path.join(os.path.dirname(__file__), "exp1_results.json"), "w"), indent=1)
