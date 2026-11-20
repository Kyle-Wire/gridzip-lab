"""Sweep grid size x context order on a slice of enwik8."""
import bz2, gzip, json, lzma, sys, time
from gridzip import learn_grid, apply_grid, map_bytes
from v3 import code_bytes

path, n = sys.argv[1], int(sys.argv[2])
cells = [int(x) for x in sys.argv[3].split(",")]
orders = [int(x) for x in sys.argv[4].split(",")]
out = sys.argv[5]
data = open(path, "rb").read()[:n]
res = dict(bytes=n, gzip=len(gzip.compress(data, 9)), bz2=len(bz2.compress(data, 9)),
           xz=len(lzma.compress(data, preset=9 | lzma.PRESET_EXTREME)), runs=[])
print({k: v for k, v in res.items() if k != "runs"}, flush=True)
t0 = time.time()
merges_all, _ = learn_grid(data, max(cells))
print(f"learned {len(merges_all)} cells in {time.time()-t0:.0f}s", flush=True)
for c in cells:
    merges = merges_all[:c]
    toks = apply_grid(data, merges)
    for o in orders:
        size = map_bytes(merges) + code_bytes(toks, 256 + c, o)
        res["runs"].append(dict(cells=c, order=o, total=size))
        print(f"cells={c:5d} order={o} total={size:12,.0f}  ({size/res['bz2']:.3f} x bzip2)", flush=True)
        json.dump(res, open(out, "w"), indent=1)
