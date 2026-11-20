"""Does 'where am I on the grid' add information on top of a strong byte model?
Same compressor, grid signal on vs off, several grid sizes, on N bytes of enwik8.
Map cost is added to every grid run. Real coded sizes (encoder output bytes)."""
import json, sys, time
from gridzip import learn_grid, map_bytes
import v4

path, n, sample, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
cells = [int(x) for x in sys.argv[5].split(",")]
data = open(path, "rb").read()[:n]
M, _ = learn_grid(data[:sample], max(cells))
res = []
def log(r):
    res.append(r); print(json.dumps(r), flush=True); json.dump(res, open(out, "w"), indent=1)
t = time.time(); nb, _ = v4.compress(data, [bytes([i]) for i in range(256)], True, False)
log(dict(n=n, cells=0, grid=False, size=nb, secs=round(time.time() - t)))
for c in cells:
    table = [bytes([i]) for i in range(256)]
    for a, b in M[:c]: table.append(table[a] + table[b])
    t = time.time(); nb, _ = v4.compress(data, table, True, True)
    log(dict(n=n, cells=c, grid=True, size=nb + map_bytes(M[:c]), map=map_bytes(M[:c]), secs=round(time.time() - t)))
