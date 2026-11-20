"""Write the grid map file lpaq1g reads: cells, coarse regions, fine regions."""
import struct, sys
import numpy as np
from gridzip import learn_grid, map_bytes
from geometry import features, pca_coords, kmeans

src, n_sample, cells, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
data = open(src, "rb").read()[:n_sample]
M, toks = learn_grid(data, cells)
table = [bytes([i]) for i in range(256)]
for a, b in M: table.append(table[a] + table[b])
m = len(table)
P = pca_coords(features(toks, m))
ra, rb = kmeans(P, 16), kmeans(P, 128, seed=1)
with open(out, "wb") as f:
    f.write(struct.pack("<i", m))
    for c in table: f.write(bytes([len(c)]) + c)
    f.write(np.asarray(ra, "<i4").tobytes()); f.write(np.asarray(rb, "<i4").tobytes())
# cost a decoder must be sent: grid (2 coords per cell) + 4 + 7 bits per cell
print(f"map cost {map_bytes(M) + m * 11 / 8:.0f}")

if len(sys.argv) > 5:   # append topics: cells clustered by co-occurrence in passages
    import v10
    NT = int(sys.argv[5])
    tl = v10.topic_map(toks, m, NT)
    with open(out, "ab") as f:
        f.write(struct.pack("<i", NT)); f.write(np.asarray(tl, "<i4").tobytes())
    import math
    print(f"topic map cost {m * math.log2(NT) / 8:.0f}")
