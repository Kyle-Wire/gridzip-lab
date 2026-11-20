# lpaq1g

[lpaq1](https://mattmahoney.net/dc/) (Matt Mahoney, 2007, GPL), plus GridZip's
contexts. This file is a derivative of lpaq1 and is distributed under the GNU
GPL (version 3 or later, as lpaq1 allows). The rest of this repository is
separate work.

With `GFLAGS` unset or `0`, lpaq1g reproduces lpaq1's output byte for byte.

| `GFLAGS` bit | Adds |
|---|---|
| 1 | grid location: the longest grid cell ending at the last byte (Aho-Corasick), the cell before it, and the cell paired with order 3 |
| 2 | spatial map: coarse and fine regions of the last three cells |
| 4 | symmetry: order-3 and order-6 contexts on case-folded, digit-folded bytes |
| 8 | topic specialists: a decaying topic vote selects the mixer weights and joins two contexts |

The map file comes from `src/write_map.py`. Its cost is counted in every
result.

```bash
g++ -O2 -o lpaq1g lpaq1g.cpp
```

```bash
GRIDMAP=map.bin GFLAGS=15 ./lpaq1g 9 input output
```

```bash
GRIDMAP=map.bin GFLAGS=15 ./lpaq1g d output restored
```
