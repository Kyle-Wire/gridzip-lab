# Champion: enwik8 in 17,623,113 B

| part | bytes |
|---|---|
| lpaq1g archive (level 9, "full8" settings) | 17,581,133 |
| map and region tree | 6,494 |
| article permutation | 18,750 |
| header side stream | 16,736 |
| **total** | **17,623,113** |

This is about 10.8% under lpaq1 on the original file (19,755,948 B at level 9). It was decoded and rebuilt to enwik8 byte for byte on 2026-11-21. The phylogeny (`docs/gridzip-phylogeny.html`) has the history, and `HANDOFF.md` section 6 has the accounting.

## Reproduce (WSL or Linux, CPU only)
```
g++ -O2 -o lpaq1g/lpaq1g lpaq1g/lpaq1g.cpp
PYTHONHASHSEED=0 python champion/build_order.py data/enwik8 champion/out    # ~15 s: order + header transform
python champion/run.py lpaq1g/lpaq1g data/enwik8 champion/out             # ~35 min: encode, decode, check, total
```
`data/enwik8` comes from `fetch_data.sh`. The only Python dependency is numpy. `champion/out/` is git-ignored: it holds the reordered and stripped files, about 300 MB.

## Files
- `build_order.py`: the article order. A greedy chain over word-similarity vectors (hashed tf-idf). Then the header transform, with a rebuild check.
- `headers.py`: cuts page, revision and contributor ids that the decoder can work out, and rebuilds them.
- `run.py`: the compressor settings, then compress, decode, compare, rebuild enwik8, and count every byte.
- `map_tree.bin`: the grid map lpaq1g reads. It holds 1,280 cells (256 bytes + 1,024 merges), their coarse and fine regions, and a 3-level region tree. It is `map.bin` plus the tree from `make_tree_map.py`. Both maps were learned from the first 5 MB of the earlier cell-bag article chain (`src/write_map.py`, `src/reorder.py`). Their cost to send is counted in the total. `make_tree_map.py` is kept as it was run, with its original paths (`work/enwik8.chain`, `~/gz_exp`), as a record of how the tree was made.

Checked on 2026-11-21: a binary built from the current `lpaq1g/lpaq1g.cpp` gives a byte-identical archive to the one the champion run used (`lpaq1g41`), on a 1 MB test with these settings.

## What is in it
| step | total | change |
|---|---|---|
| map-native stack + header transform, cell-bag article order ("full5") | 17,868,078 | |
| + articles chained by word similarity | 17,683,150 | −184,928 |
| + contexts scoped to the current article, word novelty ("full7") | 17,648,738 | −34,412 |
| + article-hit weight set, subject name ("full8") | 17,623,113 | −25,625 |

The article order and the header transform do not depend on lpaq: on full enwik8 they take 1.4% to 4.7% off zlib, bzip2, zstd, xz and plain lpaq1 (`work/exp/transfer_test.py`). The map contexts are lpaq-only; in cmix they added nothing.
