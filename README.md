# gridzip-lab

Experiments on one idea: **store data as coordinates on a grid of learned
"functions"**, and see how small it gets once every byte of the map is paid for.

Each grid cell is a function of two earlier cells: `cell(k) = cell(a) + cell(b)`.
Cells are learned by repeatedly merging the most frequent adjacent pair. This is
byte-pair encoding, the same way GPT-style models split text into tokens. A file
becomes a list of coordinates, and the coordinates are then arithmetic-coded.

## Rules every result follows

- **The map counts.** Each learned cell costs two fixed-width coordinates,
  unless the experiment says the map is pre-shared.
- **Lossless, and checked.** v1 runs a real arithmetic encoder and decoder, and
  asserts that the round trip restores the file byte for byte. v2 and v3 report
  the exact adaptive code length (the sum of −log₂ p as the model learns
  online), which that same coder reaches to within a few bytes.
- **Fixed baselines.** Every table shows gzip -9, bzip2 -9 and xz -9e on the
  same bytes.

## Versions

| Version | How coordinates are coded |
|---|---|
| v1 (`src/gridzip.py`) | by how often each cell has appeared so far |
| v2 (`src/exp1b.py`, `src/enwik.py`) | conditioned on the previous cell, falling back to v1 |
| v3 (`src/v3.py`) | conditioned on the last K cells (K up to 4), blended from long contexts down to short ones (PPM style) |
| v3.1 (`src/v31.py`) | v3 plus a learned online mixer over context lengths |
| v4 (`src/v4.py`) | a bitwise context-mixing compressor (byte contexts, words, long matches) plus **grid-location contexts**: the longest grid cell ending at the last byte, found causally with an Aho-Corasick automaton. Real encoder and decoder |
| v5 (`src/v5.py`) | v4 plus a **spatial map**: cells are placed so that cells used in similar ways share a region, and nearby regions become context |

## Results

Sizes in bytes, with the map always included. gzip, bzip2 and xz are run on
the same bytes, except the 100 MB row, where the figures come from Matt
Mahoney's Large Text Compression Benchmark (mattmahoney.net/dc/text.html).

### enwik8 (Hutter Prize data)

| Compressor | 1 MB | 10 MB | 100 MB |
|---|---|---|---|
| gzip -9 | 355,791 | 3,688,458 | 36,445,248 |
| bzip2 -9 | 281,323 | 2,916,026 | 29,008,736 |
| xz -9e | 290,732 | 2,720,256 | 24,703,772 |
| 7zip PPMd (from the benchmark) | not run | not run | 21,197,559 |
| Grid v2 | not run | not run | 29,638,696 |
| Grid v3 (512 cells, 4 cells of context) | 266,346 | 2,420,606 | not run |
| v4, bytes only | 243,008 | 2,339,828 | 22,548,741 |
| v4 + grid (1,024 cells) | 241,215 | 2,320,451 | **22,284,286** |
| v5 + learned spatial map (16 regions) | 243,270 | **2,307,970** | **22,175,978** |
| v5 + random map (control) | 244,868 | 2,320,100 | not run |

The 1 MB column for v4 and v5 differs slightly: the v5 rows add the grid's map
cost (about 2.6 KB) and use an untuned mixer.

### What the grid is worth on top of a strong byte model

| Data size | Grid saving | Learned spatial map, extra saving | Random map |
|---|---|---|---|
| 1 MB | 0.6% | 0.2% | +0.4% (worse) |
| 10 MB | 0.8% | 0.5% | 0.0% |
| 100 MB | 1.2% | 0.5% | not run |

Both signals grow with data. The random-map control shows the gain comes from
the **arrangement** of the grid, not from simply having more contexts.

### Reorganise, then compress: BWT on grid coordinates (`src/bwt.py`)

A negative result. On 1 MB of enwik8, the Burrows-Wheeler transform followed
by move-to-front and context coding gives 287,415 B on bytes. On grid
coordinates it gives 291,201 B (256 cells), 298,442 B (1,024) and 312,596 B
(4,096). The BWT and the grid exploit the same structure, so stacking them adds
nothing, and the larger alphabet makes the reordered output harder to code.
The inverse transform is verified at every grid size.

### v6: the map at two zoom levels (`src/v6.py`)

At 1 MB it is slightly worse than v5: 244,232 B (16 coarse / 128 fine regions)
against 243,270 B, because the larger map costs more and the mixer has more to
learn. At 10 MB it helps: 2,303,247 B against v5's 2,307,970 B (−0.2%).

### Which geometry of the grid? (`src/geometry.py`, `src/em_map.py`, `src/emeans_map.py`)

All rows at 10 MB, 1,024 cells, 16 regions, map cost included.

| Map | Size |
|---|---|
| k-means on PCA coordinates (v5) | **2,307,927** |
| literal coordinate grid in PCA space (3 axes × 3 bins) | 2,312,498 |
| coordinate grid, 2 axes × 4 bins | 2,314,169 |
| coordinate grid, 4 axes × 2 bins | 2,314,813 |
| Laplacian eigenmaps (manifold learning), then k-means | 2,313,945 |
| Laplacian coordinate grid, 4 axes × 2 bins | 2,315,245 |
| random regions (control) | 2,320,100 |

Clusters that follow the data's own shape beat an axis-aligned coordinate grid
and the manifold embedding. At 1 MB, EM (Gaussian mixtures, soft membership)
and e-means (arXiv 2611.09501) match k-means to within 0.02%, even though
e-means reaches 3% lower clustering error. Only 68 of 1,280 cells have a real
second region. The clustering algorithm is not the bottleneck; the features
being clustered are.

### v7: symmetry (quotient) contexts (`src/v7.py`)

Contexts are computed on each byte's orbit under a symmetry group: letter case
is swapped, and the ten digits are permuted. At 10 MB on top of v6: case only
gives 2,301,768 B (−0.06%), and case plus digits gives **2,296,319 B (−0.30%)**,
the best 10 MB result for that engine. At 1 MB, case alone hurts (+0.2%).
On all of enwik8, v7 gives **22,048,955 B**.

### Engine upgrades: v8, v9, v10

| Engine | What changed | 1 MB | 10 MB |
|---|---|---|---|
| v7 | grid + two-zoom map + symmetry | 244,058 | 2,296,319 |
| v8 | 16M-slot tables, more byte contexts (order 5 and 8, sparse, markup), 3 mixers plus a final mixer, 2 SSE stages | 237,689 | 2,201,795 |
| v8 + code-length map | cells clustered by what byte follows, minimising cross-entropy | 238,524 | 2,202,564 |
| v9 | fast and slow estimates per context; decaying mixer learning rate | 226,338 | **2,108,842** |
| v10 | v9 + topic specialists: cells clustered by which passages they co-occur in; a running vote picks the topic, which selects a 4th mixer | 225,763 | running |

The v9 and v10 1 MB figures exclude the map cost (about 4.3 KB). The
code-length map keeps 92% of the exact cell's information in 64 regions (3.47
against 3.15 bits per byte, with 5.06 for no context), yet adds nothing in v8:
the exact cell is already a context there, so the map is redundant.

### Grafted into lpaq1, a mature compressor (`lpaq1g/`)

lpaq1 is Matt Mahoney's "lite paq" (GPL). With the graft switched off, it
reproduces lpaq1 byte for byte. Sizes include the map (4,319 B, plus 640 B for
topics). Every full-size run was decoded and compared with the original.

| lpaq1 -9 plus | 1 MB | 10 MB | enwik8 (100 MB) |
|---|---|---|---|
| nothing (original) | 226,908 | 2,122,364 | 19,755,948 |
| grid location | 229,718 | 2,112,803 | 19,620,916 |
| + spatial map | 227,981 | 2,097,388 | 19,509,152 |
| + symmetry | 227,021 | 2,087,990 | **19,456,911 (−1.5%)** |
| + topic specialists (16) | not run | **2,083,294** | **19,395,627 (−1.8%)** |

At 1 MB the map costs more than it saves; from 10 MB on it pays for itself.
Context on the enwik8 ladder (Mahoney's benchmark): drt|lpaq9m 17,964,751;
paq8hp12 16,230,028; paq8px 15,849,084; cmix 14,623,723.

### Curriculum: reorder the articles (`src/reorder.py`)

cmix ships an article order for enwik9. Here our own geometry picks the
order: each of enwik8's 12,346 articles becomes a TF-IDF vector over grid
cells, reduced with SVD, and articles are arranged so neighbours are similar.
The decoder needs only the permutation, at most log2(N!) bits = 18,750 B,
counted below.

| Order | lpaq1 | Net change |
|---|---|---|
| original | 19,755,948 | baseline |
| sorted by topic (64 clusters) | 19,543,338 | −1.0% |
| nearest-neighbour chain | 19,513,523 | **−1.1%** |
| random (control) | 19,957,970 | +1.0% |
| chain + lpaq1g (grid, map, symmetry) | 19,212,282 | −2.6% (19,235,351 with map and permutation) |
| chain + lpaq1g + topic specialists | 19,170,935 | **−2.8%** (**19,194,644** with map, topics and permutation) |

The two gains stack almost exactly: 1.5% from the contexts and 1.1% from
ordering. The best row was decoded and the article order inverted, and the
result is byte-identical to enwik8. The random control shows the original order already holds useful
locality, and placing related articles side by side adds another 1%. This is
experiment 3's lesson in a compressor: what an online learner sees next to
what changes how well it specialises.

### Can position alone hold data? (`src/positional.py`)

89,619 bytes of random data were hidden purely in the **order** of 109,056
distinct words, without changing any letter, and recovered exactly. It isn't
free: the words plus the hidden data cost 346,608 B, against 315,051 B for
coding the word set efficiently and storing the payload as plain bytes. The
reverse is the useful part: when order doesn't matter, sorting reclaims it.
The same words take 438,008 B shuffled and 225,432 B sorted with set-aware
coding.

### Early results (WikiText-2 and flood/heat docs)

Sizes in KB, map included.

| Data | gzip | bzip2 | xz | v1 | v2 |
|---|---|---|---|---|---|
| English Wikipedia, 500 KB | 161 | 126 | 138 | 160 | 138 |
| Flood/heat docs, 98 KB | 37 | 32 | 34 | 41 | 35 |
| Random bytes, 500 KB | 500 | 503 | 500 | 500 | not run |
| Sparse bits, 1% ones (true minimum 40) | 55 | 48 | 49 | 40 (no grid) | 40 (no grid) |

## Experiment 3: kid to PhD (`src/exp3.py`)

A small byte-level transformer (0.9 M parameters, 4,000 steps) is trained five
ways on the same budget: broad only, narrow only, broad then narrow, mixed, and
narrow then broad. Each is scored by bits per byte on held-out flood/heat docs
and on held-out Wikipedia, which is the size a model-driven compressor would
produce. Three seeds per condition. Results are in
`results/exp3_results.jsonl` once the runs finish.

## Experiment 4: transformer on grid tokens (`src/exp4.py`)

The same small transformer reads either raw bytes or grid coordinates, with
the same steps, batch and context length in tokens. Both are scored on a
held-out 500 KB slice of enwik8, in bits per test byte. Seed 0: bytes
**2.289**, grid (1,024 cells) **2.108**, so the grid model is 7.9% better at
equal compute. More seeds are pending.

## Run it

```bash
bash fetch_data.sh
```

```bash
cd src && python sweep_v3.py ../data/enwik8 10000000 0,512,1024,2048,4096 2,3,4 ../results/out.json
```

The specialist docs set comes from a private project and is not included.
