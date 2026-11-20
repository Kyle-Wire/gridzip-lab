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
| v5 + learned spatial map (16 regions) | 243,270 | **2,307,970** | running |
| v5 + random map (control) | 244,868 | 2,320,100 | not run |

The 1 MB column for v4 and v5 differs slightly: the v5 rows add the grid's map
cost (about 2.6 KB) and use an untuned mixer.

### What the grid is worth on top of a strong byte model

| Data size | Grid saving | Learned spatial map, extra saving | Random map |
|---|---|---|---|
| 1 MB | 0.6% | 0.2% | +0.4% (worse) |
| 10 MB | 0.8% | 0.5% | 0.0% |
| 100 MB | 1.2% | running | not run |

Both signals grow with data. The random-map control shows the gain comes from
the **arrangement** of the grid, not from simply having more contexts.

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
held-out 500 KB slice of enwik8, in bits per test byte. Running.

## Run it

```bash
bash fetch_data.sh
```

```bash
cd src && python sweep_v3.py ../data/enwik8 10000000 0,512,1024,2048,4096 2,3,4 ../results/out.json
```

The specialist docs set comes from a private project and is not included.
