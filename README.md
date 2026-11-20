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

## Results so far

Sizes in KB, map included. Bold is the smallest in each row.

| Data | gzip | bzip2 | xz | v1 | v2 | v3 |
|---|---|---|---|---|---|---|
| English Wikipedia, 500 KB (WikiText-2) | 161 | **126** | 138 | 160 | 138 | not run |
| Flood/heat project docs, 98 KB | 37 | **32** | 34 | 41 | 35 | not run |
| Random bytes, 500 KB | 500 | 503 | 500 | 500 | not run | not run |
| Sparse bits, 1% ones (true minimum 40) | 55 | 48 | 49 | **40** (no grid) | **40** (no grid) | not run |
| enwik8, first 1 MB | 356 | 281 | 291 | not run | not run | **266** |
| enwik8, all 100 MB | 36,445 | **29,009** | 24,704 | not run | 29,639 | running |

The gzip, bzip2 and xz figures for all of enwik8 come from Matt Mahoney's
Large Text Compression Benchmark (mattmahoney.net/dc/text.html), where xz uses
a 1 GiB dictionary. The bold marks the best result run with a 9 MB dictionary
on the same data. Every other row was measured here.

**What we've learned:**
- Context about where you just were on the grid is worth about 14% (v1 to v2).
- Each grid has a best size. Past it, the map costs more than it saves, and the
  cells become too rare to predict well.
- Random data never shrinks (the pigeonhole principle). Sparse data is
  already at its minimum as raw bits, and a grid only makes it bigger.
- A grid learned on another field fits worse: a Wikipedia grid is 8% worse on
  the flood docs.

## Experiment 3: kid to PhD (`src/exp3.py`)

A small byte-level transformer (0.9 M parameters, 4,000 steps) is trained five
ways on the same budget: broad only, narrow only, broad then narrow, mixed, and
narrow then broad. Each is scored by bits per byte on held-out flood/heat docs
and on held-out Wikipedia, which is the size a model-driven compressor would
produce. Three seeds per condition. Results are in
`results/exp3_results.jsonl` once the runs finish.

## Run it

```bash
bash fetch_data.sh
```

```bash
cd src && python sweep_v3.py ../data/enwik8 10000000 0,512,1024,2048,4096 2,3,4 ../results/out.json
```

The specialist docs set comes from a private project and is not included.
