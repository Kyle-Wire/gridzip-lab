# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

gridzip-lab tests one idea: store data as coordinates on a grid of learned cells (byte-pair merges, `cell(k) = cell(a) + cell(b)`), arrange those cells on a spatial map so that cells used alike sit together, and use the map as context for a compressor. It started as standalone engines (`src/v1`–`v10`) and now lives mainly as contexts grafted into lpaq1 (`lpaq1g/`) and cmix (`cmix-graft/`), scored on enwik8/enwik9 (Hutter Prize data). `README.md` holds every result table; `HANDOFF.md` holds the current champion, the exact commands, and what to run next on a big machine; `docs/gridzip-phylogeny.html` is the family tree of experiments.

## Commands

Linux or WSL, g++ (C++17), Python 3.12+, numpy (numba for some engines, torch only for experiments 3 and 4). No test runner or linter is configured.

```bash
python3 -m venv .venv && . .venv/bin/activate && pip install numpy==2.5.5 numba==0.68.1
bash fetch_data.sh                                    # WikiText-2 and enwik8 into data/
g++ -O2 -o lpaq1g/lpaq1g lpaq1g/lpaq1g.cpp            # the grafted compressor
GRIDMAP=map.bin GFLAGS=15 lpaq1g/lpaq1g 9 in out      # compress; `lpaq1g d out restored` decodes
PYTHONHASHSEED=0 python champion/build_order.py data/enwik8 champion/out   # champion: article order + header transform (~15 s)
python champion/run.py lpaq1g/lpaq1g data/enwik8 champion/out              # champion: encode, decode, compare, count (~20 min)
python -m codelen.test_codelen                        # codelen self-tests (from the repo root)
python -m codelen preds.npz --verify                  # score a forecast by code length
cd src && python v8_test.py <file> <bytes> <out.json> # engine tests are scripts taking data paths, not pytest
```

Settings for lpaq1g are environment variables: `GFLAGS` (bit field of graft contexts; 0 reproduces lpaq1 byte for byte), `GRIDMAP`, `MIX2`, `MIXSEL`, `IMMUNE`, `GENOME`, `PRIMER`, `REFMODEL`. `lpaq1g/README.md` documents the bits. cmix needs ~32 GB RAM at default settings and ~2.7 h per 10 MB pass; ask Kyle before a run of that size.

## How the pieces fit

- `src/gridzip.py` learns the grid (`learn_grid`, `map_bytes`); `src/geometry.py`, `em_map.py`, `emeans_map.py`, `codelen_map.py` place cells into regions; `src/write_map.py` writes the binary map lpaq1g and cmix read; `champion/make_tree_map.py` adds the region tree.
- `src/reorder.py` and `champion/build_order.py` reorder enwik8's articles (a greedy chain by similarity); the permutation is sent and counted. `champion/headers.py` drops page, revision and contributor ids the decoder can rebuild (the side stream).
- `src/evolve.py` and `src/speciation.py` evolve context genomes for lpaq1g (`GENOME`); they need `LPAQ_DIR` and prepared slices (see each docstring).
- `src/exp3.py`, `src/exp4.py`: small byte-level transformers (curriculum order; bytes vs grid tokens).
- `codelen/`: a separate, CPU-only scorer that rates probabilistic forecasts by the bits needed to code real readings, with an integer arithmetic coder to verify the round trip. Its conventions follow the SPARC field-model scores in the separate `sparc-transformers` repo (`docs/design/field-model.md`, `src/sparc_tf/bits.py`), and its tests check agreement within 1e-5 bits.
- `results/` holds committed JSON/logs of finished runs; `work/` (git-ignored) holds local experiment harnesses (`work/exp/run_round*.py`, `run_full8.py`).

## Rules every result follows

- **The map counts.** Every size includes the map, region tree, permutation and side stream unless the experiment says the map is pre-shared.
- **Lossless, and checked.** A size is reported only after the archive is decoded and compared byte for byte (`cmp`). v2/v3 report exact adaptive code length instead, and say so.
- **Fixed baselines.** Tables show gzip -9, bzip2 -9 and xz -9e on the same bytes (100 MB rows quote Mahoney's benchmark).
- Report negative results (BWT on coordinates, dead ends in `HANDOFF.md`) as results. Compare on held-out slices before claiming a gain from a tuned setting.
- **Licences.** `lpaq1g/` and `cmix-graft/` are GPL derivatives. Everything else, including `codelen/`, is written fresh: never copy lpaq1, cmix or other GPL code into it. A Hutter Prize submission must run CPU-only within 10 GB RAM; GPUs are for research only.
- `data/docs_*.txt` (experiment 3's specialist set) is private and must not be committed.
