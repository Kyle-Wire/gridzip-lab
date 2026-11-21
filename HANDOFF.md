# Handoff: continuing gridzip-lab on your own computer

Everything here was built in a cloud session with 4 cores and 15 GB of RAM.
cmix-class work needs more than that, so it continues on your machine. This
file gives the exact steps. A fresh Claude Code session can follow it from the
top.

## Where things stand

| What | Result | Notes |
|---|---|---|
| enwik8 champion | **18,608,516 B** (−5.8% vs lpaq1) | lpaq1g, chain order, `GFLAGS=95 MIX2=1 MIXSEL=1 IMMUNE=1`; map and permutation counted; decode verified |
| enwik9 | 161,433,444 B (−1.9% vs lpaq1) | `GFLAGS=15` only; full stack not yet run on enwik9 |
| cmix graft | compiles, **untested** | `cmix-graft/`; cmix needs about 32 GB of RAM, so cut-down "cmix-lite" was killed here at 8.3 GB |
| Evolution, speciation, experiment 3 | see `results/` | finished runs are committed |
| enwik8, local PC (2026-11-20) | **18,362,762 B** (−7.1% vs lpaq1) | round-3 settings + page, revision and contributor ids rebuilt from relations (`work/exp/headers.py`); see section 6 |
| **enwik8, full5 (2026-11-21)** | **17,868,078 B** (about −9.5% vs lpaq1; −494,684 B vs 18,362,762) | map-native stack from rounds 10–31 (cell contexts and memories, cell-region weight set, cell trigram and 4-gram) + the three header fields; map, tree, permutation and side stream counted; decode verified (section 6) |

The family tree of every experiment is `docs/gridzip-phylogeny.html`.

## 1. Set up

Linux or WSL recommended. Needs g++ (C++17), Python 3.12+ and about 10 GB of
free disk.

```bash
git clone https://github.com/Kyle-Wire/gridzip-lab && cd gridzip-lab
```

On Windows, clone with `git -c core.autocrlf=false clone ...`. Windows line
endings break `patch` and the shell scripts. In a WSL without `pip`,
`python3-venv` or `unzip` (and no sudo), make the venv with
`python3 -m venv --without-pip`, install into it with any other pip's
`--python` option, and unzip with `python3 -m zipfile -e file.zip .`.

```bash
python3 -m venv .venv && . .venv/bin/activate && pip install numpy==2.5.5 numba==0.68.1
```

`torch==2.15.0` is needed only for experiments 3 and 4.

```bash
bash fetch_data.sh
```

That downloads WikiText-2 and enwik8. For enwik9:

```bash
cd data && curl -O https://mattmahoney.net/dc/enwik9.zip && unzip enwik9.zip && cd ..
```

## 2. Build the compressors

```bash
g++ -O2 -o lpaq1g/lpaq1g lpaq1g/lpaq1g.cpp
```

```bash
git clone https://github.com/byronknoll/cmix cmix && git -C cmix checkout 194af9cd133b8b741d2a53afca13ed0ca453c276
```

```bash
cp cmix-graft/gridzip.h cmix/src/contexts/ && (cd cmix && patch -p1 < ../cmix-graft/cmix-v21-gridzip.patch)
```

```bash
cd cmix && make && cd ..
```

cmix's makefile uses `clang++-17`. With a different compiler, run
`make CXX=g++` or `make CXX=clang++`.

## 3. Reproduce the enwik8 champion

```bash
cd src && python reorder.py ../data/enwik8 ../work chain && python write_map.py ../work/enwik8.chain 5000000 1024 ../work/map.bin 16 && cd ..
```

```bash
IMMUNE=1 MIX2=1 MIXSEL=1 GRIDMAP=work/map.bin GFLAGS=95 lpaq1g/lpaq1g 9 work/enwik8.chain work/champ8
```

```bash
IMMUNE=1 MIX2=1 MIXSEL=1 GRIDMAP=work/map.bin GFLAGS=95 lpaq1g/lpaq1g d work/champ8 work/check && cmp work/check work/enwik8.chain
```

`reorder.py`'s faster article features (the ones in the repo now) give a
chain order about 0.16% worse than the run that set the record. Expect a size
within about 30 KB of 18,584,807, plus 4,959 B for the map and 18,750 B for
the permutation.

## 4. The cmix experiment (the reason for moving)

Question: do our contexts still help a far stronger compressor? Run the same
10 MB slice with the graft off and on. Same binary, same everything else:
only `GRIDMAP` differs.

```bash
head -c 10000000 data/enwik8 > work/e10 && (cd src && python write_map.py ../work/e10 2000000 1024 ../work/map10.bin 16)
```

```bash
cmix/cmix -n work/e10 work/e10.cmix
```

```bash
GRIDMAP=work/map10.bin cmix/cmix -n work/e10 work/e10.cmixg
```

```bash
GRIDMAP=work/map10.bin cmix/cmix -d work/e10.cmixg work/e10.check && cmp work/e10 work/e10.check
```

Expect about 2.7 hours per 10 MB pass (Ryzen 9 9900X). `CMIX_PPMD_MB`
(default 14000) and `CMIX_PAQ` (default 11) set the memory; use the same
values on **both** runs. Add the map's cost (4,959 B at 1,024 cells and 16
topics) before comparing.

Memory cmix reserves (VmPeak), measured; a 10 MB run grows into nearly all of it:

| `CMIX_PPMD_MB` / `CMIX_PAQ` | Reserved |
|---|---|
| 14000 / 11 (default) | 32.6 GiB |
| 4000 / 10 | 17.9 GiB |
| 2000 / 10 | 16.0 GiB |
| 4000 / 9 | 15.4 GiB |
| 2000 / 9 | 13.5 GiB |
| 1000 / 9 | 12.5 GiB |

These two settings alone can't reach the Hutter Prize limit of 10 GB: about
9 GiB sits in the rest of cmix (FXCM, the hash models, a 2 GB shared map). A
prize-legal run needs those trimmed too, as fx2-cmix does.

If the graft wins by 1% or more on 10 MB and on full enwik8, a prize attempt
is worth planning: fx2-cmix is the strongest open base under the 10 GB limit.

## 4b. Put the evolved genome in the champion (quick win)

Speciation's winner beat the hand-designed contexts by 1.0% on held-out 2 MB.
It has not been run at full scale. `GENOME` replaces the flag-based contexts;
keep `GFLAGS=95` so the grid, map, topic, structure and multi-scale state is
still tracked.

```bash
G=0:3,0:7,0:f,0:3f,80:0,1:0,103:0,1c:4,19:0,70:1,840:7,800:3f,24a0:3,580:0,200:3,1290:0,4060:1
```

```bash
GENOME=$G IMMUNE=1 MIX2=1 MIXSEL=1 GRIDMAP=work/map.bin GFLAGS=95 lpaq1g/lpaq1g 9 work/enwik8.chain work/champ8g
```

```bash
GENOME=$G IMMUNE=1 MIX2=1 MIXSEL=1 GRIDMAP=work/map.bin GFLAGS=95 lpaq1g/lpaq1g d work/champ8g work/check && cmp work/check work/enwik8.chain
```

Compare with step 3's size, then add the map and the permutation. The other
champions are in `results/speciation.json`, and the single-population winner
is in `results/evolve.json`.

## 5. Other open threads

- **Evolution:** `src/evolve.py` and `src/speciation.py` need `LPAQ_DIR`
  pointing at a folder holding `lpaq1gx` (build `lpaq1g.cpp` to that name),
  the slices `e2`/`h2`/`hab_*` and `map10t16.bin`. The commands are in each
  file's docstring. Raise `WORKERS` in `evolve.py` to use all your cores.
- **enwik9 full stack:** as in step 3, with `data/enwik9` (243,425 articles;
  permutation ≤ 500,557 B). It needs about 6 GB of RAM for reordering and
  about an hour to compress.
- **GPU:** allowed for research, such as tuning many variants or training
  neural parts, but not in the submitted program (Hutter rules: CPU only,
  10 GB RAM, 70,000/T hours).

## 6. Local PC, 2026-11-20 and overnight

All on WSL (Ryzen 9 9900X, 26 GB), every size counting map, permutation and
side stream, every archive decoded and compared. Quick-experiment harnesses are
`work/exp/run_round*.py` (results in `round*.jsonl`, logs `round*.log`);
`work/exp/run_full8.py` runs full enwik8.

- **Champion, 17,868,078 B (2026-11-21):** `CONFIG=full5 TRANSFORM=hdr run_full8.py` (binary lpaq1g34 or later):
  17,826,098 B archive + map with region tree 6,494 + permutation 18,750 + side stream 16,736. Encode and decode
  444 s each, peak 1.46 GB; enwik8 rebuilt byte for byte.
- **Previous champion, 18,362,762 B:** `CONFIG=candidate TRANSFORM=hdr run_full8.py`.
  Page ids are rebuilt from the article permutation, revision ids from their
  timestamps, and contributor ids from usernames seen before; side stream
  16,736 B. The transform helps every compressor tried (10 MB: lpaq1 −0.34%,
  xz −0.32%).
- **Map-native stack, level 9, 10 MB, against the round-4 host:**
  - cell contexts (the last three cells + a byte) and immune memories keyed
    on the last four and two cells: −0.64%;
  - with co-stimulation, next cell and the tree weight set ("full3"): −0.82%;
  - plus a mixer weight set chosen by the current cell × the previous cell's
    coarse region ("full4"): −1.27%;
  - held out at 20 MB and 85 MB into the chain: −1.24% and −1.36%.

  It is ready to run as `CONFIG=full4 TRANSFORM=hdr run_full8.py` (needs
  `~/gz_exp/lpaq1g25` and `map_tree.bin`).
- **Dead ends:**
  - in lpaq: a fluid (moving) map, article order by physics, particle-layout
    regions, a reversed chain, cell keys longer than four;
  - covered by cmix: a final APM keyed on cells (a byte-keyed APM of the same
    size does nearly as well).
- **cmix:** the graft as first written loses on 10 MB once the map is
  counted (+0.27%). New graft switches:
  - GZFLAGS 2048/4096: four- or three-cell contexts;
  - GZFLAGS 8192: a mixer by cell × previous region.

  They're built in `~/gz_exp/cmixc` and queued in `work/exp/run_cmix_cell4.py`
  after the 1 MB feature test (`run_cmix_feat.py`).

## Notes

- `data/docs_*.txt` (experiment 3's specialist set) came from a private
  project and is not redistributed. Experiment 3 needs a substitute set.
- `lpaq1g/` and `cmix-graft/` are GPL, as lpaq1 and cmix are.
