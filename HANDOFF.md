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

The family tree of every experiment is `docs/gridzip-phylogeny.html`.

## 1. Set up

Linux or WSL recommended. Needs g++ (C++17), Python 3.12+ and about 10 GB of
free disk.

```bash
git clone https://github.com/Kyle-Wire/gridzip-lab && cd gridzip-lab
```

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

Expect several hours per 10 MB run and about 32 GB of RAM. To fit the Hutter
Prize limit of 10 GB, set `CMIX_PPMD_MB` (default 14000) and `CMIX_PAQ`
(default 11) on **both** runs, for example `CMIX_PPMD_MB=4000 CMIX_PAQ=10`.
Add the map's cost (4,959 B at 1,024 cells and 16 topics) before comparing.

If the graft wins by 1% or more on 10 MB and on full enwik8, a prize attempt
is worth planning: fx2-cmix is the strongest open base under the 10 GB limit.

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

## Notes

- `data/docs_*.txt` (experiment 3's specialist set) came from a private
  project and is not redistributed. Experiment 3 needs a substitute set.
- `lpaq1g/` and `cmix-graft/` are GPL, as lpaq1 and cmix are.
