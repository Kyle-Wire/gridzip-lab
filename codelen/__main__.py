"""CLI: python -m codelen FILE.npz [--verify] [--res R --lo L --hi H --eps E] [--out bits.npz]

FILE.npz (the format agreed with sparc-transformers): one entry per real reading.
  y          float64 (n,)  readings, degC
  mu, sd     float64 (n,)  Gaussian forecast (sd in K)                     } give one of these
  samples    float32 (n, S) generative forecast                            }
  q_levels   (k,) and q_values (n, k): quantile forecast                  }
  bar_borders (B+1,) or (n, B+1) and bar_probs (n, B): TabPFN bar distribution }
  block      str (n,)      independent unit (week, city-week): bits are also reported per block
  hour       int (n,)      hour index: with block, fixes the coding order (block, then hour; stable)
  sensor     int (n,)      optional: sensor id, for the persistence baseline (keyed per block and sensor)
  ctx        int (n, k)    optional: row indices of context readings (coded earlier), -1 = none, for the context-mean baseline
  meta       JSON string   model, set, units, delta, eps, lo, hi, post_hoc. delta/eps/lo/hi set the grid
                           unless --res/--eps/--lo/--hi are given.
Prints one JSON object: grid, model and baseline bits (total, per reading, per block) and the savings.
--out writes the per-reading bits (row order of FILE) for the model and each baseline.
"""
import argparse, json
import numpy as np
from . import dists, scoring as sc

ap = argparse.ArgumentParser(prog="python -m codelen")
ap.add_argument("file"); ap.add_argument("--res", type=float); ap.add_argument("--lo", type=float)
ap.add_argument("--hi", type=float); ap.add_argument("--eps", type=float)
ap.add_argument("--verify", action="store_true", help="run the arithmetic coder both ways and assert the round trip")
ap.add_argument("--bandwidth", type=float, default=None, help="kernel width for samples (default Silverman)")
ap.add_argument("--out", help="write per-reading bits here (.npz)")
a = ap.parse_args()
d = np.load(a.file, allow_pickle=False)
meta = json.loads(str(d["meta"])) if "meta" in d else {}
res = a.res if a.res is not None else meta.get("delta"); eps = a.eps if a.eps is not None else meta.get("eps", dists.EPS)
lo = a.lo if a.lo is not None else meta.get("lo"); hi = a.hi if a.hi is not None else meta.get("hi")
if None in (res, lo, hi): ap.error("grid not given: pass --res/--lo/--hi or put delta/lo/hi in meta")
g = dists.Grid(res, lo, hi, eps)

n = len(d["y"]); order = np.arange(n); by_hour = False
block = d["block"].astype(str) if "block" in d else None
if block is not None and "hour" in d: order = np.lexsort((d["hour"], block)); by_hour = True   # block, then hour
inv = np.empty(n, np.int64); inv[order] = np.arange(n)
take = lambda k: d[k][order]
y = take("y")

out = {"file": a.file, "meta": meta, "grid": {"resolution": g.res, "lo": g.lo, "hi": g.hi, "eps": g.eps, "bins": g.K,
       "max_bits_per_reading": g.max_bits()}, "n": n, "coding_order": "block, hour" if by_hour else "rows"}
runs = {}
if "mu" in d: runs["model"] = sc.Given(g, take("mu"), take("sd"))
elif "samples" in d: runs["model"] = sc.GivenSamples(g, take("samples"), a.bandwidth)
elif "q_values" in d: runs["model"] = sc.GivenQuantiles(g, d["q_levels"], take("q_values"))
elif "bar_probs" in d: runs["model"] = sc.GivenBars(g, d["bar_borders"] if d["bar_borders"].ndim == 1 else take("bar_borders"), take("bar_probs"))
if "sensor" in d:
    key = take("sensor").astype(str) if block is None else np.char.add(np.char.add(block[order], "|"), take("sensor").astype(str))
    runs["baseline_persistence"] = sc.Persistence(g, key)
if "ctx" in d:
    c = d["ctx"][order]; c = np.where(c >= 0, inv[np.maximum(c, 0)], -1)          # row indices -> coding positions
    runs["baseline_context"] = sc.ContextMean(g, c)
bits = {}
for name, pred in runs.items():
    r = sc.score(y, pred, g, a.verify, per_reading=True); b = r.pop("bits"); bits[name] = b[inv]
    if block is not None:
        r["per_block"] = {k: {"n": int(m.sum()), "total_bits": float(b[m].sum()), "bits_per_reading": float(b[m].mean())}
                          for k in np.unique(block) for m in [block[order] == k]}
    out[name] = r
if "model" in out:
    for k in [k for k in out if k.startswith("baseline_")]:
        out[f"saving_vs_{k[9:]}_bits_per_reading"] = out[k]["bits_per_reading"] - out["model"]["bits_per_reading"]
        if block is not None:
            wins = [out[k]["per_block"][b]["total_bits"] > out["model"]["per_block"][b]["total_bits"] for b in out["model"]["per_block"]]
            out[f"blocks_won_vs_{k[9:]}"] = f"{sum(wins)}/{len(wins)}"
if "ref_bits" in d and "model" in bits:                                          # a fixture: compare with its reference
    dif = np.abs(bits["model"] - d["ref_bits"])
    out["vs_ref_bits"] = {"max_abs_diff": float(dif.max()), "mean_diff": float(np.mean(bits["model"] - d["ref_bits"])),
                          "ref_total_bits": float(d["ref_bits"].sum()), "worst_row": int(dif.argmax())}
if a.out: np.savez(a.out, **bits)
print(json.dumps(out, indent=1))
