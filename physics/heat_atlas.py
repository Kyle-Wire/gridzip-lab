"""Study 4 of HeatSim's daily swing (HeatSim validation/prereg/2026-11-23_swing-4-atlas.md): which day and station
conditions explain the log ratio of the model's daily range to the station's, by code length.

python physics/heat_atlas.py --synthetic <days_dev.csv>         step 0: planted law and pure noise
python physics/heat_atlas.py <days_dev.csv> <label> <out.json>   step 2: search on r of <label> (col070 or era5)

Only the development half is read. Provenance: original; uses codelen.lawsearch and physics/runner.py.
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import runner  # noqa: E402
from codelen import lawsearch as L  # noqa: E402

G = L.Grid(0.01, -3.0, 3.0)
DAY = ["sm1", "sm2", "ef", "blh", "cloud", "wind", "dd", "sun"]          # registered, day level
STATION = ["lai_lo", "lai_hi", "cvh", "dz", "sdor", "lat"]               # registered, station level
PRESCREEN = 400


def load(path, label=None):
    d = pd.read_csv(path, dtype={"sid": str})
    assert (d["half"] == "dev").all(), "the confirmation half is not read by the search"
    need = DAY + STATION + ([f"{label}_r"] if label else [])
    return d.dropna(subset=need).reset_index(drop=True)


def library(d):
    """Standardised features (mean 0, sd 1 on these rows) and products of two day-level features."""
    z = {k: ((d[k] - d[k].mean()) / d[k].std()).values for k in DAY + STATION}
    lib = dict(z)
    for i, a in enumerate(DAY):
        for b in DAY[i + 1:]:
            lib[f"{a}*{b}"] = z[a] * z[b]
    return lib, {k: (float(d[k].mean()), float(d[k].std())) for k in DAY + STATION}


def folds(d):
    sids = sorted(d["sid"].unique(), key=lambda s: hashlib.sha256(s.encode()).hexdigest())
    grp = {s: i % 5 for i, s in enumerate(sids)}
    f = [(f"stations fold {k}", np.where(d["sid"].map(grp) != k)[0], np.where(d["sid"].map(grp) == k)[0]) for k in range(5)]
    east = d["lon"].values > -95.0
    f += [("leave out east of 95W", np.where(~east)[0], np.where(east)[0]), ("leave out west of 95W", np.where(east)[0], np.where(~east)[0])]
    return f


def synthetic(path):
    d = load(path)
    lib, _ = library(d)
    rng = np.random.default_rng(0)
    y = 0.15 - 0.6 * (d["sm1"].values - 0.25) + rng.normal(0, 0.12, len(d))
    top, n = L.search(lib, y, G, max_terms=3, keep=5, prescreen=PRESCREEN)
    slope = top[0]["fit_coef"].get("sm1", np.nan) / d["sm1"].std()            # back to per unit of soil water
    ok_a = top[0]["terms"] == ["sm1"] and abs(slope - (-0.6)) <= 0.12
    y0 = 0.1 + rng.normal(0, 0.12, len(d))
    top0, _ = L.search(lib, y0, G, max_terms=3, keep=3, prescreen=PRESCREEN)
    ok_b = top0[0]["terms"] == []
    out = {"n": len(d), "candidates": n, "planted": {"top": [t["terms"] for t in top], "slope_per_unit_sm1": float(slope), "ok": bool(ok_a)},
           "null": {"top": [t["terms"] for t in top0], "ok": bool(ok_b)}}
    print(json.dumps(out, indent=1))
    return out


def main(path, label, dest):
    d = load(path, label)
    y = d[f"{label}_r"].values.astype(float)
    lib, scale = library(d)
    res = runner.case(lib, y, G, folds(d), max_terms=3, keep=10, n_hold=5, prescreen=PRESCREEN)
    out = {k: v for k, v in res.items() if not k.startswith("_")}
    out.update(label=label, n_days=len(d), n_stations=int(d["sid"].nunique()), scale=scale,
               mean_r=float(y.mean()), sd_r=float(y.std()))
    Path(dest).write_text(json.dumps(out, indent=1, default=float))
    print(label, "n", len(d), "mean r", round(y.mean(), 3), "sd", round(y.std(), 3), "one", round(res["one_correction"]["bits_per_reading"], 4))
    for t in res["top"][:6]:
        print("  ", t["terms"], round(t["bits_per_reading"], 4), {k: round(v, 3) for k, v in t["fit_coef"].items()})
    for h in res["holdout"]:
        print("  HO", h["terms"], h["transfers"], [round(f["one_correction_bpr"] - f["formula_bpr"], 3) for f in h["folds"]])


if __name__ == "__main__":
    if sys.argv[1] == "--synthetic":
        synthetic(sys.argv[2])
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
