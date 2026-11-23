"""PROTOCOL-4 part B, W2 "after": WindSim Layer 1 errors after F5 (station z0 from its own gusts).

python physics/layer1_f5.py <WindSim checkout with docs/validation/layer1/f5> <isd-history.csv>

F5 was adopted for ordinary storms only (WindSim docs/validation.md, "F5 result"), so the "after" set is the
F5 run's ordinary-storm pairs from variant A (`sources/`) and its thunderstorm pairs from `off/` (unchanged by
F5), calibration and blind places together. Same features, library, folds and test as physics/layer1.py
(PROTOCOL-3); only the data differ. Station elevations come from NOAA's isd-history.csv (ELEV(M)) instead of
the per-station header csvs.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import layer1 as P3

OUT = P3.OUT


def meta_from_history(path):
    h = pd.read_csv(path, dtype=str)
    h["id"] = h["USAF"] + h["WBAN"]
    return {r.id: (float(r.LAT), float(r.LON), float(r["ELEV(M)"])) for _, r in h.dropna(subset=["LAT", "LON", "ELEV(M)"]).iterrows()}


def main(ws, hist):
    meta = meta_from_history(hist); f5 = Path(ws) / "docs" / "validation" / "layer1" / "f5"
    def grab(mode):
        return pd.concat([P3.load(sorted((f5 / g / mode).glob("*.json")), g, meta) for g in ("cal", "blind")])
    src, off = grab("sources"), grab("off")
    out = {"sources": {"windsim_f5_results": "33af0a7 + e41b6c5", "f5_adopted_types": ["non-thunderstorm"]}}
    th = P3.features(off[off.type == "thunderstorm"])
    q1 = P3.run_case(th, "Q1 thunderstorm, after F5 (unchanged by F5; 5 new blind places)")
    q1["data"] = th[["place", "group", "station", "storm", "y", "dirn", "hl", "month", "lnreg", "dz_elev", "storm_loo", "stn_loo"]].to_dict("list")
    out["Q1_thunderstorm"] = P3.strip(q1)
    nt = P3.features(src[src.type == "non-thunderstorm"])
    q2 = P3.run_case(nt, "Q2 ordinary storms, after F5 (own-gust z0 at the source stations)")
    q2["data"] = nt[["place", "group", "station", "storm", "y", "dirn", "hl", "month", "dz_elev", "elev", "z0"]].to_dict("list")
    out["Q2_ordinary"] = P3.strip(q2)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "layer1_f5.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main(*sys.argv[1:3])
