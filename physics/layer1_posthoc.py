"""Post hoc (PROTOCOL-3, after A1): does removing each station's persistent offset narrow the thunderstorm range,
and which side are the Phoenix and Atlanta misses on? Not part of the protocol's tests."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import layer1 as P

ws, a3, hdr = sys.argv[1:4]; meta = P.stations_meta(hdr); l1 = Path(ws) / "docs/validation/layer1"
af = pd.concat([P.load(sorted((l1 / "f2/cal").glob("*.json")), "cal", meta), P.load(sorted((l1 / "f2/blind").glob("*.json")), "blind", meta)])
out = {}
for typ in ("thunderstorm", "non-thunderstorm"):
    d = P.features(af[af.type == typ]); r = np.exp(d.y)          # observed / predicted
    cal = d.group == "cal"; lo, hi = np.percentile(r[cal], [10, 90])
    # station offset from the station's OTHER storms, shrunk toward 0 by its count (k / (k + 3))
    k = d.groupby(["place", "station"]).y.transform("count") - 1
    off = d.stn_loo * k / (k + 3); r2 = np.exp(d.y - off)
    lo2, hi2 = np.percentile(r2[cal], [10, 90])
    res = {"range": [lo, hi], "range_after_station_offset": [lo2, hi2], "places": {}}
    for p, g in d.groupby("place"):
        m = d.place == p
        inside = (r[m] >= lo) & (r[m] <= hi); inside2 = (r2[m] >= lo2) & (r2[m] <= hi2)
        res["places"][p] = {"n": int(m.sum()), "group": g.group.iloc[0], "coverage": float(inside.mean()), "coverage_after_offset": float(inside2.mean()),
                            "miss_high": float((r[m] > hi).mean()), "miss_low": float((r[m] < lo).mean()),
                            "sd_y": float(g.y.std()), "sd_after_offset": float((g.y - off[m]).std()),
                            "worst_stations": g.groupby("station").y.agg(["mean", "count"]).sort_values("mean").round(3).to_dict("index")}
    out[typ] = res
    print(typ, "range", np.round([lo, hi], 3), "after offset", np.round([lo2, hi2], 3))
    for p, v in res["places"].items():
        print(f"  {p:14s} {v['group']:5s} n{v['n']:4d} cov {v['coverage']:.2f} -> {v['coverage_after_offset']:.2f}  high {v['miss_high']:.2f} low {v['miss_low']:.2f}  sd {v['sd_y']:.3f} -> {v['sd_after_offset']:.3f}")
(P.OUT / "layer1_posthoc.json").write_text(json.dumps(out, indent=1, default=float))

# station level: is the persistent offset explained by roughness, elevation, record length or a major-airport flag?
st = {}
for f in list((l1 / "f2/cal").glob("*.json")) + list((l1 / "f2/blind").glob("*.json")):
    j = json.loads(f.read_text())
    for s in j["stations"]: st[(j["region"], s["icao"])] = s
d = P.features(af[af.type != "thunderstorm"]); rows = []
for (p, s), g in d.groupby(["place", "station"]):
    if len(g) < 5: continue
    m = st[(p, s)]; nm = m["name"].upper()
    rows.append(dict(place=p, station=s, n=len(g), mean_y=g.y.mean(), lnz0=np.log(g.z0).mean(), elev=g.elev.iloc[0],
                     years=m.get("complete_years", 0), intl=float(("INTL" in nm) or ("INTERNATIONAL" in nm)), dist=m["distance_km"]))
S = pd.DataFrame(rows)
cors = {c: float(np.corrcoef(S.mean_y, S[c])[0, 1]) for c in ("lnz0", "elev", "years", "intl", "dist")}
print("station-level, ordinary+tropical, n stations", len(S), {k: round(v, 2) for k, v in cors.items()})
print(S.groupby("intl").mean_y.agg(["mean", "count", "std"]).round(3))
out["station_level"] = {"n_stations": len(S), "correlations": cors, "by_intl": S.groupby("intl").mean_y.agg(["mean", "count"]).to_dict(), "rows": S.to_dict("list")}
(P.OUT / "layer1_posthoc.json").write_text(json.dumps(out, indent=1, default=float))
