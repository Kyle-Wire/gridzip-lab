"""Points table of the method 1 control at Pawnee (OPERA DISP-S1, pin 3612eda) for the compressor.

python m1_table.py <mintpy dir> <grid dir (export_grid.py out)> <chi repo> <insar_pipeline repo> <m2 npz> <gnss cache> <out dir>
(WSL, insar env; read-only on the project and on V2.)

The series is CHI's exported series after the reference jump check (GROUND-POC-PLAN.md A11, A17, A19:
ground_readings.corrected), and the rate is re-fitted from it where the check corrected anything.
"""
import json
import sys
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np

MP, GD, CHI, IP, M2, GC, OUT = (Path(a) for a in sys.argv[1:8])
OUT.mkdir(parents=True, exist_ok=True)
sys.path[:0] = [str(CHI / "tools" / "gate4")]
import export_grid as eg          # noqa: E402
import ground_readings as gr      # noqa: E402

sys.path.insert(0, str(IP.parent))
from insar_pipeline.validate import ground_truth as gt  # noqa: E402

rng = np.random.default_rng(20261124)
WL = 0.05546576 * 1000

# ---- the exported grid (cells inside the box) and its series
g = dict(np.load(GD / "grid.npz"))
meta = json.loads((GD / "grid.json").read_text())
z = np.load(GD / "series.npz")
ser_raw, dates = z["mm"], [str(d) for d in z["dates"]]
t = gr.decimal_years(dates)
g2, ser, check = gr.corrected(g, (ser_raw, t))
ser_c = ser[0] if check["result"] == "corrected" else ser_raw
vel_c = g2["vel"]

# ---- coherence and the cell's grid position, built exactly as export_grid builds the box
with h5py.File(MP / "velocity.h5", "r") as f:
    shape, attrs = f["velocity"].shape, dict(f.attrs)
lat_all, lon_all = eg.cell_centres(attrs, shape)
s, w, n, e = meta["box"]
inside = (lat_all >= s) & (lat_all <= n) & (lon_all >= w) & (lon_all <= e)
rows, cols = np.nonzero(inside)
tc = np.full(rows.size, np.nan, np.float32)
if (MP / "temporalCoherence.h5").exists():
    with h5py.File(MP / "temporalCoherence.h5", "r") as f:
        tc = f["temporalCoherence"][()][inside].astype(np.float32)
assert rows.size == g["lat"].size

lat, lon = g["lat"], g["lon"]
valid = g2["valid"]
inc, az = np.radians(g["inc_deg"]), np.radians(g["az_deg"])
le, ln, lu = -np.sin(inc) * np.sin(az), np.sin(inc) * np.cos(az), np.cos(inc)
jumps = (np.abs(np.diff(ser_c, axis=1)) > WL / 4).sum(axis=1)

# ---- stations: positions as data/insar_m2, NGL IGS20 daily ENU from V2's own tenv3 copies, onto this LOS
m2 = np.load(M2)
names = [str(x) for x in m2["stations"]]
k = 111320.0
near = np.zeros(lat.size, bool)
pad = lambda d, k_: (datetime.strptime(d, "%Y%m%d") + __import__("datetime").timedelta(days=k_)).strftime("%Y%m%d")  # noqa: E731
gs = gt.load_gnss_series(names, window=(pad(dates[0], -10), pad(dates[-1], 10)), cache_dir=GC, refresh=False)
s_series, s_gnss, s_enu, s_n, s_los = [], [], [], [], []
for i, nm in enumerate(names):
    la, lo = float(m2["station_lat"][i]), float(m2["station_lon"][i])
    d = np.hypot((lat - la) * k, (lon - lo) * k * np.cos(np.radians(la)))
    near |= (d <= 2000) & valid
    m = (d <= 200) & valid
    s_n.append(int(m.sum()))
    s_series.append(np.nanmean(ser_c[m], axis=0) if m.any() else np.full(len(dates), np.nan))
    v = np.array([np.nanmedian(x[m]) for x in (le, ln, lu)]) if m.any() else np.full(3, np.nan)
    v = v / np.linalg.norm(v)
    s_los.append(v)
    if nm in gs:
        me = gt.match_epochs(gs[nm], dates)
        enu = me[["e_mm", "n_mm", "u_mm"]].to_numpy(float)
        los = enu @ v
        los = los - np.nanmean(los)
    else:
        enu, los = np.full((len(dates), 3), np.nan), np.full(len(dates), np.nan)
    s_enu.append(enu)
    s_gnss.append(los)

# ---- points: every valid cell within 2 km of a station, the rest a uniform random sample, 200,000 in all
iv = np.flatnonzero(valid)
keep = near[iv]
rest = iv[~keep]
pick = np.sort(np.r_[iv[keep], rng.choice(rest, size=min(max(0, 200_000 - int(keep.sum())), rest.size), replace=False)])

f32 = lambda a: np.asarray(a, np.float32)  # noqa: E731
np.savez_compressed(
    OUT / "m1_pawnee_points.npz",
    row=rows[pick].astype(np.int16), col=cols[pick].astype(np.int16), lat=f32(lat[pick]), lon=f32(lon[pick]),
    rate_mm_yr=f32(vel_c[pick]), rate_mm_yr_engine=f32(g["vel_mm"][pick]), sigma_mm_yr=f32(g["std_mm"][pick]),
    temporal_coherence=tc[pick], annual_amplitude_mm=f32(g["amp_mm"][pick]),
    incidence_deg=f32(g["inc_deg"][pick]), azimuth_deg=f32(g["az_deg"][pick]),
    n_jump_dates=jumps[pick].astype(np.int16), within_2km_of_station=near[pick],
    series_0p1mm=np.clip(np.round(np.nan_to_num(ser_c[pick], nan=-32768 / 10) * 10), -32768, 32767).astype(np.int16),
    dates=np.array(dates), date_in_m2=np.isin(np.array(dates), m2["dates"]),
    date_in_m2_window=(np.array(dates) >= str(m2["dates"][0])) & (np.array(dates) <= str(m2["dates"][-1])), stations=np.array(names), station_lat=m2["station_lat"], station_lon=m2["station_lon"],
    station_points_200m=np.array(s_n), station_radar_mean_200m_mm=f32(s_series),
    station_gnss_los_mm=f32(s_gnss), station_gnss_enu_mm=f32(s_enu), station_los_unit_enu=f32(s_los),
)
info = {"points": int(pick.size), "within_2km": int(keep.sum()), "cells_in_box": int(lat.size),
        "valid_cells": int(valid.sum()), "dates": len(dates), "first": dates[0], "last": dates[-1],
        "grid": list(shape), "stations": dict(zip(names, s_n)), "gnss_missing": dict(gs.missing),
        "jump_check": check, "export": {k_: meta[k_] for k_ in ("route", "component", "insar_commit", "box", "files_read")},
        "series_nan_code": -32768}
(OUT / "meta.json").write_text(json.dumps(info, indent=1, default=str))
print(json.dumps({k_: info[k_] for k_ in ("points", "within_2km", "valid_cells", "dates", "first", "last", "stations")}))
print("jump check:", check.get("result"), [(s_["decimal_year"], s_.get("fitted_mm")) for s_ in check.get("steps", [])],
      "bad dates", check.get("bad_dates"))
