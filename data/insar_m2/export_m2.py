"""Points table of the V2 phase-linked run (Pawnee, Oklahoma; a public test site) for the compressor.

python -I export_m2.py <project dir> <out dir>      (WSL, insar env; read-only on the project)
"""
import csv
import glob
import json
import sys
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np

P, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261124)
M = P / "mintpy"
WL = 0.05546576 * 1000                     # mm


def rd(f, k):
    with h5py.File(f, "r") as h:
        return h[k][()], dict(h[k].attrs) | dict(h.attrs)


ts, at = rd(M / "timeseries_HRRR_demErr.h5", "timeseries")          # m, LOS, radar grid
ts = ts * 1000.0                                                     # mm
with h5py.File(M / "timeseries_HRRR_demErr.h5", "r") as h:
    dates = [d.decode() for d in h["date"][()]]
T, L, W = ts.shape
t = np.array([datetime.strptime(d, "%Y%m%d").toordinal() for d in dates], float)
ty = (t - t[0]) / 365.25
vel, _ = rd(M / "velocity.h5", "velocity")
vstd, _ = rd(M / "velocity.h5", "velocityStd")
tcoh, _ = rd(M / "temporalCoherence.h5", "temporalCoherence")
geo = {k: rd(M / "inputs/geometryRadar.h5", k)[0] for k in ("latitude", "longitude", "height", "incidenceAngle")}
cc, _ = rd(M / "maskConnComp.h5", "mask")
tri, _ = rd(M / "numTriNonzeroIntAmbiguity.h5", "mask")
ry, rx = int(at["REF_Y"]), int(at["REF_X"])

# connected components: share of interferograms where the cell is outside every component (label 0)
with h5py.File(M / "inputs/ifgramStack.h5", "r") as h:
    drop = h["dropIfgram"][()]
    n0 = np.zeros((L, W), np.int32)
    for i in np.flatnonzero(drop):
        n0 += (h["connectComponent"][i] == 0)
    n_used = int(drop.sum())
    sub = {k: int(h.attrs[k]) for k in ("SUBSET_YMIN", "SUBSET_XMIN")}

# full-resolution MiaplPy layers, averaged into each 3 x 11 look cell
with h5py.File(P / "miaplpy/maskPS.h5", "r") as h:
    ps = h["mask"][()]
    msub = {k: int(h.attrs[k]) for k in ("SUBSET_YMIN", "SUBSET_XMIN")}
ad = np.fromfile(P / "miaplpy/inverted/amp_dipersion_index", "<f4").reshape(ps.shape)
y0 = 3 * sub["SUBSET_YMIN"] - msub["SUBSET_YMIN"]
x0 = 11 * sub["SUBSET_XMIN"] - msub["SUBSET_XMIN"]
blk = lambda a: a[y0:y0 + 3 * L, x0:x0 + 11 * W].reshape(L, 3, W, 11)
ps_frac = blk(ps.astype(np.float32)).mean(axis=(1, 3))
ad_med = np.nanmedian(blk(np.where(ad > 0, ad, np.nan)).transpose(0, 2, 1, 3).reshape(L, W, 33), axis=2)

# per cell: annual amplitude (line + annual fit), jump dates (consecutive step > a quarter wavelength)
A = np.c_[np.ones(T), ty, np.sin(2 * np.pi * ty), np.cos(2 * np.pi * ty)]
coef = np.linalg.lstsq(A, ts.reshape(T, -1), rcond=None)[0]
ann = np.hypot(coef[2], coef[3]).reshape(L, W)
jumps = (np.abs(np.diff(ts, axis=0)) > WL / 4).sum(axis=0)

# distance to the reference (m, local flat earth)
lat, lon = geo["latitude"], geo["longitude"]
k = 111320.0
dref = np.hypot((lat - lat[ry, rx]) * k, (lon - lon[ry, rx]) * k * np.cos(np.radians(lat[ry, rx])))

# GNSS stations: position from their tenv3, LOS series from the V2 rerun's harness table
st = {}
for f in sorted(glob.glob(str(P / "validation/protocol_V2_V2/gnss_cache/*/*.tenv3"))):
    with open(f) as fh:
        fh.readline()
        row = fh.readline().split()
    st[Path(f).stem] = (float(row[-3]), float(row[-2]))
ep = {}
with open(P / "validation/protocol_V2_V2/validation_epochs.csv") as fh:
    for r in csv.DictReader(fh):
        ep.setdefault(r["station"], {})[r["date"]] = r

valid = np.isfinite(vel) & (vel != 0) & np.isfinite(tcoh)
iy, ix = np.nonzero(valid)
near = np.zeros(valid.shape, bool)
names = sorted(st)
s_series, s_gnss, s_harness, s_n = [], [], [], []
for nm in names:
    la, lo = st[nm]
    d = np.hypot((lat - la) * k, (lon - lo) * k * np.cos(np.radians(la)))
    near |= (d <= 2000) & valid
    m = (d <= 200) & valid
    s_n.append(int(m.sum()))
    s_series.append(np.nanmean(ts[:, m], axis=1) if m.any() else np.full(T, np.nan))
    e = ep.get(nm, {})
    s_gnss.append([float(e[x]["gnss_los_mm"]) if x in e and e[x]["gnss_los_mm"] not in ("", "nan") else np.nan for x in dates])
    s_harness.append([float(e[x]["insar_los_mm"]) if x in e and e[x]["insar_los_mm"] not in ("", "nan") else np.nan for x in dates])

keep = near[iy, ix]
rest = np.flatnonzero(~keep)
n_rest = max(0, 200_000 - int(keep.sum()))
pick = np.sort(np.r_[np.flatnonzero(keep), rng.choice(rest, size=min(n_rest, rest.size), replace=False)])
iy, ix = iy[pick], ix[pick]

# ministacks: MiaplPy sequential, 10 images each; per date its index and days to the nearest join
mini = np.arange(T) // 10
joins = t[np.r_[np.flatnonzero(np.diff(mini)) + 1]]
to_join = np.min(np.abs(t[:, None] - joins[None, :]), axis=1)

f32 = lambda a: np.asarray(a, np.float32)
np.savez_compressed(
    OUT / "v2_pawnee_points.npz",
    row=iy.astype(np.int16), col=ix.astype(np.int16), lat=f32(lat[iy, ix]), lon=f32(lon[iy, ix]),
    rate_mm_yr=f32(vel[iy, ix] * 1000), sigma_mm_yr=f32(vstd[iy, ix] * 1000), temporal_coherence=f32(tcoh[iy, ix]),
    ps_fraction=f32(ps_frac[iy, ix]), amp_dispersion_median=f32(ad_med[iy, ix]),
    dem_height_m=f32(geo["height"][iy, ix]), incidence_deg=f32(geo["incidenceAngle"][iy, ix]),
    dist_to_reference_m=f32(dref[iy, ix]), in_conncomp_mask=cc[iy, ix],
    n_ifgrams_outside_component=n0[iy, ix].astype(np.int16), n_triplets_nonzero_ambiguity=f32(tri[iy, ix]),
    n_jump_dates=jumps[iy, ix].astype(np.int16), annual_amplitude_mm=f32(ann[iy, ix]),
    within_2km_of_station=keep[pick],
    series_0p1mm=np.clip(np.round(ts[:, iy, ix].T * 10), -32767, 32767).astype(np.int16),
    dates=np.array(dates), ministack=mini.astype(np.int8), days_to_ministack_join=f32(to_join),
    stations=np.array(names), station_lat=f32([st[n][0] for n in names]), station_lon=f32([st[n][1] for n in names]),
    station_points_200m=np.array(s_n), station_radar_mean_200m_mm=f32(s_series),
    station_gnss_los_mm=f32(s_gnss), station_harness_radar_mm=f32(s_harness),
    reference_row_col=np.array([ry, rx]), reference_lat_lon=f32([lat[ry, rx], lon[ry, rx]]),
)
meta = {"points": int(iy.size), "within_2km": int(keep.sum()), "valid_cells": int(valid.sum()), "dates": T,
        "first": dates[0], "last": dates[-1], "ifgrams_used": n_used, "grid": [L, W],
        "stations": dict(zip(names, s_n)), "ps_fraction_vs_coherence_r":
        float(np.corrcoef(ps_frac[valid], tcoh[valid])[0, 1])}
(OUT / "meta.json").write_text(json.dumps(meta, indent=1))
print(json.dumps(meta))
