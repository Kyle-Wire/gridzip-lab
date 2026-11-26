"""Points table of the V2 Pawnee rerun with the bridged fix (insar_pipeline PR #22, docs/BRIDGED_PLAN.md) for
the compressor, laid out as data/insar_m2/v2_pawnee_points.npz, with the bridged labels and the stated 1-sigma.

python -I export_fix.py <project dir> <run id: V2-fix> <out dir> --vdir <validation folder> [--points-from <npz>]

WSL, insar env, read-only on the project. Run it on BR1-A first; B and C are re-unwraps of A's stack on the same
grid, and with --points-from they are read at A's cells, so the three tables line up row for row.
"""
import argparse
import csv
import glob
import json
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("project", type=Path)
ap.add_argument("run")
ap.add_argument("out", type=Path)
ap.add_argument("--points-from", type=Path)
ap.add_argument("--vdir", type=Path, required=True)
a = ap.parse_args()
P, OUT, RUN = a.project, a.out, a.run
OUT.mkdir(parents=True, exist_ok=True)
tag = RUN.lower().replace("-", "_")                 # br1_a
rng = np.random.default_rng(20261125)
M = P / "mintpy"
WL = 0.05546576 * 1000                              # mm
CLASS_EDGES = (0.4, 0.6, 0.7, 0.8, 0.9, 1.01)       # analysis.significance.CLASS_EDGES


def rd(f, k):
    with h5py.File(f, "r") as h:
        return h[k][()], dict(h[k].attrs) | dict(h.attrs)


tsf = sorted(M.glob("timeseries_*demErr.h5"))[0]
ts, at = rd(tsf, "timeseries")
ts = ts * 1000.0
with h5py.File(tsf, "r") as h:
    dates = [d.decode() for d in h["date"][()]]
T, L, W = ts.shape
t = np.array([datetime.strptime(d, "%Y%m%d").toordinal() for d in dates], float)
ty = (t - t[0]) / 365.25
vel, _ = rd(M / "velocity.h5", "velocity")
vstd, _ = rd(M / "velocity.h5", "velocityStd")
tcoh, _ = rd(M / "temporalCoherence.h5", "temporalCoherence")
geo = {k: rd(M / "inputs/geometryRadar.h5", k)[0] for k in ("latitude", "longitude", "height", "incidenceAngle")}
cc, _ = rd(M / "maskConnComp.h5", "mask")
tri = None
for name in ("numTriNonzeroIntAmbiguity.h5",):
    if (M / name).exists():
        tri, _ = rd(M / name, "mask")
ry, rx = int(at["REF_Y"]), int(at["REF_X"])

with h5py.File(M / "inputs/ifgramStack.h5", "r") as h:
    drop = h["dropIfgram"][()]
    n0 = np.zeros((L, W), np.int32)
    for i in np.flatnonzero(drop):
        n0 += (h["connectComponent"][i] == 0)
    n_used = int(drop.sum())
    sub = {k: int(h.attrs[k]) for k in ("SUBSET_YMIN", "SUBSET_XMIN")}
    inverted = [k for k in h if k.startswith("unwrapPhase")]

# the bridged label (insar_pipeline analysis/bridged.py, mintpy/bridged.h5)
bfrac, breg = np.full((L, W), np.nan, np.float32), np.full((L, W), -1, np.float32)
battrs = {}
if (M / "bridged.h5").exists():
    with h5py.File(M / "bridged.h5", "r") as h:
        bfrac, breg, battrs = h["bridged_fraction"][()], h["region"][()], {k: str(v) for k, v in h.attrs.items()}
bridged = (breg > 0) | (breg == -2)

# the stated 1-sigma: class figure, and the bridged figure where a cell is bridged (max of the two)
um = json.loads((P / "deliverables/uncertainty_model.json").read_text())
try:                     # the bridged record survives QA here (bridged_sigma.json, 8a27a57)
    um["bridged"] = json.loads((P / "deliverables/bridged_sigma.json").read_text())
except (OSError, ValueError):
    pass
cls = np.digitize(np.nan_to_num(tcoh, nan=-1), CLASS_EDGES) - 1
stated = np.full((L, W), np.nan, np.float32)
for i, row in enumerate(um.get("by_class") or []):
    v = row.get("empirical_mm_yr")
    if v is not None:
        stated[cls == i] = float(v)
for k, v in ((um.get("bridged") or {}).get("stated_by_class") or {}).items():
    m = bridged & (cls == int(k))
    stated[m] = np.fmax(stated[m], float(v))

# full-resolution MiaplPy layers, averaged into each 3 x 11 look cell
ps_frac = ad_med = np.full((L, W), np.nan, np.float32)
try:
    with h5py.File(P / "miaplpy/maskPS.h5", "r") as h:
        ps = h["mask"][()]
        msub = {k: int(h.attrs[k]) for k in ("SUBSET_YMIN", "SUBSET_XMIN")}
    ad = np.fromfile(P / "miaplpy/inverted/amp_dipersion_index", "<f4").reshape(ps.shape)
    y0 = 3 * sub["SUBSET_YMIN"] - msub["SUBSET_YMIN"]
    x0 = 11 * sub["SUBSET_XMIN"] - msub["SUBSET_XMIN"]
    blk = lambda z: z[y0:y0 + 3 * L, x0:x0 + 11 * W].reshape(L, 3, W, 11)   # noqa: E731
    ps_frac = blk(ps.astype(np.float32)).mean(axis=(1, 3))
    ad_med = np.nanmedian(blk(np.where(ad > 0, ad, np.nan)).transpose(0, 2, 1, 3).reshape(L, W, 33), axis=2)
except (OSError, ValueError, KeyError) as exc:
    print("PS layers not read:", exc)

A = np.c_[np.ones(T), ty, np.sin(2 * np.pi * ty), np.cos(2 * np.pi * ty)]
coef = np.linalg.lstsq(A, np.nan_to_num(ts.reshape(T, -1)), rcond=None)[0]
ann = np.hypot(coef[2], coef[3]).reshape(L, W)
jumps = (np.abs(np.diff(ts, axis=0)) > WL / 4).sum(axis=0)

lat, lon = geo["latitude"], geo["longitude"]
k = 111320.0
dref = np.hypot((lat - lat[ry, rx]) * k, (lon - lon[ry, rx]) * k * np.cos(np.radians(lat[ry, rx])))

vdir = a.vdir
st = {}
for f in sorted(glob.glob(str(vdir / "gnss_cache/*/*.tenv3"))):
    with open(f) as fh:
        fh.readline()
        row = fh.readline().split()
    st[Path(f).stem] = (float(row[-3]), float(row[-2]))
ep = {}
if (vdir / "validation_epochs.csv").exists():
    with open(vdir / "validation_epochs.csv") as fh:
        for r in csv.DictReader(fh):
            ep.setdefault(r["station"], {})[r["date"]] = r

valid = np.isfinite(vel) & (vel != 0) & np.isfinite(tcoh)
near = np.zeros(valid.shape, bool)
names = sorted(st)
s_series, s_gnss, s_harness, s_n, s_br = [], [], [], [], []
for nm in names:
    la, lo = st[nm]
    d = np.hypot((lat - la) * k, (lon - lo) * k * np.cos(np.radians(la)))
    near |= (d <= 2000) & valid
    m = (d <= 200) & valid
    s_n.append(int(m.sum()))
    s_br.append(int((m & bridged).sum()))
    s_series.append(np.nanmean(ts[:, m], axis=1) if m.any() else np.full(T, np.nan))
    e = ep.get(nm, {})
    num = lambda x, c: float(e[x][c]) if x in e and e[x][c] not in ("", "nan") else np.nan   # noqa: E731
    s_gnss.append([num(x, "gnss_los_mm") for x in dates])
    s_harness.append([num(x, "insar_los_mm") for x in dates])

if a.points_from:
    z = np.load(a.points_from)
    iy, ix = z["row"].astype(int), z["col"].astype(int)
    keep = z["within_2km_of_station"] if "within_2km_of_station" in z.files else near[iy, ix]
else:
    iy, ix = np.nonzero(valid)
    kk = near[iy, ix]
    rest = np.flatnonzero(~kk)
    n_rest = max(0, 200_000 - int(kk.sum()))
    pick = np.sort(np.r_[np.flatnonzero(kk), rng.choice(rest, size=min(n_rest, rest.size), replace=False)])
    iy, ix, keep = iy[pick], ix[pick], kk[pick]

mini = np.arange(T) // 10
joins = t[np.r_[np.flatnonzero(np.diff(mini)) + 1]]
to_join = np.min(np.abs(t[:, None] - joins[None, :]), axis=1) if joins.size else np.full(T, np.nan)

f32 = lambda z: np.asarray(z, np.float32)   # noqa: E731
np.savez_compressed(
    OUT / f"{tag}_points.npz",
    row=iy.astype(np.int16), col=ix.astype(np.int16), lat=f32(lat[iy, ix]), lon=f32(lon[iy, ix]),
    rate_mm_yr=f32(vel[iy, ix] * 1000), formal_sigma_mm_yr=f32(vstd[iy, ix] * 1000),
    stated_sigma_mm_yr=f32(stated[iy, ix]), temporal_coherence=f32(tcoh[iy, ix]),
    bridged_fraction=f32(bfrac[iy, ix]), bridged_region=breg[iy, ix].astype(np.int16), bridged=bridged[iy, ix],
    ps_fraction=f32(ps_frac[iy, ix]), amp_dispersion_median=f32(ad_med[iy, ix]),
    dem_height_m=f32(geo["height"][iy, ix]), incidence_deg=f32(geo["incidenceAngle"][iy, ix]),
    dist_to_reference_m=f32(dref[iy, ix]), in_conncomp_mask=cc[iy, ix],
    n_ifgrams_outside_component=n0[iy, ix].astype(np.int16),
    n_triplets_nonzero_ambiguity=f32(tri[iy, ix]) if tri is not None else f32(np.full(iy.size, np.nan)),
    n_jump_dates=jumps[iy, ix].astype(np.int16), annual_amplitude_mm=f32(ann[iy, ix]),
    within_2km_of_station=keep,
    series_0p1mm=np.clip(np.round(np.nan_to_num(ts[:, iy, ix].T, nan=-3276.7) * 10), -32767, 32767).astype(np.int16),
    dates=np.array(dates), ministack=mini.astype(np.int8), days_to_ministack_join=f32(to_join),
    stations=np.array(names), station_lat=f32([st[n][0] for n in names]), station_lon=f32([st[n][1] for n in names]),
    station_points_200m=np.array(s_n), station_bridged_points_200m=np.array(s_br),
    station_radar_mean_200m_mm=f32(s_series), station_gnss_los_mm=f32(s_gnss), station_harness_radar_mm=f32(s_harness),
    reference_row_col=np.array([ry, rx]), reference_lat_lon=f32([lat[ry, rx], lon[ry, rx]]),
)
meta = {"run": RUN, "points": int(iy.size), "within_2km": int(np.sum(keep)), "valid_cells": int(valid.sum()),
        "bridged_valid_cells": int((valid & bridged).sum()), "bridged_regions": max(0, int(np.nanmax(breg))) if np.isfinite(breg).any() else 0,
        "bridged_layer_attrs": battrs, "dates": T, "first": dates[0], "last": dates[-1], "ifgrams_used": n_used,
        "inverted_datasets": inverted, "series_file": tsf.name, "grid": [L, W], "stations": dict(zip(names, s_n)),
        "stations_bridged_points_200m": dict(zip(names, s_br)), "points_from": str(a.points_from or ""),
        "uncertainty_model_bridged": {k: v for k, v in (um.get("bridged") or {}).items() if k != "by_class"}}
(OUT / f"{tag}_meta.json").write_text(json.dumps(meta, indent=1, default=str))
print(json.dumps({k: meta[k] for k in ("run", "points", "valid_cells", "bridged_valid_cells", "bridged_regions")}))
