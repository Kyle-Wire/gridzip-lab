"""V2 (Pawnee) per point: the HRRR correction applied, the series before it, the reference cell, and snaphu components.
Same row order as data/insar_m2/v2_pawnee_points.npz. Read-only on the project.

python -I export_m2_tropo.py <project dir> <v2_pawnee_points.npz> <out dir>
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

P, PTS, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
M = P / "mintpy"
z = np.load(PTS)
r, c = z["row"].astype(int), z["col"].astype(int)
q = lambda a: np.clip(np.round(a * 10), -32767, 32767).astype(np.int16)  # noqa: E731  (mm -> 0.1 mm)


def per_point(f, key="timeseries"):
    with h5py.File(f, "r") as h:
        d = h[key]
        out = np.empty((d.shape[0], r.size), np.float32)
        for k in range(d.shape[0]):
            a = d[k]
            out[k] = a[r, c]
        return out, dict(h.attrs), [x.decode() if isinstance(x, bytes) else str(x) for x in h["date"][()]] if "date" in h else None


raw, at_raw, dates = per_point(M / "timeseries.h5")
hrrr, at_h, _ = per_point(M / "timeseries_HRRR.h5")
final, at_f, _ = per_point(M / "timeseries_HRRR_demErr.h5")
assert dates == [str(d) for d in z["dates"]]
corr_mm = (raw - hrrr) * 1000.0                   # what the HRRR step took off, mm LOS (referenced as the series)
raw_mm = raw * 1000.0
ry, rx = int(at_raw["REF_Y"]), int(at_raw["REF_X"])
refs = {f: (int(a.get("REF_Y", -1)), int(a.get("REF_X", -1)), str(a.get("REF_DATE", "")))
        for f, a in (("timeseries.h5", at_raw), ("timeseries_HRRR.h5", at_h), ("timeseries_HRRR_demErr.h5", at_f))}

# the weather model's own delay (not referenced): at the reference cell, per date, and per point
with h5py.File(M / "inputs/HRRR.h5", "r") as h:
    dl = h["timeseries"]
    hd_ref = np.array([dl[k][ry, rx] for k in range(dl.shape[0])]) * 1000.0
    hd_pt = np.stack([dl[k][r, c] for k in range(dl.shape[0])]) * 1000.0
    hd_dates = [x.decode() if isinstance(x, bytes) else str(x) for x in h["date"][()]]

# interferograms: the reference cell's unwrapped phase (the dataset inverted), and the component labels
with h5py.File(M / "inputs/ifgramStack.h5", "r") as h:
    used = h["dropIfgram"][()].astype(bool)
    pairs = np.array([[a.decode() if isinstance(a, bytes) else str(a) for a in p] for p in h["date"][()]])
    ds = json.loads((M / "unwrap_error_applied.json").read_text())["dataset"]
    n = used.size
    ref_phase = np.array([h[ds][k, ry, rx] for k in range(n)], np.float32)
    ref_phase_raw = np.array([h["unwrapPhase"][k, ry, rx] for k in range(n)], np.float32)
    ref_label = np.array([h["connectComponent"][k, ry, rx] for k in range(n)], np.int16)
    diff = np.zeros((n, r.size), bool)            # the point and the reference are not in one component
    lab0 = np.zeros((n, r.size), bool)            # the point is outside every component
    nlab = np.zeros(n, np.int16)                  # components in the interferogram
    for k in range(n):
        L = h["connectComponent"][k]
        lp = L[r, c]
        diff[k] = (lp != ref_label[k]) | (lp == 0) | (ref_label[k] == 0)
        lab0[k] = lp == 0
        nlab[k] = int(np.unique(L[L > 0]).size)

# per date: of the used interferograms touching the date, how many have the point and reference apart
di = {d: i for i, d in enumerate(dates)}
per_date = np.zeros((len(dates), r.size), np.uint8)
touch = np.zeros(len(dates), np.int16)
for k in np.flatnonzero(used):
    for d in pairs[k]:
        if d in di:
            per_date[di[d]] += diff[k]
            touch[di[d]] += 1

np.savez_compressed(OUT / "v2_tropo.npz", hrrr_correction_0p1mm=q(corr_mm.T), hrrr_delay_point_0p1mm=q(hd_pt.T),
                    hrrr_delay_reference_mm=hd_ref.astype(np.float32), hrrr_dates=np.array(hd_dates), dates=np.array(dates))
np.savez_compressed(OUT / "v2_series_before_hrrr.npz", series_0p1mm=q(raw_mm.T), dates=np.array(dates))
np.savez_compressed(OUT / "v2_components.npz", ifgram_dates=pairs, ifgram_used=used,
                    apart_from_reference_bits=np.packbits(diff.T, axis=1), outside_component_bits=np.packbits(lab0.T, axis=1),
                    n_ifgrams_apart=diff[used].sum(0).astype(np.int16),
                    per_date_ifgrams_apart=per_date.T, per_date_ifgrams_used=touch, dates=np.array(dates),
                    reference_label=ref_label, components_per_ifgram=nlab,
                    reference_unwrap_phase_rad=ref_phase, reference_unwrap_phase_raw_rad=ref_phase_raw)
info = {"points": int(r.size), "dates": len(dates), "ifgrams": int(n), "ifgrams_used": int(used.sum()),
        "reference_yx": [ry, rx], "reference_by_file": refs, "reference_pin": json.loads((M / "reference_pin.json").read_text()),
        "reference_date": (M / "reference_date.txt").read_text().strip(), "inverted_dataset": ds,
        "reference_outside_component_ifgrams": int(((ref_label == 0) & used).sum()),
        "points_apart_median_of_used": float(np.median(diff[used].sum(0))),
        "hrrr_correction_rms_mm_median": float(np.median(np.sqrt(np.nanmean(corr_mm ** 2, axis=0))))}
(OUT / "v2_tropo_components_meta.json").write_text(json.dumps(info, indent=1))
print(json.dumps(info))
