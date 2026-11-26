# Pawnee points for the compressor, method 1 (OPERA DISP-S1): the control for data/insar_m2

`m1_pawnee_points.npz` (38.5 MB) is a points table from the insar_pipeline OPERA DISP-S1 route (method 1). It is the
control for `data/insar_m2/` (method 2, V2). Both use the same box, Sentinel-1 ascending track 34, and the same
stations at Pawnee, Oklahoma, a public stable-ground test site, so positions are included.

## The run
- insar_pipeline pin `3612eda` (`~/.cache/insar-pins/3612eda`, rebuilt from git on 24 Nov 2026; the same pin as
  CHI's DISP-S1 runs), `backend: disp`, direct fit. The config is `m1_pawnee_disp_asc34.yaml`, which uses the same
  settings as CHI's Kenner and Katy wide runs.
- **The GNSS tie was off** (`gnss_tie_radius_km: 0`, Kyle, 24 Nov 2026). The product was not anchored to OK02-OK07
  or to any other station. Its reference is the pixel the backend chose automatically.
- **Weather correction was off**, as at Kenner and Katy. Method 2 (V2) had HRRR applied.
- Frame F08887 (ascending track 34, 100% of the box), 128 granules. The track-136 frames were left out.
- 128 dates from 2021-10-01 to 2025-11-27. 126 of them fall in V2's window and are dates V2 also has
  (`date_in_m2`, `date_in_m2_window`). The last two (Nov 2025) are past V2's last date.
- Run time 5 h 29 min (17:55 to 23:24 UTC, 24 Nov 2026), CPU only, no GPU. The project folder is
  `/home/kylew/insar_projects/m1_pawnee_disp_asc34`, 1.4 GB.

## The export
- The grid was exported with CHI's `tools/gate4/export_grid.py --route disp --component los --series` (Climate
  Inspection main at 2e94b91). The box is S W N E = 36.300 -96.994 36.515 -96.734. `export_grid.json` gives the
  hashes of every MintPy file read.
- **Reference jump handling (CHI GROUND-POC-PLAN A11, A17, A19)** used `ground_readings.corrected`. Reading the
  box's median series, it found one lasting step at 2024.49: -28.6 mm date to date, +6.3 mm on the box's fit.
  That is under A17's 10 mm, so it was taken off every cell at the box's size. It also found 3 bad dates (2023.74, 2024.82, 2024.85), which were dropped
  (NaN in every cell). Each cell's rate was then re-fitted from the corrected series. The details are in
  `meta.json` under `jump_check`.
- 618,286 cells are in the box. 416,936 are valid: the export's mask (temporal coherence ≥ 0.4, ≥ 50% valid
  dates, MintPy's masks), and a finite rate after the re-fit.
- 200,000 points: every valid cell within 2 km of a station (61,324), and the rest a uniform random sample
  (seed 20261124). Built by `m1_table.py`.

## Points
| field | unit | what |
|---|---|---|
| row, col | cell | position on the 30 m DISP-S1 grid as MintPy wrote it (geocoded) |
| lat, lon | degrees | cell centre |
| rate_mm_yr | mm/yr, LOS | rate re-fitted from the corrected series (line plus annual term); positive is toward the satellite |
| rate_mm_yr_engine | mm/yr, LOS | the pipeline's own `velocity`, before the jump handling |
| sigma_mm_yr | mm/yr | the pipeline's `velocityStd`, the stated 1-sigma, kept as the engine gave it |
| temporal_coherence | 0-1 | MintPy temporal coherence |
| annual_amplitude_mm | mm | the pipeline's annual amplitude |
| incidence_deg, azimuth_deg | degrees | LOS geometry (MintPy convention) |
| n_jump_dates | count | consecutive-date steps larger than a quarter wavelength (13.9 mm) in the corrected series |
| within_2km_of_station | bool | kept because it lies within 2 km of a station |
| series_0p1mm | int16, 0.1 mm | corrected LOS displacement from the first date, points x 128 dates; -32768 means no value (the dropped dates) |

## Per date
`dates` (yyyymmdd), `date_in_m2` (the same date is in data/insar_m2), and `date_in_m2_window`.

## Stations (OK02 to OK07, NGL IGS20)
- `station_lat`, `station_lon` are copied from data/insar_m2.
- `station_points_200m`: valid cells within 200 m, from 83 to 144.
- `station_radar_mean_200m_mm`: the mean corrected series of those cells (stations x dates).
- `station_gnss_los_mm`: the station's daily NGL IGS20 east, north and up, read from the same tenv3 files the V2
  harness used (copied from V2's gnss_cache, not re-downloaded). Each value is prepared by insar_pipeline's
  `validate.ground_truth.load_gnss_series` (equipment steps and outliers removed), averaged over the days within
  3 days of each date (`match_epochs`), and projected onto the median LOS unit vector of the 200 m cells.
  - **It is referenced to its own mean over the dates.** The m2 table's `station_gnss_los_mm` is the V2 harness's
    value, on that harness's datum. Before comparing the two, take each series about its own mean.
- `station_gnss_enu_mm` (stations x dates x 3) and `station_los_unit_enu` are given so that any other projection
  or reference can be made.

## A first look (not graded, for the compressor)
- The site is stable ground, yet the rate spread is wide:
  - the 5th to 95th percentiles of the engine's rate are -19.8 to +29.7 mm/yr (median +2.2);
  - median temporal coherence is 0.57;
  - 98% of points have at least one date-to-date step over 13.9 mm.
- About each series' own mean, the 200 m mean series differs from the GNSS LOS by an RMS of 11 to 17 mm at the
  six stations.
- So method 1 here is noisy, as expected for DISP-S1 on rural ground without a weather correction. The reading
  is exploratory and was not tested against any protocol.
