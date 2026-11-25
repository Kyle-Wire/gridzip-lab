# V2 Pawnee points for the compressor (method 2, phase-linked)

`v2_pawnee_points.npz` (42.8 MB) is a points table from the insar_pipeline phase-linked route (method 2: MiaplPy
sequential phase linking, then snaphu and MintPy). It comes from experiment V2, Pawnee, Oklahoma: a public
stable-ground test site (insar_pipeline `configs/experiments/README.md`), so positions are included. The run was
made on 24 Nov 2026 with insar_pipeline at 68ba324 (same tree as main 54d52f2). It has 127 dates from 2021-10-01
to 2025-11-03, on Sentinel-1 ascending track 34, with HRRR weather and DEM-error corrections. The protocol V2
verdict is FAIL, as recorded (insar_pipeline PR #20, `results/D49/V2.yaml`). Exported by Claude, read-only on
the project, in about 7 s of CPU. `meta.json` holds the counts.

## Points
- There are 200,000 cells of the MintPy radar grid (674 x 720 cells; each cell is 3 azimuth x 11 range looks).
- Every valid cell within 2 km of a GNSS station is kept: 46,714 of them (`within_2km_of_station`).
- The rest are a uniform random sample (seed 20261124) of the 475,284 valid cells.
- One row per point:

| field | unit | what |
|---|---|---|
| row, col | cell | position on the MintPy radar grid |
| lat, lon | degrees | from `geometryRadar.h5` |
| rate_mm_yr | mm/yr, LOS | MintPy `velocity`; positive is toward the satellite |
| sigma_mm_yr | mm/yr | MintPy's formal `velocityStd` (white noise). **Not the pipeline's stated 1-sigma** (corrected 25 Nov 2026; see below) |
| temporal_coherence | 0-1 | MintPy temporal coherence |
| ps_fraction | 0-1 | share of the cell's 33 full-resolution pixels in MiaplPy's PS mask (the rest are distributed scatterers) |
| amp_dispersion_median | - | median MiaplPy amplitude dispersion of those pixels |
| dem_height_m | m | DEM height |
| incidence_deg | degrees | incidence angle |
| dist_to_reference_m | m | distance to the MintPy reference cell (flat earth) |
| in_conncomp_mask | bool | MintPy `maskConnComp` |
| n_ifgrams_outside_component | count of 375 | interferograms where the cell's snaphu component label is 0 |
| n_triplets_nonzero_ambiguity | count | MintPy `numTriNonzeroIntAmbiguity` (closure triplets with a non-zero integer ambiguity) |
| n_jump_dates | count | consecutive-date steps larger than a quarter wavelength (13.9 mm) in the series |
| annual_amplitude_mm | mm | amplitude of the annual term of a line + annual fit to the series |
| series_0p1mm | int16, 0.1 mm | the LOS displacement series, points x 127 dates (divide by 10 for mm) |

The full-resolution PS layers are mapped into each look cell through the two subsets. This holds:
`full_row = 3 * (row + 362) - 666` and `full_col = 11 * (col + 1928) - 19789`. I checked the mapping by
correlating MiaplPy's full-resolution temporal coherence, block-averaged, with MintPy's: r = 0.69 at zero offset
against 0.64 or less at one block off. (`ps_fraction` itself correlates only weakly with temporal coherence, r = 0.07.)

## The stated 1-sigma (added 25 Nov 2026)
Until 25 Nov 2026 this README called `sigma_mm_yr` "the stated 1-sigma". That was my error (Claude). It is MintPy's
formal velocityStd, a white-noise figure with a median of 0.72 mm/yr, which the pipeline itself labels a floor.

The pipeline states, and the V2 harness graded, the **empirical rate repeatability by temporal-coherence class**.
It comes from the run's own `deliverables/uncertainty_model.json`, measured on far-field cells: the robust sd of
(second-half rate - first-half rate) / 2. That is the full-record rate error under random-walk-like noise, so it
already counts temporally correlated noise.
- `v2_stated_sigma.npz` holds it per point, in the same row order as `v2_pawnee_points.npz`:
  - `stated_sigma_mm_yr`, median 4.75;
  - `cell_noise_mm_yr`, the max of that and the farmland cell spread, the figure the pipeline tests a cell against.
- `v2_uncertainty_model.json` is a copy of the model's table and basis, with its sha256.
- By class (temporal coherence → mm/yr): 0.4-0.6 → 6.59, 0.6-0.7 → 4.75, 0.7-0.8 → 3.20, 0.8-0.9 → 2.26,
  0.9-1.0 → 1.75.

## Per date
- `dates` (yyyymmdd).
- `ministack`: the MiaplPy sequential ministack index, 10 images each.
- `days_to_ministack_join`: days to the nearest date where one ministack hands over to the next.

## Stations (OK02 to OK07, NGL IGS20)
- `station_lat`, `station_lon`: from each station's tenv3.
- `station_points_200m`: cells within 200 m (77 to 81).
- `station_radar_mean_200m_mm`: the mean series of those cells (stations x dates).
- `station_gnss_los_mm`: the station's own daily GNSS, projected to LOS by the V2 harness (mean of the days
  within 3 days of each date), from `validation_epochs.csv`.
- `station_harness_radar_mm`: the harness's own radar reading at the station (median of reliable cells within
  100 m), for comparison.

NaN means no value. Not included: a stated sigma per date, and the per-point ministack index. The ministack index
depends on the date only, so it is in the per-date fields.
