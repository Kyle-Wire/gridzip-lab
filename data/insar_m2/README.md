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

## The weather correction, the reference and the components (added 25 Nov 2026)
`export_m2_tropo.py` reads the V2 project folder without changing it. It ran in 13 s of CPU. Every array below is
in the same point order as `v2_pawnee_points.npz`, over its 127 dates. Values in 0.1 mm are int16 (divide by 10
for mm). Line of sight, positive toward the satellite. `v2_tropo_components_meta.json` holds the counts.

**`v2_tropo.npz`, the HRRR correction (a):**
- `hrrr_correction_0p1mm` (points x dates): what MintPy's HRRR step took off, `timeseries.h5` minus
  `timeseries_HRRR.h5`. It is on the same reference cell and date (20211130) as the series, so it is about 0 at
  the reference.
  - RMS over the dates: median 6.2 mm per point.
  - Its own rate: 5th to 95th percentile -0.55 to +0.37 mm/yr. So by itself the correction puts in or takes out
    well under 1 mm/yr.
- `hrrr_delay_point_0p1mm` (points x dates): the weather model's own LOS delay, unreferenced (`inputs/HRRR.h5`).
- `hrrr_delay_reference_mm` (dates): that delay at the reference cell.

**`v2_series_before_hrrr.npz` (b):** `series_0p1mm`, from `timeseries.h5`, before the weather correction and
before the DEM-error correction. The points table's series is `timeseries_HRRR_demErr.h5`, so the DEM-error
term is the difference between the two that is left after the HRRR correction is taken off.

**The reference cell (c):**
- It stayed the same. Radar y, x = 545, 535 (36.4951, -96.8055) in every time-series file, with reference date
  20211130 (`reference_pin.json`).
- The reference is inside a snaphu component in all 375 interferograms.
- `v2_components.npz` holds the reference cell's unwrapped phase per interferogram:
  - `reference_unwrap_phase_rad`: from the dataset that was inverted, `unwrapPhase_bridging_phaseClosure`;
  - `reference_unwrap_phase_raw_rad`: from `unwrapPhase`.
- There is no per-date reference phase, because MintPy's series is zero at the reference by construction. The
  per-interferogram values above are the per-date record, and `hrrr_delay_reference_mm` is the reference's
  weather.

**The snaphu components (d), `v2_components.npz`:**
- `ifgram_dates` (375 x 2) and `ifgram_used`: all 375 interferograms were used.
- `apart_from_reference_bits`: per point and interferogram, packed bits (`np.unpackbits(.., axis=1)[:, :375]`).
  A bit is 1 where the point and the reference are not in the same component, or either is outside every
  component.
- `outside_component_bits`: the same layout; 1 where the point's label is 0.
- `n_ifgrams_apart` (per point): median 3 of 375.
- `per_date_ifgrams_apart` (points x dates, uint8): of the interferograms that use that date, how many have the
  point apart from the reference. `per_date_ifgrams_used` (dates) counts those interferograms.
- `components_per_ifgram`: the number of components in each interferogram. `reference_label`: the reference's
  own label in each.
