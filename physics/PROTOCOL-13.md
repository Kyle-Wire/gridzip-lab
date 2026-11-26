# Protocol 13: what our own InSAR processing (method 2) is missing on stable ground, V2 Pawnee (written 24 Nov 2026, before the search ran)

> Written and run as "Protocol 10" on branch claude/project-thread-dtcz3g (24-26 Nov 2026). Renumbered 13 on 26 Nov
> because main took protocols 10-12 in the meantime; nothing else changed.

**Why.** Kyle (24 Nov, 11:57 AM ET) wants the insar_pipeline's own processing ("method 2": MiaplPy phase linking,
snaphu, MintPy) to pass, since it is finer and more flexible than OPERA DISP-S1 ("method 1"), which passes. V2 at
Pawnee, Oklahoma, a stable-ground test site, is graded FAIL in insar_pipeline PR #20 (`results/D49/V2.yaml`):
noise_benchmark 7.75 mm/yr, pair share 8/15, rate RMSE 5.7, false-movement caveats; bias, series RMSE and
uncertainty calibration pass. This protocol asks which named mechanism makes the false rates, so the fix can be
written down and tested on a fresh site.

**Data.** `data/insar_m2/v2_pawnee_points.npz` (cd9ceab; README there). 200,000 MintPy cells, 127 dates
(2021-10-01 to 2025-11-03), ministacks of 10, six GNSS stations OK02 to OK07 with 77 to 81 cells within 200 m.

**Seen before this was written:** the README and meta.json; the per-point rate percentiles (median −1.5, 5–95 %
−12.3 to +6.7 mm/yr, sd 6.4), stated sigma percentiles (5/50/95 %: 0.44/0.72/0.99 mm/yr) and the per-date
ministack index. No relation between rate and any covariate, and no station series, has been looked at.

V2 is spent data (graded). Nothing here is evidence that method 2 passes. It names suspects; the fix is tested
only on a fresh site under insar_pipeline's own plan.

## Part A: the rate field (space)

- Target: y = rate_mm_yr per point. Stable ground, so the truth is about 0 up to the GNSS frame; whatever the
  formula takes out is a candidate error. Grid 0.05 mm/yr, robust coding.
- Sample: 20,000 points drawn at random (seed 10) from the valid points (finite rate, in the connected-component
  mask). The same search is repeated on a second, disjoint 20,000 (seed 11): a term counts only if it is named on both.
- Blocks and groups: 3 km tiles on lat/lon. Bootstrap by tile, B = 20; null search with tile-blocked shuffles, N = 9.
- Folds: west half → east half, east → west, south → north, north → south.
- Library (fixed now; each is a free slope unless said):
  - quality: `coherence` (temporal coherence), `ps_share` (ps_fraction), `amp_disp`, `outside_comp`
    (n_ifgrams_outside_component / 375), `closure_amb` (n_triplets_nonzero_ambiguity, log1p), `jumps` (n_jump_dates);
  - geometry: `height` (DEM height minus its mean), `incidence`, `ref_dist` (distance to the reference cell, km),
    `ramp` (east and north, 2 columns), `ramp2` (quadratic in east and north, 3 columns);
  - season: `annual_amp` (annual amplitude, mm).
- Grammar on (physlib v2), so shapes of these inputs are searched too.

## Part B: the series (time), the mechanisms a rate field can't name

For every sampled point, the line + annual fit is compared with richer fits, each adding one named piece. The
reading is the spread (robust sd) of the fitted rate across points, and its bits saved over the plain fit
(Gaussian, summed over points, with the parameter cost 0.5 log2(127) per extra coefficient per point).

- B1 `ministack steps`: one offset per ministack join (12 steps). Tests the MiaplPy sequential-estimator join bias.
- B2 `common mode`: subtract, per date, the median over all valid points of the 3 km tile (long-wavelength
  atmosphere or orbit left in every date).
- B3 `fading`: rate fitted on each ministack alone (10 dates), then the spread of those short rates around the
  whole-record rate. A phase-linking (fading-signal) bias shows as short rates that are biased in one direction,
  more strongly for low coherence and distributed scatterers.
- B4 `late start`: rate from 2022-10 onward only (drops the first year).

## Part C: the GNSS stations

Per station: radar series (200 m mean) minus the GNSS LOS series. The rate of that difference, and how much of it the
B1, B2 and B4 corrections take out. Six stations, so this is described, not gated.

## Readings, fixed now

- **A named.** A term is named if it is in both samples' winners, transfers on at least 3 of 4 folds and is picked in
  at least 60 % of bootstrap resamples. Its share of the rate variance is reported.
- **B named.** A B piece is named if it saves bits on average per point and cuts the robust sd of the rates by at
  least 20 %.
- **Unnamed.** The neighbour test on A's remainder says how much structure is left.
- **The headline is a ledger:** the rate sd at the start; the sd after each named piece (A first, then B); what is
  left. Expected (stated now, may be wrong): B1 or B3 (the phase-linking side) carry most of it, because OPERA, whose
  phase linking re-references differently, passes the same checks.

## Synthetic check first

On the real covariates and the real date list, with noise from a plain fit's residual:
- (a) no planted term in A: nothing named in at least 4 of 5 seeds;
- (b) planted `coherence` slope worth 10 % of the rate variance: named in at least 4 of 5;
- (c) planted ministack steps of 2 mm sd in B: B1 cuts the rate sd by at least 20 % in at least 4 of 5.

## Limits, stated now

- One site, one track, one season mix. Pawnee is rural Oklahoma (farmland), so distributed scatterers dominate;
  an urban site may behave differently.
- The full-resolution PS layer is mapped by look cell, so `ps_share` is a cell share, not a pixel flag.
- No land cover or soil moisture is in the table; `annual_amp` stands in for the seasonal ground.
- Nothing here changes insar_pipeline. Any fix goes through its own plan, written first, and a fresh site.

## A1 (24 Nov 2026, after parts B and C and a first timing run, before the synthetic check and part A)

**Seen by then:** parts B and C in full (results/physics/insar_m2.json), a 3,000-point timing run of part A
(no bootstrap; quality and geometry terms cut the rate sd only from 2.88 to 2.63 mm/yr, neighbour test p = 0.01),
and exploratory looks at the station rings, the 4 km tile map, a rate variogram and the series' lag structure.

**Change:** a 20,000-point search with bootstrap and null takes about 2.5 h per sample on this CPU. Part A now
uses 6,000 points per sample (seeds 10 and 11, disjoint), bootstrap B = 10 by tile, null N = 9. The readings
are unchanged.

## A2 (24 Nov 2026, same day, before any synthetic or part A result was read)

The 6,000-point run was killed for memory (the grammar's product terms on 15 GB) before it printed anything.
Part A and the synthetic check now use 3,000 points per sample. Nothing else changes.

## A3, part D: the method 1 control (24 Nov 2026 evening, written before any m1 value was computed here)

**Seen by then:** `data/insar_m1/README.md` (71b575a), including its untested first look: m1 rate 5–95 % −19.8 to
+29.7 mm/yr, median coherence 0.57, station series 11–17 mm RMS off GNSS. No m1 map, pair or station rate has
been computed. Confounds stated now: m1 had no weather correction (m2 had HRRR), and m1's rate is re-fitted after
CHI's reference-jump handling (A11, A17, A19).

**Readings, fixed now (same windows for both: dates in both tables):**
- **D1, the corner.** Using the same regions as before, NW (east < −4 km, north > 4 km) against the rest (east > −2 km, north < 4 km), take the
  median rate of each and difference them. "Ground or common physics" if m1's NW − rest ≤ −2 mm/yr; "our
  processing" if |m1 NW − rest| < 1 mm/yr; otherwise "unclear". Also the correlation of the two methods' 1 km
  tile maps.
- **D2, stations.** Rate of each station's 200 m radar series minus its GNSS LOS, relative to OK03, both methods,
  on the shared dates. The reading is the share of the 15 pairs within 2 mm/yr.
- **D3, error bars.** Stated sigma (median) against the rate spread between cells 0.2 to 0.4 km apart.
- **D4, lag structure.** RMS change at lags of 1, 10, 40 and 120 dates, with each date's median removed.
- **D5, compressor.** `physlib.insar.INSAR_RATE` on m1's rates, for the terms its table supports (coherence,
  incidence, ramp, ramp2, annual_amp, jumps), 3,000 points, B = 10 and N = 9, as in part A. Plus a cross-method
  term `m2_tile`: method 2's 1 km tile median at each m1 cell. If it is named, the two methods share the
  structure.

## A4, part E: the corner, with processing records (25 Nov 2026, written before these files were read here)

**Seen by then:** the export's README (8e77a82). Its first look said the HRRR correction's own rate runs −0.55 to
+0.37 mm/yr (5th to 95th percentile), and the median count of interferograms apart from the reference is 3 of 375.
Also corrected by then: the stated sigma reading (FINDINGS-13 correction).

**Readings, fixed now (valid points; NW and rest regions as in D1; target = the rate):**
- E1, which step makes the corner. NW − rest, taken of the rate from (i) the series before both corrections,
  (ii) the HRRR correction alone, (iii) the DEM-error correction alone (before − final − HRRR), and (iv) the final
  series. A step "makes the corner" if its own NW − rest is at least half of the final value, with the same sign.
- E2, unwrapping. New terms: `apart` (n_ifgrams_apart / 375) and `apart_late` (the share of the point's apart
  interferograms in the second half of the record). Plus the rate made only by the interferograms where the point is
  apart: the rate of the final series minus the rate re-inverted with those interferograms dropped is not
  available here, so the proxy is the corner's NW − rest inside the points with `apart` = 0. If the corner stays
  (at least half) among points never apart, unwrapping islands are not its cause.
- E3, compressor. `INSAR_RATE` plus `hrrr_rate`, `dem_rate`, `apart` and `apart_late` on 3,000 points, with
  B = 10 and N = 9, as in part A.
