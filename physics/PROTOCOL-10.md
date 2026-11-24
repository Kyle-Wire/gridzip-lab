# Protocol 10: what our own InSAR processing (method 2) is missing on stable ground, V2 Pawnee (written 24 Nov 2026, before the search ran)

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
