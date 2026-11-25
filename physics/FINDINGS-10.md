# Findings 10: method 2 (V2 Pawnee) on stable ground (protocol 10, A1, A2; 24 Nov 2026)

> **Correction (25 Nov 2026).** The "error bars about 3× too small" reading below is wrong. The table's
> `sigma_mm_yr` was MintPy's formal white-noise `velocityStd` (median 0.72), not the pipeline's stated 1-sigma.
> The stated sigma is an empirical split-half figure from far-field cells, by coherence class (median 4.75 mm/yr;
> `data/insar_m2/v2_stated_sigma.npz`, 9d4a7da). Against the 2.3 to 2.4 mm/yr seen between neighbouring cells it
> is wide, not narrow, and V2's calibration row passes (6 of 6). The lag-structure reading (temporally correlated
> noise) stands; only the comparison with the formal floor misled. Part D's sigma medians are also formal values for
> both methods. What fails V2 is the false-rate field (the OK02 patch): pair share and the noise benchmark.

Data: `data/insar_m2/v2_pawnee_points.npz` (cd9ceab). Script: `python -m physics.insar_m2 [synth|A|B|C]`;
results in `results/physics/insar_m2.json`; figure `physics/figures/p10_method2_pawnee.png` (`physics/figures10.py`).
V2 is spent data. Nothing here is evidence that method 2 passes.

## Synthetic check
- (a) null: nothing named in 5 of 5. Pass.
- (b) planted coherence slope (10 % of variance): named in 5 of 5. Pass.
- (c) planted ministack steps: FAIL. B1 as written (12 free step columns plus a line) makes the rate badly
  determined, so the rate spread rises (cut −0.9 to −1.3) with or without a plant. **B1 cannot be read**: its
  "not named" below is a broken test, not a null.

## Part A, rate field (3,000 points per sample, A2)
- Nothing named. Both samples pick grammar shapes of coherence × incidence and closure ambiguity × east, but no
  term reaches 60 % of bootstrap resamples (best 0.5 to 0.6), and the search saves 0.094 and 0.078 bits per reading
  against 0.071 and 0.058 for shuffled data (p = 0.1 both). sd 2.88 → 2.63 and 2.75 → 2.50 mm/yr (about 9 %).
- Neighbour test: structure left (p = 0.01 both): the remainder is spatially organised, and no quality or geometry
  word names it.

## Part B, series
- B2 (3 km tile common mode per date) saves bits (164 per point) but cuts the rate spread only 9 %. Not named.
- B3: within each ministack the series rises about 9 to 11 mm/yr against the reference and steps back at joins.
  This is a scene-wide sawtooth relative to the reference cell; it barely varies by point (a sawtooth term
  per point saves −2.2 bits), so it moves the common level, not the spread.
- B4 (drop the first year): per-point rates correlate with the full-record ones (r 0.97 for a winter-only fit);
  the false rates are steady, not a first-year artefact.
- Lag structure (exploratory): the typical change between two dates grows from 8 mm (1 date apart) to 14 mm (4
  years apart). The noise is temporally correlated, so the white-noise stated sigma (median 0.72 mm/yr) is about
  3× too small: the rate spread between cells 0.3 km apart is 2.4 mm/yr, between places 10+ km apart 3.0.

## Part C, GNSS (described)
Station rates relative to OK03 (radar 200 m mean / GNSS LOS, mm/yr): OK02 −7.1 / −0.4, OK04 −1.4 / +0.1,
OK05 +0.6 / −3.6 (GNSS has gaps), OK06 −0.6 / −0.1, OK07 +0.6 / −3.4.
- OK02 sits in a NW patch of about −4 mm/yr (8 km across) that the radar shows in winter-only dates too; its GNSS
  does not share it. A radar error, mechanism not named.
- OK07's GNSS drifts about 10 to 15 mm against OK03 over the first 2.7 years, while the radar ground around it,
  out to 2 km, stays with OK03. Most likely the station's own monument, not the radar (not proven).

## Reading
Method 2's misses on stable ground are (1) error bars about 3× too narrow, from temporally correlated noise, and
(2) a spatially organised false-rate field (the NW patch), which no quality, geometry or seasonal word explains.
Masking by quality would not fix it. Next: the method 1 (OPERA DISP-S1) control on the same box and dates, with
the GNSS tie off, running on THE-COW. If method 1 shows the same patch, the cause is the ground (farmland), not the
processing.

## Part D, the method 1 control (A3; data/insar_m1 at 71b575a; `python -m physics.insar_m1`)
126 shared dates. Method 1 = OPERA DISP-S1, GNSS tie off, **no weather correction**; method 2 had HRRR.
- **Method 1 is far noisier here.** Rate spread (robust sd) is 11.4 against 2.8 mm/yr. Cells 0.2 to 0.4 km apart
  differ by 8.0 against 2.3. The typical change at 120 dates apart is 46 against 14 mm. Station pairs within
  2 mm/yr: 3 of 15 against 4 of 15. Stated sigma medians are 1.41 and 0.66.
- **D1, the corner:** by the fixed rule the verdict is "ground or common physics" (m1 NW − rest −7.4 mm/yr). But
  method 1 carries a box-wide ramp and a patchwork of ±10 mm/yr tiles, and the two methods' tile maps correlate
  only r = 0.20. The compressor's cross-method term `m2_tile` was not named (bootstrap 0.3). So method 1 does not
  clearly reproduce method 2's corner; the rule's verdict is driven by method 1's own large-scale error.
- **D5, compressor on method 1:** sd 14.5 → 9.75 (inputs east/north, annual amplitude, coherence). The winner is
  `ramp` plus seasonal-amplitude and coherence shapes. Only one of them clears 60 % of resamples
  (annual_amp > 6.3 · coherence, 0.7), and the null search's p is 0.1. Neighbour test: structure left in east × north.
  Read: a long-wavelength (troposphere-like) ramp plus seasonal-ground leakage, which a weather correction and a
  GNSS tie would target.

## Reading after part D
At this rural site our own processing (method 2) beats the product the reports use (method 1, untied and
without weather correction) on every reading. V2's failure is not evidence that method 2 is worse than method 1.
Pawnee is a hard site for both. Method 1 passes our CHI checks with a GNSS tie and reference-jump handling at
urban places. So method 2's fair test is the same: CHI's checks, method 2's error bars widened for correlated
noise (as in insar_pipeline PR #21 for PS), at a fresh place, written first.
