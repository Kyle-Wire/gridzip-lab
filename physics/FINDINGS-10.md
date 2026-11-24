# Findings 10: method 2 (V2 Pawnee) on stable ground (protocol 10, A1, A2; 24 Nov 2026)

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
