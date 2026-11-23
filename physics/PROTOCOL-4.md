# Protocol 4: the Philadelphia house patch against its wide ring, real InSAR (written 23 Nov 2026, before any fit)

Ask (Kyle, via the InSAR thread): A11 flagged a +10.4 mm jump at 2019.542 in patch − wide ring for PHL-POC-001. Is it
an isolated bad date, or part of structure we don't model?

**Seen before this was written:** the InSAR thread's description (one-date dip at 2019.51, rel −10.3; robust SD of
date-to-date changes 1.44 mm in rel and 3.8 mm in the wide ring; M4 finds no change). The series' first values were
printed to check the file's layout. Nothing has been fitted or plotted.

**Data:** `data/derived/PHL-POC-001/ground/series.json` (CHI main; copy at project files `ground-a11/`), route `disp`
(DISP-S1, line of sight, after the A11 reference correction), 254 dates 2016.75-2025.92. y = patch − wide ring (mm).
Grid 0.05 mm, −60..60.

**Vocabulary** (each one term; family parameters are fitted and paid for):
- `t`: linear rate. `step(a)`: a step at year a (2017-2025.5).
- `annual(φ)`: cos 2π(t − φ), φ fitted 0-1 yr; `semi(φ)`: the 6-month harmonic. A thermal (rowhouse expansion) signal
  would show as `annual` with its peak near late July (φ ≈ 0.55-0.6, Philadelphia's warmest weeks).
- `wdev`: the wide ring's value minus its own ±2-date running median (the common-mode, atmosphere-like swing); `wdev_lag`:
  the same at the previous date. A leak of the common-mode into rel pays here.
- `ar1`: the previous date's y (what a decoder already has).
- `spike(i)`: 1 on one date, the date i fitted over all 254 (so naming a date costs the bits to pick it from 254).
- `nr`: near ring − wide ring at the same date (does the patch move with its surroundings).

**Tests, readings fixed now:**
1. MDL search, up to 3 terms. Hold-out: first half → second half and second → first. A term counts only if it holds
   both ways.
2. Q1: if `spike` is chosen and no smooth term (annual, wdev) explains the date, it is an isolated bad date. If
   `wdev` or `annual` pays and the date's residual falls below the noise after it, the date is part of that structure.
3. Q3 (spike rule): for every date, d_i = y_i − (y_{i−1} + y_{i+1}) / 2 over a robust sd of d. Report the 2019.51
   date's z and how many dates exceed |z| > 5; under Gaussian noise none of 254 should.
4. Neighbour test (99 shuffles over time): context = previous two dates. If it pays, rel has temporal structure.
