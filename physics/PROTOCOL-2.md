# Missing physics by code length: flood, heat and ground (written 22 Nov 2026, before any of these residuals was fitted)

Same method, tool and rules as [PROTOCOL.md](PROTOCOL.md) (`codelen.lawsearch`): a correction formula counts only if
it stores the errors in fewer bits, after paying for its own description. Formulas have up to 3 terms from a fixed
vocabulary, plus an intercept, and each case has a hold-out split and a model-free neighbour test.
Kyle asked for this on 22 Nov 2026 ("the same for the floods ... and the heat model ... and the insar if that works").

**What had been seen before writing this:** the summary rows of FloodSim PV-SS10 and PV-SS13 (`RESULTS.md`, F15, F21);
the HeatSim station-check summary for PHL-POC-001 (bias and MAE per station and window, HeatSim 0.4.1 runs); the
ground synthetic scenario definitions in CHI `tools/gate4/synthetic.py` and the test-run summary. No per-point residual
of these cases had been fitted or plotted.

## H. HeatSim, Philadelphia, open-site air temperature with the station correction off
- **Data:** HeatSim 0.6.1 runs `phl_noanchor_2025-07-{08,25,28}.nc` (Kyle's PC, CHI `runs/heatsim/`), CDS forcing.
  Variable `open_site_temperature` at the cell of each land ISD station in the domain (Wings Field, PHL, NE PHL),
  against the station's hourly TMP (NOAA ISD on AWS `noaa-global-hourly-pds`). Reports within each UTC hour are
  averaged. Quantity: air temperature (open site, 2 m). Not independent: ERA5 assimilates these stations.
- **Target:** r = station − model (K). Grid 0.1 K, −20 to 20, eps 1e-3.
- **Vocabulary:** `sun` = max(sin solar elevation, 0); `sun_lag(τ)` = sun τ hours earlier, τ fitted 0-6 h (heat
  stored and released late); `night` indicator; `hsin`, `hcos` (24 h harmonics, a shape with no physics name);
  `wind` and `1/wind` from the station's own WND; `sun/wind`; `dew` (station dewpoint, °C); `cloud` (sky cover
  fraction from GF1, or GA1 when missing; 0.5 where neither is reported); `Tm` (the model's own temperature);
  `stn_PHL`, `stn_NE` (station indicators; with three stations, land cover cannot be told from station identity).
- **Hold-out:** leave one day out (3 ways) and leave one station out (3 ways). A formula transfers if it beats
  "one correction" on every held-out fold.
- **Neighbour test:** each hour's error guessed from the same hour at the other stations, against 199 shuffles of
  each station's errors over time.

## F. FloodSim, inlet exchange
- **F1, the Sheffield rig (PV-SS10, Rubinato et al. 2017 S1, 15 steady tests; measurements):** r = ln(measured /
  FloodSim) for the exchange into the sewer, FloodSim = 0.5 · πD · √(2g) · h^1.5 with h = the table depth, as
  `score.mjs` computes. Vocabulary: `lnh` (a different exponent), `lnQ1` (approach flow, i.e. the velocity at the
  inlet), `datum(δ)` = 1.5 ln((h + δ)/h) with coefficient fixed at 1 and δ fitted −5 to +10 mm (the depth datum),
  `Fr` (approach Froude number from Q1, the rig's 1 m width and h). Grid 0.005, −3 to 3. Hold-out: odd tests
  against even tests, both ways. With 15 points, "one factor" winning is a likely and acceptable outcome.
- **F2, a grate on a grade (PV-SS13, 9 cases; the reference is HEC-22, not a measurement):** r = ln(FloodSim share /
  HEC-22 share). Vocabulary: `V` (gutter velocity), `T` (spread), `lnQ` (flow reaching the grate), `Eo`
  (HEC-22's frontal share), `SL` (grade), `d` (depth in the grate cell). The cause is already known (F21: the
  exchange uses depth only, HEC-22 uses velocity), so this is a test of whether the search names it: **V** or a
  velocity-carrying term must rank first. Hold-out: leave one rain rate out (3 ways).

## G. Ground chain, the four synthetic towns (truth known)
- **Data:** CHI `tools/gate4/synthetic.py` `town()` at its default seeds, CHI `e4b98ae`. Truth is the planted
  line-of-sight rate. The chain's output is the straight-line velocity of each cell.
- **Target:** r = truth − chain velocity (mm/yr), valid cells, a fixed random sample of 3,000 per scenario (seed 0).
  Grid 0.01 mm/yr, −40 to 40.
- **Vocabulary (from the chain's output only, never the truth):** `vel`, `std`, `amp`; `late` = rate over the last
  third of each series minus `vel`; `accel` = rate over the second half minus rate over the first half; `dist`
  (distance to the house, km); `near` (within 250 m of the house).
- **Expected beforehand:** in `sinkhole`, a term in `late` or `accel` (the chain's straight line misses an onset in
  2022); in `quiet` and `bad-reference`, one correction (the reference offset); in `regional`, nothing beyond a factor.
- **Hold-out:** cells west of x = 0 against cells east, both ways. Neighbour test: 4 nearest earlier cells within
  100 m, against 199 shuffles.

## Limits (fixed now)
None of this validates a model. Heat has 3 days at 3 stations, and the errors are strongly autocorrelated over
hours, so the effective sample is small. F2 and G score against references, not measurements. A term found here is
a lead for each repository's own protocol and is never applied to a report.
