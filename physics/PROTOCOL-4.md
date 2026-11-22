# Protocol 4: the physics compressor (a curated term library and sparse search), validated, then run on every model

Written 22 Nov 2026 before any real model's errors were scored with the new library (`physlib/`), and before the
full synthetic validation ran. Kyle asked for it on 22 Nov 2026: a curated, per-domain library of known physical
terms, a search over sparse combinations of them with a bit penalty against overfitting, a report of which terms
explain each model's errors and what is left, validated on planted laws with a false-discovery check on noise,
then run against WindSim, HeatSim, FloodSim and the ground (InSAR) chain after the day's fixes.

**What had been seen before this was written.** Everything in PROTOCOL 1 to 3 and their findings (the same real
cases, scored there with hand-picked per-case vocabularies). During development, one seed of each synthetic
scenario below was run to check the code runs (all named what was planted except: roughness ×3 was named as the
gust-factor term, its declared alias; the direction part of "station + direction" was missed; the 3 mm datum was
named as the weir exponent). No real case has been run with `physlib`.

## The tool
- `physlib/engine.py`: two-part code length, exactly as `codelen.mdl` / `codelen.lawsearch` (structure: Elias-delta
  token count + log2(library size + 2) per token; parameters rounded at their best step; residuals binned under a
  Gaussian with codelen's escape). Terms may carry several columns (a harmonic pair, station offsets) and one
  fitted family parameter. Search: a beam (width 6, up to 5 terms) ranked by an approximate code length, the 10 best
  formulas plus "no term" re-scored exactly, then backward elimination. The winner's residual is re-coded with the
  arithmetic coder and decoded (round trip).
- `physlib/terms.py`: the library, fixed now: wind 19 terms, heat 14, flood 13, ground series 9, ground rate 7.
  Terms that say the same thing to first order share a group, and a formula holds at most one per group
  (`log_law_z0` and `gust_factor`; sunlight and sunlight-lag; grade and ln grade; the two fetch forms).
  `diurnal_shape` is marked not-physics: when it wins, the library could not name the shape.
- "One correction" (intercept only) is the baseline every formula must beat.

## A. Synthetic validation (`python -m physlib.validate`)
Each domain has realistic inputs and three planted laws (`physlib/synth.py`), 20 seeds each, noise at the stated sd.
Pure-noise cases: 100 seeds per domain of independent noise, and 50 of AR(1) noise (ρ = 0.8). A sweep scales the
first law of each domain from 0.1× to 1×, 10 seeds each.

**Pass rules (fixed now; never changed after the run):**
1. Each planted law, at its stated size: the correction the winner fits reproduces the planted signal (R² ≥ 0.9)
   in at least 80 % of seeds ("law explained"), and the winner names exactly the planted terms (a declared alias
   counting as a match) in at least 70 % ("named exactly").
2. Independent noise: in each domain, at most 5 % of runs name any term (false-discovery rate ≤ 5 %).
3. AR(1) noise: reported, no gate. The code assumes independent readings, so correlated noise is expected to raise
   false discoveries. That is why every real case below also needs a block hold-out to transfer.

A failed rule is reported as a result. Any change to the tool after seeing part A is an amendment here, dated,
saying what was seen, and made before part B runs.

## B. Real cases, before and after the 22 Nov fixes
Same quantities, grids and hold-out folds as PROTOCOL 1 to 3, now with the domain library in place of the per-case lists.
- **W1 WindSim AIJ Case C** (WindSim 3c780ff points, LBM D/10 blind and RANS D/10): r = ln(measured/model).
  Inputs: `lp` = ln(model/U_app), `fetch` = (x + 0.5)/D, `along`, `cross`, `hc`, `corner`. Hold-out: layouts
  (leave one of 0D, 1D, 2D out). No WindSim street fix landed today, so there is no "after" here.
- **W2 WindSim Layer 1** (PR 10, after F2): thunderstorm and ordinary pairs, r = ln(observed/predicted). Inputs:
  `stn_loo`, `dz` (elevation relative to the neighbours), `dirn`, `hour`, `month`, `z0` with z = 10 m. Hold-out: leave
  one place out. "After" is the F5 per-station z0 run if its result files exist when part B runs; otherwise reported
  as pending. Never use same-storm terms (PROTOCOL-3 A1).
- **H HeatSim 0.6.1, Philadelphia, 5 hot days, station correction off, open-site air temperature at 3 ISD stations**:
  before = `fixes2_phl_noanchor_stamp_*`, after = `fixes2_phl_noanchor_hour_ending_*` (the sunlight-timing fix,
  HeatSim M1). r = station − model (K). Inputs: sun (NOAA solar position at each station), model temperature and
  its hourly slope, station wind, cloud, dewpoint, station tree and built fractions from the run's cell, station id.
  Hold-out: leave one day out (5) and one station out (3).
- **F1 FloodSim Sheffield rig (PV-SS10)**: before Cw 0.5 (the PROTOCOL-2 data), after Cw 0.32 (FloodSim 8f525e9).
  After is calibration (0.32 was fitted to these tests), so it is reported as a check that one factor is gone, not
  as evidence. **F2 grate on a grade (PV-SS13)**: before depth capture (17 Nov run), after `capture: 'hec22'`
  default (8f525e9, rerun on the CPU here). Hold-outs as PROTOCOL-2.
- **G ground chain, four synthetic towns (truth known)**: before = the straight-line chain (PROTOCOL-2), after = CHI
  A11 (reference step check, after-onset rate near the house). Two levels: (i) per cell, r = truth − chain rate,
  with the ground-rate library; (ii) per series, the network mean of (series − the chain's per-cell fit) and a
  fixed sample of 200 cells (seed 0, all within 300 m of the house plus a random rest), with the ground series
  library. Hold-out for (i): west against east. For (ii): first half of the dates against the second.

**What is reported for each case:** the winning terms with their fitted values, the bits ladder (one correction,
then each winning term in order of what it saves), whether it transfers on every hold-out fold, the coder round
trip, the error sd before and after, and the remainder. The remainder is said in two ways: the sd left, and where a
neighbour context exists (W1 points, H same-hour stations, G cells), the model-free neighbour test on the remainder
(bits gained by guessing each remaining error from its neighbours, against 99 shuffles): structure left that the
library cannot name.

## Limits (fixed now)
A winning term is a lead for that model's own protocol, never a fix applied to a report. Heat errors are strongly
correlated in time, so its effective sample is far below 360. F2 and G score against references, not measurements.
The library only names what is in it: a gap with no term reads as a remainder, never as a discovery.
