# BR1 Las Vegas Valley points for the compressor (method 2, phase-linked, with the bridged labels)

These tables come from experiment BR1 (insar_pipeline `docs/BRIDGED_PLAN.md` A1 and A2, protocol
`review_2026-11-25/protocol_BR1.yaml` version 2, branch `bridged-components`, PR #22). It is the fresh-site test of
the bridged-cells rule. Las Vegas Valley is a public test site, and its box was drawn from public facts only: the 7
NGL stations within 15 km of the city centre, plus 0.03°. So positions are included.

Each table is exported by `export_lv.py`, read-only on the project. The layout follows `data/insar_m2/v2_pawnee_points.npz`.

| file | run | what changed |
|---|---|---|
| `br1_a_points.npz`, `br1_a_meta.json` | BR1-A | as now: snaphu CONNCOMPTHRESH 300, MintPy bridging+phase_closure |
| `br1_b_points.npz`, `br1_b_meta.json` | BR1-B | CONNCOMPTHRESH 1000, re-unwrapped from A's phase linking |
| `br1_c_points.npz` (when C is done) | BR1-C | no bridging (phase_closure only), re-unwrapped from A's phase linking |

B and C are read at A's cells (`--points-from`), so all three tables line up row for row.

**Run A:**
- Sentinel-1 descending track 173, IW2, 3 bursts a date.
- 148 dates, 2021-10-04 to 2025-11-12. 20251013 was dropped because the catalogue lists it twice.
- 438 interferograms, with HRRR weather and DEM-error corrections.
- Grid 713 × 754 radar cells (3 azimuth × 11 range looks), 537,526 valid cells.
- 200,000 points: every valid cell within 2 km of a station, plus a uniform random sample of the rest (seed
  20261125).

## Fields (the m2 fields, plus these)

| field | unit | what |
|---|---|---|
| formal_sigma_mm_yr | mm/yr | MintPy's formal `velocityStd` (white noise). Not the stated 1-sigma |
| stated_sigma_mm_yr | mm/yr | the pipeline's stated 1-sigma: the empirical split-half figure by coherence class (`deliverables/uncertainty_model.json`), and the bridged figure where a cell is bridged and one was measured |
| bridged_fraction | 0-1 | share of the kept interferograms in which the cell is outside the reference pixel's snaphu component (label differs, or either is 0) |
| bridged_region | int | 0 = the reference's component; 1..k = bridged regions by size; -2 = bridged but under 50 cells; -1 = no data |
| bridged | bool | bridged_fraction ≥ 0.5 (region > 0 or -2) |
| station_bridged_points_200m | count | bridged cells within 200 m of each station |

`series_0p1mm` is int16 in 0.1 mm. **-32767 means no value**, where the m2 table had NaN.

## What run A showed (recorded, not a finding of this folder)
- **Bridged cells:** 38 of 537,526 valid cells (0.007%), with no numbered region and none within 200 m of a
  station. Las Vegas has essentially no ground that only MintPy's bridging joins on.
  - So the bridged rule changed nothing here. This test shows the rule does no harm. It does not show that it
    helps, as the plan said a site without bridged ground would.
- **Graded rows** (identical with bridged cells left out and with every cell counted):

| row | result | grade |
|---|---|---|
| rate_rmse_dd | 1.1 mm/yr (95% CI 0.47 to 1.59) | pass |
| pair_share | 21 of 21 pairs within 5 mm/yr | pass |
| rate_bias_loo | +0.04 mm/yr | pass |
| uncertainty_calibration | 7 of 7 inside | pass |
| noise_benchmark | worst 68th percentile 2.05 mm/yr | pass |
| CHI A12 pair (NVNL, NVCA) | 1.54 mm/yr difference against 1.64 tolerance | pass |

## What run B showed (recorded, not a finding of this folder)
- CONNCOMPTHRESH 1000 reached snaphu: every pair re-unwrapped, because its config had changed. But the components were the same at every one of the 200,000 points (identical `bridged_fraction`, `n_ifgrams_outside_component`, `in_conncomp_mask`), and the bridged cells were still 302 in the whole grid. Rates differ from A by at most 0.003 mm/yr.
- Every graded row passes, both with bridged cells left out and with every cell counted, and the figures match A's. CHI A12 pair: 1.54 against 1.64, pass.
