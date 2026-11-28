# Protocol 16: name the LBM's Case D error from the district's geometry (WindSim F3 steps 8 and 9) (written 28 Nov 2026, before any step 9 run)

**Why.** WindSim F3 step 8 ran the LBM on AIJ Case D (a tower among low-rise blocks) for the first time. It reads low,
×0.711 overall, and the error grows with distance from the tower: ln(model/measured) falls by 0.78 per 0.1 m (r −0.70),
from ×0.85 within 0.1 m to ×0.63 beyond. Case C, a cube array resolved with 16 cells per block height, reads ×1.25
high. D's low-rise blocks are 25 mm tall, which is 4 cells at 6.25 mm. Kyle's standing rule: extend the compressor with
the domain's formulas and run it against the model's errors.

**Seen before this was written:** step 8's D points (`WindSim docs/validation/f3_step8/D/points.csv`), including the
distance reading above, and the step 9 protocol. No step 9 value exists yet.

**Data.**
- Step 8 D, raw (no level factor): 66 points at z = 5 mm, wind 0°. Folded about y = 0, as the points come in mirror pairs.
- Each step 9 D run (D9-W walls, D9-16 finer grid), when it lands.
- The district's height map: from the case STL (`BD_caseD.stl`, Zenodo 10.5281/zenodo.15402668, CC BY 4.0),
  voxelised on THE-COW at 3.125 mm and committed as `results/physics/caseD_heights.npz`.
- Target: y = ln(model / measured) per point, on grid 0.01.

**Terms (fixed now).** Each one is computed from geometry only, never from model output:
- `d_tower`: horizontal distance from the tower's axis, in tower widths B.
- `wake_side`: cos of the angle from +x (downstream is +1, upstream −1).
- `shelter_up`: the height of the tallest low-rise block within 2 block heights upwind of the point (along −x), over H_low.
- `gap`: the open width across the flow at the point (the distance between the facing block faces in y), over B.
- `in_canyon`: 1 if blocks stand on both sides of the point within 1 B across the flow, else 0.
- `tower_view`: 1 if the straight line from the point to the tower crosses no block, else 0.
- `cells_low`: cells per low-rise block height (4.0 at step 8, 6.4 at D9-16). This is a run term. It is only used in
  the pooled search below.

**Searches and readings (fixed now):**
- `LF.search` with `T.WIND` and these terms, grammar on, the same guards as protocols 9 and 14. Folds: x < 0 against
  x > 0, and near (d_tower < 1.6) against far.
- **S8, step 8 alone.** Does a geometry term carry the distance pattern? A term **pays** with inclusion ≥ 0.8, a
  per-input bootstrap ≥ 60 % and transfer on both fold pairs, at null p ≤ 0.05.
  - First a synthetic check: pure noise at the residual sd must name no term in ≥ 9 of 10 seeds (seeds 1600 to 1609).
    If that fails, S8 is reported as "untrustworthy".
- **S9, each step 9 run.** The same search on that run's y. Does the paying term (if any) keep its sign, and does its
  coefficient shrink by half or more?
- **Pooled.** Steps 8 and 9 stacked, with `cells_low` and `walls` as run terms. Does `cells_low` or `walls` pay, or
  does the geometry term's coefficient change with them (an interaction in the grammar)?
- Each reading is reported as named, not named, or untrustworthy. Nothing here changes WindSim's scoring.

**Limits.**
- 33 independent points after the fold. Only one term or two at most can pay.
- The terms come from geometry chosen after seeing the distance pattern, so S8 can only name a word for that pattern.
  It cannot be counted as a test.
- Runs differ by one setting at a time, but there is one seed each.

Outputs:
- code: `physics/district_fields.py`;
- figure: `physics/figures/district/`;
- results: `results/physics/district_fields.json`;
- write-up: `physics/FINDINGS-16.md`.
