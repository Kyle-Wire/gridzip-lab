# Findings 16: Case D's miss and the district's geometry (protocol 16, S8 part; 28 Nov 2026)

**Run:** `python physics/district_fields.py --runs f3s8_D`. Data: WindSim F3 step 8 Case D, raw
(`results/physics/f3s8_D_points.csv`, copied from WindSim `docs/validation/f3_step8/D/points.csv`), and the height map
`results/physics/caseD_heights.npz`, built on THE-COW from the AIJ STL. The map has one tower 0.25 m tall and low blocks
all 0.025 m tall. Mirror pairs are folded, giving 33 points. Figure: `physics/figures/district/district_f3s8_D.png`.

**Synthetic check: passes.** On pure noise at y's sd, no geometry term was named in 10 of 10 seeds (seeds 1600-1609).

**S8: the winner is `d_tower` alone, but it does not pay under the protocol's rule.**
- The winner is ln(model/measured) = −0.34 − 0.48 · (d / B). It cuts the sd from 0.31 to 0.23, saving 0.014 bits per
  reading, at null p 0.05.
- On the rule's other tests: inclusion 0.99, bootstrap 70 %.
- Folds: it transfers upstream ↔ downstream both ways and near → far, but **not far → near**. The rule needs both
  directions of both pairs, so `d_tower` does **not pay**.
- The guards add that a squared form, d², cannot be ruled out (+2.8 bits).
- None of the mechanism words has any weight: `shelter_up`, `in_canyon` and `tower_view` have inclusion 0, `wake_side`
  0 and `gap` 0.17.

**Reading.** As fixed in the protocol, S8 is reported as **not named**. The pattern is real and repeats across the
upstream and downstream halves. But the compressor finds no local geometric cause for it. Sheltering by an upwind
block, sitting in a street canyon, and having a clear line to the tower all explain nothing beyond distance itself. A
miss that grows smoothly with distance, with no local geometric cause, is what one would expect if the whole low-rise
canopy flow is too slow and the tower's downwash props up the points near it. That favours a canopy-wide cause, such as
under-resolved low blocks (4 cells tall) or missing wall shear, over a local one. Step 9's D9-W and D9-16 runs test
exactly those two. As the protocol says, this cannot be counted as a test: the terms were chosen after seeing the
distance pattern.

**Still to run under this protocol:** S9 on D9-W (queued on THE-COW after flood's Lumber run) and on D9-16 (waiting on
Kyle), then the pooled search.

## S9 on D9-W (walls on) and the pooled search (28 Nov 2026)

`python physics/district_fields.py --runs f3s8_D,f3s9_DW`. The D9-W points are WindSim `f3_step9/D_walls/points.csv`, at 7389651.
- **D9-W alone: not named.** The synthetic check passes (10 of 10 silent). The search keeps no term: the winner is the
  intercept only, at null p 1.0, although the distance pattern is the same as step 8 by eye (WindSim's step 9 reading:
  slope −0.77 against −0.78). On 33 points, distance falls just short of paying for its own description in this run.
  So S9's question, whether the term keeps its sign and shrinks by half, cannot be answered from the search. The direct
  reading in WindSim's step 9 result is that it does not shrink.
- **Pooled (steps 8 and 9 stacked, 66 rows): not named.** The winner is a grammar shape, max(0, d_tower − 1.72) · lp:
  the miss grows beyond 1.7 tower widths, where the model's own speeds are low. `d_tower` has inclusion 1.0 and bootstrap
  100 %, but it fails the fold rule. `walls` has inclusion 0. `cells_low` cannot be read yet, because both runs have 4
  cells per block height. That needs D9-16.
- **Reading.** The wall model changes neither the pattern nor the level. What remains is a miss that grows with
  distance, worst in the far streets where the model is slow. That is consistent with an under-resolved low-rise canopy,
  which D9-16 tests.

## S9 on D9-16 (16 cells per tower) and the three-run pool (28 Nov 2026)

`python physics/district_fields.py --runs f3s8_D,f3s9_DW,f3s9_D16`. The D9-16 points are WindSim `f3_step9/D_b16/points.csv`, at ad43f9a.
- **D9-16 alone: not named.** The synthetic check is silent in 10 of 10 seeds, and the winner is the intercept alone
  (level −0.134 in ln, ×0.874). The distance slope has halved, by WindSim's direct reading, and the search no longer
  finds it worth a term.
- **Pooled (99 rows): the compressor writes the grid effect itself.** The winner is
  ln(model/measured) = −0.28 + 0.042 · gap − 0.22 · d_tower² / cells_low.
  - The miss grows with the square of the distance from the tower and **shrinks in proportion to the cells per
    low-block height**.
  - `cells_low` has inclusion 1.0 and bootstrap 67 %; `d_tower` has inclusion 1.0 and bootstrap 100 %. Null p is 0.05,
    and `walls` stays at 0.
  - It transfers upstream ↔ downstream and far → near, but **not near → far**, so under the fixed rule nothing
    **pays**. The reading stays "not named".
  - `gap`'s bootstrap is 53 %, and it does not pay.
- **Reading.** Every guard short of the fold rule points the same way as WindSim's step 9 reading. The low read away
  from the tower is a resolution error of the low-rise canopy, and it falls as the blocks get more cells. The fold that
  fails is the extrapolation from near points to far ones, which is expected: the near points barely carry the effect.
