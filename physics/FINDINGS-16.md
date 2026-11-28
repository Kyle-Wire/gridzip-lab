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
