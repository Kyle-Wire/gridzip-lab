# Surprise score: protocol (written 22 Nov 2026, before the synthetic street was run)

## What the score is
`surprise/core.py`. For each house, the code length (bits) of its modelled readings given its neighbours, minus what a
typical draw from the same model would cost. The neighbour model is a mixture: "like neighbour j" kernels (distance
weights, sigma 15 m within 30 m), "like itself in the last storm" (from the second reading on), and one wide fallback,
coded storm by storm with the chain rule and a 5 % fixed share. Readings are coded with codelen's integer tables, and on
the real block the whole code is arithmetic-coded and decoded (round trip checked). "Surprise given height" restricts
the neighbours to those at about the same ring ground height (tau fitted by code length).

Settings fixed while developing on the South Street block (so South Street is calibration, not a test): coding scale
log10(depth + 3 cm) at 0.02, radius 30 m, expert priors (self 0.3, wide 0.02), fixed share 0.05, kernel width chosen
from 9 values by shortest code.

## Test 1: synthetic street with planted answers (benchmarks/emulator/street.mjs in FloodSim)
16 lots, ordinary random variation per lot (seed 5), three planted lots:
- lot 4: back-yard low spot behind the house (-0.22 m)
- lot 9: the next lot east fills its back yard by 0.30 m and sheds onto lot 9
- lot 13: the whole lot sits 0.25 m lower

Readings: 90th-percentile max depth in the 2 m ring around each house, 2-, 10- and 100-yr 1-hr storms, CPU solver.
Neighbours here are the lots within 40 m (2.5 lots each side), since the street has one row of 16.

Expected answer and pass rule, fixed now:
- P1: lots 4 and 9 are both in the top 3 by surprise. Pass if both are.
- P2: lot 13 is in the top 4 by surprise, and its surprise given height is at most half its surprise (its height explains
  at least half). Pass if both hold.
- P3: no unplanted lot scores above the lowest-scoring of lots 4 and 9. Reported either way.
A fail is reported as a result, and the settings are not changed after seeing it.

## Test 2: the real block, as a check of the machinery
- Round trip: decode equals the readings (pass/fail).
- Location carries information: the 95th percentile of surprise on the real block exceeds that with readings shuffled
  between houses within each storm, and surprise is more correlated between neighbours than after shuffling. Reported.
- Planted wetter houses (+2, +5, +10, +20 cm in every storm, 30 houses, 3 repeats): the median extra bits and the share
  ranked above the block's 95th percentile. Reported as a sensitivity curve, with no pass rule (no prior expectation was
  written for it before the first look).
