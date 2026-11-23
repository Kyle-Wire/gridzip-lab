# Surprise score: findings (22–23 Nov 2026)

Protocol: `PROTOCOL.md` (committed 49429b0 before the synthetic street ran). Code: `core.py` (the score),
`flood_block.py` (real block), `street_test.py` (synthetic street), `figures_block.py`, `figures_street.py`.
Figures in `figures/`. Climate Home Inspection's proposal for the report is `SPEC-SURPRISE.md` there.

## What the score is
Bits to describe a house's modelled readings given its neighbours, minus what a typical house would cost. A mixture of
"like neighbour j" kernels, "like itself last storm" and a wide fallback, coded storm by storm (chain rule, 5 % fixed
share), with codelen's integer tables. 0 means as its neighbours suggest; each extra bit halves how often a typical house
would look that odd. On null data (independent noise, 400 units, 4 readings) the 95th / 99th percentiles are 3.6 / 5.9
bits, against 3.9 / 6.7 for an ideal Gaussian score, so the score is close to calibrated.

Two early versions failed and were replaced before any test was scored: a single Gaussian per house (neighbour spread too
noisy, typical houses at -3 bits and a heavy tail) and a pure kernel mixture without the self and wide experts (a house
unlike every neighbour paid the 14-bit escape in every storm).

## Test 2: the real South Street block (calibration, settings were tuned here)
1,062 L&I footprints, FloodSim runs of 16 Nov 2026, six design storms (10/100/500-yr, 1 h and 24 h; 1-h storms
averaged over 4–5 run variants, whose spread is the model noise).
- Round trip: the whole block (22,356 bits) was arithmetic-coded and decoded exactly. **Pass.**
- Location carries information: real 95th percentile 12.4 bits, shuffled 5.3; neighbour correlation of surprise 0.16
  real vs 0.04 shuffled. Odd houses cluster. **As expected.**
- Height explains almost nothing here: median and 95th percentile barely move given height (-0.30 / 12.3 bits). On this
  flat rowhouse block the oddities are not about ground height.
- Planted wetter houses (30 houses, 3 repeats): median extra bits +0.9 (+2 cm), +1.6 (+5 cm), +2.0 (+10 cm), +5.0
  (+20 cm); share then ranked in the block's top 5 %: 7 %, 3 %, 2 %, 8 %. **The score is a weak detector of small
  local changes on this block**, because the block is already very varied (many neighbours within 30 m already have
  wet rings, so "like that neighbour" is cheap).
- PHL-POC-001: -4.8 bits, 10th percentile: one of the most typical houses on the block, about 0.5–1 cm wetter than its
  neighbours' weighted ring depth. (Printed only in the private run; the committed JSON leaves the homes out.)

## Test 1: synthetic street with planted answers (the test, fixed before running)
FloodSim CPU solver, 16 lots, 2-, 10- and 100-yr 1-hr storms (`benchmarks/emulator/street.mjs`).

| Planted | Surprise (rank of 16) | Ring depth, 100-yr | What the model did |
|---|---|---|---|
| lot 4, back-yard low spot | **11.8 bits (1st)** | 16.4 cm vs 4–9 cm around | a pond against the back wall |
| lot 9, neighbour's fill | -0.2 (8th) | 6.9 cm, like its neighbours | barely any extra water reached lot 9 |
| lot 13, whole lot 25 cm lower | -0.2 (7th) | 8.5 cm, like lots 10–11 | the street water rose with it; no concentration |

- P1 (lots 4 and 9 in the top 3): **fail**. P2 (lot 13 top 4, half explained by height): **fail**. P3 (no unplanted
  lot above lots 4 and 9): **fail** (lots 0, 1, 2, 6 and 12 score 2–5.5 bits; lot 12 is a dry lot beside the lower one).
- Diagnosis (after the result, so not a test): two of the three plants did not change the modelled water at the house.
  Lot 9's back-yard 90th-percentile depth is 2.3 cm (100-yr) vs 1.2–4.3 cm on the others; the 0.30 m fill on lot 10 slopes
  toward lot 10's own house, and only its 1.5 m ramp sheds west. The score can only flag what the flood model shows.
- The one plant the model expressed was found clearly, first by 6 bits.

## What this means
- The score finds houses whose modelled water breaks the local pattern. It is not a detector of site features the flood
  model does not turn into water at the house. That is said in the report wording (SPEC-SURPRISE).
- On a varied block it needs a big local difference (about 20 cm) to stand out. It will be sharper on blocks of similar
  houses (suburban streets) than on dense rowhouse blocks.
- Next: a fairer street test with the fill shedding into the neighbour (protocol first), the same score on HeatSim's
  street PET and on ground rates, and a "surprise over time" version for change detection (a house whose score jumps
  between runs).
