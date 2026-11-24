# Findings FIRE-1 (running record; rules in PROTOCOL-FIRE.md)

## 1. Planted-mechanism checks (amendment A4), 23 Nov 2026

FireSim `5f00fe0`, `python physics/fire_planted.py`, result `results/physics/fire_planted.json`, chart
`physics/figures/fire_planted.png`. Twelve made-up landscapes (seeds 20-31), about 1,500 arrival-time readings per
case, noise sd 0.05 in log time.

| switched off | winner | holds in hold-out both ways | verdict under A4 |
|---|---|---|---|
| slope | `slope_run` (5.25 to 4.97 bits per reading) | yes | **PASS** |
| crown fire | `canopy`, `shelter`, `slope_run` (7.44 to 5.88) | yes | **FAIL**: the crown-fire term is not chosen; a slope term is |
| canopy sheltering | `wind_run`, `shelter`, `canopy`, `crown_fuel` (6.38 to 5.94) | no (one way only) | **FAIL**: sheltering is named, but with crown-fire and shape terms beside it |
| nothing (twin) | no term | - | **PASS**: silent on a model with nothing missing |

What this allows on real fires, and what it does not:
- A missing slope response can be named from arrival times, and the search stays silent when nothing is missing.
- A miss under forest canopy is detected (it is the first term chosen in both canopy cases, and worth 1.3 and 0.2
  bits per reading), but the search **cannot tell crown fire from wind sheltering**, and in the crown case it also
  credits slope, which was not touched. On real fires a canopy term may be read only as "something under canopy",
  never as crown fire or sheltering specifically, and a slope term found together with a canopy term is not read.
- Why, as far as these runs show: the three canopy terms are averages along a straight line and are strongly
  correlated on these landscapes (low canopy is a subset of any canopy; sheltering is cover times a constant), and
  crown fire acts through a threshold on the surface fire, which a path average cannot express. A sharper vocabulary
  (a crown term built from the candidate's own surface intensity against the canopy base height) is the candidate
  fix; it needs a new amendment and fresh seeds before it is scored.
- Embers and house-to-house spread have not had this check; they get it when they exist.

## Camp 2018, first real-fire runs (23 Nov 2026; calibration fire, nothing here is validated)

Protocol A5 and A6. FireSim `scripts/run_fire.py`, results and maps in FireSim `results/camp_2018/`. One seed each,
ember and house numbers are the uncalibrated starting guesses. Building record: USGS points, 29,948, of which 16,493
burned. "Bits" is per building with one seed (probabilities 0 or 1, clipped to 0.005 and 0.995); the fire's own
burned share would cost 0.99 bits, so every run below is still worse than knowing only that share.

| run | wind | burned at 21 h (mapped 376 km2) | Jaccard 9.45 h | Jaccard 21 h | AUC | bits | houses lit | wall time |
|---|---|---|---|---|---|---|---|---|
| M1 spread only | stations (A5) | 61 | 0.01 | 0.16 | 0.51 | 4.10 | 918 | 4 min |
| M3 embers + houses | stations (A5) | 332 | 0.03 | 0.26 | 0.52 | 3.96 | 4,552 | 37 min |
| M1 | HRRR, gust share 0 | 124 | 0.12 | 0.31 | 0.53 | 3.99 | 1,514 | 7 min |
| M1 | HRRR, gust share 0.5 | 158 | 0.28 | 0.35 | 0.55 | 3.79 | 2,795 | 10 min |
| M3 | HRRR, gust share 0 | 381 | 0.33 | 0.52 | 0.68 | 2.56 | 12,586 | 53 min |
| M3 | HRRR, gust share 0.5 | 922 | 0.32 | 0.35 | 0.69 | 2.31 | 19,166 | 131 min |

Jaccard is against the NASA FEDS satellite perimeters. The GeoMAC outline stamped 3.35 h is left out (A6).

- Station wind fails: the fire runs south-southwest and far too slowly. HRRR wind turns it west through Paradise.
  HRRR against the canyon station (Jarbo Gap): speed ratio 0.91, direction +1 degree.
- Spread only stops at the east edge of Paradise: the town is non-burnable fuel and houses light only by flame.
  Embers and house-to-house spread carry it through. Share of house ignitions with gust share 0: flame 68 %,
  neighbour 30 %, ember 2 %.
- Gust share 0 matches area and overlap at 21 h but is slow at 9 h, sits too far south, and misses the north-west
  lobe toward Chico. Gust share 0.5 is right at 3 h and then burns 2.5 times the mapped area.
- A6's rule for the gust share (lower bits per building) picks 0.5, by 2.31 to 2.56. That rule was badly made: with
  one seed every probability is 0 or 1, so the score rewards burning more houses whatever else burns. It is recorded
  here as picked, and it is not carried forward: amendment A7 replaces it before any further run.
- To inspect: straight-edged wedge shapes in the south-west of the gust 0.5 run (spot fires spreading as separate
  patches; the 600 m wind lattice or the solver's band may set the edges).
