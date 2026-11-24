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
- Wedge shapes inspected (FireSim results/camp_2018/camp_M3_24h_hrrr_g05_s1/sw_corner_fuel.png): they are fields
  bounded by LANDFIRE road lines, one 30 m non-burnable cell wide. Surface spread cannot cross a cell, so each field
  burns only when an ember lands in it. Neither the wind lattice nor the band makes them. Amendment A8 adds a
  road-crossing rule to test.

### After the incremental band (24 Nov 2026, FireSim 12f61ec; calibration, uncalibrated priors, gust share 0)

| run | seeds | wall time | km2 at 21 h (mapped 376) | Jaccard at 21 h | bits per building (own rate 0.99) | AUC | A7 score |
|---|---|---|---|---|---|---|---|
| camp_M3_24h_hrrr_g0_s1_v2 | 1 | 25 min | 345 | 0.54 | 2.59 | 0.673 | 65.0 |
| camp_M3_24h_hrrr_g0_ens | 3 | 81 min | seed 0: 345 (whole run 375, 393, 463) | 0.54 (seed 0 only) | 1.80 | 0.735 | 61.3 |

- The band change moved the one-seed result a little (381 to 345 km2 at 21 h, Jaccard 0.52 to 0.54): spot fires
  are random, and one cell of difference changes every later draw.
- Three seeds already lower the bits (2.59 to 1.80) and lift the AUC, because house probabilities are no longer
  0 or 1. Still worse than the fire's own base rate (0.99 bits).
- The A7 score is nearly all area (59.5 of 61.3). A8 replaces it with two ratios.
- Still wrong: slow in the first 9 h, too far south, no north-west lobe toward Chico, about 3,250 standing
  houses near Oroville modelled as lost.

### Road rule R1 on Camp (24 Nov 2026, FireSim 25b1403, one seed each, gust share 0; calibration)

63,237 thin road cells in the Camp box. Mapped area: 285 (incident) and 331 (satellite) km2 at 9.45 h, 376 at 21.45 h.

| run | L_cross | km2 at 3.35 h | km2 at 9.45 h | km2 at 21.45 h | Jaccard at 21.45 h | bits per building | AUC | A8 score (2.0 = no skill) |
|---|---|---|---|---|---|---|---|---|
| ..._g0_s1_v2 | rule off | | 153 | 345 | 0.54 | 2.59 | 0.673 | |
| ..._road0_s1 | 0 m | 56 | 208 | 442 | 0.51 | 2.29 | 0.705 | 3.55 |
| ..._road1_s1 | 1 m | 67 | 206 | 453 | 0.49 | 2.42 | 0.690 | 3.74 |
| ..._road2_s1 | 2 m | 74 | 210 | 444 | 0.48 | 2.46 | 0.687 | 3.82 |

- Roads were part of the slow start but not most of it: area at 9.45 h rises from 153 to about 208 km2, still
  well short of the mapped 285 to 331.
- The rule adds overrun in the south by 21 h (about 445 against 376 km2 mapped), so overlap falls a little.
- L_cross 0, 1 and 2 m are indistinguishable at one seed (seeds alone moved the area from 375 to 463 km2).
- Every A8 score is above 2.0: no skill yet against the fire's own rates. This is with uncalibrated priors.

### Search stage 1 and wind aloft for embers (24 Nov 2026; Camp, first 10 h, 2 seeds, road rule on unless said; calibration)

Scored at the satellite outline at 9.45 h (331 km2 mapped; incident outline 285 km2). A8 score 2.0 = no skill.
Area ratio below 1 means the burned-area forecast beats the mapped burned share of the box.

| setting (FireSim results/camp_2018/search1/) | km2 at 10 h (two seeds) | Jaccard at 9.45 h | area ratio | AUC | A8 score |
|---|---|---|---|---|---|
| g0_off_d150 (road rule off) | 159, 159 | 0.31, 0.30 | 1.28 | 0.577 | 4.79 |
| g0_road_d150 | 216, 219 | 0.30, 0.29 | 1.43 | 0.581 | 4.87 |
| g0_road_d400 (ember median 400 m) | 403, 420 | 0.35, 0.34 | 1.61 | 0.607 | 4.60 |
| g025_road_d150 (gust share 0.25) | 421, 394 | 0.33, 0.33 | 1.63 | 0.632 | 4.44 |
| g0_road_d150_E80 (embers on 80 m wind) | 308, 326 | 0.35, 0.37 | 1.40 | 0.609 | 4.49 |
| g0_road_d150_E850 (embers on 850 hPa wind) | 254, 289 | 0.49, 0.46 | 0.90 | 0.675 | 3.41 |

- Stage 1 as planned did not fit its time cap (2 of 12 settings); the rest above were run after a fix that made
  runs 3.3 times faster with identical output (the fields were being rebuilt on every ember lookup).
- More embers, further throws or more gust all burn more area but in the wrong place: overlap stays near 0.3.
- Transects of the model's inputs (no outcome data) showed the HRRR 10 m wind dying over Paradise and reversing
  west of it, while the 850 hPa wind holds 10 to 19 m/s toward Paradise and Chico. The 80 m wind follows the surface.
- With embers on the 850 hPa wind (A10, written before the run) the overlap at 9.45 h rises to 0.46 to 0.49 with
  less area burned, and the burned-area forecast beats its no-skill rate for the first time. Buildings are still
  about 2.5 times worse than the fire's own rate at 10 h (the record is not timed; A9).
- One fire, and the idea came from looking at it. A10's check stands: it must also help Tubbs or it is dropped.

### Camp stage 2 and the Tubbs cross-check of wind aloft (24 Nov 2026; gust share 0, road rule on, ember median 150 m; calibration)

| run (FireSim results/) | seeds | burned km2 | mapped km2 | overlap at last outline | area ratio | bits per building (own rate) | AUC | A8 score |
|---|---|---|---|---|---|---|---|---|
| camp_2018/..._g0_ens (road off, E0, earlier) | 3 | 375, 393, 463 | 376 | 0.54 (seed 0) | | 1.80 (0.99) | 0.735 | |
| camp_2018/..._road0_E850_ens3 | 3 | 525, 569, 568 (490 at 21.45 h, seed 0) | 376 | 0.61, 0.51, 0.57 | 0.92 | 1.64 (0.99) | 0.736 | 2.57 |
| tubbs_2017/tubbs_M3_44h_g0_road0_E0 | 2 | 238, 237 | 125 | 0.44, 0.44 | 1.50 | 1.62 (0.63) | 0.753 | 4.08 |
| tubbs_2017/tubbs_M3_44h_g0_road0_E850 | 2 | 268, 338 | 125 | 0.39, 0.32 | 1.91 | 2.87 (0.63) | 0.657 | 6.47 |

- Camp: embers on the 850 hPa wind hold up over 24 h and 3 seeds. Best Camp result so far, still short of no
  skill on buildings (1.64 against 0.99 bits) and overburning by 21 h (490 against 376 km2).
- Tubbs: the same change makes every score worse. By A10's rule E850 is dropped.
- Tubbs overruns with either wind (1.9 times the area, 2.2 times the houses with E0): the model reaches the
  south-west tip of the outline within about 11 h and then keeps spreading for 30 h with no suppression.
- House-to-house spread is 44 to 51 % of modelled house losses on Tubbs, about 32 % on Camp.
- First Tubbs runs were first tests of the whole chain on a second fire: the fire runs the right way (south-west
  to Santa Rosa's north edge) with nothing tuned to Tubbs.

### A11: embers lofted by intense fire (option EP), and the Tubbs first night (24 Nov 2026; calibration)

| run | seeds | burned km2 | overlap at last outline | area ratio | bits per building (own rate) | AUC | A8 score |
|---|---|---|---|---|---|---|---|
| tubbs E0 (rerun, identical) | 2 | 238, 237 | 0.44, 0.44 | 1.50 | 1.62 (0.63) | 0.753 | 4.08 |
| tubbs EP, I_loft 10 MW/m | 2 | 276, 319 | 0.39, 0.34 | 1.91 | 2.97 | 0.635 | 6.63 |
| tubbs EP, I_loft 30 MW/m | 2 | 390, 260 | 0.27, 0.40 | 1.93 | 2.95 | 0.643 | 6.61 |
| camp EP, I_loft 10 MW/m | 3 | 442, 480, 497 | 0.57, 0.50, 0.58 | 0.94 | 1.58 (0.99) | 0.752 | 2.53 |

- EP fails A11's rule: it makes Tubbs much worse than E0 at both I_loft values. It is dropped, with E850 and E80.
  Embers stay on the 20-ft wind. Camp's missed run west of Paradise is carried as a stated limit of the HRRR
  surface wind until wind comes from WindSim (candidate M4). Camp EP at 30 MW/m was not run (nothing to decide).
- Tubbs first night, E0, seed 0: 96 % of 594 satellite detections at 5.5 h lie on or within 190 m of modelled
  fire, with 128 km2 modelled burned at that time; 87 to 100 % at the four later passes. The mapped area at 43 h
  is 125 km2. So the model's first night covers the detections, and its 110 km2 of later growth (to 234 km2) is
  the overrun. The share measures misses only: it cannot show overrun.
- Tubbs houses, E0: 11,619 modelled destroyed against 5,185 recorded.

### A12: house-to-house spread scaled by wind (rule H1) (24 Nov 2026; embers on 20-ft wind, road rule on; calibration)

| run | seeds | burned km2 | overlap at last outline | houses modelled / recorded | bits per building (own rate) | AUC | A8 score |
|---|---|---|---|---|---|---|---|
| tubbs, H1 off | 2 | 238, 237 | 0.44, 0.44 | 11,619 / 5,185 | 1.62 (0.63) | 0.753 | 4.08 |
| tubbs, H1 on | 2 | 260, 242 | 0.40, 0.43 | 7,105 / 5,185 | 1.23 | 0.731 | 3.43 |
| camp, H1 off | 3 | 479, 513, 508 | 0.51, 0.51, 0.48 | 14,172 / 16,493 | 1.63 (0.99) | 0.753 | 2.81 |
| camp, H1 on | 3 | 461, 483, 512 | 0.50, 0.47, 0.47 | 12,002 / 16,493 | 2.03 | 0.709 | 3.27 |

- As predicted in A12: Tubbs improves (A8 4.08 to 3.43, modelled losses 2.2 to 1.4 times the record), Camp
  worsens (2.81 to 3.27, AUC 0.753 to 0.709), because the HRRR wind over Paradise is weak while the town burned.
  By A12's rule H1 is not kept.
- Five changes tested on 24 Nov. Kept: none yet as a fitted choice; the road rule is carried as an option
  (small gain in early spread). Dropped by their own rules: E80, E850, EP, H1. Each wind-dependent change helped
  one fire and hurt the other. The common cause is the wind input (HRRR 3 km surface wind over Paradise), not
  the rules themselves; they may be retested once the wind is better, under a new amendment.
- Standing with no fitted numbers (priors, road rule on, E0): Camp A8 2.81, AUC 0.753; Tubbs A8 4.08, AUC 0.753.
  Both worse than the no-skill 2.0. Nothing passes.

### A13: the HRRR wind against every station around each fire (24 Nov 2026; weather only, no fire outcome read)

FireSim `scripts/wind_check.py`, tables in FireSim `results/wind_check/` (commit 13502ca). Stations inside each
fire's box plus 0.75 degrees, hour by hour over the run window, HRRR analysis at the nearest lattice point.

| fire | stations | all: median HRRR / station speed | windy stations (observed mean 5 m/s or more) | their median speed ratio | their median gust ratio | ratio against height above HRRR terrain |
|---|---|---|---|---|---|---|
| Camp, 24 h | 46 (11 airport, 35 RAWS) | 1.46 | 13 | 0.79 | 0.93 | r = -0.03 |
| Tubbs, 44 h | 38 (12 airport, 26 RAWS; one dead anemometer left out) | 0.92 | 7 | 0.75 | 0.86 | r = -0.12 |

- Speed, by A13's rule: the windy-station median is under 0.85 on both fires, so the surface wind is multiplied
  by 1 / mean(0.79, 0.75) = 1.30 in fire runs from here on. HRRR's 10 m wind is about a quarter too slow where the
  wind was really strong: Jarbo Gap 9.8 m/s observed against 7.5, Saddleback 11.2 against 5.9, Pike County
  8.2 against 4.7, Chalk 5.8 against 1.5.
- Gusts, by the rule: 0.93 and 0.86 are inside 0.85 to 1.15 on both fires, so the HRRR gust is taken as sound and
  the gust share stays a number to fit (A7).
- Terrain, by the rule: no link between the miss and station height above the model's terrain, so no height
  correction is designed.
- What the rule does not capture, reported as seen: at sheltered stations HRRR is far too fast (Camp: 16 stations
  with a ratio above 2, all with observed means under 3.5 m/s), so one multiplier makes those worse. The analysis
  wind is too smooth: too slow on exposed ridges and gaps, too fast in sheltered canyons and forest. That is a
  terrain-and-roughness pattern, which is what the WindSim terrain run (A15) is for. Direction errors are 10 to
  30 degrees at windy stations and 50 to 170 at calm ones.
- Tubbs OVYC1 reported 0.0 m/s every hour and is left out; Camp ORRC1 has no wind readings.

