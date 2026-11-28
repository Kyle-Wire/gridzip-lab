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

### A14: the compressor on the Camp and Tubbs errors (24 Nov 2026; calibration fires; E0 runs, seed 0)

`python physics/fire_real.py Camp=<table> Tubbs=<table>`, result `results/physics/fire_real.json`; tables from
FireSim `results/error_tables/` (13502ca). 26 fire formulas in `physlib/terms.py` (list `FIRE`); 23 had inputs on
both fires (no weather-model terrain height in either table; house terms only in T3).

| table | readings (Camp, Tubbs) | each fire's level | spread before, after the winning formula | winning formula on the two fires pooled | bits saved per reading (largest of 19 shuffled nulls) | transfers both ways |
|---|---|---|---|---|---|---|
| T1 arrival, ln(observed / modelled) | 3,252, 1,360 | Camp -0.37, Tubbs +1.34 | 1.00, 0.41 | wind along the run, clock (bookkeeping), grass share, canopy cover, slope | 1.07 (0.07) | no |
| T2 reach (+1 model only, -1 record only) | 6,391, 2,802 | Camp +0.03, Tubbs +0.37 | 0.69, 0.50 | ln wind, wind aloft over surface, distance, shrub share, slope | 0.37 (0.06) | no |
| T3 houses, observed minus modelled share | 1,202, 962 blocks | Camp -0.02, Tubbs -0.26 | 0.41, 0.31 | contrast, timber share, canopy cover, night share, wind aloft | 0.34 (0.08) | one way only |

- **No finding under A14's rules.** No pooled formula transfers both ways: fitted on one fire, it codes the
  other worse than a plain offset (T1 and T2), or helps in one direction only (T3). And the noise-only twin failed
  on all three tables: on shuffled errors that keep their 3 km structure, the search still names two to five terms.
  So with two fires and errors this smooth in space, a named term is not yet evidence. The planted wind term was
  not recovered on T1 and the planted calm term not on T3; lee descent, embers and house gaps were recovered
  wherever they were planted.
- Levels (plain facts, not search results): on Camp the record shows fire at a point at 0.69 times the modelled
  time on average (the model is late, as the outlines said). On Tubbs the observed clock reads 3.8 times the
  modelled time, but that clock is an upper bound read on outlines from hour 35 on, so it overstates; Tubbs'
  overrun is the cleaner number: +0.37 on T2, and modelled house losses 26 points of share above the record.
- Leads, each from one fire or one table, driving nothing by themselves:
  1. Houses, both fires separately: `contrast` wins with a coefficient of -0.9 (Camp) and -0.6 (Tubbs), picked in
     30 of 30 bootstrap resamples. The model's block losses are too all-or-nothing: blocks it calls lost are lost
     less often, blocks it calls safe are lost more often. Part of this is the ensemble itself (2 or 3 seeds give
     shares of 0, 1/3, 2/3, 1). First step: at least 8 seeds, as A7 already asks, then see what contrast is left.
  2. Camp alone, arrival and reach: ln wind and wind-aloft-over-surface enter with nearly equal coefficients
     (3.7 and 3.5 in T1; 1.6 and 2.1 in T2), which together is a term in the 850 hPa wind: Camp's timing error
     follows the wind aloft, not the surface wind the model burned on. Spread on Camp falls from 1.00 to 0.28.
     This agrees with the transects and with A13 (surface analysis a quarter slow where the wind was strong). It
     does not hold on Tubbs.
  3. Tubbs alone, reach: overrun goes with calm hours (coefficient -2.8 on the calm share with the sign of fewer
     misses, more overrun in the long quiet tail) and lee ground; consistent with growth after the wind dropped,
     where the record has suppression and the model has none.
- Structure left after the winning formulas is large on all three tables (a point's leftover error is 54 to 70 %
  predictable from its neighbours in input space): the errors are smooth patches, not point noise.
- What this changes: the compressor cannot separate causes on two fires. The next calibration fires (Marshall,
  East Troublesome) matter more for it than more terms. Until then the levels and the two-fire tests of each
  change stay the deciding evidence.

### A13's speed scale in the fire runs, and A15's first station score (24 Nov 2026 evening; calibration)

Fire runs with the surface wind times 1.30 (A13's rule), otherwise the E0 runs (priors, road rule on, embers on
the 20-ft wind). FireSim results `camp_M3_24h_hrrr_g0_ws130_road0_E0_ens3`, `tubbs_M3_44h_g0_ws130_road0_E0`.

| run | seeds | burned km2 (mapped) | overlap at last outline | Camp at 3.35 h, model / mapped km2 | houses modelled / recorded | bits per building (own rate) | AUC | A8 score |
|---|---|---|---|---|---|---|---|---|
| camp, wind x1.00 | 3 | 479, 513, 508 (376) | 0.51, 0.51, 0.48 | 56 / 222 | 14,172 / 16,493 | 1.62 (0.99) | 0.753 | 2.81 |
| camp, wind x1.30 | 3 | 821, 805, 852 (376) | 0.42, 0.43, 0.40 | 182 / 222 | 17,888 / 16,493 | 1.56 (0.99) | 0.737 | 3.48 |
| tubbs, wind x1.00 | 2 | 238, 237 (125) | 0.44, 0.44 | | 11,619 / 5,185 | 1.62 (0.63) | 0.753 | 4.08 |
| tubbs, wind x1.30 | 2 | 296, 314 (125) | 0.36, 0.34 | | 15,730 / 5,185 | 2.55 (0.63) | 0.684 | 6.01 |

- The station-true wind level makes both fires score worse. It fixes what it should (Camp's first three hours:
  182 km2 modelled against 222 mapped, up from 56) and then feeds the error that was already the largest: fire
  that keeps running where the real ones stopped. Camp ends at 2.1 times the mapped area, Tubbs at 2.4 times.
- So the wind level was not the main fault. With the wind nearer the truth at exposed places, what is left is
  (a) the same multiplier also speeding up sheltered ground, where A13 showed the analysis is already 2 to 8
  times too fast, and (b) no end to spread: no suppression, no burnout, in the model by design (section 7).
- A13's rule stands as written (the scale is the station-based level); the numbers above are the standing with
  it. Nothing passes.

WindSim terrain wind W3 (A15), three bearings (30, 45, 60 degrees), each 2.30 million cells, worst
non-orthogonality 45, no divergence, no limiter; station table FireSim `results/wind_check/camp_2018_w3.txt`.

| class | station-hours (stations) | W3 RMSE, median ratio | W3 from 80 m over a clearing | HRRR 10 m |
|---|---|---|---|---|
| all | 28 (2) | 8.2 m/s, 2.09 | 6.4 m/s, 1.63 | 4.3 m/s, 0.92 |
| windward | 4 (1) | 7.6, 3.51 | 7.9, 3.60 | 1.5, 0.94 |
| lee | 24 (1) | 8.3, 1.70 | 6.1, 1.32 | 4.6, 0.87 |

- By A15's rule: only 2 stations lie inside the terrain runs (Jarbo Gap and Cohasset), so the check is too thin
  to decide and W3 is not used in any fire run.
- What the two stations show anyway: W3 is about twice too fast. Its level comes from the HRRR 80 m wind over
  the high ground upwind, carried down by a neutral log law; real air that night was stable and did not mix that
  momentum down everywhere. The prediction in A15 (good on ridges, wrong in the lee) was too kind: the level is
  wrong at both stations. The pattern (which ground is exposed, which sheltered) is untested, because two
  stations cannot test a pattern.


## A16 part 2: the overrun, searched with flank, backing, late-arrival and burn-duration terms (24 Nov 2026)

Run: `python physics/fire_real.py Camp=... Tubbs=... --a16` on the A14 tables (E0 runs, plain HRRR wind), terms and
table fixed in A16 and its note before the run. Output `results/physics/fire_real_a16.json`.

**No finding passes.** On both tables no formula transfers from one fire to the other, and the noise-only twin
names terms (5 on each table), so the search still over-reads these two fires. What follows are leads.

| table | readings | level (Camp, Tubbs) | pooled winner | bits saved per reading | best null | transfers |
|---|---|---|---|---|---|---|
| T2 reach | 9,193 | +0.03, +0.37 | wind_speed, wind_aloft, fuel_shrub, fuel_developed, flank_spread | 0.36 | 0.06 | no |
| T2o (points the model burned) | 7,519 | 0.34, 0.48 | wind_aloft, calm_share, night_share, fuel_shrub, burn_duration | 0.46 | 0.08 | no |

- The level itself is the main fact: 34 % of the ground the model burned on Camp and 48 % on Tubbs is outside
  the record.
- Late spread (predicted): the late group is named in every one of the 30 block-bootstrap draws on T2o, and its
  planted check passes. `burn_duration` is -0.019 per hour: ground the model reached in the last hours is about
  0.4 more likely to be overrun than ground it reached at the start of a 20-hour window. Tubbs alone names
  `late_arrival` (+0.84). This is the one overrun lead both fires share, and it says the modelled edge keeps
  running late in the window when the real one had slowed or was held.
- Flanks (predicted, not borne out): `flank_spread` is named pooled on T2 (bootstrap share 0.47, under the 0.5
  rule) with a negative sign, and on Tubbs alone with a positive one (+1.95). Camp alone names `backing_spread`
  with opposite signs on the two tables. No consistent flank or backing story.
- Camp alone: canopy cover on the line is the strongest single term (+3.5 on T2, +4.8 on T2o, 0.9 bits saved
  per reading). The Camp overrun sits in forest, which fits a fire that in the record was a wind-driven run
  through the towns and not a spread through timber on all sides.
- What it means for the model: a term that pays is a lead, not a correction. The lead to carry to Marshall and
  East Troublesome under its own amendment is a spread that weakens with hours since the wind event (night
  recovery, crews), and for Camp the timber and canopy spread rate under wind.

## A17: the free-burning window (25 Nov 2026)

Kyle's question: "can we test fires up until intervention? would that help?" Rule and score fixed in A17 before the
cut was computed. Numbers are from each run's `result.json` and the A14 tables; chart
`/mnt/project-files/wildfire/11_free_burning_window.png`.

The cut (end of the wind event at A13's windy stations): **Tubbs hour 13** (peak 14.7 m/s, under half from hour
13 on). **Camp: no cut inside the 24 hours** (peak 9.9 m/s, still 4.9 at hour 23). The prediction of 8 to 12 hours
was right for Tubbs and wrong for Camp.

| Camp, hours | mapped km² | model km², HRRR wind | overlap | model km², wind x1.30 | overlap |
|---|---|---|---|---|---|
| 3.35 | 222 | 56 | 0.10 | 182 | 0.25 |
| 9.45 | 331 | 208 | 0.30 | 379 | 0.36 |
| 21.45 | 376 | 442 | 0.51 | 778 | 0.42 |

| Tubbs, hours | model km² (HRRR wind) | share outside the final outline | detections on the modelled fire |
|---|---|---|---|
| 5.5 | 128 | 0.18 | 0.96 |
| 15.2 (first pass after the cut) | 171 | 0.31 | 0.87 |
| 43.0 | 234 (mapped 125) | 0.48 | |

- It helps, and it changes the reading of A13's wind scale. Up to 9.45 h the station-corrected wind is the better
  Camp run at both mapped times (overlap 0.10 to 0.25, and 0.30 to 0.36; area 0.82 and 1.15 of mapped against 0.25
  and 0.63). It only loses at 21 h, where the real fire grew 14 % in twelve hours and the model doubled.
- Tubbs at 5.5 h is already the size of the final fire with 18 % of it outside the final outline; every later
  hour adds mostly ground the real fire never took. No outline exists before hour 35, so there is no early overlap
  number for Tubbs, only this bound.
- Camp's slowdown is not the weather: the wind at the windy stations was still two thirds of its peak when the
  mapped fire all but stopped growing. That leaves crews, burnout and the ground it had reached (valley farmland,
  the lake, earlier burn scars) as the causes, none of which the model has. The weather cut cannot stand in for
  an intervention time on Camp; the incident status reports are the next source to open.
- Limits: three mapped times on Camp and none early on Tubbs; the 3.35 h Camp outline has a doubtful time stamp
  (noted in A5); houses cannot be scored by time. No gate changes.

## A16 part 1: the wide terrain wind and the contrast form W3c against stations (25 Nov 2026)

Three terrain-only WindSim runs (30, 45, 60 degrees; 144 km turned domain, 600 m cells, 2,073,600 cells each,
worst non-orthogonality 38, about 42 minutes each on 8 cores; residuals for p near 2e-3 at the end, so not fully
converged). Table: FireSim `results/wind_check/camp_2018_w3c_wide.txt` (`6d69a04`). 21 stations inside, 17 with
hours of 3 m/s or more, 158 station-hours in the 12 of 24 hours when the reference bearing was in range.

| | speed RMSE (m/s) | median ratio | bias | direction error |
|---|---|---|---|---|
| HRRR 10 m | 3.72 | 0.87 | -1.07 | 23° |
| W3c (HRRR level, terrain pattern) | 4.81 | 1.03 | +0.64 | 24° |
| W3 (terrain run's own level) | 3.89 | 1.15 | +0.48 | 58° |
| W3 over a clearing (decides nothing) | 3.23 | 1.08 | +0.44 | 58° |

Windward (13 stations): W3c 4.43, HRRR 3.19. Lee (5): W3c 5.81, HRRR 4.99.

**By A16's rule W3c does not earn a fire test, and by A15's rule neither does W3.** The prediction (W3c better at
sheltered stations) was wrong.
- W3c takes out HRRR's low bias (ratio 0.87 to 1.03) but adds scatter: at ridge stations it lifts an HRRR wind
  that was already about right (CDEC1 observed 3.9, HRRR 9.1, W3c 16.7; CBXC1 9.7, 9.1, 14.5; JBGC1 6.8, 7.7,
  12.3). At 600 m cells with a 3 km box the pattern is mostly "ridges faster", and HRRR's 3 km terrain already
  has the big ridges, so the speed-up is counted twice.
- With enough stations, W3's level is no longer twice too fast (1.15, against 2.09 on the two stations of the
  first check), but its direction is far worse than HRRR's (58° against 23°): one inflow bearing for a 144 km
  domain cannot follow the turning of the real flow.
- The clearing table has the lowest RMSE of all (3.23 against HRRR's 3.72, 13 % lower) but the same poor
  direction. It was fixed beforehand as deciding nothing, so it is a lead: the terrain run's 80 m wind brought
  down over open-ground roughness, with HRRR's direction kept, is the form to pre-register next, on Tubbs's
  stations as well as Camp's so that it is not tuned to one fire.
- Code fault found by the run: `box_mean` read the column spacing from distinct x values, which jitter by under
  a metre on a warped mesh; fixed in FireSim `2eef747` (nearest-neighbour spacing), tests pass.

## Correction, 25 Nov 2026 (A19): Camp's NASA FEDS outline times were read wrongly

FEDS outline times are half-day labels (night passes, afternoon passes), not clock times. The two Camp outlines
scored at 9.45 h and 21.45 h were really mapped at about 18.5 h and about 30 h. What this touches above:
- Every Camp A8 score, and A12's and A13's verdicts on Camp, compared the model at 21.45 h with a fire mapped
  about 8 hours later. The mapped area at the model's time was smaller than 376 km², so the Camp overrun is
  larger than stated everywhere above, and the x1.30 run's 778 km² is more than twice the real area, not twice.
- A17's Camp table: the rows at 3.35 h and 9.45 h against GeoMAC (222 and 285 km²) stand; the FEDS figures 331
  and 376 km² belong to about 18.5 h and 30 h. "The real fire grew 14 % in twelve hours" is still about right but
  for 18.5 h to 30 h. At 9.45 h the x1.30 run (379 km²) is 1.33 times GeoMAC's 285, not 1.15.
- A14, A16 and A18 tables: Camp's "record" in T2 and T2o is the outline of about 30 h, so Camp's overrun level
  (0.34 in T2o) is understated, and Camp's T1 arrival times at the FEDS outlines are about 9 hours early.
  Tubbs is not affected.
The Camp runs are being rescored with satellite pass times; the corrected numbers will be set beside the old ones
here, and the Camp error table rebuilt before the compressor reads Camp again.

## A18: the slow-edge stopping term on Camp and Tubbs (25 Nov 2026)

Run `fire_real.py --a16 --a18`, table T2o, output `results/physics/fire_real_a18.json`; chart
`/mnt/project-files/wildfire/12_overrun_by_speed_and_time.png`. Camp's table still carries the outline-time fault
of A19, so Camp's part is provisional until its table is rebuilt.

**No finding.** `slow_edge` is found when planted, and is named in the pooled winner (+0.24 for edges slower than
about 2.9 m/min), but it is chosen in only 13 % of bootstrap draws, does not transfer either way, and neither fire
alone names it. The prediction that it would take over from `burn_duration` was wrong: the late group is still
named in every draw.

Plain counts behind it (share of model-burned points the record does not show burned, by the model's own spread
rate there):

| m/min | under 1 | 1 to 2 | 2 to 5 | 5 to 10 | 10 to 20 | over 20 |
|---|---|---|---|---|---|---|
| Tubbs | 0.90 | 0.83 | 0.59 | 0.40 | 0.20 | 0.20 |
| Camp | 0.25 | 0.37 | 0.48 | 0.44 | 0.20 | 0.20 |

- Fast ground (over 10 m/min) is right four times in five on both fires: the wind-driven run is the part the
  model does well.
- On Tubbs, nearly everything the model burned slowly is ground the real fire never took. On Camp the slow ground
  is no worse than the middle. So "slow edges go out" would fix Tubbs and do little for Camp as tabled; it is a
  one-fire lead until Camp's table is corrected and Marshall and East Troublesome are in.
- By time both fires agree: the later the model reaches ground, the more likely the real fire never did (Tubbs
  from 23 % in the first sixth of the window to over 90 % after the first third; Camp from 18 % to 54 %).

## A19 Camp rescore and A18 Marshall: three fires now show the same fault (25 Nov 2026)

Runs on THE-COW (FireSim `18cc962` Camp, `8e137ea` and `77211e4` Marshall; E0 priors, 3 seeds unless said). Chart
`/mnt/project-files/wildfire/13_camp_marshall_area_over_time.png`. Old Camp `result.json` files are kept;
`result_a19.json` sits beside each.

Camp with corrected outline times (model km², 3-seed mean; x1.30 is seed 0 only):

| hours | mapped km² | E0, HRRR wind | ratio | overlap | E0, wind x1.30 | ratio |
|---|---|---|---|---|---|---|
| 3.35 (GeoMAC, doubtful time) | 222 | 76 | 0.34 | 0.10 | 182 | 0.82 |
| 6.95 (FEDS, was dropped) | 166 | 179 | 1.08 | 0.28 | 310 | 1.87 |
| 9.45 (GeoMAC) | 285 | 214 | 0.75 | 0.24 | 379 | 1.33 |
| 20.08 (FEDS) | 331 | 440 | 1.33 | 0.44 | 747 | 2.26 |

- Camp A8 at the corrected fit outline: E0 3.02 (was 2.81), x1.30 4.02 (was 3.48).
- **A17's reading of Camp is withdrawn in part.** "Early, the x1.30 wind is the better Camp run" rested on the
  3.35 h GeoMAC outline. The FEDS outline at 6.95 h is smaller than that outline (166 against 222 km²), which a
  growing fire cannot be, so the 3.35 h time stamp is wrong and that row is no longer used. At 6.95 h the plain
  HRRR run has the area right (1.08) and x1.30 is nearly twice too large. On Camp the wind scale does not help
  at any trustworthy time.

Marshall (ran unheld; ignition point from a landmark description, A18 note; outlines at 2.95 h and 15.27 h):

| hours | mapped km² | E0, HRRR wind | ratio | overlap | wind x1.30 | ratio | overlap |
|---|---|---|---|---|---|---|---|
| 2.95 | 19.9 | 13.5 | 0.68 | 0.46 | 18.7 | 0.94 | 0.41 |
| end (14 h) against 15.27 h | 25.3 | 59.5 | 2.35 | 0.27 | 78.8 | 3.1 | 0.21 |

- Both A18 predictions hold: x1.30 puts the area during the run inside 0.5 to 1.5 (0.94), and the model keeps
  growing after the wind eases (from 8.4 h to the end, +97 % and +92 %, against +27 % mapped over twelve hours).
- Overlap at 2.95 h (0.41 to 0.48) is at the published Marshall range (0.42 to 0.51) with no tuning, on area
  alone and one early outline; it is not a pass of anything.
- The model's growth speeds up between hours 8 and 11, when the observed wind is falling. That needs a look
  (wind turning the flank into unburned grass, or the HRRR wind not easing) before any stopping rule is written.
- HRRR wind at Marshall's 14 windy stations: median ratio 0.91, so A13's x1.30 is not supported there; foothill
  stations are 2 to 3 times over. A13's scale was set on two California fires and does not carry to Colorado.
- Houses (USGS CC0 points, stored, not scored under A18): 1,196 destroyed in the record; model 2,570 (E0) and
  4,209 (x1.30).

Across Camp, Tubbs and Marshall: while the wind drives the fire the modelled area is within about a third of the
mapped one; afterwards the model keeps growing (Camp 1.33 at 20 h, Tubbs about 1.9 at 43 h, Marshall 2.35 at
14 h) while the real fires all but stop, with 0 to 5 % official containment on Camp and Tubbs and none possible
on Marshall. The missing physics is what ends a run, and on Marshall crews cannot be the explanation.

## A18 East Troublesome: the model misses the run (25 Nov 2026)

Runs FireSim `f2a7f86` (E0) and `920d039` (x1.30), started from the 21 Oct 21:10 UTC outline (123.0 km²), 43 h,
3 seeds; chart `/mnt/project-files/wildfire/14_east_troublesome_area_over_time.png`.

| hours | mapped km² | E0 | ratio | overlap | x1.30 | ratio | overlap |
|---|---|---|---|---|---|---|---|
| 12.3 | 578 | 166 | 0.29 | 0.27 | 234 | 0.40 | 0.37 |
| 23.7 | 724 | 221 | 0.30 | 0.26 | 367 | 0.51 | 0.40 |
| 36.8 | 741 | 269 | 0.36 | 0.28 | 499 | 0.67 | 0.44 |

- A18's prediction (area ratio 0.5 to 1.5 during the run with x1.30) fails: 0.40 at 12.3 h. A8 3.49 and 3.69.
- The real fire took 455 km² in its first twelve hours, overnight, while the windy stations read 5 to 8 m/s. The
  model adds 43 km² (E0) or 109 km² (x1.30) in those hours and does its growing at hours 19 to 27, under the
  surface wind peak. The run was not in the 10 m wind the model is given.
- By A17's rule the wind does not ease until hour 40, so all three outlines are inside the free-burning window
  and there is no outline after it: this fire says nothing about stopping.
- Houses (stored, not scored): 618 destroyed in the record; model 74 (E0) and 596 (x1.30), with no skill in
  which ones (AUC 0.50 and 0.56).
- Candidates for what the model lacks, each a term to search before anything is changed: crown fire in
  beetle-killed lodgepole (LANDFIRE 2016 fuels predate much of the kill and carry no dead-canopy load), long
  spotting under the wind aloft, and the wind at ridge height. The model is the opposite of wrong here compared
  with the other three fires: in grass and oak under a wind event it runs right and then fails to stop; in
  timber it fails to run.

## A20: the search on four fires, leave one fire out (25 Nov 2026)

Outputs `results/physics/fire_real_a20_T2.json` and `fire_real_a20_T2o.json`. Tables: Camp (corrected, 6,145
points), Tubbs (2,802), Marshall (230), East Troublesome (7,400).

**No finding.** On both tables the pooled formula fails to beat the one-correction baseline on every one of the
four held-out fires (0 of 4; three were needed), and the noise-only twin again names terms (five on T2, four on
T2o).

| table | level per fire (Camp, Tubbs, Marshall, East Troublesome) | pooled winner | bits saved | best null |
|---|---|---|---|---|
| T2 reach | +0.07, +0.37, -0.37, -0.71 | wind_run, distance, fuel_timber, canopy_cover, slope_any | 0.29 | 0.05 |
| T2o | 0.41, 0.48, 0.11, 0.34 | wind_speed, wind_aloft, canopy_cover, humidity, burn_duration | 0.30 | 0.04 |

- The levels carry most of what is known: Tubbs overruns, East Troublesome misses (71 % of its reach points are
  ground only the real fire took), and Camp's corrected overrun share is 0.41 (it was 0.34 on the wrong times).
  Marshall's low 0.11 on T2o is against its final outline at 15 h, inside which most of the model's early burn
  lies; its overrun is in area, later.
- Leads that survive resampling (block-bootstrap share): on T2o the late group 1.0 (third search running),
  moisture 0.93 (humid hours go with overrun: the model burns on through the night recovery), canopy cover 0.83;
  on T2 distance 1.0, canopy cover 0.93 (canopy on the line goes with overrun, +3.1), fuel 0.77 (timber on the
  line goes with misses, -1.5). The canopy and timber signs pull opposite ways, which is the Camp and East
  Troublesome contrast in one line: the model burns too far through sheltered canopy at a creep and not far
  enough through timber that crowned.
- `slow_edge` and the flank terms are not in either winner.
- What the search itself shows: with four fires the unit of evidence is the fire, and four is still too few for
  a formula to carry from three to the fourth; and a noise twin that names terms on every table since A14 means
  the term cost is too low for errors this patchy. That is a fault in how the search is set for fire, to be
  fixed (a cost raised until the twin is silent, set on the twin alone) before the next reading, under its own
  amendment.

## A21 and A22 (term): a search that is silent on noise, and the held-edge term (25 Nov 2026)

Outputs `results/physics/fire_real_a21.json` and `fire_real_a22.json`.

A21, thinning chosen on noise alone (five noise-only twins per step; count that named any term, T2 / T2o):
300 m 5 / 5, 600 m 5 / 4, 900 m 4 / 3, **1.5 km 0 / 0**. From here the fire search runs on the 1.5 km lattice
(659 reach points and 351 model-burned points over four fires). The prediction (step 5 or 10) held. At this
setting the noise twin passes for the first time since A14, and every planted term but the flank term on T2o is
still found.

The real search at that setting: **still no finding** (T2: 0 of 4 folds transfer; T2o: 1 of 4, Marshall, which
has 5 points). What is left standing after the noise is taken out, by block-bootstrap share:
- T2: distance 0.97, canopy cover 0.83 (+3.6: canopy on the line goes with overrun), fuel 0.60 (timber -1.6:
  goes with misses).
- T2o: canopy cover 0.93 (+4.9), fuel 0.77 (timber -2.4), moisture 0.70 (humid hours with overrun), late 0.57.
- Gone at this setting: the wind terms, flank, slow_edge. Those were patch fits.

A22's term `held_edge(L)` (model flame length at arrival under 1.22 or 2.44 m): found when planted, **not named**
on the real errors (the prediction was wrong). The plain counts say why (share of model-burned points outside
the record, by the model's flame length there):

| | under 1.22 m | 1.22 to 2.44 m | over 2.44 m |
|---|---|---|---|
| Tubbs | 0.62 | 0.56 | 0.31 |
| Camp | 0.26 | 0.28 | 0.52 |
| East Troublesome | 0.25 | 0.56 | 0.41 |
| Marshall (130 points) | 0.12 | 0.00 | 0.11 |

Tubbs's overrun is weak edge, as a held-edge rule assumes. Camp's is the opposite: its wrong ground is where the
model has tall flames, which with the canopy lead means the model drives a strong fire through forest the real
fire did not take. One stopping rule will not serve both. The H2 model runs (A22) are in progress and will be
read by their fixed rule; this term result says to expect a gain on Tubbs and none on Camp.

## A23 (term): the dry-crown term (25 Nov 2026)

Output `results/physics/fire_real_a23.json` (T2, 1.5 km lattice, four fires, tables `*_E0_a23.npz`).

`dry_crown(m)` (a cell whose model surface intensity would start a crown fire at foliar moisture m, 10 to 100 %,
but does not at the model's fixed 100 %): found when planted, **not named** on the real errors. The winner is
unchanged from A21 (wind run, distance, timber -1.6, canopy cover +3.6) and 0 of 4 folds transfer. The
prediction (named, with a negative sign, on East Troublesome) was wrong.

The plain counts say why. Share of reach points the model missed on East Troublesome, by crown class:

| no canopy | crowns at 100 % | crowns only if drier | never crowns |
|---|---|---|---|
| 0.65 | 0.80 | 0.80 | 0.83 |

The miss is about 0.8 in every class, cells the model already crowns included. So the timber miss is not a
crown-initiation threshold, and lowering foliar moisture in the model would not close it. By A23's own wording,
the next terms are the wind aloft over the fire and spotting distance. On Camp the cells that crown at 100 %
are overrun 0.48 of the time (missed 0.13), the same picture as the flame-length table above: the model's
strong fire in forest goes too far there.

## A22 (model runs): the held-edge rule H2 is not kept (26 Nov 2026)

Runs in FireSim `results/*/..._E0_H2_122*` and `..._E0_H2_244*` (M3, E0 priors, HRRR as analysed, 3 seeds, Tubbs 2).

| | E0 | H2(1.22 m) | H2(2.44 m) |
|---|---|---|---|
| Area ratio at last mapped time: Camp (20.08 h) | 1.33 | 0.54 | 0.00 |
| Tubbs (42.98 h) | 1.86 | 1.03 | 0.69 |
| Marshall (end state) | 2.35 | 1.35 | 0.56 |
| East Troublesome (36.83 h) | 0.36 | 0.20 | 0.18 |
| First mapped time: Camp overlap (3.35 h) | 0.106 | 0.011 | 0 |
| Marshall overlap (2.95 h) | 0.461 | 0.385 | 0.375 |
| Tubbs outside-share (5.5 h) | 0.18 | 0.13 / 0.16 | 0.18 / 0.15 |
| A8: Camp | 3.02 | 4.35 | 5.92 |
| Tubbs | 4.08 | 2.75 | 2.79 |
| Marshall | 2.45 | 1.94 | 2.07 |
| East Troublesome | 3.49 | 3.66 | 3.53 |

By the rule fixed in A22, **neither setting is kept**.
- H2(1.22) fails three ways: Camp's ratio moves away from 1 (1.33 to 0.54), the first-time overlap drops by more
  than 0.02 on Marshall (0.076) and Camp (0.095), and Camp's A8 is worse by 1.33.
- H2(2.44) never lets Camp start: the head-fire flame at the ignition point is under 2.44 m in the first hour
  and the rule is on from minute 0, so all three seeds stay at one 30 m cell. Under 1.22 one Camp seed goes out
  at ignition and two escape by embers (263 and 333 km²).
- The prediction was wrong in size: H2(1.22) trimmed Tubbs and Marshall more than predicted (1.03 and 1.35, not
  above 1.4), and both settings cut Camp hard instead of doing little.

What it does show: on Tubbs, H2(1.22) lands the area at 123, 125 and 129 km² against 108, 122 and 125 mapped,
with overlap up at every outline and modelled houses destroyed 5,674 against 5,185 observed (E0: 11,619). On
Marshall the end state improves (overlap 0.27 to 0.33) but the first three hours get worse, when the real fire
was a wind-driven grass run that nobody held. This agrees with the term result above: Tubbs's overrun is weak
edge late in the fire; Camp's and Marshall's early runs are not. A rule that holds weak edges from the first
minute is the wrong shape. The lead for a new amendment is a rule that is off while the fire is in its
wind-driven run and comes on only afterwards (and never at the ignition), with "afterwards" defined from
forecast-available weather, not from the fire record.

## A25 part 1 (term): the Cruz crown fire rate (26 Nov 2026)

Output `results/physics/fire_real_a25.json` (T2 and T2o, 1.5 km lattice, four fires, A23 tables).

`cruz_crown` and `cruz_active` (Cruz, Alexander and Wakimoto 2005): found when planted on both tables, **not
named** on the real errors (block-bootstrap share 0.03 on T2, 0 on T2o). T2's winner is unchanged (wind run,
distance, timber -1.6, canopy cover +3.6; 0 of 4 folds transfer); T2o's is canopy cover alone (+1.4, share 0.93),
with fuel 0.70, moisture 0.67 and late arrival 0.53 behind it. The noise twin is silent on both.

The prediction was half right: not named on T2, as expected, because with window-mean wind and moisture East
Troublesome has the lowest Cruz rate of the four (21 m/min) and so the term cannot carry its miss. It was not
named on T2o through Camp either. What this says: a crown-rate formula fed the window-mean weather does not sort
the errors, so the test of the formula is the model run, where it sees each hour's wind and moisture (East
Troublesome's first 12 hours had 1-h moisture 3.9 to 6 %, where the Cruz rate in a timber cell at 5 m/s wind and
CBD 0.08 is 43 m/min against the model's 4 m/min). That run, setting C1, was fixed in A25 part 2 before this
result and is queued on THE-COW after InSAR and wind.

## A24 (term and model runs): the held edge after the wind event, H3 (26 Nov 2026)

Term (tables `*_E0_a24.npz`, which also carry `u700` and `u10_marr`): `held_after_event(L)` found when planted,
noise twin silent, **not named** on T2o (winner: timber, canopy cover, humidity, burn duration). The plain count
is sharper: on Tubbs, ground the model reaches after the wind event ends is wrong ground 94 to 100 % of the time
(22 points) against 34 % before the cut, whatever the flame length; on East Troublesome it is mostly right ground.

Model runs (E0 priors, HRRR as analysed, 3 seeds, Tubbs 2; East Troublesome H3(2.44) still to run):

| | E0 | H2(1.22) | H3(1.22) | H3(2.44) |
|---|---|---|---|---|
| Marshall overlap at 2.95 h | 0.461 | 0.385 | 0.461 | 0.461 |
| Marshall area ratio, end | 2.35 | 1.35 | 1.62 | 0.96 |
| Marshall overlap, end | 0.27 | 0.33 | 0.35 | 0.49 |
| Marshall A8 | 2.45 | 1.94 | 1.94 | 1.88 |
| Tubbs ratio at 42.98 h | 1.86 | 1.03 | 1.50 | 1.40 |
| Tubbs overlap at 42.98 h | 0.44 | 0.47 | 0.52 | 0.55 |
| Tubbs 5.5 h outside-share | 0.18 | 0.13 / 0.16 | 0.18 / 0.16 | 0.18 / 0.16 |
| Tubbs A8 | 4.08 | 2.75 | 3.39 | 3.28 |
| East Troublesome ratio at 36.8 h / A8 | 0.36 / 3.49 | 0.20 / 3.66 | 0.34 / 3.51 | to run |

By A24's fixed rule H3 is **not kept**: Camp's HRRR wind never falls to half its peak inside its window, so H3 is
E0 there and Camp's ratio cannot move toward 1, as A24 said before the runs. On the other two it passes every
other clause (Marshall and Tubbs ratios move toward 1, first-time overlap unchanged, A8 better on both).
Predictions: Tubbs right (1.4 to 1.7), Marshall right at 2.44 (0.96) and over at 1.22 (1.62 against 1.0 to 1.5).

What it shows: H3(2.44) is the best Marshall run so far on every score (end ratio 0.96, overlap 0.49, houses
1,750 against 1,196 observed with AUC 0.85), and it keeps the early grass run whole. Tubbs's overlap is the best
yet at every outline, but the model has 176 km² burned by the cut at hour 17, above every mapped outline, so most
of Tubbs's overrun is in the wind run itself and no stopping rule after it can remove it. Camp's overrun is the
same kind. The open fault on the overrun fires is now during the run (Camp, Tubbs), and on East Troublesome it is
the run that does not happen; A25's crown rate (C1, queued) is the next test of the second, and it will likely
make the first worse, which A25 already says.

A24 completed (27 Nov 2026). East Troublesome H3(2.44): ratio 0.29 / 0.30 / 0.34 at 12.3 / 23.7 / 36.8 h, overlap
0.272 / 0.256 / 0.266, A8 3.51 (1.82 + 1.69), houses 68 of 618 (AUC 0.50); identical to H3(1.22) to the tenth of a
km² until hour 28, then flat at 251 km². The reading above stands. The plain counts behind the term result:

| share of model-burned points the record does not show burned | before the cut | after, flame < 1.22 m | after, 1.22 to 2.44 m | after, > 2.44 m |
|---|---|---|---|---|
| Tubbs | 0.34 (74) | 0.94 (17) | 1.00 (4) | 1.00 (1) |
| East Troublesome | 0.38 (53) | 0.17 (6) | 0.00 (2) | 0.00 (1) |
| Camp / Marshall | 0.42 (188) / 0.00 (5) | none | none | none |

Before each cut every H3 run is identical to E0 (Marshall and East Troublesome hour by hour, Tubbs seed 0), so the
rule changes nothing it should not.

## A25 part 2 (model runs): the Cruz crown rate, setting C1 (28 Nov 2026)

Runs (FireSim fd54baa, 1dae85c, 219c0df, 4195173): E0 priors, HRRR as analysed, 3 seeds (Tubbs 2). East
Troublesome finished seed 0 only inside the 115-minute cap (A25 note 2), so its row is one seed, at all three
outlines.

| | E0 | C1 |
|---|---|---|
| East Troublesome ratio / overlap at 12.3 h | 0.29 / 0.272 | 1.19 / 0.411 |
| at 23.7 h | 0.30 / 0.256 | 1.47 / 0.395 |
| at 36.8 h | 0.36 / 0.277 | 1.70 / 0.401 |
| East Troublesome A8 / houses (618 observed) | 3.49 / 74 | 12.67 / 4,127 |
| Camp ratio at 20.08 h / A8 | 1.33 / 3.02 | 5.09 / 9.06 |
| Tubbs ratio at 42.98 h / A8 | 1.86 / 4.08 | 7.49 / 13.93 |
| Marshall ratio at end / A8 | 2.35 / 2.45 | 2.43 / 2.58 |

By A25's fixed rule C1 is **not kept**: East Troublesome's ratio moves toward 1 at 12.3 h and 23.7 h but away at
36.8 h (0.64 off to 0.70 off) and its A8 is far worse; Camp and Tubbs are far worse on every clause. Predictions:
East Troublesome's 12.3 h ratio was predicted 0.4 to 0.7 and came out 1.19 (too low a guess); Camp overrunning
more was right.

What it shows: the timber fire now runs. At 12.3 h C1 has the best East Troublesome shape yet (overlap 0.41
against 0.27), so the miss there was the crown rate, as A23's class table pointed to. But the same rate in every
canopy cell is far too strong elsewhere and later: 99.5 % of torching cell-hours were outside Cruz's fitted range
(1-h moisture under 7 %, canopy bulk density under 0.12, wind over 32 km/h), Tubbs's oak woodland is not the
conifer stands the formula was fitted on, and embers rose 40 to 60 times on East Troublesome and about 3 times on
Tubbs, which drives the house counts. Leads for a new amendment, compressor first: the crown rate only where the
stand is conifer (LANDFIRE existing vegetation type), held to the fitted range at its edges; and the ember source
tied to the crown fraction rather than to intensity alone.
