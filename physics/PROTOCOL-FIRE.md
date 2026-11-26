# Protocol FIRE-1: scoring wildfire models at the house, real fires (written 23 Nov 2026, before any data was opened)

Ask (Kyle, 23 Nov 2026): research wildfire models, test them with the compressor, and build a frontier property-level
model for the house report. He agreed the scope the same day: run anywhere in the US, prove on fires that reached
homes, report (1) chance and direction of fire reaching the block, (2) ember exposure at the house, (3) flame length
and heat at the walls, (4) what the owner can change; route B (fast spread + ember and house ignition first, then
GPU fire physics coupled to WindSim with house-to-house spread).

**Seen before this was written:** catalogue pages for CAL FIRE DINS, WFIGS perimeters, LANDFIRE FBFM40, NIST Camp
Fire observations and USFS Wildfire Risk to Communities (what each holds and its licence line); the ELMFIRE software
page (EPL-2.0); the abstract of Qin 2025 (UMD); Szasdi-Bardales et al. 2025 (IJWF WF24027: WRF-Fire into SWUIFT on
the Marshall Fire; ignited structures −13 % in the best case, Jaccard 0.42-0.51, Sørensen 0.59-0.67); the abstract of
Young et al. 2025 (Fire Technology 61:7). No perimeter, damage point, fuel raster or model output has been opened,
plotted or scored. A download of the data listed below was started on Kyle's PC with the instruction to run no
analysis of outcomes.

## 1. Fires and the split (fixed now)

| set | fire | year | state | why |
|---|---|---|---|---|
| calibration | Camp | 2018 | CA | ember-driven, NIST timed observations |
| calibration | Tubbs | 2017 | CA | wind-driven run into a city grid |
| calibration | Marshall | 2021 | CO | grass fire into suburbs; the published score to beat |
| calibration | East Troublesome | 2020 | CO | timber crown run, sparse homes |
| blind | Almeda | 2020 | OR | urban corridor, house-to-house |
| blind | Lahaina | 2023 | HI | grass into a dense town, outside the West's fuels |
| blind | Eaton | 2025 | CA | ember storm into a street grid |
| blind | Palisades | 2025 | CA | chaparral into hillside homes |
| blind | Chimney Tops 2 | 2016 | TN | Southeast hardwood and pine |
| blind | Smokehouse Creek | 2024 | TX | plains grass, few structures |

A fire is dropped only if it has no licence-clean damage record or no pre-fire fuels; a dropped fire is not replaced
and no fire changes set. Each drop is listed in the results with its reason. The blind set needs at least four fires
with damage records, at least one of them outside California, or the gate cannot be scored and the result is
"not scored".

**Seal.** Blind-fire damage records and perimeters sit in the data folder but are not read by any code until the
candidate's code and parameters are committed and that commit is written to `physics/fire/FREEZE.json`. The scoring
script refuses a blind fire if the candidate's commit differs from the freeze. One freeze per candidate; a second
freeze of the same candidate is a new candidate and is reported as such, with the first result kept.

## 2. What a candidate is given, and what it may not see

Given: ignition point and time (from the incident record); fuels, canopy and terrain from the LANDFIRE version that
precedes the fire; building footprints; weather. Weather is a hindcast (station observations, and HRRR analysis where
used), and every result says so: this scores fire behaviour given the weather, not forecast skill.

Not given: the final or any progression perimeter, satellite fire detections, damage records, post-fire fuels or
imagery, suppression actions. A candidate that reads any of these is void.

## 3. Candidates

- B0 base rate: every building gets the calibration fires' destroyed fraction.
- B1 federal risk map: USFS Wildfire Risk to Communities "Risk to Potential Structures" at the building, mapped to a
  probability by a one-variable logistic fit on calibration fires.
- B2 distance to wildland fuel: metres from the building to the nearest burnable non-urban fuel cell, same mapping.
- M1 spread only: Rothermel surface spread with crown fire, level set, our own code (and ELMFIRE as a reference
  implementation where it runs); a building is "reached" if fire arrives at its cell.
- M2 = M1 + ember transport and spot ignition, + structure ignition from flame, radiation and embers.
- M3 = M2 + structure-to-structure spread.
- M4 = M3 with wind and fire coupled at street scale (GPU).
Each candidate outputs, per building, a probability of destruction and (M1-M4) an arrival time, from an ensemble
over stated input uncertainty. A candidate with no ensemble states a single probability map and is scored the same.

## 4. Scores (fixed now)

Population A: buildings with a damage record. "Destroyed" is damage class above 50 %; everything else is "survived".
Population B: A plus every footprint within 2 km outside the final perimeter with no record, counted as survived
(tests whether a model stops the fire). The perimeter is used here only to define who is scored, the same for every
candidate.

1. **Bits per building** (primary): −log2 of the probability the candidate gave to what happened, probabilities
   clipped to [0.005, 0.995], via `codelen`. Reported per fire and pooled, next to B0, B1, B2.
2. **AUC** of destroyed against survived, per fire, populations A and B.
3. **Area overlap**: Jaccard and Sørensen of modelled against final burned area at 30 m, at the time of the final
   perimeter, to compare with published scores.
4. **Arrival time**, where timed observations exist (NIST Camp Fire points; dated progression perimeters): bits per
   reading for log(observed / modelled hours since ignition) on a 0.05 grid under a Gaussian with calibration-fit sd,
   and median absolute error in hours.
5. **Count**: ignited structures, modelled over observed.

## 5. Gate for "passes" (blind fires only, fixed now)

A candidate passes if all hold on the blind set:
- G1 pooled bits per building (population B) below each of B0, B1 and B2, and below each of them on at least
  two-thirds of the scored blind fires;
- G2 median AUC (population A) at least 0.05 above B1's and at least 0.70;
- G3 median Jaccard at least 0.51 (the top of the published Marshall range);
- G4 structure count within a factor of 1.5 on at least two-thirds of scored blind fires;
- G5 calibration: of buildings given 20-80 % the observed destroyed fraction lies within 15 points of the mean
  stated probability, pooled.
Failing any is a fail, reported as a result. The gate is not changed after a blind result is seen.

## 6. The compressor's part

For each candidate on calibration fires (and, after its freeze, on blind fires as a separate table):
- y1 = log(observed / modelled arrival time) per timed observation; y2 = observed − predicted destroyed fraction per
  300 m block with at least 5 scored buildings. Grid 0.05 for both.
- Vocabulary, one term each (families pay for their parameter): wind speed at 10 m; wind-aligned slope; fuel
  model class; dead fuel moisture from station RH and temperature; canopy cover; `ember(L)`: fuel upwind within
  distance L (L fitted 50 m-5 km) weighted exp(−x/L); building density within 100 m; `gap(d)`: share of
  neighbours closer than d (d fitted 3-30 m); roof and siding class where the record has them; hours since ignition;
  night; distance to the nearest road wider than 15 m (firebreak); per-fire offset.
- Reading rules: MDL search up to 4 terms through `physlib.engine`; a term counts only if it holds in hold-out both
  ways across fires (fit on half the calibration fires, test on the other half, and the reverse). The neighbour test
  (99 shuffles) says how much unnamed structure is left. A term that pays names physics the candidate lacks; that
  finding may drive the next candidate, which then needs its own freeze before any blind fire.
- Planted checks first: on synthetic fires from the candidate itself with one mechanism switched off (embers; slope;
  house-to-house), the search must name the missing mechanism and name nothing on an unaltered twin. If it cannot,
  its findings on real fires are not read.

## 7. Limits stated now

Calibration fires are spent for testing. Damage records outside California differ in classes and completeness; each
fire's record is described beside its score. Suppression and structure defence are not modelled and not seen, so
survival near engines reads as model error. Footprints may postdate a fire; the release date is reported.

## Amendments

**A1 (23 Nov 2026, after the data inventory, before any data file was opened by scoring code).** Seen: the list of
what was downloaded for each fire and each source's licence text; no perimeter, damage point or raster was read,
plotted or scored.
- Damage records that are openly licensed exist only for the four California fires (CAL FIRE DINS: Camp 23,624
  points, Tubbs 6,023, Eaton 18,428, Palisades 12,137). Boulder County and FEMA (Marshall) state no licence, Maui
  County (Lahaina) is not open, and Almeda, East Troublesome, Chimney Tops 2 and Smokehouse Creek have none. As
  section 1 stands, the blind set has two fires with damage records, both in California, so gates G1, G2, G4 and G5
  read **"not scored"** until at least two more blind fires, one outside California, have open damage records. The
  rule is not relaxed. A search for open records (FEMA, NIST, DesignSafe and state sources) is added; any found for
  a fire already in section 1 is used in that fire's existing set.
- G3 (area overlap) and the arrival-time score need only perimeters and are scored on every fire that has pre-fire
  fuels. Chimney Tops 2 has none yet (LANDFIRE 1.4.0 is archived, on request) and is dropped from G3 until it does.
- Fuels: Almeda and East Troublesome use LANDFIRE 2016 (the 2019 version is no longer served; 2020 includes the
  fires). Terrain is LANDFIRE 2020 for all fires (terrain does not carry a burn scar).
- B1 (federal risk map) is built on LANDFIRE 2020 fuels, so for fires of 2016-2020 it has seen the burn scar. It is
  scored as a baseline only on fires of 2021 or later (Eaton, Palisades; Marshall, Lahaina and Smokehouse Creek for
  area), and that is stated beside every B1 number. B0 and B2 are unaffected.
- Footprints are Overture release 2026-11-18.0, after every fire, so buildings destroyed and not rebuilt may be
  missing. For population A the damage record's own point is the building; footprints are joined to it where one
  exists within 15 m, and the share without a footprint is reported per fire.
- Ignition times marked APPROXIMATE in the data (Camp, Tubbs, Chimney Tops 2, Lahaina, Marshall) are replaced from the
  incident record or the official investigation, with the source cited, before any run of that fire.
- Reference code: ELMFIRE changed licence on 18 Jun 2026 (AGPL-3.0 with Commons Clause). Only the tag `2025.0717`
  (EPL-2.0) may be run as a reference, and no ELMFIRE code enters FireSim either way.
- DINS points carry addresses, parcel numbers and assessed values. They are used for scoring only; no per-address
  output is published.

**A2 (23 Nov 2026, before the file below was downloaded or opened).** Seen: the USGS catalogue page for "Wildfire
building destruction and exposure in the conterminous United States, 2002-2022" (released 16 Dec 2025,
doi.org/10.5066/P1QX6UXD, CC0 1.0): building points classed burned or unburned from pre- and post-fire aerial imagery,
inside each perimeter and a 2.4 km buffer, for the 362 fires with at least 10 buildings destroyed. Whether it holds
the section 1 fires is not yet known.
- For a section 1 fire with no DINS record, this release is its damage record: "burned" is destroyed, "unburned" is
  survived, and its 2.4 km buffer stands in for population B's 2 km. Where both this release and DINS exist (Camp,
  Tubbs), DINS scores and the agreement between the two is reported.
- A1's "not scored" for G1, G2, G4, G5 lifts only if this gives at least two more blind fires, one outside
  California, as A1 says. Lahaina (Hawaii) and the 2024-2025 fires are outside the release.
- The release's other fires are not added to either set by this amendment. Using them (a wider blind test across
  many fires) needs its own amendment naming the split rule before their outcomes are read.

**A3 (23 Nov 2026, after the USGS release was downloaded; its points were not opened).** Seen: the release's file
list, attribute names, licence text (CC0-1.0, use constraints "None"), and its per-fire summary table of exposed and
destroyed building counts for all 362 fires, which includes the section 1 fires Marshall, Almeda, East Troublesome,
Chimney Tops 2, Camp and Tubbs. So the destroyed totals of those fires are known before any run; G4 (count) is
read with that said beside it. No building point, perimeter or raster has been opened.
- With A2, blind fires with a damage record are Almeda, Eaton, Palisades and Chimney Tops 2. A1's condition (four,
  one outside California) is met only when Chimney Tops 2 has pre-fire fuels (LANDFIRE 1.4.0, on request). Until
  then the house-level gates stay "not scored" on three fires, and the three are reported.
- `expos_type` "Direct" (inside the perimeter or within 100 m) is population A for these fires; "Direct" plus
  "Indirect" is population B.
- Ignition times now come from official records for Camp, Tubbs, Lahaina and Chimney Tops 2 (the last is the time
  the fire was reported). Marshall's is still approximate; no Marshall run until an official time is cited.

**A4 (23 Nov 2026, before any planted run).** The planted checks of section 6, made concrete for the mechanisms
FireSim has today (commit pinned in the result). Embers and house-to-house spread get the same check when they exist.
- Truth: FireSim on 12 made-up landscapes (`firesim.synth`, seeds 20-31, never used before), 6 m/s wind, 10 m cells,
  90 min. Candidates: the same with one mechanism off: (a) slope (slope set to zero), (b) crown fire, (c) canopy
  sheltering of the wind. Twin: nothing off.
- Readings: 150 random points per landscape that both runs burn after minute 5. y = ln(truth arrival / candidate
  arrival) + noise N(0, 0.05), grid 0.01. The noise stands in for timing error in real observations and gives the
  twin something to code.
- Vocabulary, each averaged along the straight line from ignition to the point (the fire's time is a sum along its
  path): `slope_run` = tan(slope) squared times the cosine between upslope and the line (the slope factor's form);
  `slope_any` = tan squared; `crown_fuel` = share of the line under canopy with base height 3 m or less;
  `canopy` = share under any canopy; `shelter` = canopy cover times crown ratio / 3 (the sheltering fill);
  `grass`, `shrub`, `timber` = share in each fuel group; `wind_run` = cosine between wind and the line;
  `ln_dist` = log distance. Slope terms are one group, canopy terms (`crown_fuel`, `canopy`, `shelter`) are scored
  separately because telling crown fire from sheltering is part of the test.
- Pass, per candidate: the winning formula contains the planted mechanism's term (`slope_run` or `slope_any` for
  (a); `crown_fuel` for (b); `shelter` for (c)), it holds in hold-out both ways (seeds 20-25 against 26-31), and no
  term of another mechanism's group is in the winner. Twin passes if the winner has no term. A fail is reported as a
  limit of what the compressor can tell apart from arrival times; it then says which real-fire findings can be read.

**A5 (23 Nov 2026, before any real-fire run).** How a real fire is set up for candidates M1 to M3, fixed for the
first calibration fire (Camp) and used for every fire after it unless a later amendment says otherwise. What had been
seen when this was written: the file layout, grids, field names and label counts of the Camp inputs (no model output,
no overlay of any input on the perimeter or the damage record); the synthetic and closed-case results in FireSim
`sandbox/REPORT.md`; one end-to-end run of the runner on a made-up fire.
- **Grid.** LANDFIRE 30 m fuels and canopy; solver on 10 m cells (each fuel cell split 3 x 3), because 30 m solver
  cells failed the synthetic landscapes. Slope and upslope direction are computed from the LANDFIRE elevation layer
  on its own grid. Wind bearings are turned from true north to grid north (about 15 degrees in northern California).
- **Ignition.** One point at the place and time of the official record, radius 15 m. Camp: the origin in NIST TN
  2135 Appendix F (39.81461, -121.434004), 2018-11-08T14:33Z (CAL FIRE). The NIST timeline is otherwise an outcome
  record: only its origin row is read by a candidate.
- **Window.** The first 24 hours. Perimeter overlap is scored at every mapped perimeter time inside the window
  (GeoMAC incident perimeters; NASA FEDS satellite perimeters reported apart, they are 375 m estimates), not against
  the final perimeter, because suppression is not modelled and the final perimeter is days later. Buildings are
  scored on the whole record, which for a fire that kept destroying buildings after 24 hours counts against the
  candidate; the share of the record inside the 24-hour modelled or mapped area is reported beside it.
- **Weather.** Hourly means of the ASOS and RAWS stations in the fire folder. Wind: inverse-distance (power 2) mean
  of the station wind vectors at each place, ASOS 10 m wind divided by 1.15 to 20 ft, RAWS as read; gusts not used.
  Dead fuel moisture, one value per hour for the whole fire, from the inverse-distance mean temperature and humidity
  at the ignition point: EMC by Bradshaw et al. 1983 (INT-169) Appendix C, 1-hour = 1.03 EMC, 10-hour = the RAWS
  fuel stick where read, else 1.28 EMC, 100-hour = 1-hour + 2 points (the spacing of the Scott and Burgan scenarios).
  No fuel-level sun correction (INT-169 Table 6).
- **Live moisture.** Not measured. Per fire, one of the Scott and Burgan live scenarios, stated in the fire's config
  with the reason, chosen from the season and recent rain before any run of that fire. Camp: L1 (herbaceous 30 %,
  woody 60 %; fully cured: November, after a rainless autumn). Foliar moisture 100 %.
- **Buildings.** Positions from the USGS CC0 building points (the pre-fire inventory; amendment A2) for fires in
  that release, each a 150 m2 footprint. This replaces the Overture footprints of A1 for those fires: Overture's
  release is from 2026, so buildings destroyed and not rebuilt are missing from it. The runner reads positions for
  the model and keeps the burned/unburned field apart for scoring. For California fires the DINS record is scored
  as well, matched to the nearest building point within 30 m. For a blind fire this file is opened only after the
  freeze, as section 1 requires.
- **Ensemble.** M1 has no random part beyond house ignition by flame. M2/M3: seeds 0 to N-1, N stated per result;
  a building's probability is the share of runs in which it lights. Wind speed and direction error enter the
  ensemble only under a later amendment (the synthetic landscapes showed burned area changing up to 74 times
  between 4 and 8 m/s, so this matters and is stated as a limit of the first runs).
- **Calibration.** The ember and house numbers in FireSim `embers.PRIORS` and `structures.PRIORS` are starting
  guesses. The first Camp runs use them as they stand and are reported as such. Fitting them is done on calibration
  fires only, by bits per building, and the fitted values and the runs tried are listed in the findings.
- **Runs.** FireSim `scripts/run_fire.py` with `fires/<fire>.json`, on THE-COW; every result pins the FireSim
  commit and refuses nothing silently (a dirty tree is recorded in the result).

**A6 (23 Nov 2026, after the first Camp runs, before any HRRR value was read).** What had been seen: the three
first Camp runs under A5 (FireSim f8a466b, `results/camp_2018/`), which failed. Spread only: 91 km2 against 376 mapped
at 21 h, overlap (Jaccard) 0.16, building AUC 0.51. With embers and houses, one seed, uncalibrated: 332 km2, Jaccard
0.26, AUC 0.52, 3.96 bits per building against 0.99 for the fire's own burned share. The modelled fire ran
south-southwest and was far too slow in the first 9 hours (10 km2 against about 300); the mapped fire ran west.
Reading: the wind input. The inverse-distance mean blends the canyon station (Jarbo Gap, wind from the north-north-
east) with valley stations (from the north-north-west), uses 10-minute or 2-minute means of 3 to 6 m/s, and no gusts.
- **Wind source W2.** NOAA HRRR analysis (hour f00, 3 km, US government work, public domain), 10 m wind components
  and surface gust, fetched from NOAA's open archive, turned from the HRRR grid to true north to the fire grid,
  interpolated bilinearly, divided by 1.15 to 20 ft. This replaces station wind W1 for candidates M1 to M3; W1
  results stay on record. Dead fuel moisture stays as in A5 (stations).
- **Gust share g.** Effective wind = mean + g x (gust - mean), direction unchanged. g is a calibrated number: 0 and
  0.5 are run on Camp, and the one with the lower bits per building (ties: higher Jaccard at the 21 h satellite
  perimeter) is carried to the other calibration fires, where it is checked, not refitted, before any freeze.
- **Check of the source, reported with every W2 result:** HRRR 10 m wind at the stations against the station
  readings over the window (mean speed ratio and mean direction difference per station). This is a description, not
  a gate.
- **Perimeters.** The GeoMAC outline stamped 2018-11-08 17:54Z carries 222 km2 at 3.35 h after ignition, which the
  9.45 h outlines (285 and 331 km2) and the incident record make implausible as a mapping time. It is reported but
  left out of any pooled perimeter score, here and for any outline whose stamped area exceeds the next later
  outline's area times 0.75 while being under 4 hours old. The rule is stated now and applied by time and area only.
- **Blind fires.** W2 is weather, allowed by section 2 (hindcast, stated). HRRR starts in 2014 (v1) and is used at
  whatever version the archive holds for the date; a fire with no HRRR hour falls back to W1 and says so.
- A terrain-following wind at finer scale than 3 km (WindSim, candidate M4's part) is not in W2.

**A7 (24 Nov 2026, after the six Camp runs in FINDINGS-FIRE.md, before any further real-fire run).** What had been
seen: those six runs. A6's rule for the gust share picked 0.5 on bits per building while that run burned 2.5 times
the mapped area; with one seed the rule only counts houses. It is replaced, and the pick is not used.
- Calibrated numbers (gust share, and the ember and house numbers in FireSim `PRIORS`) are fitted together on the
  calibration fires only, with ensembles of at least 8 seeds so that building probabilities are not 0 or 1.
- Fitting score, lower is better: bits per building (population B) plus bits for the burned area, the latter being
  the cost of coding each 30 m cell inside the box as burned or not at the latest satellite or incident perimeter
  inside the window, under the ensemble's burn share for that cell (clipped as buildings are), divided by the number
  of buildings so that the two parts are on one scale. A run that burns the whole map pays for it.
- The fitted values, every setting tried and its score are listed in the findings. A setting is fitted on Camp and
  Tubbs and then checked, unchanged, on Marshall and East Troublesome before any freeze.
- Shorter or cropped runs may be used to search (stated per result); the reported calibration scores are full runs.

**A8 (24 Nov 2026, after the Camp reruns and 3-seed ensemble in FINDINGS-FIRE.md, before any run that uses a road
rule or the new fitting score).** What had been seen: the Camp run with gust share 0 and the incremental band
(one seed: 345 km2 at 21 h, Jaccard 0.54, AUC 0.67, 2.59 bits per building, A7 score 65.0) and its 3-seed ensemble
(375, 393 and 463 km2; 1.80 bits, AUC 0.735, A7 score 61.3). Two things in them prompt this amendment.
- Roads. LANDFIRE marks roads as lines of non-burnable cells one 30 m cell wide (codes 91 and 99). FireSim's
  surface spread cannot cross a non-burnable cell, so a grass field ringed by roads burns only if an ember lands in
  it (straight-edged patches in the gust 0.5 run; the map `sw_corner_fuel.png`). A two-lane road is far narrower
  than 30 m, and real grass fires in wind cross such roads. Rule R1, an option that is off unless asked for:
  - A thin cell is a 30 m cell of code 91 or 99 with burnable cells on both sides along at least one of the four
    axes through it (north-south, east-west, the two diagonals). Blocks two or more cells thick along all four
    axes (towns, lakes, rock, farm fields of code 93) are untouched, as are water (98) and snow (92).
  - A thin cell takes the fuel model of the lower-numbered of the two burnable cells of the first such axis
    (order as listed). In each hour it burns as that fuel if the flame length that fuel would have there is at
    least L_cross, and is non-burnable in that hour otherwise. Embers cannot start a fire on a thin cell.
  - L_cross is one more fitted number under A7, tried at 0 m (roads never stop fire), 1 m and 2 m, against R1
    off. These values are unsourced priors; no firebreak study has been opened for them.
  - Sandbox case S14 (FireSim `sandbox/CASES.md`), written before it runs, fixes what R1 must do on a made-up
    strip. With R1 off every earlier result is bit-identical.
- Fitting score. In the A7 score the area part was 59.5 against 1.80 for buildings, so buildings hardly count,
  as A7 said would happen; the size of the imbalance (33 to 1) is now known. From here the fitting score is the
  sum of two ratios, each 1.0 for a forecast with no skill: bits per building over the bits of this fire's own
  destroyed share, plus bits per 30 m cell (burned or not at the A7 perimeter) over the bits of the mapped
  burned share of the box. Lower is better. The A7 score is still printed beside it. This is a fitting score for
  calibration fires only; section 4 scores and section 5 gates are unchanged.
- Perimeter overlap is reported for every seed (it was seed 0 only).
- Tubbs is run with the same setup as Camp (A5, A6 wind W2): ignition from the official record, stated in its
  config with the source, 24 h window.

**A9 (24 Nov 2026, after the three Camp road-rule runs, before any search run or Tubbs run).** What had been seen:
Camp, one seed each, gust share 0, R1 at L_cross 0, 1 and 2 m (FINDINGS-FIRE.md): 479 to 485 km2 burned against
345 to 463 with R1 off, area at 9.45 h up from 153 to about 208 km2 (mapped 285 to 331), overlap at 21 h 0.48 to
0.51 against 0.54, AUC 0.69 to 0.71. The three L_cross values differ by less than seeds do. Also seen: the Tubbs
data folder has no mapped outline before 35 h after ignition and no satellite outline set.
- R1 is not decided by those runs. L_cross 1 and 2 m are dropped (no measurable difference from 0); the search
  carries R1 off and R1 with L_cross 0.
- Search stage 1, Camp only, cheap runs to find where to look: 10 h window, 10 m cells, 2 seeds per setting,
  ranked on the area ratio of A8 at the satellite outline at 9.45 h alone (the building record is not timed, so
  buildings cannot be scored at 10 h; their numbers are printed, not ranked). Settings: gust share 0, 0.25, 0.5;
  ember median throw 150 m (prior) and 400 m; R1 off and on. Twelve settings. Runs that do not fit the time cap
  are listed as not run.
- Stage 2: the best three stage-1 settings run 24 h with 3 seeds and are ranked on the A8 score. Stage 3 is A7's
  fit proper (8 or more seeds, Camp and Tubbs). Every setting tried is listed with its score.
- Tubbs: 44 h window so that the first three incident outlines (35.25, 40.9 and 43.0 h) fall inside; the A7
  outline is the last of these. The satellite fire detections (VIIRS 375 m points) in the first 24 h are shown
  on the map as a check of the first night and are not scored. An outline's edit time is not a fire time; only
  the progression file's times are used.

**A10 (24 Nov 2026, after search stage 1 and the Camp transects, before any run with wind aloft).** What had been
seen: stage 1 fitted only two settings in its time cap (FINDINGS-FIRE.md): at 10 h the road rule lifts burned area
from 159 to 217 km2 but the area ratio gets worse (1.28 to 1.43), because the added area is in the wrong place.
The modelled fire heads south-west in a wide blob; the mapped fire ran west in a long narrow run through Paradise
to the edge of Chico. Tables of the model's own inputs along the line from the ignition through Paradise toward
Chico (no outcome data) show why: the HRRR 10 m wind at the ignition blows about 38 degrees to the left of that
line; over Paradise it falls to 1 to 3 m/s by hour 4; west of Paradise it blows back toward the fire (the valley
wind). Timber litter under canopy spreads at about 1 m/min there. So in the model nothing carries fire west.
- Reasoning written before any test: lofted embers ride the wind at plume height, not the 20-ft wind. In a
  downslope windstorm the flow aloft keeps its direction and speed over ground where the surface wind is weak or
  reversed. Ember transport (throw distance and direction) may therefore use wind aloft; surface spread keeps the
  20-ft wind of A6.
- Options, all from the same HRRR analysis file (no new source): E0 embers on the 20-ft wind (as now); E80 on the
  wind at 80 m above ground; E850 on the 850 hPa wind. The winds aloft enter as they are (no gust share, no
  height adjustment). This is one more choice fitted under A7/A9, on Camp and Tubbs, and it must hold unchanged on
  Marshall and East Troublesome before a freeze.
- The throw-distance law is unchanged (distance scales with wind speed), so a faster wind aloft also throws
  further; the ember median distance is fitted together with this choice, starting from the prior 150 m.
- First look, Camp 10 h, 2 seeds, gust share 0, road rule on at L_cross 0: E0 (already run), E80, E850, each at
  ember median 150 m. Ranked on the area ratio at the satellite outline at 9.45 h as in A9 stage 1.
- A stated risk: this is a second fix found by looking at Camp. If it helps Camp and not Tubbs it is dropped.

**A11 (24 Nov 2026, after Camp stage 2 and the Tubbs cross-check of A10, before any run of the option below).**
What had been seen (FINDINGS-FIRE.md): Camp 24 h, 3 seeds, embers on the 850 hPa wind: overlap at 21.45 h 0.51 to
0.61, cells beat no skill, 1.64 bits per building, AUC 0.736, A8 score 2.57. Tubbs 44 h, 2 seeds: with embers on
the 20-ft wind 1.9 times the mapped area, overlap 0.44, AUC 0.753, A8 4.08; with embers on the 850 hPa wind
worse on every score (overlap 0.32 to 0.39, AUC 0.657, A8 6.47), the extra fire going north-east behind the
ignition over the later, calm hours.
- By A10's own rule, E850 as tested is dropped: it helps Camp and hurts Tubbs. E80 is dropped too (little gain).
- A narrower hypothesis, written before any test: how high an ember rides depends on how hard the fire under it
  burns. A weak fire's embers stay in the surface wind; only an intense fire lofts them into the flow aloft. In
  E850 every ember, even from a creeping fire in the calm second day of Tubbs, was carried by the 850 hPa wind.
  Option EP: the ember's wind is the vector blend (1 - w) x 20-ft wind + w x 850 hPa wind, with
  w = min(1, I / I_loft)^(2/3), I the fireline intensity of the cell (or house) that throws it. The 2/3 power is
  the usual scaling of plume rise with heat release, taken as a form only; I_loft is an unsourced prior, tried at
  10 and 30 MW/m.
- Test, same as A10 and decided the same way: Camp 24 h 3 seeds and Tubbs 44 h 2 seeds, against E0 on each. EP
  is kept only if it does not make Tubbs worse than E0 on the A8 score (beyond seed spread) and still improves
  Camp. Otherwise embers stay on the 20-ft wind and the miss on Camp is carried as a stated limit until the wind
  comes from WindSim (candidate M4).
- Stated risk: this is the third wind fix found on these two fires. It is the last one tried before the freeze.
- Tubbs overrun, noted not fixed: with E0 the model burns 1.9 times the mapped area by 43 h and 2.2 times the
  houses. Most of the mapped fire burned in the first night; the model keeps spreading for 30 more hours with no
  suppression, which it does not model (section 7). Tubbs is therefore also scored on the first night: the
  share of satellite fire detections (VIIRS 375 m, clipped to the final Tubbs outline plus 500 m so that the
  other fires of that night are left out) at the 10:15Z pass on 9 Oct (5.5 h) that fall on cells the model has
  burned by then, and the modelled area at that time. Reported beside the outlines, not part of the A8 score.

**A12 (24 Nov 2026, after A11's results and a table of modelled house ignition times, before any run of the rule
below).** What had been seen: Tubbs E0, 2 seeds: 11,619 houses modelled destroyed against 5,185 recorded. By
modelled ignition time, 94 % of the recorded-destroyed houses the model lights are lit in the first 9 h, while the
20-ft wind where the houses burned is 8 to 12 m/s; after hour 8 that wind is 1 to 3 m/s for the rest of the run,
and in that calm the model lights 4,646 houses that stood against 231 that burned. (It also lights about 6,200
houses that stood during the first 9 h; this amendment does not address those.) Camp: 58 % of modelled losses
fall at 12 to 24 h, when the model's wind over Paradise is 7 to 9 m/s; in hours 3 to 10 it is 0.5 to 5 m/s there.
- Reasoning, written before any test: in FireSim a burning house lights its neighbours and throws brands at a
  rate that does not depend on wind. Fire spreading house to house through a town is a wind-driven event; in
  calm air separation distances and crews hold it. Rule H1, off by default: the neighbour strength k_s and the
  number of brands a burning house throws are both multiplied by min(1, U / U_ref)^2, U the 20-ft wind at the
  burning house at its ignition time. U_ref = 10 m/s. Both numbers are unsourced priors; one setting is tried.
- Test: Tubbs 44 h 2 seeds and Camp 24 h 3 seeds, embers on the 20-ft wind, road rule on, gust share 0, against
  the same without H1 (for Camp that pairing is run fresh, since the earlier E0 ensemble had the road rule off).
  H1 is kept if the A8 score on Tubbs improves by more than the seed spread and Camp's does not worsen by more
  than its seed spread. Predicted before running: Tubbs improves; Camp may worsen, because the HRRR wind over
  Paradise is weak in the hours when the town burned.
- This does not model suppression. The limit in section 7 stands.

**A13 (24 Nov 2026, after A12's result, before any station table is read).** What had been seen: four wind-dependent
changes (E80, E850, EP, H1) each helped one calibration fire and hurt the other. The modelled wind over Paradise
is 0.5 to 5 m/s in the hours the town burned. The only station comparisons seen so far are the five or six
stations in each fire's config, as 24 h means (HRRR over station speed 0.80 to 1.94 on Camp, 0.64 to 2.13 on
Tubbs).
- Wind is now checked against stations alone, with no fire outcome in the loop. For each calibration fire's
  window, every ASOS and RAWS station inside the fire box plus 0.75 degrees is compared hour by hour with the
  HRRR analysis at the nearest lattice point: 10 m speed and direction, surface gust against the station's peak
  gust, and the 80 m and 850 hPa speeds for reference (FireSim `scripts/wind_check.py`).
- Rules fixed before the tables are read, the same for every fire:
  - Speed scale. If at stations whose observed mean is 5 m/s or more the median HRRR-to-station speed ratio is
    below 0.85 on both Camp and Tubbs, the surface wind is multiplied by 1 / (the mean of the two medians) in all
    fire runs. If it is 0.85 to 1.15 on both, no scale. Anything else: no scale, and the disagreement is reported.
  - Gusts. If the median HRRR-to-station gust ratio at those stations is 0.85 to 1.15 on both fires, the HRRR
    gust is taken as sound and the gust share stays a fitted number (A7); otherwise the gust share is fixed at 0.
  - Terrain. If the speed ratio falls with station height above the HRRR terrain (correlation below -0.5 on
    both fires, at least 8 stations each), a height correction is designed under a new amendment; otherwise not.
- Stations are few and sit where stations sit (airports, ridges). The check can show a bias; it cannot show the
  wind in a place with no station, such as Paradise. That limit stays stated.

**A14 (24 Nov 2026, after A12's results, before any error table of a real fire is built or read).** Kyle's
directive of this date: "when we run into issues, we should always try to update compressor to include formulas for
the specific domain, and then run it against the model." What had been seen: the result files of the runs listed in
the findings (areas, overlaps, building bits, detection shares), the Camp wind transects, and that four
wind-dependent changes each helped one fire and hurt the other. No per-point error table of a real fire exists yet.
Section 6 is carried out now on the two calibration fires that have runs, with these additions.
- Model runs read: the priors with the road rule on and embers on the 20-ft wind (E0), seed 0, Camp 24 h and
  Tubbs 44 h. Tables come from FireSim `scripts/error_table.py`; the search from `physics/fire_real.py` here.
- Three error tables per fire, on a 300 m lattice of the 30 m grid (water left out):
  - T1 arrival: points that both the record and the model show burned. Observed arrival = the time of the first
    outline or satellite pass (a detection within 190 m) that shows the point burned, so it is an upper bound read
    on a coarse clock. y = ln(observed / modelled arrival), modelled time floored at 10 min. Grid 0.05.
  - T2 reach: points that the record or the model shows burned by the last outline inside the run. y = +1 where
    only the model burned (overrun), -1 where only the record did (missed), 0 where both. Grid 0.25.
  - T3 houses: 300 m blocks with at least 5 scored buildings. y = observed minus modelled destroyed share. Grid
    0.05. Blocks carry no coordinates, addresses or parcel facts.
- Fire formulas added to the term library (`physlib/terms.py`, list `FIRE`), each read at the point over the hours
  before it burned, or along the straight line from the ignition:
  - wind: ln of the 10 m speed (Rothermel's wind factor is a power of speed, so a wrong exponent or a speed bias
    shows as a slope on ln U); ln(850 hPa speed / 10 m speed) (strong wind aloft that the surface analysis does not
    bring down); lee descent, the fall of the ground along the 850 hPa wind over 1 km and over 5 km (downslope
    windstorm); ground height above the weather model's own terrain; the share of the ignition-to-point line that
    runs with the wind.
  - calm creep and burnout: the share of hours with wind under 2 m/s; the share of hours at night; ln hours since
    ignition; hours of daylight after the first 9 h (when crews work); distance to the nearest road line.
  - fuel and terrain along the line: shares of grass, shrub, timber, non-burnable and developed cells; canopy
    cover; slope along the run (tan squared times cos) and slope of any aspect.
  - moisture: dead 1-h fuel moisture and relative humidity over the hours before the point burned.
  - embers: `ember(L)`, canopy fuel upwind within L weighted exp(-x/L), L from 100 m to 3 km.
  - houses (T3 only): buildings within 100 m; `gap(d)`, the share of buildings whose nearest neighbour is closer
    than d, d from 3 to 30 m; the block's modelled share itself (contrast: the model over- or under-states its own
    range); the modelled share lit by embers, flame and neighbours.
  - a per-fire offset in every table.
- Reading rules, fixed now. The search is `physlib.guards.analyse` with up to 4 terms. Points 300 m apart are not
  independent, so resampling and shuffles go by 3 km blocks. A term is a finding only if all four hold: it is in
  the winning formula on the two fires pooled; it transfers both ways (fitted on Camp it shortens the code on
  Tubbs, and the reverse); the block bootstrap (30 resamples) picks its group at least half the time; and the
  shuffled-null search (19 shuffles) never saves as many bits. A term that wins on one fire only is reported as
  that fire's lead and drives nothing.
- Planted check on the real tables, before the real answers are read: for each of five terms (ln wind, lee
  descent, calm share, ember(L), gap(d)) a made-up error equal to that term times a known coefficient plus
  shuffled noise of the real error's spread is searched; the term's group must be named, and a noise-only twin must
  name nothing. This tests whether these columns can be told apart on these fires; it is weaker than A4, which
  switched a mechanism off inside the model. A term group that fails its planted check is not read on the real
  errors.
- What a finding may do: name the next change to FireSim or to its wind input, which then needs its own amendment
  and the usual Camp and Tubbs test. No coefficient fitted by the search goes into FireSim. Both fires are
  calibration fires, so nothing here is a validation result.
- Stated limits: the observed clock is coarse (two to four outlines or passes per day), so T1 cannot see errors
  shorter than the gap between passes; T2 treats the record as complete, though outlines can include unburned
  islands and passes miss cool or hidden fire; two fires give one transfer test in each direction and no more.
- A14 note (24 Nov 2026, same day, after the code ran on a made-up fire only, before any real table exists). How
  the code carries out the rules above, where the text left room: hourly inputs are averaged over the hours before
  the observed arrival in T1, and over the whole window to the last outline in T2 and T3, because a point that only
  the model burned has no observed time. In T1 the term ln(observed hours) holds the observed time that the error
  itself is built from, so it is entered as bookkeeping and can never be a finding; T2 and T3 use ln(distance from
  the ignition) for the same question. Each fire's mean error is taken out before the search and reported as that
  fire's level (the per-fire offset). The search runs with the engine's standard limit of 5 terms, not 4, grammar
  off. Noise for the planted checks and the null is the real error with its 3 km block means dealt to other blocks
  and the remainders shuffled inside each block; a planted term is given half the real error's spread. Camp's
  house labels are the USGS building-loss points, as in every Camp run. HRRR relative humidity is not used; the
  humidity and moisture inputs are the station values the model itself read.

**A15 (24 Nov 2026, before A13's station tables are read and before any WindSim run for a fire).** Kyle, same day:
"if we can incorporate windsim into it at all, that'd be pretty cool". What had been seen: the Camp transects
(the 10 m analysis wind is 0.5 to 5 m/s over Paradise while the 850 hPa wind holds 10 to 19 m/s), and the WindSim
thread's notes on what its terrain solver needs at this size. No WindSim run on fire terrain exists.
- What is built: wind option W3 for FireSim, `--wind windsim`. WindSim's steady RANS solver is run on terrain alone
  (LANDFIRE elevation, roughness from the fuel and canopy layers, no buildings) over a Camp sub-domain in UTM 10N
  that holds the ignition, Concow, Paradise and the canyon mouths: longitude -121.70 to -121.38, latitude 39.68 to
  39.88. Cells about 150 m, first cell 15 to 25 m tall, inflow from 30, 45 and 60 degrees (the range of the 850 hPa
  wind in the first ten hours). CPU only, each run capped at 2 hours.
- How the wind is made: each run gives a speed-up map S (10 m speed over the run's own 80 m speed averaged on the
  upwind third of the domain) and a turning map. For fire hour h the reference is the HRRR 80 m wind averaged over
  the same upwind third; the 10 m wind at a cell is S times the reference speed, turned by the turning map, with S
  and turning interpolated between the three bearings. Hours whose reference bearing is outside 15 to 75 degrees,
  and cells outside the sub-domain, keep the HRRR 10 m wind. Nothing is fitted.
- Scored on stations before any fire run reads it, with A13's station set inside the sub-domain and hours with an
  observed speed of 3 m/s or more: RMSE of speed, median ratio and mean direction error, for W3 and for the plain
  HRRR 10 m wind, pooled and split into windward and lee stations. A station is lee if the ground 5 km upwind
  (along that hour's reference bearing) is more than 100 m higher than the station.
- Rule: W3 goes into a Camp fire test only if its pooled speed RMSE is at least 10 % lower than HRRR's and it is
  not more than 10 % worse in either class. With fewer than 4 stations in the sub-domain the check is reported as
  too thin to decide, and W3 is not used. A Camp fire test with W3, if earned, is E0 priors against the same run on
  HRRR wind, 3 seeds, read on the A8 score; Tubbs needs its own domain and the same station rule before W3 counts
  for anything.
- Predicted before running: W3 beats HRRR at ridge and windward stations and fails in the lee, because a neutral
  steady solver has no stable layer and so no downslope windstorm or gap jet. If so, that failure is the finding,
  and the next step is a stability term, not a tuned speed-up.
- A15 note (24 Nov 2026, same day, code written and tested on made-up runs only; no WindSim terrain run exists).
  Cells are 200 m, not 150 m, with a 20 m cell on the ground, so that one bearing fits the 2-hour cap on 8 cores
  (about 2 million cells); the domain is a 52 km square turned with the wind, flattened over its outer 4 km. The
  HRRR reference for an hour is taken over the upwind third of the middle (45 degree) run for every bearing. A
  reference bearing up to 15 degrees outside 30 to 60 uses the nearest run with the hour's own bearing. W3 ignores
  the gust share. Code: WindSim `d72ebc3` (terrain-only cases, `projects/camp_fire_terrain.yml`), FireSim
  `firesim/windsim.py` and `run_fire.py --windsim`.
- A15 note 2 (24 Nov 2026, after the WindSim thread's review of the code, still before any terrain run). Changes
  at WindSim `c31a54b`: terrain and roughness are fetched over the square that holds the turned domain, and the
  domain's edges follow a plane from the inlet side's median ground level to the outlet side's, instead of one
  level all round, so a wind off high ground does not meet a made-up ramp. Stated limits added: the inflow is a
  neutral log law with no Coriolis turning and no capping inversion (the inversion is what makes a downslope jet);
  roughness in forest is up to 2 m with no displacement height, so the model's 10 m wind there is near canopy top
  while RAWS masts stand in clearings. The rule is read on the 10 m wind as modelled; a second table, reported
  beside it and deciding nothing, brings the model's 80 m wind down to 10 m over a clearing roughness of 0.1 m with
  the log law. The 10 m speed along the line from the inlet to Paradise is drawn for each run so that an edge
  artefact would show. checkMesh's worst non-orthogonality is reported; above 70 the run is repeated with two
  non-orthogonal correctors before it is read.


### A16 (24 Nov 2026, 8:50 PM ET): a wider terrain wind check, a contrast form of W3, and the overrun

Seen before writing this: A15's first station score (2 stations inside the 28 km focus square, W3 about twice too
fast at both, too thin by A15's rule, W3 not used) and the x1.30 reruns (Camp A8 2.81 to 3.48, Tubbs 4.08 to 6.01;
both fires burn two to two and a half times the mapped area). Both are in `FINDINGS-FIRE.md`. No blind fire has
been read.

Part 1, the wind. Weather only; no fire run reads it until its rule is met.
- W3c ("contrast form"): the 10 m wind is HRRR's own 10 m wind times S / mean(S), where S is the terrain run's
  10 m speed over its 80 m reference as in A15 and mean(S) is S averaged over the 3 km HRRR cell around the point
  (a 3 km box mean on the terrain grid). The level stays HRRR's; the terrain run supplies only the pattern inside
  each HRRR cell. Direction is the terrain run's turning added to HRRR's 10 m direction. W3 as in A15 is scored
  beside it and stays unused.
- Domain: one Camp terrain project wide enough for the A13 stations: a 120 km focus square on the same centre,
  12 km margin, 500 m cells, 20 m first cell, top 12 km, the same three bearings (30, 45, 60 degrees). Each run
  must fit A12's 2-hour CPU cap on 8 cores; if the mesh passes 2.5 million cells the focus square is cut to
  100 km, and that is said in the result.
- Score, same hours and 3 m/s floor as A15: speed RMSE, median ratio and direction error for HRRR, W3 and W3c,
  pooled and split windward and lee.
- Rule: W3c earns a Camp fire test only with 8 or more stations inside, pooled speed RMSE at least 10 % under
  HRRR's, and not more than 10 % worse in either class. Fewer than 8 stations is "too thin" and nothing changes.
  If earned, the test is A15's (E0, 3 seeds, A8 score against the same run on HRRR wind) with the x1.30 scale off,
  because the scale and the contrast would otherwise correct the same thing twice.
- Predicted before running: W3c beats HRRR at sheltered and valley stations (where A13 found HRRR 2 to 8 times too
  fast) and gains little on ridges; W3 stays about twice too fast everywhere.

Part 2, the overrun. Under the directive of 24 Nov the next step is terms, then the search, before any setting.
- New inputs in the T2 table (FireSim `scripts/error_table.py`), all from the run and the public record, none
  from the blind fires: `flank` (sine of the angle between the hour-weighted wind and the line from the ignition
  to the point, 0 on the head axis, 1 on the flank), `back` (1 if the point is upwind of the ignition), `late`
  (share of the window already gone when the model arrives), and `edge_age` (hours between the model's arrival and
  the end of the window).
- New terms in `physlib.terms.FIRE`: `flank_spread`, `backing_spread`, `late_arrival`, `burn_duration`, each with
  its formula and plain meaning. A planted check on each runs before the real errors are read, and the A14 rules
  for a finding stand unchanged (two-way transfer, bootstrap share, block nulls, noise twin).
- Predicted before running: the overrun sits on the flanks and behind the ignition and late in the window (the
  model's flank and backing fire never go out, and nobody holds a line), so `flank_spread` and `late_arrival` are
  named on both fires. If neither transfers, the overrun is not a shape the model's own inputs explain, and the
  next step is the two extra calibration fires, not a spread limit.
- Nothing fitted goes into FireSim from this. A term that pays becomes a spread rule only through its own
  amendment, tested on Marshall and East Troublesome before any blind fire.
- A16 note (24 Nov 2026, 9:05 PM ET, before the search is run). The four overrun inputs need no new run: `flank`
  and `back` follow from the hour-weighted wind cosine already in the A14 tables (`run_all`: flank is
  sqrt(1 - run_all^2), back is 1 where run_all < 0), and `late` and `edge_age` from the model's arrival time and
  the window's end. `late` and `edge_age` exist only where the model burned, and on the full T2 table "the model
  never arrived" would by itself mark every missed point, so they are searched on a new table T2o: the points the
  model burned by the last outline, y = 1 where the record does not show them burned and 0 where it does, each
  fire's mean removed. `flank_spread` and `backing_spread` are searched on T2 and T2o. T1 and T3 are not rerun.
  The tables are the A14 ones (E0 runs on the plain HRRR wind). Output: `results/physics/fire_real_a16.json`.
- A16 note 2 (24 Nov 2026, 9:35 PM ET, before any wide terrain run). Mesh arithmetic done first: a 120 km focus
  square needs a 180 km turned domain, about 4.7 million cells at 500 m, and the 100 km cut still gives about
  3 million. So the project (`projects/camp_fire_terrain_wide.yml` in WindSim) is a 96 km focus square in a 144 km
  turned domain with 600 m cells (about 2 million cells), DEM at 120 m, top 12 km. A station counts when it is
  inside a run's inner square (4 km in from its edge) for all three bearings. The 3 km box mean is 5 cells.
  Code: FireSim `firesim/windsim.py` (`contrast=True`, `box_mean`) and `scripts/wind_check_w3.py` (column `w3c`).
- A16 note 3 (24 Nov 2026, 10:00 PM ET, no terrain run yet). The wide domain's relief is 3,153 m (valley floor to
  Lassen), and WindSim refuses a top under 4 times the relief, so the top is 13 km, not 12. Nothing else changes.

### A17 (25 Nov 2026, 12:05 AM ET): the free-burning window (Kyle: "can we test fires up until intervention?")

Seen before writing this: every run's area and overlap at each mapped time (they are in each `result.json`: Camp
at 3.35, 9.45 and 21.45 h; Tubbs satellite passes from 5.5 h and outlines from 35 h), and A16 part 2's lead that
the model's late ground is its wrong ground. Not yet seen: the hour-by-hour station wind series, and how much of
Tubbs's modelled fire at 5.5 h lies inside the final outline.

- Why: the model has no crews, no burnout and no end to a wind event's run except the weather it is given. A
  score at the last outline mixes two questions: does the fire spread right while nothing holds it, and does it
  stop when the real one did. They need different fixes, so they are read apart.
- The cut, from weather only: the end of the wind event. Stations are A13's windy class (window mean of 5 m/s or
  more) in `results/wind_check/<fire>.csv`. Take the hourly mean of their observed speed and its 3-hour running
  mean; the event ends at the first hour after the peak from which that running mean stays under half its peak to
  the end of the window. If it never does, the cut is the end of the window. Reported beside it and deciding
  nothing: the first containment figure in the public incident status reports, once opened.
- Free-burning score: at each mapped time up to and including the first one after the cut, the area ratio
  (model / mapped) and the overlap (Jaccard). For a fire with no outline that early (Tubbs), the share of the
  ground the model has burned at each satellite pass before the cut that lies outside the final outline (certain
  overrun), and the share of that pass's detections on the modelled fire. Houses have no time in the record, so
  the house score stays whole-window.
- What it decides: nothing at the gate (section 5 stands, read at the last outline on the blind fires). From
  here each calibration run is reported with both readings, and a spread or wind setting is judged on the
  free-burning score while a stopping rule (night recovery, crews, burnout), which needs its own amendment, is
  judged on the growth after the cut.
- Predicted before computing the cut: both events end 8 to 12 hours after ignition; before it the x1.30 wind is
  the better of the two on Camp, and after it the modelled fires roughly double while the mapped ones grow under
  a fifth.

### A18 (25 Nov 2026, 9:20 AM ET): containment, fires that ran unheld, and a stopping term that needs no crew record

Kyle, 25 Nov: "I think the fire should have the containment times ... containment times aren't going to be
included in the model though either. So are there wildfires that weren't contained that we can try to replicate?"

What the record says (pages opened today; news figures, to be replaced by the daily incident status reports once
their licence is confirmed, which the dataset's page did not show):
- Camp: "five percent contained" at 70,000 acres in a report updated 9 Nov 2018 2:28 PM EST, about 29 hours
  after ignition (PBS NewsHour), and still 5 percent at 90,000 acres that evening (ABC30). No containment figure
  falls inside the 24-hour window.
- Tubbs: 0 % on 9 and 10 October, 10 % on 11 October (Wikipedia's table, whose prose differs by a day in places).
  The 44-hour window ends on 10 October: 0 % throughout.
- So inside the windows already scored, both fires were officially uncontained. Containment percent cannot be
  the stopping rule for them: the mapped fires slowed with no line around them. Crews defending houses and
  holding flanks are not in that figure.

Fires that ran unheld, both already calibration fires in section 1, are the pure-physics test:
- Marshall, 30 Dec 2021: "two reported initial ignition points, occurring at 18:08 and 19:00 UTC" (Juliano et
  al. 2023, Geophysical Research Letters); "wind gusts over 50 m/s made it impossible to get the flames under
  control", gusts "finally weakening to less than 25 m/s by 1930 LT" and "the real break" a cold front "around
  0100 LT 31 December"; "burned through 24 km2" (Benjamin et al. 2023, Weather and Forecasting). These published
  ignition times lift A5's hold on Marshall runs. Window: 18:08 UTC 30 Dec to 08:00 UTC 31 Dec (14 hours),
  scored on area at every mapped time and at the end; houses are not scored (no licence stated for the damage
  record, section 1 note).
- East Troublesome, the run of 21 to 22 October 2020: "growing from 30,000 acres to 170,000 acres in about 24
  hours" (NASA Earth Observatory); "you couldn't have stopped it" (Grand Lake fire marshal, CPR News); snow on 25
  October slowed it. The run starts from the last outline mapped before the wind (the model is lit over that
  whole outline, a new start mode) and is scored at each outline to the first one after the wind eases by A17's
  rule.
- Both use E0 priors, HRRR wind as analysed and with A13's x1.30, 3 seeds, and are read on A17's free-burning
  score and on what happens after the wind eases. Predicted: with x1.30 the area ratio during the run is between
  0.5 and 1.5 on both; after the wind eases the model keeps growing on Marshall (no snow, no crews in it) by more
  than a quarter.

The stopping term (compressor first, on the Camp and Tubbs tables in hand, before any model change):
- Input `ros`: the model's own spread rate where it reached a point, from the arrival times on the 300 m lattice
  (300 m over the largest arrival-time difference to a lattice neighbour the fire came from, in m/min).
- Term family `slow_edge(v)`: 1 where `ros` is under v, for v = 0.5, 1, 2, 5, 10 m/min, on table T2o. Meaning: a
  slowly moving edge goes out on its own or is held by whoever is there, and the model lets it creep on for
  hours. It needs no crew record, so it could be a rule in a forecast. Planted check first; A14's rules for a
  finding stand. Run: `fire_real.py --a16 --a18`, T2o only.
- Predicted: named on both fires with a positive sign (slow modelled edges are the overrun), best near 1 to 2
  m/min, and it takes over from `burn_duration`. If it transfers both ways it becomes a stopping rule proposal
  (an edge slower than v for N hours goes out), tested unchanged on Marshall and East Troublesome.

### A19 (25 Nov 2026, 10:10 AM ET): NASA FEDS outline times are half-day labels, not times

Found while setting up Marshall: the FEDS outline labelled 30 Dec 2021 12:00 already shows 19.85 km², but the
fire started at 18:08 UTC. FEDS `t` is 00:00 for a date's night satellite passes and 12:00 for its afternoon
passes, and the runner has read it as a clock time.
- From now on a FEDS outline's time is the last VIIRS pass over the fire in its half-day (from the fire's own
  detections file), written to a derived outline file with the original label kept and the rule noted in the
  data manifest. Where no pass time can be found the outline is not scored.
- Camp's past readings are affected and are corrected, not reinterpreted: its FEDS outlines were scored at 9.45 h
  and 21.45 h and are really about 18.5 h and about 30 h (outside the 24-hour window); the afternoon outline of
  8 Nov (about 6 h) was wrongly dropped. Every Camp score against a FEDS outline so far (A8 scores, A12, A13's
  x1.30 rerun, A17's table) compared the model with an outline mapped hours later, which flatters the model's
  area. Tubbs uses GeoMAC outlines and is not affected. The GeoMAC times stand.
- Camp is rescored with corrected times on the existing runs (same seeds, no new physics): the last outline in
  the window becomes the fit outline. Old and new scores are both kept in the findings. No rule or gate changes;
  the Camp window stays 24 hours.
- A18 note (25 Nov 2026, 10:45 AM ET, before any Marshall run). No source that could be opened gives a
  coordinate for either Marshall ignition (the paper's are in supporting text behind a bot check; the Sheriff's
  release of 8 Jun 2023 gives a street address for the first and a description for the second: "just south of the
  Marshall Mesa Trailhead"). The run is lit at one place: 200 m south of the trailhead's OpenStreetMap position,
  at 18:08 UTC, with a 150 m radius to cover the stated uncertainty and the "several 100s of meters" between the
  two starts. It is labelled in the fire's file as derived from a landmark description, not surveyed. The private
  address is not geocoded. If the paper's supporting text is opened later, the point is replaced and the run
  repeated, both kept.

### A20 (25 Nov 2026, 1:15 PM ET): the search on four fires

Seen before writing this: the area tables of Camp (rescored), Marshall and East Troublesome (in FINDINGS-FIRE and
the runs' result files): three fires overrun after the wind, East Troublesome is far too slow in its first night
(model 0.29 of the mapped area at 12 h). The new error tables have not been read.
- Tables: Camp from `camp_2018_E0_a19.npz` (corrected outline times; replaces the A14 Camp table from here on),
  Tubbs as before, `marshall_2021_E0.npz`, `east_troublesome_2020_E0.npz` (distance and line columns measured
  from the start outline's edge; points inside it left out). All E0 runs on the HRRR wind as analysed.
- Search: T2 (reach) and T2o (points the model burned), with every term of A14, A16 and A18 that all four fires
  have inputs for. Each fire's mean is removed as before.
- Hold-out: leave one fire out, four ways. "Transfers" now means the formula fitted on three fires beats the
  one-correction baseline on the fourth, and a finding needs it on at least three of the four folds with the
  other rules unchanged (bootstrap share 0.5, every block null beaten, planted check, noise twin naming nothing).
  Requiring three of four is set now because Marshall's table is small (230 points) and one fold may be noise.
- Predicted: on T2 the winner is led by a term that separates East Troublesome's misses (timber or canopy on the
  line, wind aloft) from the others' overruns, and it does not transfer to or from East Troublesome; on T2o the
  late group is still named in most draws and holds on at least three folds.

### A21 (25 Nov 2026, 4:30 PM ET): a search setting under which noise names nothing

Seen: every fire search since A14 has had its noise-only twin name terms, on every table. The cause is in the
accounting: the code length counts each 300 m lattice point as an independent reading, but fire errors come in
patches kilometres across, so a term that fits one patch is paid for as if it explained hundreds of readings.
- Fix, set on noise alone: the lattice is thinned to every s-th point in both directions (s = 1, 2, 3, 5, 10,
  that is 300 m to 3 km spacing) and, for each s, five noise-only twins (block-shuffled as in A14, seeds 0 to 4)
  are searched. The setting is the smallest s at which all five twins name nothing. The real errors are not
  read while s is chosen. If no s silences the twins, that is reported and no real search is run.
- The real search then runs once on the thinned tables of the four fires (A20's tables, T2 and T2o, the same
  terms, leave one fire out, three of four folds to transfer), with the planted checks repeated on the thinned
  tables. A planted term that can no longer be found at that s is reported as out of reach, and a real term in
  its group cannot be a finding or a null.
- Predicted: s is 5 or 10; at that spacing the late group is still named on T2o, and little else is.

### A22 (25 Nov 2026, 5:15 PM ET): the held-edge rule H2, a stopping rule that needs no crew record

Seen: the three overrun fires and their hourly weather as the runs used it (humidity never above 37 % in the
Tubbs, Camp and Marshall windows, so dead fuel stays far below its moisture of extinction and night recovery
cannot stop these modelled fires); A18's and A20's leads (slow and late ground is the wrong ground). A21's
calibrated search has not been read.
- Source for the numbers (opened): NWCG Fire Behavior Field Reference Guide, PMS 437, "Interpreting expected
  surface fire behavior": flame length 0 to 4 ft "Direct, Hand"; 4 to 8 ft "Direct, Equip"; 8 ft and more
  "Indirect". Nothing is fitted.
- Compressor first: the error tables gain `flame_arr`, the model's head-fire flame length (m) at each
  model-burned point in the hour it arrived. Term `held_edge(L)`: 1 where `flame_arr` is under L, L = 1.22 m and
  2.44 m (4 and 8 ft), on T2o at A21's thinning, four fires, leave one out. A planted check first.
- The rule H2(L) in the model: in any hour a cell whose head-fire flame length is under L does not carry fire
  (its spread rate is zero for that hour; it can burn in a later hour if the flame there rises). It stands for
  an edge weak enough to be worked directly, wherever people are fighting the fire, and it uses nothing a
  forecast would not have. It is on from the first minute. Embers can still light such a cell, as now.
- Test: Camp, Tubbs, Marshall and East Troublesome, E0 priors, HRRR wind as analysed, 3 seeds, H2(1.22) and
  H2(2.44), against the E0 runs in hand. Read on area ratio and overlap at every mapped time and A8.
- Rule, fixed now: H2(L) is kept as a candidate setting only if, against E0, the area ratio at the last mapped
  time moves toward 1 on each of Camp, Tubbs and Marshall, the overlap at each fire's first mapped time (or
  Tubbs's 5.5 h outside-share) is not worse by more than 0.02, and A8 improves on at least two of the three with
  none worse by more than 0.1. East Troublesome is reported and cannot veto, since its fault is the opposite
  one. If both L pass, the smaller is kept. A kept rule is still only a candidate: it goes to the freeze with
  the rest and is judged on the blind fires.
- Predicted: H2(1.22) trims Tubbs and Marshall but not enough (ratio at the end still above 1.4); H2(2.44)
  brings Marshall under 1.5 and Tubbs under 1.5, does little on Camp, and makes East Troublesome slightly worse.
  `held_edge(2.44)` is named on T2o with a positive sign.

### A23 (25 Nov 2026, 6:00 PM ET): why the timber fire does not run, as a term first

Seen: East Troublesome's area table (model 0.29 of mapped at 12 h), its flame lengths at arrival (median 0.71 m:
the model has a surface fire creeping under the canopy), and A21's surviving leads (timber on the line goes with
misses). Not seen: any crown-fire quantity at the missed points.
- Idea, from the model's own crown rule (Van Wagner's initiation intensity as in Scott and Reinhardt 2001,
  RMRS-RP-29 eq. 11, already in `firesim/crown.py`): the surface intensity needed to light the canopy rises
  steeply with foliar moisture, which the model fixes at 100 %. Beetle-killed lodgepole holds dead needles at a
  small fraction of that. At 10 % the threshold is about a ninth of its value at 100 %.
- New table columns (all four fires, every lattice point, so missed points have them too): `fli_max`, the
  largest head-fire surface intensity the model gives the point in any hour of the window (kW/m), and the
  canopy base height and bulk density there.
- Term family `dry_crown(m)`: 1 where there is canopy and `fli_max` exceeds the initiation intensity at foliar
  moisture m, for m = 10, 30, 60, 100 %. On T2 at A21's thinning, four fires, leave one out; planted check first.
- Predicted: named on T2 with a negative sign (ground that would crown with dry foliage is ground the model
  missed), best at 30 % or under, carried by East Troublesome, and so not transferring. If named, the lead is a
  dead-canopy map (the Forest Service's public insect and disease survey) feeding foliar moisture, proposed
  under its own amendment. If not named, the miss is not a crown-initiation threshold and the next terms are the
  wind aloft over the fire and spotting distance.

### A25 (26 Nov 2026, 10:40 AM ET): the crown fire rate, as a term first and then as a model setting

Seen: A23's class table (East Troublesome's missed share is about 0.8 even in cells the model already crowns, so
the miss is in how fast crown fire runs once it starts, not whether it starts); East Troublesome's E0 weather as
the runs used it (1-h dead fuel 3.9 to 6 % over the first 16 hours, humidity 15 to 28 %); a preview of the Cruz
rate below on the A23 tables (median active rate with the window-mean wind and the tables' moisture: Camp 60,
Tubbs 34, Marshall 87, East Troublesome 21 m/min). A24 (the device's run-phase held edge) is separate and not read.
- Source (opened 26 Nov): Cruz, Alexander and Wakimoto 2005, "Development and testing of models for predicting
  crown fire rate of spread in conifer forest stands", Can. J. For. Res. 35: 1626-1639. Active rate eq. 3,
  CROS_A = 11.02 U10^0.90 CBD^0.19 exp(-0.17 EFFM) (m/min; U10 the 10 m open wind in km/h, CBD kg/m3, EFFM the
  1-h dead fuel moisture in %); criterion for active crowning eqs. 1-2, CAC = CROS_A / (3.0 / CBD); passive rate
  eq. 4, CROS_P = CROS_A exp(-CAC). Fitted on 5 to 32 km/h, CBD 0.12 to 0.48, EFFM 7 to 12 %, level ground;
  "not applicable for zero wind speed". The PDF's equation text is garbled in extraction; the form is checked
  against the paper's own worked numbers (137 m/min at 50 km/h, CBD 0.2, EFFM 4; within 3 %). The model now
  uses Rothermel 1991 (3.34 x fuel model 10 at 40 % of the open wind) with Scott and Reinhardt's crown fraction.
- Part 1, compressor first. Term `cruz_crown`: b · [canopy] · ln(CROS / 1 m/min), CROS = CROS_A where CAC >= 1
  and CROS_P otherwise, from each point's own window-mean HRRR 10 m wind, its LANDFIRE CBD, and the fire's
  window-mean 1-h moisture from the run's weather (not the tables' `moist`, which averages the three dead
  classes). Term `cruz_active` (family over CAC 0.5, 1, 2): [CAC >= c]. Group "crownrate". T2 and T2o, A21's
  thinning, four fires, leave one out, planted check first.
- Predicted: not named on T2. In the window means East Troublesome has the lowest Cruz rate of the four, so the
  term cannot carry its miss; if anything it is named with a positive sign on T2o through Camp (fast crown
  ground the model overran). That would say the formula's inputs, not the formula, are where East Troublesome's
  miss lives.
- Part 2, the model setting C1 (pre-registered now; it runs whatever part 1 says, because it is a published
  alternative and not a fit): in each hour and canopy cell the crown rate is Cruz's instead of Rothermel 1991.
  Crown initiation stays Van Wagner's; once a cell torches, it is active at CROS_A when CAC >= 1 and passive at
  the larger of the surface rate and CROS_P otherwise. Intensity keeps the Scott and Reinhardt canopy fraction
  (cfb) for its canopy load, so flame and ember behaviour change only through the rate. Inputs out of Cruz's
  fitted range are run as they are and counted (share of torching cell-hours out of range), a label, not a switch.
  Runs: the four fires, E0 priors, HRRR as analysed, 3 seeds (Tubbs 2), CPU on THE-COW after InSAR and wind.
- Rule, fixed now: C1 is kept as a candidate only if, against E0, East Troublesome's area ratio moves toward 1
  at every mapped time and its A8 improves, and on Camp, Tubbs and Marshall A8 is not worse by more than 0.1 and
  the area ratio at the last mapped time does not move away from 1 by more than 0.1. If East Troublesome passes
  and the others fail, C1 is not kept alone, and is re-tested only together with a stopping rule that passed its
  own test (A24 or later), under a new amendment.
- Predicted: East Troublesome's ratio at 12.3 h rises from 0.29 to between 0.4 and 0.7; Camp overruns more
  (last ratio above 1.33), so C1 alone fails the rule.
