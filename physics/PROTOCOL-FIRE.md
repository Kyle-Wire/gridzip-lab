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
