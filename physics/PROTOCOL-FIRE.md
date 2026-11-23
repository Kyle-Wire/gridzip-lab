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
