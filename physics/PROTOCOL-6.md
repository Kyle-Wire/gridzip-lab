# Protocol 6: storm physics for WindSim's thunderstorm gust misses (written 24 Nov 2026, before any storm weather was read)

Ask (Kyle, 23 Nov, via the WindSim thread): the thunderstorm gusts look noisy. Does the compressor consider what makes a
thunderstorm gust different from an ordinary one: the cold air behind the gust front, the temperature exchange across
the front, the storm's structure, the stability of the air?

**What the library has now** (protocol 3, `FINDINGS-3.md`): station siting (each station's own mean miss), station height,
wind direction, hour of day, month, roughness, distance to neighbours. None of it is weather inside the storm.
Protocol 3 found the thunderstorm misses are mostly per-station offsets, and nothing else paid.

**Seen before this was written:** protocol 3's results; the F2 pair files' fields (`observed_ms`, `predicted_ms`,
`storm_time`, `station`); one station-year header from the ISD bucket. No temperature, dewpoint or pressure of any
storm has been read.

**Data.** The 903 thunderstorm pairs of protocol 3 (WindSim `docs/validation/layer1/f2/{cal,blind}/*.json`, 13 places,
704 calibration and 199 blind). Target y = ln(observed / predicted), grid 0.005, as protocol 3. For each pair, the
station's own ISD hourly record for that year (NOAA, AWS Open Data, public domain), the same rows Layer 1 reads.

**Finding the gust in the record.** Parse the station-year with WindSim's `layer1.isd.parse_year`. Within ±3 h of
`storm_time` (the chain's thunderstorm window), take the row whose gust is within 0.3 m/s of `observed_ms`; if
several, the one nearest `storm_time`; if none, the largest gust in the window, flagged. Report the match rate.
Gust time t_g = that row's gust time.

**Storm terms, from the station's own record (step 1a)**, all from FM-15/FM-16 rows with quality codes 0, 1, 4, 5, 9 accepted:
- `cold_pool`: sqrt(ΔT⁺ / T₀), with ΔT = T₀ − min T over (t_g, t_g + 90 min], T₀ = max T over [t_g − 120 min, t_g]
  (kelvin). The density-current speed scales as sqrt(g h Δθ/θ), so the gust's excess over the chain should grow with it.
- `pressure_jump`: Δp⁺ in hPa, max station pressure (MA1, else SLP) over (t_g, t_g + 90 min] minus the last value at or
  before t_g − 30 min. Cold, dense outflow raises the surface pressure.
- `dew_drop`: Td at or before t_g − 30 min minus min Td over (t_g, t_g + 90 min], positive part (rain-cooled air).
- `dryness`: T₀ − Td₀ before the storm (K). A dry layer under the cloud evaporates rain and drives strong downdrafts
  (a proxy for DCAPE); dry microbursts hit hardest.
- `gust_front`: 1 if ΔT ≥ 3 K and Δp ≥ 1 hPa, else 0 (a gust front passed the station).
- `heat`: T₀ itself (°C, centred): hot afternoons and convective strength.
These join protocol 3's library for thunderstorms (`station_bias`, `elevation`, `direction`, `stability`, `season`),
so a storm term counts only if it pays beyond the station offsets.

**Search and readings, fixed now:**
- MDL beam search, up to 4 terms. Folds: leave one place out (13), and calibration → blind. A storm term **contributes**
  only if MDL chooses it and it transfers on calibration → blind and on at least 10 of 13 place folds.
- Report each storm term's single-term saving and its coefficient's sign. Sign expected from the physics: positive for
  `cold_pool`, `pressure_jump`, `dew_drop`, `dryness` and `gust_front` (a stronger outflow gives a bigger gust than the chain).
- **Range reading:** the thunderstorm range is the 10th-90th percentile of observed/predicted. If `gust_front` or another
  storm term contributes, compute the range separately for each half of that term (split at the calibration median, or
  by the tag) on the calibration pairs, apply it to the blind pairs, and report blind coverage against the 70-90 % mark,
  next to the current single range (67 % at the third blind test, `docs/validation.md`).
- Missing values (no pressure, no temperature after t_g): the term is 0 after centring for that pair, the count is reported,
  and a pair missing all storm terms is dropped from the search but counted.

**Synthetic check first** (real covariates, real per-station offsets of the size protocol 3 found, noise sd from the
protocol 3 bias-only fit): 10 seeds each of (a) station offsets only → no storm term named in 9 of 10; (b) + 0.10 per
unit sd of `cold_pool`; (c) + 0.10 · `gust_front`. Power is reported as protocol 5 A1 did: the smallest sd found in
8 of 10, over 0.02 to 0.10.

**Step 1b (ERA5: storm motion, CAPE), only by amendment:** ERA5 needs either a CDS key or the NSF-NCAR copy on AWS. If one
is reachable, the terms (CAPE, 0-6 km shear, 700 hPa wind speed along the gust direction) are added by an amendment
written before any ERA5 value is read. **Step 2 (NEXRAD, GOES)** only if a step 1 term contributes.

**Limits, stated now:** 903 gusts at 97 stations; one station per pair; hourly or special reports, so a fast temperature
drop can fall between rows; the chain's thunderstorm tag comes from the METAR text. Nothing here changes WindSim; a
finding goes there as a lead under its own protocol.
