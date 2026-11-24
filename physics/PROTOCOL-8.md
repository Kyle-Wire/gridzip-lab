# Protocol 8: the physics compressor on electric grid load (written 24 Nov 2026, before any demand value is fitted)

Ask (Kyle, 24 Nov, via the coordinator): could the compressor work for grid loads, and could it beat what grid
operators already have? Research only; nothing here goes on the website. Runs in the cloud container, not on THE-COW.

**Question.** A load forecaster's errors are not noise: they follow weather, calendar and behaviour the model leaves
out. Can physlib v2 (PROTOCOL-7 engine, unchanged) look at a standard load model's errors, name what is missing in
plain physics, stay quiet when nothing is, and does adding what it names actually improve the forecast on years it
never saw, against an expert's hand-built model and a machine-learning baseline?

**Data.** GEFCom2017 "defined data" track: ISO New England hourly demand (MW), dry-bulb and dew-point temperature (°F)
and a holiday flag for 8 load zones (ME, NH, VT, CT, RI, SEMASS, WCMASS, NEMASSBOST) and the system total, March 2003
to April 2017. Source: ISO-NE's public SMD hourly files, as packaged in `camroach87/gefcom2017data`
(`data/gefcom.rda`, fetched 24 Nov 2026 by `python -m physlib.gridload fetch`). Not committed.

**Seen before this was written:** the file's column names and its first five rows (CT, 1 March 2003, hours 0 to 4;
demand 3,157 to 3,386 MW). One timing run of the engine on pure Gaussian noise with the 2011-12 TOTAL weather inputs
(n = 17,536 hours, 15 s; no demand used). Nothing else. From memory (papers not opened for this work), the
load-forecasting literature reports that Hong's "vanilla" benchmark misses recency (lagged and moving-average
temperature), humidity and holidays; that is my prior, written here so a rediscovery can't be dressed up as a surprise.

## Code (`physlib/gridload.py`)
- **Vanilla benchmark (M0):** Hong's GEFCom2012 multiple-regression benchmark, OLS on MW:
  Trend + Month + Weekday×Hour + T + T² + T³ + T×Month + T²×Month + T³×Month + T×Hour + T²×Hour + T³×Hour.
- **Load term library (14 terms):** recency (max(0, T̄24 − 65 °F); max(0, 50 °F − T̄24); yesterday's max over 80 °F),
  3-hour thermal lag, heat-wave day count, humidity (max(0, T − 70)·RH; dew point over 60 °F), holiday, day next to a
  holiday, lighting after sunset (06-22 h), behind-the-meter solar (sun × years after 2010), weekend weather
  sensitivity, curved trend, and a shape-only daily harmonic. Grammar (PROTOCOL-7) on eight raw inputs: T, T_d, T̄24,
  yesterday's max, hour, day of year, trend, sun.
- **Residual:** r = ln(measured / vanilla) per hour, coded AR(1) (hours are serially correlated); grid 0.002.
  Every hour is read (step 1).

## A. Planted laws on real weather (synthetic outcome)
Inputs: TOTAL 2011-12 weather and calendar. Outcome: a planted law plus AR(1) noise (sd 0.03, ρ 0.9 per hour). Signal
sd = 0.5 × noise sd, split evenly over the planted terms. P1 recency (max(0, T̄24 − 65)); P2 humidity
(max(0, T − 70)·RH); P3 holiday dip; P4 behind-the-meter solar; P5 P1 + P2 + P3 together; P6 a shape the library
lacks, max(0, Tmax_yesterday − median) × hour (grammar must find it). 10 seeds each; N: pure AR(1) noise, 20 seeds.
- **Rule A:** each of P1-P6 is named (a winning term of the planted term's group; for P6 any grammar term) in ≥ 8 of
  10 seeds.
- **Rule A-null:** on pure noise it names nothing in ≥ 19 of 20 seeds.

## B. What does the vanilla benchmark miss? (exploration, real demand)
Vanilla fitted on 2004-2010; its errors on 2011-12 go to the compressor (every zone; TOTAL also gets a 9-shuffle
null and a 10-resample bootstrap). The winning formula is fitted on 2011-12 and must code 2013 shorter than a single
correction does.
- **Rule B:** on TOTAL the formula is not empty and transfers to 2013.
- Reading (not a gate): in how many of the 9 series does it name recency, humidity or holiday; anything else it names.

## C. Head-to-head on sealed years 2014-2016 (read once, after B is pushed)
Every model is fitted on 2004-2013 and scored on every hour of 2014-2016 with the actual weather (an ex-post
hindcast, as in GEFCom; operators forecast with forecast weather, so these errors are not comparable with an ISO's
published day-ahead accuracy).
- M0 vanilla. M1 vanilla + the terms part B named for that series (as extra regression columns, nothing else
  changed). M2 my expert model, fixed now: vanilla + the T, T², T³-by-hour blocks for T̄24 and for 1, 2 and 3-hour
  lags + RH, RH·T, RH·T² by month + holiday and day-next-to-holiday by hour. M3 gradient boosting
  (scikit-learn HistGradientBoosting on log MW, 18 inputs, 800 trees, rate 0.05, 63 leaves, seed 0, no tuning).
- Metric: MAPE over all hours; also over the top 1 % of hours (the peaks utilities plan for), and per year.
- **Rule C1:** M1 beats M0 on TOTAL and on at least 6 of the 8 zones.
- **Rule C2:** on TOTAL, M1 closes at least half of the gap M2 opens (gain M1 ≥ 0.5 × gain M2).
- Reading: where M1 and M2 stand against M3. No rule: the compressor is a diagnostic, not a forecaster.

Pass = A, A-null, B, C1, C2. Results in `results/gridload/`, written up in `physics/FINDINGS-8.md` with every
failure reported. If part A fails, part C still runs, and the write-up says the diagnosis is not trustworthy.
