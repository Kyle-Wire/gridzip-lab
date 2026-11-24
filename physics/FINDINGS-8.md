# Findings 8: the physics compressor on electric grid load (ISO New England, 24 Nov 2026)

Protocol: `PROTOCOL-8.md` (pushed at 5f2a5b5 before any demand was fitted; amendment A1 after parts A to C were read).
Code `physlib/gridload.py`, figures `physlib/figures8.py` -> `physics/figures/gridload/`, results `results/gridload/`.
Data: GEFCom2017 ISO-NE hourly demand, dry bulb, dew point, holidays (8 zones + total), not committed.

**Verdict: pass (A, A-null, B, C1, C2, E on the total).** The engine is trustworthy on load data and fixes a
standard model by itself. It does not beat an expert or machine learning on its own. As an auditor of strong
models it helps the system total on fresh months, but region by region the help is mixed.

## A. Planted laws on real weather
| Planted | Named (rule 8/10) | Extra terms |
| --- | --- | --- |
| P1 recency, P2 humid heat, P3 holiday, P4 rooftop solar, P5 three at once, P6 grammar-only shape | 10/10 each | 0 |
| Pure AR(1) noise (rule 19/20 silent) | 20/20 silent | 0 |

## B. What Hong's vanilla benchmark misses (errors 2011-12, about 17,500 hours per series)
Every series named five terms, all transferring to 2013 except SE Mass. Total: recency (T̄24 hinge) plus humidity
x sun, temperature x season and evening humidity shapes; 0.076 bits per hour saved against 0.016 at most in 9
shuffled-data searches (p = 0.1, the floor for 9 shuffles); residual sd 0.0146 -> 0.0137. Bootstrap: individual
shapes are unstable (the top one is chosen in 80 % of resamples, recency in 20 %), so read the inputs (yesterday's
heat, dew point, sun, season), not the exact formula. Named across the 9 series: recency or lag in 9, humidity (dew point) in 9,
sun or daylight in 9, holidays in 0.

## C. Head-to-head on sealed 2014-2016 (MAPE %, all hours / top 1 % hours)
| Series | M0 vanilla | M1 vanilla + found | M2 expert | M3 boosted trees |
| --- | --- | --- | --- | --- |
| TOTAL | 3.88 / 4.25 | 3.28 / 2.39 | 3.07 / 2.21 | 3.18 / 2.16 |
C1: M1 < M0 in TOTAL and 8 of 8 zones (ME and CT by under 0.05). C2: M1 closes 74 % of the expert gain on TOTAL
(5 added columns against the expert's 372). Per zone the share ranges from 1 % (CT) to 92 % (VT). Boosted trees beat
M1 in 8 of 9 series (not SE Mass.). Full table in `results/gridload/C_sealed.json`.

## D/E (A1). Auditing the strong models; fresh test January to April 2017
Expert errors still carry named structure: total 0.095 bits/hour (null max 0.016) with recency, 3-hour lag, dew
point and rooftop solar. Boosted-tree errors carry almost none on the total (0.001 bits/hour, below the shuffled
maximum 0.069, so this is a null audit); zones save 0.003 to 0.028.
| Fresh Jan-Apr 2017, MAPE % | M2 | M2 + audit | M3 | M3 + audit |
| --- | --- | --- | --- | --- |
| TOTAL | 4.02 | 3.70 | 5.87 | 5.48 |
| Better in (of 9 series) | | 3 | | 7 |
Rule E on the total passes. The expert's additions hurt in 6 zones, worst SE Mass. (+2.2 points), whose audit had
failed the 2013 transfer; the one boosted-tree audit that failed transfer (Boston) is also one of its two losses.
Post hoc on spent 2014-16: M2 + audit better in 4 of 9, M3 + audit in 8 of 9. No control was run for "any extra
features help boosted trees", so part of M3's gain may not be specific to these terms. The fresh window is 2,876
hours of winter and spring only.

## Disclosures
- Two crashes before any outcome was printed (null-search block format; NaN lags in the expert design) and one in E
  (a grammar term 1/trend is infinite at the start of the fit years; such terms are dropped, NH lost one term).
- Hindcast with actual weather; not comparable with operators' day-ahead accuracy. Data end in April 2017, before
  data-centre growth and most rooftop solar. The expert model is mine, from memory of the literature.
