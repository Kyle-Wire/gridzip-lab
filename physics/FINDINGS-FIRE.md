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
