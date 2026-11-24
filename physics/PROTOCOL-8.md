# Protocol 8: name the LBM's lateral pattern on AIJ Case C from the model's own near-ground flow (written 24 Nov 2026, before the field was read)

**Why.** FINDINGS-7 (D) found that the LBM's Case C error depends on distance from the centreline (`lateral`),
interacting with distance into the array, with inclusion odds ≈ 1. The shape was unstable, and the library has no
physics word for it. FINDINGS-5 left an unnamed point pattern that repeats from run to run (r 0.91). Protocol 1 and
FINDINGS-5 named the missing words: flow direction, recirculation, local turbulence. WindSim's F3 step 5 run
(`docs/validation.md`, "F3 step 5") wrote the model's mean field below 0.1 m. This protocol turns that field into
terms and asks whether one of them takes the place of `lateral`.

**Seen before this was written:**
- the run's headline numbers: hit rate 42.5 %, p50 1.296, r 0.937, overall level 1.289;
- the probe-check verdict (probe factor 16 % of the level, height slope 8 %, "ruled out");
- the file name and size of `fields_1D.npz`.

No value of the field, and no per-point model value of this run, has been read.

**Data.**
- WindSim F3 step 5: Case C 1D, D/32, walls off, z0 0.0003, L 0.5 D, Cs 0.12.
- `fields_1D.npz` holds the mean u, v, w, their squares and the mean |u| per cell below 0.1 m, plus the solid mask, dx and origin.
- The run's 120 points (model) and the AIJ measurements, with the mirror fold applied as WindSim scores them.
- Target: y = ln(model / measured) per point, on grid 0.01, as in FINDINGS-4 and FINDINGS-7.

**Field terms, computed at each point by trilinear interpolation at z = 0.02 m (fixed now):**
- `turn`: |sin θ|, where θ is the angle between the local mean horizontal vector and the approach direction (+x).
- `reverse`: 1 if the local mean u < 0 (recirculation at the point), else 0.
- `ti_local`: sqrt((u'² + v'² + w'²) / 3) / |mean horizontal velocity|, with the variances taken from the stored squares.
- `downwash`: mean w at z = 0.05 m over U_app (negative means air coming down into the street).
- `side_flux`: mean v · sign(y) at the point over U_app (positive means flow pushed outward, away from the centreline).
- `street_flux`: the run's speed-up through the point's street cross-section. This is the streamwise flux below
  0.1 m, integrated across that street (between the facing block faces, at the point's x for along-wind streets and at
  its y for cross streets). It is divided by the same street's open-ground flux upstream of the array (x < −2 D).
  It is the continuity (blockage) word.

These join FINDINGS-7's Case C inputs: `lp`, `fetch`, `x_pos`, `lateral`, `hc`, `corner`, `along`, `cross`.

**Search and readings, fixed now:**
- `physlib.guards.analyse(T.WIND + field terms, ctx, y, Grid(0.01, −4, 4), grammar=True, bootstrap_B=30, null_N=19)`.
- Folds: leave one quadrant out (4, split on the sign of x and y after the fold), plus upstream half → downstream half.
- Run twice: (A) FINDINGS-7 inputs only; (B) with the field terms added.
- **Q1, does the lateral pattern recur in this run?** Yes if `lateral` (alone or in any grammar shape) has inclusion
  ≥ 0.8 in A and transfers on at least 3 of the 4 quadrant folds.
- **Q2, does a field term name it?** Yes if, in B, a field term has inclusion ≥ 0.8, is picked in ≥ 60 % of bootstrap
  resamples and transfers on ≥ 3 of 4 quadrant folds, and `lateral`'s inclusion falls below 0.5. If a field term
  pays but `lateral` stays, it is reported as "adds to" and does not name the pattern.
- **Level ledger (described, not gated):** Q(x) is the streamwise flux below 0.1 m across the whole section, over
  the upstream Q₀, together with the plan-mean w at z = 0.1 m over the array. Q above 1 inside the array means air
  drawn down from above the half-height, which is D1's roof-height-mixing lead. Report Q at the array's front, middle
  and back.
- Each field term's single-term saving and coefficient sign are reported. Expected signs: `reverse`, `ti_local` and
  `turn` negative (the scalar probe and a thin wake read low where the model's flow is weak and turning) or
  unknown; `street_flux` and `downwash` (−w) positive.

**Synthetic check first** uses the real covariates with noise sd taken from run A's residual, 10 seeds each:
- (a) no planted term: no field term named in 9 of 10;
- (b) planted `side_flux` at 0.02 to 0.10 sd;
- (c) planted `street_flux` at 0.02 to 0.10 sd.

Power is the smallest sd named in 8 of 10, as in protocol 5 A1.

**Limits, stated now.**
- One run, one layout (1D), 120 points, so the folds are spatial, not across layouts.
- The covariances u'v' and so on are not stored, so `ti_local` is the isotropic form.
- Field terms come from the same model whose error is explained, so a field term that names the pattern says where
  the model's flow differs, not that the tunnel's flow is different.
- Nothing here changes WindSim's scoring. Any fix goes through WindSim's own protocol and a held-out case.
