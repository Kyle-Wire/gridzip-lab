# Protocol 10: does the flow over the roofs name the LBM's Case C pattern or its level? (written 24 Nov 2026, before either F3 step 6 field was read)

**Why.** Protocol 9 (FINDINGS-9) found that the LBM's lateral pattern on Case C 1D recurs. No near-ground flow
term names it. `downwash` × `lateral` adds to it, and the ledger shows the street flux below 0.1 m at 0.60 to 0.69 of
upstream. The step 5 level, ×1.29, is set by nothing below half the block height. WindSim F3 step 6
(`docs/validation.md` in WindSim, "F3 step 6") reruns Case C 1D at D/16 in two sections:
- R0, the standard section, blockage 5.6 %;
- R1, a 3.6 × 2.0 m section, blockage 1.7 %.

Both keep the mean fields to z = 0.45 m (2.25 D), now with u'w', u'v' and v'w'. Kyle's standing rule (24 Nov): when
a model fails a test, extend the compressor with that domain's formulas and run it against the model's errors. This
protocol is that step for the roof-level shear layer and the section.

**Seen before this was written:**
- the step 6 protocol and its readings script (WindSim 1705cba);
- that the CUDA self-test passed, that R0 started at 19:20 UTC, and R0's estimated finish time.

No R0 or R1 value has been read in this thread: no level, no point value and no field value.

**Data.**
- Fields: `f3s6_R0_fields.npz` and `f3s6_R1_fields.npz`. Each holds the mean u, v, w, the second moments
  (uu, vv, ww, uv, uw, vw), the mean |u| per cell to z = 0.45 m, and the solid mask, dx and origin.
- Points: each run's 120 points with the AIJ measurements, mirror fold as WindSim scores them.
- Targets on grid 0.01:
  - y0 = ln(model R0 / measured);
  - y1 = ln(model R1 / measured);
  - dy = y0 − y1, the per-point change from the wider section.

**New roof terms, computed per point by trilinear interpolation on the column above the point (fixed now).** U_app(z)
is mean u at the upstream column x = −0.8, y = 0. U_ref is U_app at z = D.
- `roof_speed`: mean u at z = 1.25 D over U_app(1.25 D). This is D1's "flow above the array runs fast", per point.
- `roof_stress`: −u'w' at z = D over U_ref². It is the momentum flux down into the street.
- `roof_k`: k = (u'² + v'² + w'²) / 2 at z = D, over U_ref².
- `roof_shear`: (u(1.25 D) − u(0.75 D)) / (0.5 D) × D / U_ref. It is the strength of the shear layer.
- `coupling`: mean u at 0.02 m over mean u at 1.25 D. It says how much of the roof flow reaches the probe.
- `ti_full`: the turbulence intensity at 0.02 m from the full stress tensor, sqrt(2k/3) / |U_h|. It replaces protocol
  9's isotropic form.

Where a column point falls inside a block (z ≤ D over a roof), the term is taken at the nearest street cell at the
same height, and the point is listed. These join protocol 9's field terms and FINDINGS-7's Case C inputs.

**Search, as in protocol 9:**
- `physlib.guards.analyse(..., grammar=True, bootstrap_B=30, null_N=19)`;
- folds: leave one quadrant out (4), plus upstream half → downstream half.

**Questions and rules (fixed now):**
- **Q1, the point pattern (R0, y0).** First, (A) FINDINGS-7 inputs plus protocol 9's terms. Then (B) with the roof
  terms added.
  - A roof term *names* the lateral pattern if, in B, all four hold: inclusion ≥ 0.8; picked in ≥ 60 % of
    per-input bootstrap resamples (`boot_inputs`, fixed now, which mends protocol 9's disclosed deviation); transfers
    on ≥ 3 of 4 quadrant folds; and `lateral`'s inclusion falls below 0.5.
  - If a roof term pays but `lateral` stays, it *adds to* the pattern.
  - Q1 is first checked on y0: does `lateral` recur at D/16 in R0 (inclusion ≥ 0.8 in A, transfer on ≥ 3 folds)?
    If not, Q1 is reported as "pattern absent at D/16", and B is still run.
- **Q2, what the wider section changes (dy).** Is dy a uniform shift (blockage acts on the level only), or does it
  carry structure?
  - Search dy with the same library and terms, and also each term's R0 − R1 difference.
  - If nothing is named (the null's 95th percentile), report it as "uniform shift" with the mean and sd of dy.
  - If a term is named under Q1's rules, report it as the blockage's footprint.
- **Q3, the level (described, not gated).** The step 5 ledger repeated per run: Q(x) below 0.1 m over Q₀, at the
  array's front, middle and back. Added to it:
  - the array-mean `roof_speed`;
  - the array-mean `roof_stress`;
  - the plan-mean w at z = D.

  Read them next to WindSim's blockage share s. This reads whether the section's effect goes through the air over the
  roofs (roof_speed falls from R0 to R1 together with the level) or through the street itself.
- Expected signs: `roof_speed`, `roof_stress` and `coupling` positive (more momentum brought down means the probe
  reads fast); `roof_k` and `ti_full` unknown.

**Synthetic check first,** on R0's real covariates with noise sd from run A's residual, 10 seeds each:
- (a) nothing planted: no roof term named in 9 of 10;
- (b) planted `roof_stress` at 0.02 to 0.10 sd;
- (c) planted `roof_speed` at 0.02 to 0.10 sd.

Power is the smallest sd named in 8 of 10. If (a) fails, Q1 and Q2 are reported as untrustworthy.

**Limits, stated now.**
- Two runs of one layout (1D), at D/16 only. The run-to-run noise at D/16 is unmeasured, so a small dy structure
  could be noise. dy is read with that caveat.
- The roof terms come from the model whose error is explained. A named term says where the model's flow differs. It
  does not show the tunnel's flow differs, because the AIJ data for Case C have no roof-level profiles over the
  streets (a check against D1 in Case B is a later step).
- Nothing here changes WindSim's scoring. A fix goes through WindSim's protocol and a held-out case.

Outputs:
- code: `physics/roof_fields.py`;
- figures: `physics/figures10.py` → `physics/figures/roof/`;
- results: `results/physics/roof_fields.json`;
- write-up: `physics/FINDINGS-10.md`, every failure reported.
