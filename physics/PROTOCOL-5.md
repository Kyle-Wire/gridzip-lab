# Protocol 5: where WindSim's 1.3× street excess on AIJ Case C comes from (written 23 Nov 2026, before any fit)

Ask (Kyle, via the WindSim thread): both F3 candidates (inflow length 0.5 D, Smagorinsky Cs 0.17) left the LBM's
~1.3× excess on Case C 1D at D/32 where it was. Before trying more fixes, use the physics compressor to rank the
suspects that remain.

**Seen before this was written:** WindSim `docs/validation.md` (F3 steps, D1 table, F3 results at bf9e352); the
median model/measured ratio of each committed 1D run (D/16 1.41, D/16 tunnel z0 1.33, D/20 1.45, D/24 1.38,
D/32 1.30, D/32 walls 1.32, D/32 L0.5D 1.30, D/32 Cs 0.17 1.30); the model's upstream column (x = −0.8 m) at
z = 0.02 m and z = D in the three runs that wrote profiles (2.44 to 2.47 and 3.82 to 3.85 m/s). The points' x and y
values were listed. No residual has been fitted, mapped or correlated in this work.

**Data** (WindSim bf9e352, `docs/validation/`): `lbm_case_c_1D_resolution/{d16,d16_z0tunnel,d20,d24,d32_z0tunnel,
d32_walls,d32_L0.5D,d32_cs0.17}/points_1D.csv` (120 points each, z = 0.02 m), `profiles_1D.json` of the three
D/32 runs that have one, `rans_case_c_d20/results_points.csv` (RANS, layout 1D). Target r = ln(measured / model),
grid 0.01, as protocol 1. The measured approach speed at z = 0.02 m is 2.43 m/s (protocol 1).

## Suspects, each with the signature it would leave
| | Suspect | Signature | Test |
|---|---|---|---|
| S1 | Approach flow too fast near the ground (the inflow fit) | the model's upstream column above 2.43 m/s at 0.02 m | ledger: ln(U_model,up(0.02) / 2.43) |
| S2 | The 0.02 m evaluation height (between cell centres at D/32) | an error ≈ local slope d lnU / d ln z × the height error | ledger bound: the slope from the in-array columns, for a height error of half a cell (3.1 mm) |
| S3 | Domain blockage and the fixed-velocity top (the model's 1.8 × 1.2 m section; 5.6 % frontal blockage) | the flow over the array is fast against the model's own upstream column | ledger: over-array excess at z = 1.25 to 2 D from the D1 columns; term `lateral` (|y| / D: outer streets sit nearer the side boundaries) |
| S4 | A solver setting (resolution, ground z0, building walls, inflow length, Cs) | the level moves with that setting across runs | MDL on all 8 runs pooled (960 rows): `point_pattern` (the same point's mean r in the other runs) plus one term per setting: `res` (ln cells per D), `z0_tunnel`, `walls`, `L05`, `cs17`. Hold-out: leave one run out |
| S5 | A street-scale mechanism | r varies across the array in a nameable way | MDL on the mean r of the four D/32 runs per point (120): the physlib wind library on this geometry (`contrast`, `fetch_decay`, `fetch_linear`, `canyon_along`, `canyon_cross`, `corner`, `contrast_along`, `fetch_decay_along`) plus `lateral` and `near_wall` (distance to the nearest block face, in D). Hold-out: fit on y < 0, code y > 0, and the reverse |
| S6 | Something shared with RANS (the measurements, the probe, the geometry) | the point pattern (level removed) is the same in both engines | correlation of LBM D/32-mean and RANS D/20 point patterns, against the correlation the mirror-pair noise floor allows |
| S7 | Unnamed structure | neighbours predict the remainder | protocol 1's neighbour test on S5's remainder, 199 shuffles |

**Readings, fixed now:**
- A ledger suspect (S1 to S3) is **ruled out** if its share of ln 1.30 is under 10 %, **a part** if 10 to 50 %,
  **a main cause** if over 50 %. Shares are reported in log terms and sum with the rest as "made inside the array".
- A term (S3 `lateral`, S4, S5) **contributes** only if MDL chooses it and it transfers on every hold-out fold.
- S4: if `point_pattern` is chosen with coefficient near 1 and no setting term pays, the pattern and level are
  fixed by something every run shares, not by any setting tried.
- S6: a correlation above half of what the noise floor allows is called shared.
- The ranking lists the suspects by share (ledgers) and by bits saved (terms), with the ones ruled out at the end.

**Synthetic check first** (real geometry and the D/32 mean model speeds; noise at the measured mirror-pair sd,
0.063 in ln): 10 seeds each of (a) a uniform factor only, (b) uniform + 0.08 · `lateral`, (c) uniform + 0.10 ·
`near_wall` decay, and on the 8-run pool (d) a point pattern plus a resolution shift of 0.06 per doubling. Passes
if the planted terms are named in at least 8 of 10 and (a) names no term in at least 9 of 10. If it fails, the tool
is fixed and re-run before any real fit.

**Limits, stated now:** Case C 1D is calibration for every run here, so nothing found is validated; a fix found
this way goes through WindSim's own protocol and a held-out case. One layout, one direction, one height. The 8
runs share the inflow, domain, scoring and measurements, so S4 cannot separate those; the ledgers do that.

## A1 (23 Nov 2026, after the synthetic check, before any real fit)
**Seen:** the synthetic check. (a) uniform only: no term named in 10 of 10; (d) pool, point pattern + resolution:
10 of 10. (b) lateral 0.08: 5 of 10; (c) near_wall 0.10: 0 of 10, with no false term in any seed. The cause is the
planted size, not the tool: no point lies on a block face (the nearest is 0.25 D away), so (c)'s signal sd is
0.012 in ln, a fifth of the noise; (b)'s is 0.035. The check as written could not pass.
**Change:** the check now measures power instead. For `lateral`, `near_wall` (L = 0.25 D), `fetch_decay` (L = 1 D)
and `canyon_along`, the planted signal is scaled to sd 0.01, 0.015, 0.02, 0.03, 0.04 and 0.06 (ln), 10 seeds each,
and the **minimum detectable effect** (MDE) is the smallest sd named in at least 8 of 10. (a) and (d) stay as pass
tests. On the real data, "does not pay" for a term is read as "smaller than its MDE", reported with it. Nothing
else changes.
