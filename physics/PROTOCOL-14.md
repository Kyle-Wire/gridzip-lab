# Protocol 14: the LBM's Case C errors with an inflow that carries shear stress (WindSim F3 step 7), with a guard for roof_shear (written 26 Nov 2026, before any step 7 field or score was read)

> Written and run as "Protocol 13" on branch claude/project-thread-fnbdyl (26 Nov 2026; protocol a5fd5b9, script 9983bea).
> Renumbered 14 on 26 Nov because main took protocol 13 (InSAR method 2) in the meantime; nothing else changed.

**Why.** In protocol 10 (FINDINGS-10) the synthetic check failed. On pure noise built around run A's fitted
structure, `roof_shear` was named in 2 of 10 runs: it stood in for the `ti_local` · `lp` shape A had already found.
WindSim F3 step 7 (`docs/validation.md`, WindSim e36f445) reruns Case C 1D at D/16 with three changes:
- C1: an inflow that carries −u'w' = u*²(1 − z/δ);
- W1: the same in the wide section;
- N1: a noise run, step 6's R0 with seed 1.

Kyle's standing rule: extend the compressor with the domain's formulas and run it against the model's errors. This
protocol does that with a guard, and asks what the inflow change does point by point.

**Seen before this was written:**
- the step 7 protocol and code;
- that the self-test passed and that N1 started at 20:37 ET on 25 Nov.

No N1, T1, C1 or W1 value has been read: no level, no point and no field.

**Data.**
- `f3s7_C1_fields.npz` and `f3s7_W1_fields.npz`, with the points CSVs, from WindSim F3 step 7.
- Step 6's R0 (`f3s6_R0_*`) and N1's points.
- Targets on grid 0.01:
  - yC = ln(C1 / measured);
  - dC = ln(R0) − ln(C1), the per-point effect of the inflow stress;
  - dN = ln(R0) − ln(N1), the per-point noise.

**Terms.** These are protocol 9's field terms and protocol 10's roof terms (`physics/roof_fields.py`), unchanged and
computed on C1's field, plus FINDINGS-7's Case C inputs.

**The guard (fixed now).** Roof terms are judged on what is left after run A's structure, not alongside it.
- Run A: FINDINGS-7 inputs plus protocol 9's terms, searched on yC.
- Run B searches the residual yC − μA (A's fitted values) with the roof terms and grammar.
- The synthetic check is made the same way: pure noise at A's residual sd, searched with the roof library. It must
  name no roof term in ≥ 9 of 10 seeds (seeds 1300 to 1309).
  - If it fails again, the roof terms are dropped from WindSim's diagnosis for this project, and that is reported.
  - The power check plants `roof_stress` and `roof_speed` at 0.02 to 0.10 sd on that residual (seeds 1400+).

**Questions and rules (fixed now):**
- **Q0, noise.** s_pt = sd(dN). Every per-point reading below is compared with it.
- **Q1, the point pattern after the inflow fix.**
  - Does `lateral` recur in run A on yC? It recurs with inclusion ≥ 0.8 and transfer on ≥ 3 of 4 quadrant folds.
  - Does a roof term pay in run B? It pays with inclusion ≥ 0.8, per-input bootstrap ≥ 60 % and ≥ 3 of 4 folds,
    with the guard passed.
- **Q2, what the inflow stress changes.** Is dC uniform (a level shift) or structured?
  - It counts as structured only if sd(dC) > 1.5 × s_pt and the search on dC (all terms plus their C1 − R0
    differences) names a term under Q1's rules at null p ≤ 0.05.
  - Otherwise it is reported as a level shift of mean(dC).
- **Q3, the ledger (described).** Protocol 10's ledger for C1 and W1:
  - Q(x) below 0.1 m;
  - the array-mean roof speed and roof stress;
  - the plan-mean w at D;
  - upstream −u'w' / U_ref² at 0.5, 1, 1.25, 1.5 and 2 D, against the target u*²(1 − z/δ) / U_ref².

  This shows whether the stress reached the array and where the street speed changed.

**Limits.**
- One layout, D/16, and a single noise pair, so s_pt is one sample.
- The terms come from the model being explained.
- The target stress shape is assumed, not measured, for Case C.
- Nothing here changes WindSim's scoring.

Outputs:
- code: `physics/inflow_fields.py`, which reuses `roof_fields.py` and `lateral_fields.py`;
- figure: `physics/figures14.py` → `physics/figures/inflow/`;
- results: `results/physics/inflow_fields.json`;
- write-up: `physics/FINDINGS-14.md`.
