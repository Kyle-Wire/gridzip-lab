# Findings 14: the LBM's Case C errors with a stressed inflow (WindSim F3 step 7, 26 Nov 2026)

> Written and run as "Protocol 13" on branch claude/project-thread-fnbdyl (26 Nov 2026; protocol a5fd5b9, script 9983bea).
> Renumbered 14 on 26 Nov because main took protocol 13 (InSAR method 2) in the meantime; nothing else changed.

Protocol: `PROTOCOL-14.md` (a5fd5b9). The script `physics/inflow_fields.py` (9983bea) was written before any step 7 value
was read. Figure: `physics/figures14.py` → `physics/figures/inflow/inflow_01_step7.png`. Results:
`results/physics/inflow_fields.json`.

**Verdict.** The incoming wind's turbulence is not the cause of the ×1.25 level. With the guard in place the compressor
is now trustworthy: its synthetic check passes, silent in 10 of 10 noise runs. After the pattern that run A already
explains, it finds no roof-level term in C1's errors. Step 6's "no stress in the approach flow" lead was a reading
artefact, because the stress was read at x = −0.8, inside the inlet's adjustment zone.

## Q0, noise
- The repeat run N1 has the same level as R0: |Δ ln| = 0.0002.
- Per point, s_pt = 0.051. That is the floor for any per-point comparison at D/16.

## Guarded synthetic check (on the residual after run A, sd 0.102)
- Nothing planted: no roof term named in 10 of 10 runs. **Pass.**
- Power: `roof_stress` and `roof_speed` are both found in 10 of 10 runs at 0.06 sd, and in 2 or 3 of 10 at 0.04.
  Protocol 10 had no power for `roof_stress` up to 0.10. Searching the residual fixed that as well.

## Q1, C1's point pattern
- Run A: `[ti_local > 0.45] · lp`. It takes the sd from 0.172 to 0.102 and transfers on 4 of 4 folds.
- `lateral` inclusion is 0.62, under the 0.8 rule, so it does not recur at D/16. The same was true of R0 in protocol 10.
- Run B searches the residual with the roof terms. **It finds nothing**: the winner is empty, and no roof term reaches
  inclusion 0.03. Once the model's own near-ground turbulence pattern is accounted for, the flow over the roofs adds
  nothing a term of 0.06 sd or more could carry.

## Q2, what the inflow stress changes per point
- dC = ln R0 − ln C1: mean +0.009, sd 0.049. That is under s_pt (0.051), so it reads as **a level shift of +0.009
  in ln, within noise per point**.
- The search names `side_flux` and `d_coupling` at null p 0.05, but the sd rule fails, so these are not read.

## Q3, the ledger
| | R0 | C1 | W1 |
| --- | --- | --- | --- |
| Level | ×1.254 | ×1.243 | ×1.283 |
| Roof speed over array / upstream | 1.089 | 1.089 | 1.090 |
| −u'w' at D over the array / U_ref² | 0.0090 | 0.0098 | 0.0097 |
| Q below 0.1 m / Q₀ (front, middle, back) | 0.79, 0.67, 0.67 | 0.79, 0.66, 0.67 | 0.90, 0.86, 0.87 |

**The reading correction (figure, panel 2).** Along the centreline at roof height, −u'w' is negative near the inlet in
both runs, passes zero at about x = −0.8, and reaches 0.035 to 0.043 m²/s² by x = −0.6. The target is 0.027. The
imposed inlet correlation changes this curve only slightly. Both runs build the boundary layer's stress over the
fetch, and by the first blocks it is already at or above target. Step 6's ledger and approach check both used
x = −0.8. That point is where the inlet is still adjusting, so the "no stress" lead was an artefact of where it was
read. The 5 % approach mismatch between the sections may be the same artefact. A future blockage comparison should
use x ≥ −0.65.

## What it leaves
- Every suspect tested so far moves the level by under 2.5 % (ln 0.025). The one exception is the floor z0, which
  moved it 8 % at ×4.
- Dividing by the level, R0 hits 81 % of points, so the model's pattern is right and its level is not.
- The compressor has nothing more to name at the point level.
- The open physics is a bulk momentum budget: the array removes too little momentum overall. That is where any next
  physics term must act.

## Disclosures
- The `step7_readings.py` WindSim readings script was written after the runs, but applies the protocol's numbers
  exactly. Unlike step 6, it was not committed before the runs.
- The stress-along-x figure panel is post hoc, made after the ledger showed C1 ≈ R0 at x = −0.8.
