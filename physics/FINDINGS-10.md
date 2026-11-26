# Findings 10: the flow over the roofs and the wider section, read by the compressor (WindSim F3 step 6, 25 Nov 2026)

Protocol: `PROTOCOL-10.md` (4080861, pushed before either field was read). Code `physics/roof_fields.py`, figure
`physics/figures10.py` → `physics/figures/roof/roof_01_step6.png`, results `results/physics/roof_fields.json`.
Data: WindSim F3 step 6 (bf85eca, CUDA, Case C 1D, D/16, L 0.5 D, z0 0.0003, Cs 0.12). R0 is the standard section;
R1 is 3.6 × 2.0 m. Fields `f3s6_R{0,1}_fields.npz`, points `f3s6_R{0,1}_points.csv`.

**Verdict: no trustworthy answer.** The synthetic check failed. With nothing planted, a roof term was named in 2 of 10
runs, against a limit of 1, and both false picks were `roof_shear`, the very term the real search then named. By the
protocol, Q1 and Q2 are reported as untrustworthy. WindSim's own approach check also failed (5.1 % RMS against 3 %),
so neither suspect, blockage or the roof shear layer, is settled by this step.

## Synthetic check (on R0's real covariates, noise sd 0.109)
- (a) Nothing planted: no roof term in 8 of 10 (rule: 9). **Fail.** The false picks were `[roof_shear > 0.34] · lp`
  and `roof_shear`.
- Power: `roof_speed` is found in 8 of 10 at sd 0.10, and not reliably below that. `roof_stress` is never found
  reliably up to 0.10 (at most 2 of 10). The library cannot see a roof-stress effect of the size that matters here.

## Q1, the point pattern (R0, untrustworthy)
- At D/16 the lateral pattern does not recur: `lateral` inclusion is 0.34 (rule 0.8). Run A's winner is
  `[ti_local > 0.41] · lp`, which takes the sd from 0.203 to 0.109, transfers on 4 of 4 quadrant folds, and sits at
  the null floor (p 0.05).
- Run B's winner is `e^(−roof_shear/0.52) · roof_k`: sd 0.102, inclusion 1.00 for both inputs, bootstrap 83 % and
  77 %, and it transfers on 4 of 4 folds. On the protocol's words this "names" the pattern. Given the false-alarm test,
  it is not read as a finding.

## Q2, what the wider section changes (untrustworthy)
- dy = ln(R0) − ln(R1) per point: mean −0.012, sd 0.079. The level hardly moves, but the points rearrange (figure
  panel 2). The along-wind streets read higher in the wider section and the cross streets lower, in a smooth and
  symmetric pattern.
- The search names `d_coupling` (the change in near-ground over roof-level speed) plus `roof_shear` and `fetch`. It
  takes the sd from 0.078 to 0.021, at null p 0.05, transferring on all folds. This is close to circular, because
  coupling's numerator is the model's own near-ground speed at the point. The run-to-run noise at D/16 is unmeasured,
  so even the pattern's size cannot yet be called a section effect.

## Q3, the level (described)
| | R0 | R1 |
| --- | --- | --- |
| Level (model / measured) | ×1.254 | ×1.269 |
| Roof speed over the array / upstream at 1.25 D | 1.089 | 1.091 |
| −u'w' at D over the array / U_ref² | 0.0090 | 0.0098 |
| Plan-mean w at D / U_app | +0.065 | +0.069 |
| Upstream −u'w' / U_ref² at 0.5 D, D, 1.25 D | 0.0013, 0.0007, −0.0001 | 0.0022, 0.0034, 0.0002 |

The flow over the roofs runs about 9 % fast in both sections, and the level does not follow the section. One lead
comes out of the ledger, not out of the compressor. The approach flow carries almost no turbulent shear stress above
half the block height in R0, and a different profile of it in R1 (figure panel 1). That is not the near-constant stress
layer of a developed boundary layer. The inflow, not the section, is the next thing to pin down.

## Next (WindSim, a new protocol before any run)
1. Run-to-run noise at D/16: R0 repeated with a different averaging window (about 17 min GPU), so that dy and Δ
   level can be judged.
2. An inflow check against the AIJ Case C approach profiles (U and turbulence), and R1 rerun with its inflow
   matched to R0's within 3 %.
3. In the compressor: the roof terms stay in the library, but `roof_shear` needs a guard. A term that is named on pure
   noise 2 times in 10 is too willing to stand in for `lp` and `ti_local` structure.

## Disclosures
- The script ran search B and the Q2 search even though the synthetic check had failed. The protocol said to run
  them and report them as untrustworthy, which is what is done here. Protocol 9's script had skipped B in that case.
- Field values were first read by this script. The step 6 readings (WindSim 1705cba) had been seen before.
