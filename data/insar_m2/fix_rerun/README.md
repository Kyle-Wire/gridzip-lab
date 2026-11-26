# V2 Pawnee rerun with the bridged fix (insar_pipeline PR #22)

`v2_fix_points.npz` is V2's run (Pawnee, ascending track 34, method 2) with the bridged rule of PR #22
(`docs/BRIDGED_PLAN.md`), read at the same cells as `../v2_pawnee_points.npz` (`--points-from`), so the two line up
row for row. Its fields are the `data/insar_lv` fields (see `../../insar_lv/README.md`), and `export_fix.py` made it.

**How it was run** (26 Nov 2026, CPU): `reunwrap_linked` from MintPy onward on V2's existing unwrapped stack. No
unwrap config had changed, so the stack wasn't re-unwrapped, and MintPy found its inversion up to date and kept it. So
the rates are V2's own. What is new is the bridged layer, the bridged stated 1-sigma and the grading with bridged cells
left out. The grading before the fix is kept in the project as `validation/protocol_V2_V2_before_fix`.

## What it showed (recorded, not a finding of this folder)
- **Bridged cells:** 423,350 of the grid, with 58,541 of the 475,284 valid cells (12.3%) in 2 numbered regions.
  Their stated 1-sigma is 1.75 to 7.97 mm/yr by coherence class.
- **The OK02 corner:** all 77 cells within 200 m of OK02 are bridged. Its rate is unchanged (MintPy kept the
  inversion). Under the rule OK02 is untestable and is left out of the graded rows. Before the fix, OK02's five pairs
  were in error by -2.7, -4.9, -5.1, -9.3 and -10.5 mm/yr.

| row (V2 protocol) | before the fix (= every cell counted) | after: bridged left out (graded) |
|---|---|---|
| rate_rmse_dd | 5.70 mm/yr (CI 3.06-7.87), marginal | 4.83 mm/yr (CI 1.82-6.67), pass |
| pair_share | 8 of 15 (53%), fail | 6 of 10 (60%), fail |
| rate_bias_loo | +0.13, pass | +0.33, pass |
| uncertainty_calibration | 6 of 6, pass | 5 of 5, pass |
| noise_benchmark | worst 7.75 mm/yr, fail | 7.75 mm/yr, fail (nothing excluded, by the protocol) |
| route | fail | fail |

The pairs that still fail involve OK05 (57 epochs) and OK07: OK04-OK05 -7.9, OK05-OK06 +6.7, OK04-OK07 -6.6 and
OK03-OK05 -5.5. These stations are not bridged.
