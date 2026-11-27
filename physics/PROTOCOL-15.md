# Protocol 15: what flood depth still misses after the thin-film fix (FloodSim, tuning towns)

*Written 27 Nov 2026, before any fit below ran. Seen before writing: the A13 scores of both tuning towns (CHI
`docs/validation/walls-tuning/*-A13.json`): Lumber River marks 0.2 to 1.0 m low with the fix and the A11 edge, Highland
Falls good except marks 44719/44724 (+2.6/+2.7 m, a scorer reading). Also seen: the spill heights from the Lumber
river to each mark on the 600 m ground crop (CHI `runs-summary/walls/hf-grids/`, computed while checking the crop
loads): 15.7 to 16.3 m NAVD88, all below the held edge (16.62 m). No fit had been made. Test towns untouched.*

**Why.** Kyle, 27 Nov 2026: flood depth is the next focus; run the compressor again "if there's anything that we're
missing", then plan. No model runs.

**Data.** One row per wet scored mark, A13 runs (thin-film fix on): Lumber (edge 16.62 m) and Highland Falls. Context
from CHI `tools/gate2/walls_terms.py` (A11/A9 runs of the same marks; the context is terrain and distance, run-free
except `h_model`, which is taken from the A13 run's own level minus ground). New this round:
- `storage_fill`: b · d_ch/(d_ch + L), L 30 to 1000 m. A steady run cut short leaves the far floodplain below the
  river's level.
- `spill_rise`: b · (spill height from the river − the mark's ground). A ridge the water must climb (road, railway)
  with culverts the model lacks. Lumber only (no crop for Highland Falls' other marks); left out where blank.

**Rules (fixed now).**
1. Protocol 11's bar stands: the formal search runs at 20 or more wet scored marks. There are 15 (9 Lumber, 6 Highland
   Falls), so this round is **exploratory**: whatever wins is a lead for a planned run, never a finding.
2. Exploratory fit: physlib `guards.analyse` over FLOOD_WALL, folds leave-one-mark-out within Lumber, null N 19,
   bootstrap 30. A term is a lead only if it saves bits with null p ≤ 0.05 and bootstrap inclusion ≥ 0.6.
3. Highland Falls' 44719/44724 are reported apart (a scorer reading, CHI A12/A14) and the fit is also run without them.
4. Synthetic check first: a made-up town of 25 marks with `storage_fill` planted (−0.8 · d/(d + 150 m), noise 0.15 m)
   must name it; the null must name nothing.
5. Physics checks that need no fit, reported beside it: (a) spill rise above the river's model level (rules out
   blocked paths if every spill is below it); (b) the model's level at the marks and at the Nichols gauge against the
   held edge (a river lower than its own downstream edge means the reach was still filling).

**Expected (fixed now).** If the run was cut short, `storage_fill` with a negative slope leads and `spill_rise` does
not pay (spills sit below the river). If neither leads, the miss is a whole-town offset (inflow volume or edge level)
the marks cannot separate from `town_offset`.
