# Protocol 11: what FloodSim's misses at real high-water marks follow (flood wall check, tuning towns only)

*Written 24 Nov 2026, before any FloodSim run of a wall-check town was scored. Nothing has been read from the Highland
Falls run, which was still on the GPU at that time. No town under the new rule (CHI `FLOOD-WALLS-PLAN.md` A4) was placed.*

**Why.** Kyle, 24 Nov 2026: "when we run into issues, we should always try to update compressor to include formulas for
the specific domain, and then run it against the model." Flood's launch test is the wall check in Climate Home
Inspection (`docs/validation/FLOOD-WALLS-PLAN.md`, D65). Before any fix is hand-tuned on the tuning towns, this protocol
asks which physical term the model's errors at the marks follow. FloodSim has missed real marks by 0.61 to 0.85 m RMSE
(`FloodSim/docs/ACCEPTANCE.md`, all Ida or Charlotte). The street sandbox's inlet terms (`physlib.terms.FLOOD`) don't
apply here, so a new term set is added.

**Data.**
- One row per scored mark from CHI `tools/gate2/walls_terms.py`, run on THE-COW on each tuning-town run, and on
  spent cases when their marks are in the wall-check CSV form.
- **Never on a test town before its one scored run.** After that run, a test town's file may be read once, for the
  results page only. Nothing read from it goes back into FloodSim; the plan's section 5 "one look" rule stands.
- The target is rule D's error: model minus mark, in metres. That is a level for NAVD88 marks and a depth for
  height-above-ground marks.
- Marks the model left dry are not residuals. They are tallied apart by term value (median and range of HAND, distance
  to channel, distance to edge, building share and slope).

**Terms** (`physlib.terms.FLOOD_WALL`, fixed now; context from the run and its inputs only):

| Term | Formula | What it would mean |
|---|---|---|
| `volume_scale` | b · h | too much or too little water (rain total, losses) or roughness scaling every depth |
| `depth_power` | b · ln(h + 0.1) | shallow and deep water off by different factors |
| `floodplain_height` | b · HAND | the channel carries too little (no bathymetry under lidar water) or too much |
| `channel_distance` | b · ln(1 + d_ch / 10 m) | overbank exchange or floodplain storage wrong |
| `edge_loss` | b · exp(−d_edge / L), L 30 to 600 m | water leaves or piles up near the model edge (FloodSim's edge-loss flag, D54) |
| `blockage` | b · building share within 30 m | buildings change water level more than the model's walls do |
| `ground_slope` | b · S | friction or timing on steep ground |
| `velocity_head_wall` | b · V² / 2g | run-up on upstream faces (only if the run saved peak speed) |
| `datum_kind` | b · [height-above-ground mark] | the datum rule or the lidar ground reading is off |
| `on_building` | b · [mark on a building] | wall treatment, run-up or shelter at walls |
| `town_offset` | one offset per town (shape only) | a whole town off by one amount: rain product or inflow |

`floodplain_height` and `channel_distance` are one group (only one may enter).

**Search and readings, fixed now:**
- `physlib.guards.analyse(FLOOD_WALL terms present, ctx, y, Grid(0.01 m, −5, 5), bootstrap_B=30, null_N=19)`.
- Folds: leave one town out when there are two towns or more, else leave one waterbody out.
- The search runs once 20 or more wet scored marks are pooled over the tuning towns and spent cases. With fewer, only
  the residuals and the dry tally are charted, and the page says the search wasn't run for lack of power.
- **Named:** a term that wins, saves bits over the null search (p ≤ 0.05) and keeps its sign in every fold. A winning
  `town_offset` alone means a town-level cause (rain or inflow), and the next step is that town's rain product.
- **The lead goes to the tuning list** in the plan's section 4 as a fix to try, ranked by bits saved. A fix is kept only
  if the tuning towns' rule D numbers improve. The compressor's result is a lead, never a pass.
- **Charts:** each tuning town's residual against the winning term, and the dry tally.

**Validation before use** (`python physics/flood_walls.py synth`, 24 Nov 2026, 3 made-up towns, 75 marks, noise
0.15 m). Results are in `results/physics/flood_walls_synth.json`:
- planted edge loss (−0.4 · exp(−d_edge / 120 m)) → `edge_loss`, alone;
- planted volume scale (+0.25 · h) → `volume_scale`, alone;
- noise only → nothing named (null p = 1.0).

**Limits.** Marks are few: Highland Falls has 10, and a new tuning town is still to be placed. Power is low, and a
null here means "no term is big enough to see with these marks", not "no missing physics". Context comes from
the 1 m run grids, so the terms can't see sub-grid drains or culverts. HAND here is ground minus the 5th-percentile ground
within 300 m, not a flow-routed HAND.
