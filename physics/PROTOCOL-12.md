# Protocol 12: depth audit. Teaching the compressor the physics of depth on planted cases (FloodSim)

*Written 25 Nov 2026, before the full case set ran. Seen before writing: one quick smoke run of family C, four cases
(W 3 and 6 m, dx 1 and 2 m, S 0.01), used only to check that the harness runs. In it, the model matched 2D normal depth
at its own rasterised width, and was below the wall-friction answer at W 3 m. Families B and R had not run.*

**Why.** Kyle, 25 Nov 2026: "something's just missing with the flood model, or something's just wrong... if we're getting
an error for depth, then we should start planting like synthetic tests that are focused specifically on depth. Have our
compressor understand the physics behind depth." The real-storm wall check (CHI `FLOOD-WALLS-PLAN.md`) shows three depth
misses that protocol 11 could not name: too few marks.
- **Highland Falls marks 44719/44724:** 2.7 m too high in a narrow walled stretch.
- **Lumber River:** 0.6 to 1.3 m too low with a free downstream edge.
- **Hillside film:** wall depths near 0.16 m.

Here each miss gets a family of cases where the right depth is known exactly, so the compressor can learn the terms on
clean data before it is pointed at the field.

**Cases** (FloodSim `benchmarks/physics/depth/cases.mjs`, the fp64 CPU solver, same physics as the GPU):
- **C, walled channel.** Steady flow down a rectangular channel between walls. The walls are mask-0 cells, as buildings
  and walled banks are.
  - Settings: true width W ∈ {2.5, 3, 4, 5, 6, 8, 12, 20} m; cell size dx ∈ {0.5, 1, 2, 4} m; (S, n, Q) ∈ {(0.01, 0.04, 6),
    (0.003, 0.035, 15), (0.02, 0.05, 3)}.
  - A cell is channel when its centre is inside the channel, which is how FloodSim rasterises a footprint.
  - Answers: Manning normal depth with walls (R = A/P, the field truth), and without walls (R = h, the 2D ideal).
- **B, backwater.** Steady flow in a 20 m channel at dx 2 m, gentle slopes.
  - The downstream edge is held at a level above normal depth, or left as a free outfall.
  - Five (S, n, q, rise) settings.
  - Answer: the gradually-varied-flow profile up from the held level (direct step).
- **R, rain sheet flow.** Steady rain on a uniform plane.
  - Settings: dx ∈ {1, 2, 4} m; five (S, n, rain) settings.
  - Answer: kinematic-wave steady depth h(x) = (i·x·n/√S)^(3/5).

**Expected answers (fixed now):**
- **C.** The model's depth equals the 2D ideal at its own rasterised width Wm within 3 % wherever there are 2 or more cells
  across. Against the field truth it is off by two named causes:
  - rasterised width, ln(h_model/h_true) ≈ −0.6·ln(Wm/W);
  - missing wall friction, ≈ −0.4·ln(1 + 2h/W).
  With one cell across or less, a further grid term may appear.
- **B.**
  - Held edge: within 5 % of the profile everywhere past the inlet reach.
  - Free edge: low near the edge, the error fading upstream over a length that scales with the backwater length
    Lb = h_down/S.
- **R.** Within 15 % of the kinematic depth for x ≥ 20 m at every dx. Any steady bias is reported.

**Compressor** (`physlib.terms.FLOOD_DEPTH`, new; target y = ln(h_model / h_true) on grid 0.005):

| Term | Formula | What it would mean |
|---|---|---|
| `raster_width` | −0.6 · ln(Wm/W) (fixed) | the channel the model sees is wider or narrower than the real one |
| `wall_friction` | b · ln(1 + 2h/W) | walls drag on the flow and the 2D model has no wall friction |
| `cells_across` | b · 1/(Wm/dx) | too few cells across the channel |
| `edge_drawdown` | b · exp(−x_e / (λ·Lb)), λ 0.05 to 2 | a free edge drains the river within λ backwater lengths of it |
| `froude` | b · Fr | the miss grows with flow speed |
| `sheet_scale` | b · ln(x) | the sheet-flow depth grows along the slope with the wrong power |
| `slope_term` | b · ln S | the miss tracks slope |
| `cell_size` | b · ln dx | the miss tracks cell size |
| `roughness` | b · ln n | the miss tracks roughness |

Family by family, then pooled: `physlib.guards.analyse(..., bootstrap_B=30, null_N=19)`, leaving one setting out per
fold. **Instrument check:** the compressor must name `raster_width` and `wall_friction` on C, and `edge_drawdown` on free-edge B.
It must name nothing on held-edge B. If it fails these, it isn't fit to point at the field.

**Then the audit.** Each named term is turned into a prediction for the real misses, without fitting to them:
- Highland Falls: W, Wm and dx at the walled stretch, read from the run and lidar.
- Lumber: the distance of the marks to the edge, over Lb.
- Hillside film: the kinematic depth for the local slope and rain.

The predicted and actual misses go on one chart. A term that explains a real miss becomes a FloodSim fix candidate, kept only
if it is rule-based and improves the tuning towns.

**Limits.** Steady, prismatic, uniform roughness. Wall friction as Manning on R = A/P is the textbook form; real walls
vary. The compressor's result is a lead, never a pass.

**A1 (25 Nov 2026, after the first 13 family C cases, before the rest).** *Seen:* with a free outfall at the east end, W 2.5 m
and S 0.003 gave 2.03 m against a 2D ideal of 2.24 m at the model's own width. The backwater length there (h/S, ~750 m) is
far longer than the 160 m reach, so the free edge was drawing the measured section down. That is family B's question.
Family C now holds its east edge at the model-width normal depth (`cases.mjs`, `--free-edge` keeps the old run). Nothing
else in C changes. The 13 free-edge cases are kept as a log (`results/physics/flood_depth_C_free_partial.log`).

**A2 (25 Nov 2026, after the full C search; a reading of the result, not a rule change).** Under the rule as written
the C instrument check **fails**. The search named `raster_width` (coefficient 1.00 on the fixed −0.6·ln(Wm/W)), but
for wall drag it picked the grammar form 1/W·√h plus two small terms, not `wall_friction`. The physics is confirmed. On
the true wall depth, ln(1 + 2h_true/W) takes −0.399 (residual sd 0.001), the textbook 2/5. But `wall_friction` is
built on the model's depth (the only depth a field audit has), and in narrow deep channels that is half the true depth,
which bends the relation. A term defined from the model's depth by inverting the Manning relation is a new term and
would need its own planted check before use. Until then, C's wall-drag lead is read from the planted cases directly.
