# Findings, protocol 15: flood depth after the thin-film fix (27 Nov 2026)

Exploratory. 15 wet scored marks (protocol 11 needs 20), so nothing here is a finding of a cause, only leads for
planned runs. Results `results/physics/flood_walls_a13*.json`, figure `physics/figures/p15-compressor.png`.

1. **Synthetic check (rule 4): failed at 25 marks.** A planted 0.8 m far-floodplain shortfall was not named at n = 25.
   At 50 and 100 marks it was named in 3 of 3 tries (as `channel_distance`, its same-group lookalike), null p 0.05
   (`flood_walls_a13_power.json`). So the compressor needs about 50 wet marks to name an effect this size. Fifteen
   cannot.
2. **Fit.** Lumber alone (10 marks) and both towns without Highland Falls 44719/44724 (14): nothing pays beyond a
   whole-town offset (Lumber −0.58 m). With the pair, a grammar term (velocity × nearness to the brook) takes the
   two +2.7 m marks: that is the known scorer reading (CHI A12/A14), not physics.
3. **Physics checks (rule 5).**
   - (a) The ridge water must cross from the Lumber to each mark is 15.7 to 16.3 m, below the held edge (16.62 m). No
     road or railway cuts the town off, and `spill_rise` has nothing to explain (as expected).
   - (b) The model's level at every mark (16.27 to 16.58 m) and at the Nichols gauge (16.55 m) is below its own held
     downstream edge (16.62 m). A river lower than its own downstream edge was still filling when the run stopped:
     A8's steady run is 8 h fill plus 4 h hold through about 18 km² of swamp. The far marks (320 to 370 m from the
     channel) miss most (−0.82 to −0.99 m), the near ones least (−0.22 to −0.64 m): what a cut-short fill looks like,
     though 10 marks cannot prove it (point 1).
   - (c) Even filled to the edge, the town would sit at 16.62 m, against marks of 16.79 to 17.25 m and an NWS crest at
     Nichols of 16.89 m. The A11 edge is itself 0.27 m below the one gauge inside the box.
4. **Leads for the plan (each a run with its rule written first):**
   - L1: run Lumber's steady flow until the level at the marks stops rising (a convergence stop, e.g. < 1 cm/h),
     not a fixed 12 h. First check, CPU only, from the A13 run's saved frames: were the marks still rising at 12 h?
   - L2: the downstream edge: A11's even fall over 83 km puts it 0.27 m under the Nichols crest. A test town with a
     gauge inside its box can check its edge against that gauge without touching the marks.
   - L3: more marks before any further compressor search. The three test towns cannot be used (one look). New
     tuning marks would come from spent STN events (reserves San Lorenzo and Neuse need an amendment to become
     tuning towns).
