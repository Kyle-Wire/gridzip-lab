"""Test 1 of surprise/PROTOCOL.md: the synthetic street with three planted lots (FloodSim benchmarks/emulator/street.mjs).
Usage: python surprise/street_test.py <street.json> <out.json>"""
import json, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from surprise import core

d = json.load(open(sys.argv[1])); L = d["lots"]; P = d["planted"]
y = np.log10(np.array(d["ring"]).T + 0.03)                      # lots x storms
xy = np.column_stack([np.arange(L) * d["W"] + d["W"] / 2, np.full(L, 23.0)])
r = core.score(y, xy, 40.0, 0.02, -2.0, 1.0, descriptor=np.array(d["ringGround"]), verify=True)
s, sg = r["surprise"], r["surprise_given"]; order = list(np.argsort(-s))
rank = {q: order.index(q) + 1 for q in range(L)}
a, b, c = P["low_spot"], P["neighbour_fill"], P["low_lot"]
p1 = rank[a] <= 3 and rank[b] <= 3
p2 = rank[c] <= 4 and sg[c] <= 0.5 * s[c]
p3 = all(s[q] <= min(s[a], s[b]) for q in range(L) if q not in (a, b, c))
out = {"surprise": s.tolist(), "surprise_given_height": sg.tolist(), "rank": rank, "h": r["h"], "tau": r["tau"], "roundtrip_ok": r["roundtrip_ok"],
       "ring_cm": (np.array(d["ring"]).T * 100).round(2).tolist(), "P1": p1, "P2": p2, "P3": p3}
json.dump(out, open(sys.argv[2], "w"), indent=1)
for q in range(L): print(f"lot {q:2d} rank {rank[q]:2d} surprise {s[q]:6.2f} given height {sg[q]:6.2f} ring cm {np.round(np.array(d['ring'])[:, q]*100, 1)}")
print("P1", p1, "P2", p2, "P3", p3)
