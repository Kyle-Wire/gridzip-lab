"""Step 2 of the champion: compress, decode, compare, rebuild enwik8, and count every byte.

Compresses OUT/enwik8.word.hdr (from build_order.py) with lpaq1g at level 9 and the "full5" settings below, decodes
it, checks the decoded file byte for byte, rebuilds enwik8 from it (headers.rebuild with the order and side stream)
and checks that too. Prints one JSON line with the total:
  archive + map and region tree 6,494 B + permutation 18,750 B + header side stream
The map cost: 4,959 B for the grid map as counted in every result since round 1 (HANDOFF.md), plus 1,535 B for
the region tree's labels and split order (make_tree_map.py). The decoder is lpaq1g built from lpaq1g/lpaq1g.cpp.
Encode and decode take about 9 minutes each on a Ryzen 9 9900X (WSL), peak 1.5 GB.
Usage (repo root):  python champion/run.py path/to/lpaq1g data/enwik8 champion/out"""
import json, math, os, subprocess, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import headers

HERE = os.path.dirname(os.path.abspath(__file__))
MAP = os.path.join(HERE, "map_tree.bin")
MAP_COST = 4959 + 1535
FULL5 = {"GFLAGS": str(29399 | 8388608 | 1048576 | 67108864 | 134217728), "MIX2": "1", "MIXSEL": "1", "IMMUNE": "1",
         "MXBIAS": "1", "MXSETS": str(2 | 128 | 8192), "IMMTHR": "1", "IMMUNE2": "4", "IMMUNE3": "1", "IMMCONS": "2",
         "CTXMODE": "2"}
# GFLAGS bits: 29399 = the round-4 host (grid cells, spatial map, previous word, fields, page structure, ...),
# 8388608 cell contexts, 1048576 next cell, 67108864 cell trigram guess, 134217728 cell 4-gram guess.
# MXSETS 2|128|8192: mixer weight sets by byte context, region tree, and cell x previous cell's coarse region.
# IMMTHR/IMMUNE2/IMMUNE3/IMMCONS: immune memories keyed on the last four and two cells, with co-stimulation.


def main(exe, enwik8, out):
    src, side_p, order_p = (os.path.join(out, f) for f in ("enwik8.word.hdr", "enwik8.word.hdr.side", "order.wordchain.npy"))
    arc, back = os.path.join(out, "champion.lpq"), os.path.join(out, "champion.back")
    env = {**{k: v for k, v in os.environ.items() if k not in FULL5}, **FULL5, "GRIDMAP": MAP}
    t = time.time(); subprocess.run([exe, "9", src, arc], env=env, check=True, stdout=subprocess.DEVNULL); enc = time.time() - t
    t = time.time(); subprocess.run([exe, "d", arc, back], env=env, check=True, stdout=subprocess.DEVNULL); dec = time.time() - t
    decoded = open(back, "rb").read(); ok = decoded == open(src, "rb").read()
    order = np.load(order_p); side = open(side_p, "rb").read()
    rebuilt = ok and headers.rebuild(decoded, order, side) == open(enwik8, "rb").read()
    os.remove(back)
    perm = math.ceil(math.lgamma(len(order) + 1) / math.log(2) / 8)
    size = os.path.getsize(arc)
    print(json.dumps({"archive": size, "map_and_tree": MAP_COST, "permutation": perm, "side_stream": len(side),
                      "total": size + MAP_COST + perm + len(side), "decoded_ok": ok, "enwik8_rebuilt_ok": bool(rebuilt),
                      "enc_s": round(enc), "dec_s": round(dec)}))


if __name__ == "__main__":
    main(*sys.argv[1:4])
