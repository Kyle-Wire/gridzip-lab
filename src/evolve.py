"""Survival of the fittest for compressor designs.

An organism is a genome: a list of genes, each gene one context for lpaq1g
("fmask:bmask", see lpaq1g.cpp). Fitness = compressed size of a 2 MB training
slice (smaller is fitter). Each generation keeps the elite, then fills the
population with offspring by tournament selection, crossover (genes drawn from
both parents) and mutation (flip a byte, flip a feature, add, drop or
duplicate a gene). The best genome of each generation is also scored on a
held-out 2 MB slice it never evolved on, to watch for overfitting.

Seeds: lpaq1's own contexts, our hand-designed best, and random genomes.
"""
import json, os, random, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

LP = os.environ["LPAQ_DIR"]
EXE = f"{LP}/lpaq1gx"
TRAIN, HOLD = f"{LP}/e2", f"{LP}/h2"
MAP = f"{LP}/map10t16.bin"
POP, GENS, ELITE, MAXG, WORKERS = 24, int(os.environ.get("GENS", 20)), 4, 17, 2
FEATS = 13
rng = random.Random(int(os.environ.get("SEED", 0)))

LPAQ = [(0, 0x3), (0, 0x7), (0, 0xF), (0, 0x3F), (0x80, 0)]
HAND = LPAQ + [(0x1, 0), (0x3, 0), (0x1, 0x7), (0x1C, 0), (0x19, 0), (0x70, 0x1),
               (0x800, 0x7), (0x800, 0x3F), (0x400, 0x3), (0x480, 0), (0x200, 0x3), (0x1280, 0)]


def gstr(g):
    return ",".join(f"{f:x}:{b:x}" for f, b in g)


def fitness(g, data):
    with tempfile.NamedTemporaryFile(dir=LP, delete=True) as tmp:
        env = dict(os.environ, GENOME=gstr(g), MIX2="1", MIXSEL="1", GRIDMAP=MAP, GFLAGS="31")
        subprocess.run([EXE, "9", data, tmp.name], env=env, stdout=subprocess.DEVNULL, check=True)
        return os.path.getsize(tmp.name)


def rand_gene():
    f = 0
    for i in range(FEATS):
        if rng.random() < 0.12: f |= 1 << i
    b = rng.getrandbits(8) & rng.choice([0x3, 0x7, 0xF, 0x3F, 0xFF, 0x5, 0xA])
    return (f, b) if (f or b) else (0, 0x3)


def mutate(g):
    g = list(g)
    for _ in range(rng.choice([1, 1, 2, 3])):
        op = rng.random()
        if op < 0.3 and g:
            i = rng.randrange(len(g)); f, b = g[i]; g[i] = (f, b ^ (1 << rng.randrange(8)))
        elif op < 0.6 and g:
            i = rng.randrange(len(g)); f, b = g[i]; g[i] = (f ^ (1 << rng.randrange(FEATS)), b)
        elif op < 0.75 and len(g) < MAXG:
            g.append(rand_gene())
        elif op < 0.9 and len(g) > 3:
            g.pop(rng.randrange(len(g)))
        elif g and len(g) < MAXG:
            f, b = rng.choice(g); g.append((f ^ (1 << rng.randrange(FEATS)), b))
    g = [x if (x[0] or x[1]) else (0, 0x1) for x in g]
    return list(dict.fromkeys(g))[:MAXG]


def crossover(a, b):
    pool = list(dict.fromkeys(a + b))
    k = rng.randint(min(len(a), len(b)), min(MAXG, max(len(a), len(b))))
    keep = [x for x in pool if rng.random() < 0.6]
    return (keep + [x for x in pool if x not in keep])[:k]


def evaluate(pop):
    with ThreadPoolExecutor(WORKERS) as ex:
        return list(ex.map(lambda g: fitness(g, TRAIN), pop))


if __name__ == "__main__":
    out = sys.argv[1]
    pop = [LPAQ, HAND] + [mutate(HAND) for _ in range(8)] + \
          [[rand_gene() for _ in range(rng.randint(5, 12))] for _ in range(POP - 10)]
    fit = evaluate(pop)
    base = dict(lpaq_like=fit[0], hand=fit[1],
                lpaq_like_hold=fitness(LPAQ, HOLD), hand_hold=fitness(HAND, HOLD))
    log = dict(base=base, gens=[])
    print(json.dumps(base), flush=True)
    for gen in range(GENS):
        order = sorted(range(len(pop)), key=lambda i: fit[i])
        pop = [pop[i] for i in order]; fit = [fit[i] for i in order]
        best = pop[0]
        rec = dict(gen=gen, best=fit[0], median=sorted(fit)[len(fit) // 2],
                   best_hold=fitness(best, HOLD), genes=len(best), genome=gstr(best))
        log["gens"].append(rec); json.dump(log, open(out, "w"), indent=1)
        print(json.dumps({k: v for k, v in rec.items() if k != "genome"}), flush=True)
        def pick():
            c = rng.sample(range(len(pop)), 3); return pop[min(c, key=lambda i: fit[i])]
        kids = []
        while len(kids) < POP - ELITE:
            child = crossover(pick(), pick()) if rng.random() < 0.6 else list(pick())
            kids.append(mutate(child))
        kf = evaluate(kids)
        pop = pop[:ELITE] + kids; fit = fit[:ELITE] + kf
    print("final", gstr(pop[min(range(len(pop)), key=lambda i: fit[i])]), flush=True)
