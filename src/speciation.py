"""Speciation: an island model of evolution.

Three populations evolve apart, each in its own habitat (a 1 MB slice of
enwik8 that is prose-heavy, markup-heavy or table-heavy). Every 3 generations
the best organism migrates to the next island (gene flow). At the end:
  * a specialisation matrix: each island's champion scored in every habitat
    (a strong diagonal means the species adapted to their niches)
  * a hybrid: genes from all champions merged greedily, keeping a gene only if
    it helps on the general training slice
All scored on a held-out general 2 MB slice as well.
"""
import json, os, random, sys
os.environ.setdefault("FEATS", "15")
import evolve as E

LP = E.LP
HAB = {k: f"{LP}/hab_{k}" for k in ["prose", "markup", "table"]}
POP, GENS, MIG = 12, int(os.environ.get("SGENS", 10)), 3
rng = E.rng


def evolve_island(pop, fit, habitat):
    order = sorted(range(len(pop)), key=lambda i: fit[i])
    pop = [pop[i] for i in order]; fit = [fit[i] for i in order]
    def pick():
        c = rng.sample(range(len(pop)), 3); return pop[min(c, key=lambda i: fit[i])]
    kids = [E.mutate(E.crossover(pick(), pick()) if rng.random() < 0.6 else list(pick())) for _ in range(POP - 3)]
    with E.ThreadPoolExecutor(E.WORKERS) as ex:
        kf = list(ex.map(lambda g: E.fitness(g, habitat), kids))
    return pop[:3] + kids, fit[:3] + kf


if __name__ == "__main__":
    out = sys.argv[1]
    names = list(HAB)
    islands = {}
    for k in names:
        pop = [E.HAND] + [E.mutate(E.HAND) for _ in range(5)] + [[E.rand_gene() for _ in range(rng.randint(5, 12))] for _ in range(POP - 6)]
        with E.ThreadPoolExecutor(E.WORKERS) as ex:
            fit = list(ex.map(lambda g: E.fitness(g, HAB[k]), pop))
        islands[k] = [pop, fit]
    log = dict(gens=[])
    for gen in range(GENS):
        for k in names:
            islands[k][0], islands[k][1] = evolve_island(islands[k][0], islands[k][1], HAB[k])
        rec = {k: min(islands[k][1]) for k in names}; rec["gen"] = gen
        log["gens"].append(rec); print(json.dumps(rec), flush=True)
        if gen % MIG == MIG - 1:  # migration around the ring
            bests = {k: islands[k][0][min(range(POP), key=lambda i: islands[k][1][i])] for k in names}
            for a, b in zip(names, names[1:] + names[:1]):
                islands[b][0][-1] = list(bests[a]); islands[b][1][-1] = E.fitness(bests[a], HAB[b])
        json.dump(log, open(out, "w"), indent=1)
    champs = {k: islands[k][0][min(range(POP), key=lambda i: islands[k][1][i])] for k in names}
    matrix = {k: {h: E.fitness(g, HAB[h]) for h in names} for k, g in champs.items()}
    matrix["hand"] = {h: E.fitness(E.HAND, HAB[h]) for h in names}
    # hybrid: greedy gene merge, judged on the general training slice
    best_k = min(names, key=lambda k: E.fitness(champs[k], E.TRAIN))
    hybrid = list(champs[best_k]); hf = E.fitness(hybrid, E.TRAIN)
    pool = [g for k in names if k != best_k for g in champs[k] if g not in hybrid]
    for g in pool:
        trial = (hybrid + [g]) if len(hybrid) < E.MAXG else None
        if trial is None:
            # swap with the gene whose removal hurts least
            trial = hybrid[:-1] + [g]
        f = E.fitness(trial, E.TRAIN)
        if f < hf: hybrid, hf = trial, f
    log.update(matrix=matrix, champs={k: E.gstr(g) for k, g in champs.items()}, hybrid=E.gstr(hybrid),
               hybrid_train=hf, hybrid_hold=E.fitness(hybrid, E.HOLD),
               hand_train=E.fitness(E.HAND, E.TRAIN), hand_hold=E.fitness(E.HAND, E.HOLD),
               champs_hold={k: E.fitness(g, E.HOLD) for k, g in champs.items()})
    json.dump(log, open(out, "w"), indent=1)
    print(json.dumps({k: v for k, v in log.items() if k not in ("gens", "champs")}, indent=1), flush=True)
