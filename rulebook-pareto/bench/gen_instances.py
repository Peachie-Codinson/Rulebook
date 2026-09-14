#!/usr/bin/env python3
"""Generate the instance suite: graphs, rulebooks and queries.

Graph file    : "n m N" header, then m lines "u v c_0 .. c_{N-1}".
Rulebook file : same layout as the reference repo's rules files.
Query file    : one "s t" pair per line.
"""
import argparse
import os
import random

# Rulebooks used in the study. Each is (num_rules, equivalence classes,
# strict-priority relations "a is more important than b").
#
# chain / top / diamond mirror the hierarchies in the RA*pex paper; flat is the
# pure-Pareto control; wide is the case both new solvers should struggle with.
RULEBOOKS = {
    # r0 > r1 > r2 > r3 -- a total order, so purely lexicographic
    "chain4": (4, [], [(0, 1), (1, 2), (2, 3)]),
    # r0 > r1, r2, r3 -- the paper's 4-objective BAY hierarchy (r4 > r1,r2,r3)
    "top1": (4, [], [(0, 1), (0, 2), (0, 3)]),
    # r0 ~ r1 on top, r2 below r0, r3 below r1 -- the paper's Fig. 3 preorder
    "diamond4": (4, [[0, 1]], [(0, 2), (1, 3)]),
    # r0 > r1 > r2,r3 -- two peelable layers then an incomparable pair
    "chain2top": (4, [], [(0, 1), (1, 2), (1, 3)]),
    # no relations at all: plain Pareto dominance
    "flat4": (4, [], []),
    # 3-objective version of the paper's (S2) setting: r0 > r1, r0 > r2
    "top1x3": (3, [], [(0, 1), (0, 2)]),
    "flat3": (3, [], []),
}


def write_rulebook(path, name, eps):
    n, classes, rels = RULEBOOKS[name]
    covered = {r for c in classes for r in c}
    all_classes = [list(c) for c in classes] + [[r] for r in range(n) if r not in covered]
    with open(path, "w") as f:
        f.write(f"{n}\n")
        f.write(" ".join(str(eps) for _ in range(n)) + "\n")
        f.write(f"{len(all_classes)}\n")
        for c in all_classes:
            f.write(f"{len(c)} " + " ".join(map(str, c)) + "\n")
        f.write(f"{len(rels)}\n")
        for a, b in rels:
            f.write(f"{a} {b}\n")
    return n


def gen_random_graph(n, avg_degree, num_rules, cmax, rng):
    """Directed random graph with a Hamiltonian backbone so every query is
    reachable, plus random chords. Costs are drawn independently per objective,
    which is the standard way to get a large Pareto frontier."""
    edges = []
    perm = list(range(n))
    rng.shuffle(perm)

    def cost():
        return [rng.randint(1, cmax) for _ in range(num_rules)]

    for i in range(n):
        u, v = perm[i], perm[(i + 1) % n]
        edges.append((u, v, cost()))
    extra = max(0, int(n * avg_degree) - n)
    seen = {(u, v) for u, v, _ in edges}
    while len(edges) - n < extra:
        u = rng.randrange(n)
        v = rng.randrange(n)
        if u == v or (u, v) in seen:
            continue
        seen.add((u, v))
        edges.append((u, v, cost()))
    return edges


def gen_layered_graph(layers, width, num_rules, cmax, rng):
    """Layered DAG: `layers` ranks of `width` vertices, complete bipartite
    between consecutive ranks, plus a single source and sink. Path length is
    fixed at `layers+1`, so the frontier grows with depth in a controlled way --
    the standard hard construction for multi-objective shortest paths."""
    n = layers * width + 2
    src, dst = n - 2, n - 1
    edges = []

    def vid(l, i):
        return l * width + i

    def cost():
        return [rng.randint(1, cmax) for _ in range(num_rules)]

    for i in range(width):
        edges.append((src, vid(0, i), cost()))
        edges.append((vid(layers - 1, i), dst, cost()))
    for l in range(layers - 1):
        for i in range(width):
            for j in range(width):
                edges.append((vid(l, i), vid(l + 1, j), cost()))
    return n, edges, src, dst


def gen_grid_graph(k, num_rules, cmax, rng):
    """4-connected k x k grid, edges in both directions -- the robot-navigation
    shape that motivates the rulebook formalism."""
    edges = []

    def vid(r, c):
        return r * k + c

    for r in range(k):
        for c in range(k):
            for dr, dc in ((0, 1), (1, 0), (0, -1), (-1, 0)):
                nr, nc = r + dr, c + dc
                if 0 <= nr < k and 0 <= nc < k:
                    edges.append((vid(r, c), vid(nr, nc),
                                  [rng.randint(1, cmax) for _ in range(num_rules)]))
    return edges


def write_graph(path, n, edges, num_rules):
    with open(path, "w") as f:
        f.write(f"{n} {len(edges)} {num_rules}\n")
        for u, v, c in edges:
            f.write(f"{u} {v} " + " ".join(map(str, c)) + "\n")


def write_queries(path, n, count, rng, fixed=None):
    with open(path, "w") as f:
        if fixed:
            for s, t in fixed:
                f.write(f"{s} {t}\n")
            return
        for _ in range(count):
            s = rng.randrange(n)
            t = rng.randrange(n)
            while t == s:
                t = rng.randrange(n)
            f.write(f"{s} {t}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="instances")
    ap.add_argument("--seed", type=int, default=20260914)
    ap.add_argument("--queries", type=int, default=5)
    ap.add_argument("--cmax", type=int, default=100)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    rng = random.Random(args.seed)

    # Rulebooks (eps is overridden on the command line per run, so 0 here).
    for name in RULEBOOKS:
        write_rulebook(os.path.join(args.out, f"rb_{name}.txt"), name, 0)

    manifest = []

    # Random graphs: short paths, small frontiers -- the easy end of the ladder.
    for num_rules in (3, 4):
        for n in (200, 800, 3200):
            name = f"rand_n{n}_N{num_rules}"
            edges = gen_random_graph(n, 4, num_rules, args.cmax, rng)
            write_graph(os.path.join(args.out, f"{name}.graph"), n, edges, num_rules)
            write_queries(os.path.join(args.out, f"{name}.query"), n, args.queries, rng)
            manifest.append((name, n, len(edges), num_rules, "random"))

    # Grids: the navigation-shaped case, and a clean difficulty ladder --
    # frontier size grows steeply with k because path length grows with k.
    for num_rules in (3, 4):
        for k in (6, 8, 10, 12, 14, 16, 20, 25, 30):
            n = k * k
            name = f"grid_k{k:02d}_N{num_rules}"
            edges = gen_grid_graph(k, num_rules, args.cmax, rng)
            write_graph(os.path.join(args.out, f"{name}.graph"), n, edges, num_rules)
            write_queries(os.path.join(args.out, f"{name}.query"), n, 0, rng,
                          fixed=[(0, n - 1), (k - 1, n - k), (0, n - k)])
            manifest.append((name, n, len(edges), num_rules, "grid"))

    # Layered DAGs: frontier size controlled directly by depth.
    for num_rules in (3, 4):
        for layers, width in ((6, 6), (10, 6), (14, 6), (20, 6), (26, 6), (34, 6)):
            name = f"layer_L{layers:02d}W{width}_N{num_rules}"
            n, edges, src, dst = gen_layered_graph(layers, width, num_rules,
                                                   args.cmax, rng)
            write_graph(os.path.join(args.out, f"{name}.graph"), n, edges, num_rules)
            write_queries(os.path.join(args.out, f"{name}.query"), n, 0, rng,
                          fixed=[(src, dst)])
            manifest.append((name, n, len(edges), num_rules, "layered"))

    # Tiny graphs reserved for brute-force validation.
    for i in range(12):
        n = rng.randint(6, 11)
        num_rules = rng.choice([3, 4])
        name = f"tiny_{i}_N{num_rules}"
        edges = gen_random_graph(n, 3, num_rules, 9, rng)
        write_graph(os.path.join(args.out, f"{name}.graph"), n, edges, num_rules)
        write_queries(os.path.join(args.out, f"{name}.query"), n, 3, rng)
        manifest.append((name, n, len(edges), num_rules, "tiny"))

    with open(os.path.join(args.out, "manifest.tsv"), "w") as f:
        f.write("name\tn\tm\tnum_rules\tfamily\n")
        for row in manifest:
            f.write("\t".join(map(str, row)) + "\n")
    print(f"wrote {len(manifest)} graphs and {len(RULEBOOKS)} rulebooks to {args.out}/")


if __name__ == "__main__":
    main()
