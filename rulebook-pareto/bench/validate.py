#!/usr/bin/env python3
"""Cross-validate every solver against brute-force path enumeration on the tiny
graphs, and check the two properties the new solvers are supposed to have:

  exact      == brute force, exactly
  topolex    subset of the exact frontier (soundness), possibly smaller
  peel-exact == exact frontier (peeling is claimed to be lossless)
  rapex(0)   covers the exact frontier under rule-dominance
"""
import glob
import itertools
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.environ.get("RBSEARCH_BIN", os.path.join(ROOT, "bin", "rbsearch"))


def run(graph, rules, query, alg, timeout=60, eps=None, extra=()):
    cmd = [BIN, "--graph", graph, "--rules", rules, "--query", query,
           "--alg", alg, "--timeout", str(timeout)]
    if eps is not None:
        cmd += ["--eps", str(eps)]
    cmd += list(extra)
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout + 60)
    if out.returncode != 0:
        raise RuntimeError(f"{alg} failed: {out.stderr[:400]}")
    return json.loads(out.stdout)


def parse_rulebook(path):
    """Returns (num_rules, list of equivalence classes, strict relations)."""
    toks = open(path).read().split()
    i = 0
    n = int(toks[i]); i += 1
    i += n  # eps
    k = int(toks[i]); i += 1
    classes = []
    for _ in range(k):
        c = int(toks[i]); i += 1
        classes.append([int(toks[i + j]) for j in range(c)])
        i += c
    p = int(toks[i]); i += 1
    rels = []
    for _ in range(p):
        rels.append((int(toks[i]), int(toks[i + 1])))
        i += 2
    return n, classes, rels


class Rulebook:
    """Python mirror of the C++ dominance test, so validation is independent of
    the implementation under test."""

    def __init__(self, path):
        n, classes, rels = parse_rulebook(path)
        self.n = n
        self.cls_of = {}
        for ci, c in enumerate(classes):
            for r in c:
                self.cls_of[r] = ci
        self.classes = classes
        # strictly-above relation on rules, transitively closed
        above = {r: set() for r in range(n)}
        edges = set()
        for a, b in rels:
            for x in classes[self.cls_of[a]]:
                for y in classes[self.cls_of[b]]:
                    edges.add((x, y))
        changed = True
        while changed:
            changed = False
            for a, b in list(edges):
                for c in range(n):
                    if (b, c) in edges and (a, c) not in edges:
                        edges.add((a, c))
                        changed = True
        for a, b in edges:
            if self.cls_of[a] != self.cls_of[b]:
                above[b].add(a)
        self.above = above

    def dominates(self, u, v, eps=None):
        """u weakly rule-dominates v (Def. 8 / Def. 5)."""
        def worse(j):
            lim = v[j] * (1 + eps[j]) if eps else v[j]
            return u[j] > lim + 1e-9

        def better(j):
            lim = v[j] * (1 + eps[j]) if eps else v[j]
            return u[j] < lim - 1e-9

        for j in range(self.n):
            if worse(j) and not any(better(i) for i in self.above[j]):
                return False
        return True

    def strictly(self, u, v):
        return self.dominates(u, v) and not self.dominates(v, u)

    def filter_optimal(self, vs):
        vs = sorted({tuple(v) for v in vs})
        return [v for v in vs if not any(self.strictly(w, v) for w in vs if w != v)]


def main():
    inst = os.path.join(ROOT, "instances")
    tinies = sorted(glob.glob(os.path.join(inst, "tiny_*.graph")))
    rulebooks = sorted(glob.glob(os.path.join(inst, "rb_*.txt")))
    failures = []
    checks = 0
    stats = {"topolex_subset": 0, "topolex_complete": 0, "topolex_total": 0}

    for gpath in tinies:
        base = gpath[:-6]
        qpath = base + ".query"
        num_rules = int(base.split("_N")[-1])
        for rpath in rulebooks:
            rb = Rulebook(rpath)
            if rb.n != num_rules:
                continue
            runs = {a: run(gpath, rpath, qpath, a)
                    for a in ("brute", "exact", "topolex", "seed-exact",
                              "peel-exact", "peel-seed-exact", "peel-rapex",
                              "rapex", "rapex-nodr")}
            runs["rapex0"] = run(gpath, rpath, qpath, "rapex", eps=0.0)

            for qi in range(len(runs["brute"]["results"])):
                checks += 1
                tag = f"{os.path.basename(base)} {os.path.basename(rpath)} q{qi}"
                truth = {tuple(c) for c in runs["brute"]["results"][qi]["costs"]}
                got = {a: {tuple(c) for c in runs[a]["results"][qi]["costs"]}
                       for a in runs}

                # Sanity: brute force should already be a rulebook-optimal set.
                if set(map(tuple, rb.filter_optimal(truth))) != truth:
                    failures.append(f"{tag}: brute output not self-consistent")

                if got["exact"] != truth:
                    failures.append(
                        f"{tag}: exact != brute  exact={sorted(got['exact'])} "
                        f"brute={sorted(truth)}")

                for a in ("peel-exact", "seed-exact", "peel-seed-exact"):
                    if got[a] != truth:
                        failures.append(
                            f"{tag}: {a} != brute  got={sorted(got[a])} "
                            f"brute={sorted(truth)}")

                # topolex must be sound (subset), completeness is not claimed.
                stats["topolex_total"] += 1
                if got["topolex"] <= truth:
                    stats["topolex_subset"] += 1
                else:
                    failures.append(
                        f"{tag}: topolex returned a non-optimal solution "
                        f"{sorted(got['topolex'] - truth)}")
                if got["topolex"] == truth:
                    stats["topolex_complete"] += 1

                # RA*pex at eps=0 must cover the frontier under rule-dominance.
                for a in ("rapex0", "rapex", "rapex-nodr", "peel-rapex"):
                    for v in truth:
                        if not any(rb.dominates(list(u), list(v)) for u in got[a]):
                            failures.append(f"{tag}: {a} fails to cover {v}")
                            break

    print(f"checked {checks} (instance, rulebook, query) triples")
    print(f"topolex sound (subset of frontier): "
          f"{stats['topolex_subset']}/{stats['topolex_total']}")
    print(f"topolex complete (== frontier):     "
          f"{stats['topolex_complete']}/{stats['topolex_total']}")
    if failures:
        print(f"\n{len(failures)} FAILURES:")
        for f in failures[:40]:
            print("  " + f)
        sys.exit(1)
    print("\nall invariants hold")


if __name__ == "__main__":
    main()
