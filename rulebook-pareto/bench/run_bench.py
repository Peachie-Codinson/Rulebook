#!/usr/bin/env python3
"""Run every solver on every (graph, rulebook, query) under a hard wall-clock
cap and record one row per run.

The cap is enforced twice: the binary stops its own search at --timeout and
reports timed_out=true, and the subprocess is killed at a slightly longer wall
limit in case a solver blows up outside the checked loop. Either way the run is
recorded as capped, and analyze.py drops capped cases from the comparison.
"""
import argparse
import csv
import glob
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin", "rbsearch")

# alg name -> extra CLI args. eps is filled in per sweep value.
ALGS = [
    ("exact", []),
    ("rapex-nodr", []),
    ("rapex", []),
    ("topolex", []),
    ("seed-exact", []),
    ("peel-only", []),
    ("peel-exact", []),
    ("peel-seed-exact", []),
    ("peel-rapex", []),
]


def rulebook_size(path):
    return int(open(path).read().split()[0])


def run_one(graph, rules, query, alg, timeout, eps, extra, num_queries=1):
    """The binary applies `timeout` to each query in the file and keeps going,
    so the subprocess budget has to cover every query, not just one. A tighter
    kill would throw away the fast queries that share the file with a slow one."""
    cmd = [BIN, "--graph", graph, "--rules", rules, "--query", query,
           "--alg", alg, "--timeout", str(timeout), "--eps", str(eps)] + extra
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=timeout * num_queries * 1.5 + 60)
    except subprocess.TimeoutExpired:
        return None, time.time() - t0, "killed"
    if p.returncode != 0:
        return None, time.time() - t0, f"error: {p.stderr.strip()[:200]}"
    try:
        return json.loads(p.stdout), time.time() - t0, ""
    except json.JSONDecodeError:
        return None, time.time() - t0, "unparseable output"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--instances", default=os.path.join(ROOT, "instances"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "runs.csv"))
    ap.add_argument("--timeout", type=float, default=600.0,
                    help="per-query wall-clock cap in seconds (default 600 = 10 min)")
    ap.add_argument("--eps", default="0,0.01",
                    help="comma-separated approximation factors to sweep")
    ap.add_argument("--families", default="random,grid")
    ap.add_argument("--skip-tiny", action="store_true", default=True)
    ap.add_argument("--only-alg", default="")
    ap.add_argument("--max-ext", type=int, default=0)
    ap.add_argument("--num-rules", default="",
                    help="restrict to graphs with this many objectives")
    args = ap.parse_args()

    eps_values = [float(x) for x in args.eps.split(",")]
    families = set(args.families.split(","))

    graphs = []
    for line in open(os.path.join(args.instances, "manifest.tsv")).read().splitlines()[1:]:
        name, n, m, num_rules, fam = line.split("\t")
        if fam not in families:
            continue
        if args.num_rules and num_rules != args.num_rules:
            continue
        graphs.append((name, int(n), int(m), int(num_rules), fam))

    rulebooks = sorted(glob.glob(os.path.join(args.instances, "rb_*.txt")))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    fields = ["graph", "family", "n", "m", "num_rules", "rulebook", "alg", "eps",
              "s", "t", "runtime", "timed_out", "killed", "num_solutions",
              "expansions", "generations", "dijkstras", "extensions",
              "residual_edges", "peeled_rules", "costs", "note"]
    f = open(args.out, "w", newline="")
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()

    # Difficulty inside a family is monotone in graph size, so once a solver
    # exceeds the cap at one size it will exceed it at every larger size. We
    # record those as capped without running them -- analyze.py drops capped
    # cases anyway, so this only saves wall time, it does not change the
    # surviving comparison.
    abandoned = set()   # (family, rulebook, alg, eps, num_rules)

    def capped_row(name, fam, n, m, num_rules, rname, alg, eps, wall, note):
        return {"graph": name, "family": fam, "n": n, "m": m,
                "num_rules": num_rules, "rulebook": rname, "alg": alg,
                "eps": eps, "s": -1, "t": -1, "runtime": wall,
                "timed_out": True, "killed": True, "num_solutions": 0,
                "expansions": 0, "generations": 0, "dijkstras": 0,
                "extensions": 0, "residual_edges": 0, "peeled_rules": 0,
                "costs": "[]", "note": note}

    graphs.sort(key=lambda g: (g[4], g[3], g[1], g[2]))

    total = 0
    for name, n, m, num_rules, fam in graphs:
        gpath = os.path.join(args.instances, f"{name}.graph")
        qpath = os.path.join(args.instances, f"{name}.query")
        for rpath in rulebooks:
            if rulebook_size(rpath) != num_rules:
                continue
            rname = os.path.basename(rpath)[3:-4]
            for eps in eps_values:
                for alg, extra in ALGS:
                    if args.only_alg and alg != args.only_alg:
                        continue
                    # eps only changes the approximate solvers; skip duplicates.
                    if eps != eps_values[0] and alg in (
                            "exact", "topolex", "seed-exact", "peel-only",
                            "peel-exact", "peel-seed-exact"):
                        continue
                    key = (fam, rname, alg, eps, num_rules)
                    if key in abandoned:
                        w.writerow(capped_row(name, fam, n, m, num_rules, rname,
                                              alg, eps, args.timeout,
                                              "skipped: capped at a smaller size"))
                        total += 1
                        print(f"  {name:18s} {rname:10s} {alg:12s} eps={eps} "
                              f"-> skipped (capped earlier)", flush=True)
                        continue

                    ex = list(extra)
                    if alg == "topolex" and args.max_ext:
                        ex += ["--max-ext", str(args.max_ext)]

                    nq = sum(1 for line in open(qpath) if line.strip())
                    data, wall, err = run_one(gpath, rpath, qpath, alg,
                                              args.timeout, eps, ex, nq)
                    if data is None:
                        w.writerow(capped_row(name, fam, n, m, num_rules, rname,
                                              alg, eps, wall, err))
                        abandoned.add(key)
                        total += 1
                        print(f"  {name:18s} {rname:10s} {alg:12s} eps={eps} "
                              f"-> {err}", flush=True)
                        continue
                    for r in data["results"]:
                        w.writerow({
                            "graph": name, "family": fam, "n": n, "m": m,
                            "num_rules": num_rules, "rulebook": rname,
                            "alg": alg, "eps": eps, "s": r["s"], "t": r["t"],
                            "runtime": r["runtime"],
                            "timed_out": r["timed_out"], "killed": False,
                            "num_solutions": r["num_solutions"],
                            "expansions": r["expansions"],
                            "generations": r["generations"],
                            "dijkstras": r["dijkstras"],
                            "extensions": r["extensions"],
                            "residual_edges": r["residual_edges"],
                            "peeled_rules": r["peeled_rules"],
                            "costs": json.dumps(r["costs"]), "note": ""})
                        total += 1
                    f.flush()
                    slowest = max((r["runtime"] for r in data["results"]), default=0)
                    to = sum(1 for r in data["results"] if r["timed_out"])
                    if to:
                        abandoned.add(key)
                    print(f"  {name:18s} {rname:10s} {alg:12s} eps={eps} "
                          f"max={slowest:8.3f}s timeouts={to}", flush=True)

    f.close()
    print(f"\nwrote {total} rows to {args.out}")


if __name__ == "__main__":
    main()
