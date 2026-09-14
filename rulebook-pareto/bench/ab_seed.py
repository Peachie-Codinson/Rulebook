#!/usr/bin/env python3
"""Contention-free A/B: does warm-starting the exact search with TopoLex's
certified-optimal solutions actually make it faster?

The main sweep runs several solvers in parallel, which is fine for reach and
correctness but not for a head-to-head on runtime: `exact` and `seed-exact` were
measured under different amounts of CPU contention. This script runs the two
back to back, single-threaded, alternating order within each case so that any
drift (thermal, cache, noisy neighbour) hits both arms equally, and repeats each
measurement so the spread is visible rather than assumed.

Only cases the exact search already finishes inside the cap are eligible, since
a ratio needs both arms to produce a number.
"""
import argparse
import csv
import json
import os
import statistics
import subprocess
import sys

csv.field_size_limit(1 << 30)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin", "rbsearch")


def run(graph, rules, s, t, alg, timeout, tmpdir):
    qf = os.path.join(tmpdir, "ab_query.txt")
    with open(qf, "w") as f:
        f.write(f"{s} {t}\n")
    cmd = [BIN, "--graph", graph, "--rules", rules, "--query", qf,
           "--alg", alg, "--timeout", str(timeout), "--eps", "0"]
    p = subprocess.run(cmd, capture_output=True, text=True,
                       timeout=timeout * 1.5 + 60)
    if p.returncode != 0:
        return None
    r = json.loads(p.stdout)["results"][0]
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True,
                    help="completed sweep CSVs, used to pick eligible cases")
    ap.add_argument("--instances", default=os.path.join(ROOT, "instances"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "ab_seed.csv"))
    ap.add_argument("--timeout", type=float, default=600.0)
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--min-runtime", type=float, default=0.05,
                    help="skip cases the exact search already solves instantly")
    ap.add_argument("--max-cases", type=int, default=40)
    ap.add_argument("--tmpdir", default="/tmp")
    args = ap.parse_args()

    # Eligible cases: exact finished inside the cap and took long enough that a
    # ratio means something.
    cases = []
    for path in args.runs:
        for r in csv.DictReader(open(path)):
            if r["alg"] != "exact" or r["s"] == "-1":
                continue
            if r["timed_out"].lower() == "true" or r["killed"].lower() == "true":
                continue
            if float(r["runtime"]) < args.min_runtime:
                continue
            cases.append((float(r["runtime"]), r["graph"], r["rulebook"],
                          r["s"], r["t"]))
    cases.sort(reverse=True)
    cases = cases[:args.max_cases]
    print(f"{len(cases)} eligible cases "
          f"(exact runtime {cases[-1][0]:.3f}s .. {cases[0][0]:.3f}s)"
          if cases else "no eligible cases")
    if not cases:
        return

    rows = []
    for i, (prev, graph, rbname, s, t) in enumerate(cases):
        gpath = os.path.join(args.instances, f"{graph}.graph")
        rpath = os.path.join(args.instances, f"rb_{rbname}.txt")
        ex, sd = [], []
        ok = True
        for rep in range(args.repeats):
            # Alternate which arm goes first so order effects cancel.
            order = ["exact", "seed-exact"] if rep % 2 == 0 else \
                    ["seed-exact", "exact"]
            got = {}
            for alg in order:
                r = run(gpath, rpath, s, t, alg, args.timeout, args.tmpdir)
                if r is None or r["timed_out"]:
                    ok = False
                    break
                got[alg] = r
            if not ok:
                break
            # Both arms are exact, so they must agree on the frontier. If they
            # ever disagree the speed comparison is meaningless.
            a = sorted(map(tuple, got["exact"]["costs"]))
            b = sorted(map(tuple, got["seed-exact"]["costs"]))
            if a != b:
                print(f"  MISMATCH {graph} {rbname} {s}->{t}: "
                      f"exact={len(a)} seed-exact={len(b)}")
                ok = False
                break
            ex.append(got["exact"]["runtime"])
            sd.append(got["seed-exact"]["runtime"])
        if not ok or not ex:
            continue
        me, ms = statistics.median(ex), statistics.median(sd)
        rows.append({
            "graph": graph, "rulebook": rbname, "s": s, "t": t,
            "frontier": len(got["exact"]["costs"]),
            "exact_median": me, "seed_median": ms,
            "speedup": me / ms if ms > 0 else float("inf"),
            "exact_spread": max(ex) - min(ex),
            "seed_spread": max(sd) - min(sd),
            "exact_expansions": got["exact"]["expansions"],
            "seed_expansions": got["seed-exact"]["expansions"],
        })
        print(f"  [{i+1}/{len(cases)}] {graph:18s} {rbname:10s} "
              f"{s}->{t}  exact={me:8.3f}s  seed={ms:8.3f}s  "
              f"x{me/ms if ms else 0:5.2f}", flush=True)

    if not rows:
        print("no usable pairs")
        return
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    sp = [r["speedup"] for r in rows]
    print(f"\n{len(rows)} paired cases")
    print(f"speedup  median {statistics.median(sp):.2f}x   "
          f"min {min(sp):.2f}x   max {max(sp):.2f}x")
    faster = sum(1 for x in sp if x > 1.05)
    slower = sum(1 for x in sp if x < 0.95)
    print(f"seed-exact faster in {faster}/{len(rows)}, "
          f"slower in {slower}/{len(rows)}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
