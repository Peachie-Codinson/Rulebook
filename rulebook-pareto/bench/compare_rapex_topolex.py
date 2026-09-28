#!/usr/bin/env python3
"""Head-to-head: RA*pex vs multi-level Dijkstra (TopoLex).

Just the two algorithms, on one dataset, with no exact ground truth required --
which matters, because on a road network the size of BAY the exact
rulebook-optimal frontier is not computable inside any reasonable budget.

The trick that makes a ground-truth-free comparison meaningful: every solution
TopoLex returns is *provably* rulebook-optimal (lex-minimal under some linear
extension of the rulebook cannot be strictly rule-dominated -- see the soundness
lemma in the README). So TopoLex's output is itself a certificate:

  * an RA*pex solution strictly rule-dominated by a TopoLex solution is
    PROVABLY non-optimal, no ground truth needed;
  * a TopoLex solution strictly rule-dominated by an RA*pex solution would
    contradict the lemma, so it is reported as a soundness violation and should
    never occur.

Because TopoLex's points are certified optimal, the reverse direction checks
RA*pex's own guarantee: Theorem 1 of the RA*pex paper says every rulebook-optimal
solution is eps-rule-dominated by something in the returned set, so measuring the
eps RA*pex actually needs to cover TopoLex's points is a live test of that
theorem on data far too large to compute a frontier for. It should come out at
or below the eps RA*pex was run with.

Note on what is NOT reported: the plain (eps = 0) coverage fractions between the
two sets are structurally zero and carry no information. Weak rule-dominance is
antisymmetric on distinct vectors -- if u <~ v and v <~ u then u == v -- so two
distinct *optimal* points can never weakly dominate one another, and exact
coverage can only ever come from solutions both sets happen to share. That is
already reported as `shared`. The eps-coverage columns are the informative ones.

Reported per query:
  runtime, |solutions|, and the Dijkstra/expansion counts
  shared              solutions both returned
  rapex_provably_bad  RA*pex solutions strictly dominated by a TopoLex solution,
                      i.e. proven non-optimal with no ground truth
  topolex_violations  TopoLex solutions strictly dominated by RA*pex (expect 0;
                      anything else contradicts the soundness lemma)
  eps_rapex_covers_topolex   eps RA*pex needs to cover TopoLex's certified
                             points -- a check on RA*pex's Theorem 1, should be
                             <= the eps it was run with
  eps_topolex_covers_rapex   eps TopoLex needs to cover what RA*pex found --
                             the price of TopoLex's incompleteness
"""
import argparse
import csv
import json
import os
import statistics
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate import Rulebook  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin", "rbsearch")


def run(alg, graph_args, rules, query, timeout, eps, extra=()):
    cmd = [BIN] + graph_args + ["--rules", rules, "--query", query,
                                "--alg", alg, "--timeout", str(timeout),
                                "--eps", str(eps)] + list(extra)
    p = subprocess.run(cmd, capture_output=True, text=True,
                       timeout=timeout * 400 + 120)
    if p.returncode != 0:
        raise SystemExit(f"{alg} failed:\n{p.stderr[:2000]}")
    return json.loads(p.stdout)


def min_eps_to_dominate(rb, u, v):
    cands = {0.0}
    for j in range(rb.n):
        if u[j] > v[j] and v[j] > 0:
            cands.add(u[j] / v[j] - 1.0)
    for c in sorted(cands):
        if rb.dominates(u, v, eps=[c + 1e-12] * rb.n):
            return c
    return None


def eps_to_cover(rb, cover_set, target_set):
    """Smallest uniform eps at which cover_set eps-rule-dominates all of
    target_set. None if no finite eps suffices."""
    worst = 0.0
    for v in target_set:
        best = None
        for u in cover_set:
            e = min_eps_to_dominate(rb, list(u), list(v))
            if e is not None and (best is None or e < best):
                best = e
                if best == 0.0:
                    break
        if best is None:
            return None
        worst = max(worst, best)
    return worst


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--graph", help="single-file instance")
    g.add_argument("--gr", help="comma-separated DIMACS .gr files, one per objective")
    ap.add_argument("--rules", required=True)
    ap.add_argument("--query", required=True)
    ap.add_argument("--eps", type=float, default=0.01,
                    help="approximation factor for RA*pex (TopoLex ignores it)")
    ap.add_argument("--timeout", type=float, default=600.0,
                    help="per-query cap in seconds")
    ap.add_argument("--max-ext", type=int, default=0,
                    help="cap TopoLex's linear extensions (0 = all)")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "rapex_vs_topolex.csv"))
    args = ap.parse_args()

    graph_args = ["--gr", args.gr] if args.gr else ["--graph", args.graph]
    rb = Rulebook(args.rules)

    print(f"rulebook : {args.rules} ({rb.n} rules)")
    print(f"queries  : {args.query}")
    print(f"eps      : {args.eps} (RA*pex only)\n")

    extra = ["--max-ext", str(args.max_ext)] if args.max_ext else []
    print("running rapex ...", flush=True)
    R = run("rapex", graph_args, args.rules, args.query, args.timeout, args.eps)
    print("running topolex ...", flush=True)
    T = run("topolex", graph_args, args.rules, args.query, args.timeout, 0.0, extra)

    print(f"\ngraph: n={R['n']} m={R['m']} objectives={R['N']} "
          f"quotient classes={R['classes']}\n")

    rows = []
    for rr, tr in zip(R["results"], T["results"]):
        assert (rr["s"], rr["t"]) == (tr["s"], tr["t"])
        rs = {tuple(c) for c in rr["costs"]}
        ts = {tuple(c) for c in tr["costs"]}

        # TopoLex output is certified optimal, so anything it strictly dominates
        # is certified NOT optimal -- a ground-truth-free unsoundness witness.
        bad = [u for u in rs if any(rb.strictly(list(w), list(u)) for w in ts)]
        viol = [u for u in ts if any(rb.strictly(list(w), list(u)) for w in rs)]

        e_t_covers_r = eps_to_cover(rb, ts, rs) if ts and rs else None
        e_r_covers_t = eps_to_cover(rb, rs, ts) if ts and rs else None

        row = {
            "s": rr["s"], "t": rr["t"],
            "rapex_s": rr["runtime"], "topolex_s": tr["runtime"],
            "speedup": (rr["runtime"] / tr["runtime"]) if tr["runtime"] > 0 else "",
            "rapex_n": len(rs), "topolex_n": len(ts),
            "shared": len(rs & ts),
            "rapex_provably_bad": len(bad),
            "topolex_violations": len(viol),
            "eps_rapex_covers_topolex": (
                "inf" if e_r_covers_t is None else round(e_r_covers_t, 6)),
            "eps_topolex_covers_rapex": (
                "inf" if e_t_covers_r is None else round(e_t_covers_r, 4)),
            "rapex_expansions": rr["expansions"],
            "topolex_dijkstras": tr["dijkstras"],
            "topolex_extensions": tr["extensions"],
            "rapex_timeout": rr["timed_out"], "topolex_timeout": tr["timed_out"],
        }
        rows.append(row)

        flag = ""
        if rr["timed_out"] or tr["timed_out"]:
            flag += "  [TIMEOUT]"
        if viol:
            flag += f"  [!! {len(viol)} TOPOLEX SOUNDNESS VIOLATIONS !!]"
        if e_r_covers_t is not None and e_r_covers_t > args.eps + 1e-9:
            flag += (f"  [!! rapex needed eps={e_r_covers_t:.4f} > {args.eps} "
                     f"to cover certified points !!]")
        print(f"{rr['s']:>8} -> {rr['t']:<8}  "
              f"rapex {rr['runtime']:9.3f}s n={len(rs):<6d}  "
              f"topolex {tr['runtime']:8.3f}s n={len(ts):<4d}  "
              f"shared={len(rs & ts):<4d} rapex_bad={len(bad):<4d}{flag}",
              flush=True)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    ok = [r for r in rows if not r["rapex_timeout"] and not r["topolex_timeout"]]
    print(f"\n{'='*70}")
    print(f"{len(rows)} queries, {len(ok)} finished inside the {args.timeout:.0f}s cap")
    if ok:
        sp = [r["speedup"] for r in ok if r["speedup"] != ""]
        print(f"topolex speedup over rapex : median {statistics.median(sp):.1f}x  "
              f"(min {min(sp):.1f}x, max {max(sp):.1f}x)")
        print(f"solution set size          : rapex median "
              f"{statistics.median(r['rapex_n'] for r in ok):.0f}, "
              f"topolex median {statistics.median(r['topolex_n'] for r in ok):.0f}")
        et = [r["eps_topolex_covers_rapex"] for r in ok
              if r["eps_topolex_covers_rapex"] != "inf"]
        if et:
            print(f"eps TopoLex needs to cover RA*pex's set : median "
                  f"{statistics.median(et):.3f}  (max {max(et):.3f})")
        er = [r["eps_rapex_covers_topolex"] for r in ok
              if r["eps_rapex_covers_topolex"] != "inf"]
        if er:
            worst = max(er)
            verdict = "within guarantee" if worst <= args.eps + 1e-9 else \
                      "EXCEEDS the eps it was run with <-- INVESTIGATE"
            print(f"eps RA*pex needs to cover certified points: max "
                  f"{worst:.6f} vs eps={args.eps}  ({verdict})")
        nbad = sum(r["rapex_provably_bad"] for r in ok)
        qbad = sum(1 for r in ok if r["rapex_provably_bad"])
        print(f"RA*pex solutions PROVEN non-optimal by TopoLex: "
              f"{nbad} across {qbad}/{len(ok)} queries")
        nviol = sum(r["topolex_violations"] for r in ok)
        print(f"TopoLex soundness violations: {nviol} "
              f"{'(as expected)' if nviol == 0 else '<-- INVESTIGATE'}")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
