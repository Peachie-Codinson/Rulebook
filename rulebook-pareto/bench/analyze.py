#!/usr/bin/env python3
"""Turn runs.csv into the comparison the advisor asked for.

Cases that hit the 10-minute cap are removed: a (graph, rulebook, query, eps)
case is dropped entirely -- for every solver -- as soon as any solver on it
timed out or was killed, so the surviving comparison is always like-for-like.

For each surviving case the exact frontier is the reference, and each solver is
scored on:
  runtime        seconds of search (heuristic preprocessing excluded)
  size           number of solutions returned
  sound          every returned cost is genuinely rulebook-optimal
  recall         fraction of the exact frontier actually returned
  covered        fraction of the exact frontier rule-dominated by the result
  eps_achieved   smallest uniform eps at which the result eps-rule-dominates
                 the whole exact frontier
"""
import argparse
import collections
import csv
import json
import os
import sys

# Frontiers can hold thousands of cost vectors, so the costs column blows past
# the default CSV field limit.
csv.field_size_limit(1 << 30)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate import Rulebook  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def min_eps_to_dominate(rb, u, v):
    """Smallest uniform eps >= 0 with u eps-rule-dominating v, or None."""
    cands = {0.0}
    for j in range(rb.n):
        if u[j] > v[j] and v[j] > 0:
            cands.add(u[j] / v[j] - 1.0)
    for c in sorted(cands):
        if rb.dominates(u, v, eps=[c + 1e-12] * rb.n):
            return c
    return None


def score(rb, got, truth, sample_cap=200000):
    """Score a solver's output against the exact frontier.

    Frontiers reach thousands of vectors, so the naive |got| x |truth| sweep is
    the bottleneck. Three shortcuts keep it near-linear in the common cases:
      - a returned vector that is itself in the frontier is trivially sound,
        covers itself, and needs eps = 0;
      - the frontier is an antichain, so anything drawn from it is sound;
      - only frontier points the solver did *not* return need the full sweep.
    """
    tset = {tuple(x) for x in truth}
    gset = {tuple(x) for x in got}

    recall = len(gset & tset) / len(tset) if tset else 1.0

    # Sound: every returned vector is rulebook-optimal. Vectors drawn from the
    # frontier are sound for free; only novel ones need checking.
    sound = True
    for gvec in (list(x) for x in gset - tset):
        if any(rb.strictly(list(w), gvec) for w in tset):
            sound = False
            break

    missing = [list(x) for x in tset - gset]
    got_l = [list(x) for x in gset]

    # Points the solver returned cover themselves at eps = 0, so only the
    # missing ones can lower `covered` or raise `eps*`.
    if not missing:
        return sound, recall, 1.0, 0.0
    if not got_l:
        return sound, recall, 0.0, float("inf")

    num_missing = len(missing)
    work = num_missing * len(got_l)
    if work > sample_cap and num_missing > 1:
        step = max(1, work // sample_cap)
        missing = missing[::step]

    covered_hits = 0
    worst = 0.0
    for v in missing:
        best = None
        dominated = False
        for u in got_l:
            if rb.dominates(u, v):
                dominated = True
            e = min_eps_to_dominate(rb, u, v)
            if e is not None and (best is None or e < best):
                best = e
                if best == 0.0:
                    break
        covered_hits += dominated
        if best is None:
            worst = float("inf")
            break
        worst = max(worst, best)

    # Returned points cover themselves; the rest is the (possibly sampled)
    # coverage rate over the points the solver missed.
    rate = covered_hits / len(missing)
    covered = (len(tset) - num_missing + rate * num_missing) / len(tset)
    return sound, recall, covered, worst


def fmt(x, nd=3):
    if x == float("inf"):
        return "inf"
    return f"{x:.{nd}f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=os.path.join(ROOT, "results", "runs.csv"))
    ap.add_argument("--instances", default=os.path.join(ROOT, "instances"))
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "report.md"))
    ap.add_argument("--csv-out", default=os.path.join(ROOT, "results", "scored.csv"))
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.runs)))
    for r in rows:
        r["timed_out"] = r["timed_out"].lower() == "true"
        r["killed"] = r["killed"].lower() == "true"
        r["runtime"] = float(r["runtime"])
        r["num_solutions"] = int(r["num_solutions"])

    # Case key. eps is part of it so approximate runs are compared against the
    # exact frontier of the same instance but not mixed across eps.
    def case(r):
        return (r["graph"], r["rulebook"], r["s"], r["t"])

    # --- Step 1: drop every case in which any solver exceeded the cap --------
    capped = set()
    all_cases = set()
    for r in rows:
        if r["s"] == "-1":
            # whole-file failure: drop every case of that graph+rulebook
            capped.add((r["graph"], r["rulebook"]))
            continue
        all_cases.add(case(r))
        if r["timed_out"] or r["killed"]:
            capped.add(case(r))

    def is_capped(r):
        return case(r) in capped or (r["graph"], r["rulebook"]) in capped

    kept = [r for r in rows if r["s"] != "-1" and not is_capped(r)]
    dropped_cases = {c for c in all_cases
                     if c in capped or (c[0], c[1]) in capped}

    print(f"cases total   : {len(all_cases)}")
    print(f"cases dropped : {len(dropped_cases)} (some solver exceeded the cap)")
    print(f"cases kept    : {len(all_cases) - len(dropped_cases)}")

    # Which solver caused each drop, for the report.
    blame = collections.Counter()
    for r in rows:
        if r["timed_out"] or r["killed"]:
            blame[r["alg"]] += 1

    # --- Step 2: score every kept run against the exact frontier -------------
    rbcache = {}

    def get_rb(name):
        if name not in rbcache:
            rbcache[name] = Rulebook(os.path.join(args.instances, f"rb_{name}.txt"))
        return rbcache[name]

    truth = {}
    for r in kept:
        if r["alg"] == "exact":
            truth[case(r)] = json.loads(r["costs"])

    scored = []
    for r in kept:
        c = case(r)
        if c not in truth:
            continue
        rb = get_rb(r["rulebook"])
        got = json.loads(r["costs"])
        sound, recall, covered, eps_ach = score(rb, got, truth[c])
        scored.append({
            "graph": r["graph"], "family": r["family"], "n": int(r["n"]),
            "m": int(r["m"]), "num_rules": int(r["num_rules"]),
            "rulebook": r["rulebook"], "alg": r["alg"], "eps": float(r["eps"]),
            "s": r["s"], "t": r["t"], "runtime": r["runtime"],
            "size": r["num_solutions"], "frontier": len(truth[c]),
            "expansions": int(r["expansions"]),
            "dijkstras": int(r["dijkstras"]),
            "residual_edges": int(r["residual_edges"]),
            "peeled_rules": int(r["peeled_rules"]),
            "sound": sound, "recall": recall, "covered": covered,
            "eps_achieved": eps_ach,
        })

    with open(args.csv_out, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(scored[0].keys()))
        wr.writeheader()
        wr.writerows(scored)

    # --- Step 3: report ------------------------------------------------------
    out = []
    out.append("# Rulebook Pareto-set: angle comparison\n")
    out.append(f"- cases total: **{len(all_cases)}**, dropped for exceeding the "
               f"10-minute cap: **{len(dropped_cases)}**, kept: "
               f"**{len(all_cases) - len(dropped_cases)}**\n")
    if blame:
        out.append("- runs that hit the cap, by solver: " +
                   ", ".join(f"`{a}` {n}" for a, n in blame.most_common()) + "\n")
    out.append("\n> `peel-only` performs the Dijkstra reduction and then stops "
               "without solving the residual, so it returns a solution set only "
               "when the rulebook is fully peelable. Its row is a measurement of "
               "how much work peeling removes, not of a solver; read it together "
               "with the peeling table below.\n")

    def table(rows_, keys, header):
        out.append("| " + " | ".join(header) + " |")
        out.append("|" + "|".join("---" for _ in header) + "|")
        for row in rows_:
            out.append("| " + " | ".join(str(row[k]) for k in keys) + " |")
        out.append("")

    # Per-rulebook summary, averaged over instances.
    out.append("\n## By rulebook shape\n")
    for rbname in sorted({s["rulebook"] for s in scored}):
        sub = [s for s in scored if s["rulebook"] == rbname]
        out.append(f"\n### `{rbname}`  ({len(sub)//max(1,len(set(s['alg']+str(s['eps']) for s in sub)))} cases)\n")
        agg = collections.defaultdict(list)
        for s in sub:
            agg[(s["alg"], s["eps"])].append(s)
        rows_ = []
        for (alg, eps), items in sorted(agg.items()):
            n = len(items)
            mean_rt = sum(i["runtime"] for i in items) / n
            max_rt = max(i["runtime"] for i in items)
            speedup = None
            ex = [i for i in sub if i["alg"] == "exact"]
            if ex:
                ex_rt = sum(i["runtime"] for i in ex) / len(ex)
                speedup = ex_rt / mean_rt if mean_rt > 0 else float("inf")
            rows_.append({
                "alg": f"`{alg}`", "eps": eps, "mean s": fmt(mean_rt, 4),
                "max s": fmt(max_rt, 4),
                "speedup": fmt(speedup, 1) if speedup else "-",
                "size": fmt(sum(i["size"] for i in items) / n, 1),
                "|P_R|": fmt(sum(i["frontier"] for i in items) / n, 1),
                "sound": f"{sum(i['sound'] for i in items)}/{n}",
                "recall": fmt(sum(i["recall"] for i in items) / n),
                "covered": fmt(sum(i["covered"] for i in items) / n),
                "eps*": fmt(max(i["eps_achieved"] for i in items)),
            })
        table(rows_, ["alg", "eps", "mean s", "max s", "speedup", "size",
                      "|P_R|", "sound", "recall", "covered", "eps*"],
              ["alg", "eps", "mean s", "max s", "speedup", "size",
               "\\|P_R\\|", "sound", "recall", "covered", "eps* (worst)"])

    # How much of the problem peeling removes, per rulebook shape.
    out.append("\n## Peeling effectiveness\n")
    out.append("`peeled` is how many rules were discharged by two Dijkstras each; "
               "`edges left` is the fraction of the graph that survived into the "
               "residual search.\n")
    prows = []
    for rbname in sorted({s_["rulebook"] for s_ in scored}):
        sub = [s_ for s_ in scored if s_["rulebook"] == rbname
               and s_["alg"] == "peel-only"]
        if not sub:
            continue
        n = len(sub)
        nrules = sub[0]["num_rules"]
        prows.append({
            "rulebook": f"`{rbname}`",
            "rules": nrules,
            "peeled": fmt(sum(i["peeled_rules"] for i in sub) / n, 2),
            "edges left": fmt(sum(i["residual_edges"] / i["m"] for i in sub) / n, 4),
            "fully solved": f"{sum(1 for i in sub if i['recall'] == 1.0)}/{n}",
        })
    table(prows, ["rulebook", "rules", "peeled", "edges left", "fully solved"],
          ["rulebook", "rules", "rules peeled", "edges left (frac)",
           "fully solved by peeling alone"])

    # Scaling table: runtime vs graph size for the hardest rulebook.
    out.append("\n## Scaling (runtime in seconds, mean over queries)\n")
    for rbname in sorted({s["rulebook"] for s in scored}):
        sub = [s for s in scored if s["rulebook"] == rbname]
        graphs = sorted({(s["graph"], s["n"], s["m"]) for s in sub},
                        key=lambda x: (x[1], x[2]))
        algs = sorted({(s["alg"], s["eps"]) for s in sub})
        out.append(f"\n### `{rbname}`\n")
        out.append("| graph | n | m | " +
                   " | ".join(f"`{a}`" + (f" e={e}" if e else "") for a, e in algs) + " |")
        out.append("|" + "|".join("---" for _ in range(3 + len(algs))) + "|")
        for gname, n, m in graphs:
            cells = []
            for a, e in algs:
                items = [s for s in sub if s["graph"] == gname and s["alg"] == a
                         and s["eps"] == e]
                cells.append(fmt(sum(i["runtime"] for i in items) / len(items), 4)
                             if items else "-")
            out.append(f"| {gname} | {n} | {m} | " + " | ".join(cells) + " |")
        out.append("")

    open(args.out, "w").write("\n".join(out))
    print(f"\nwrote {args.out} and {args.csv_out}")


if __name__ == "__main__":
    main()
