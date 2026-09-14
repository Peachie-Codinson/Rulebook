# Rulebook

Work on multi-objective search under **rulebooks** — preorders over objectives
that generalise both Pareto dominance and lexicographic priority.

## Projects

### [`rulebook-pareto/`](rulebook-pareto/)

Exploratory testbed for the question *"what are the fastest ways to get the
**entire** rulebook-optimal solution set, and how do they trade off against each
other?"*

Compares five angles against a brute-force-validated ground truth, under a
10-minute per-case wall-clock cap:

| angle | guarantee |
|---|---|
| `exact` | ground truth — label setting with rule-dominance pruning |
| `rapex` | RA\*pex — complete ε-coverage, members need not be optimal |
| `topolex` | multi-phase Dijkstra over linear extensions — every member optimal, coverage not guaranteed |
| `peel` | Dijkstra reduction of the globally-dominant rule prefix — exact |
| `seed-exact` | `topolex` as a warm start for the exact search |

See [`rulebook-pareto/README.md`](rulebook-pareto/README.md) for the soundness
proofs behind `topolex` and `peel`, the metrics, and how to reproduce the runs.

## Background

- Slutsky, Yershov, Wongpiromsarn, Frazzoli — *Hierarchical Multiobjective
  Shortest Path Problems* (regular cost monoids; the iterated Dijkstra
  propagation that `topolex` and `peel` build on).
- Muhammetkulyyev, Salzman, Wongpiromsarn — *Approximate Multi-Objective Search
  Under Rulebooks* (ε-rule-dominance and RA\*pex). Reference implementation:
  <https://github.com/Infus3d/Rulebook_approximation>.
