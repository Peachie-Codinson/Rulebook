# Fast routes to the full rulebook-optimal set — exploratory testbed

Exploratory work on the question: *given a rulebook (a preorder over objectives),
what are the fastest ways to get the **entire** rulebook-optimal solution set,
and how do they trade off against each other?*

Prior work this builds on:

- **Slutsky, Yershov, Wongpiromsarn, Frazzoli — "Hierarchical Multiobjective
  Shortest Path Problems."** Introduces *regular* cost monoids and **Algorithm 2**,
  an iterated Dijkstra propagation that solves the *lexicographic* (totally
  ordered) multicost problem in polynomial time by repeatedly shrinking the graph
  to its optimal subgraph, one cost coordinate at a time.
- **Muhammetkulyyev, Salzman, Wongpiromsarn — "Approximate Multi-Objective Search
  Under Rulebooks."** Introduces ε-rule-dominance and **RA\*pex**, which returns a
  compact ε-approximate rulebook-optimal set, up to two orders of magnitude faster
  than exact rulebook synthesis. Reference implementation:
  <https://github.com/Infus3d/Rulebook_approximation>.

The angle the advisor asked about — *"topological ordering, multi-phase Dijkstra
top to bottom, compared against the approximate algorithm"* — is implemented here
as `topolex`, and it turns out to have a clean soundness guarantee (§2). A second
new angle, `peel`, fell out of the same analysis and is **exact**, not approximate
(§3).

---

## 0. Problem recap

A rulebook is `R = <R, ≲>`: rules `r_1..r_N`, each a nonnegative cost on paths,
plus a preorder saying which rules matter more. Equivalent rules collapse into
equivalence classes; the classes form a DAG.

`u ≲_R v` ("u weakly rule-dominates v") iff for every rule `r_j` where `u_j > v_j`
there is a strictly more important rule `r_i > r_j` with `u_i < v_i`. The goal is
`P_R`, the set of solutions not strictly rule-dominated by any other. `P_R`
generalises both the Pareto front (no relations at all) and the lexicographic
optimum (a total order), and can be exponentially large.

---

## 1. The angles

| solver | what it is | exact? | cost |
|---|---|---|---|
| `exact` | label-setting search, rule-dominance pruning, lexicographic queue | yes — ground truth | exponential |
| `rapex` / `rapex-nodr` | RA\*pex, with / without dimensionality reduction | ε-covering | exponential, much smaller constant |
| `topolex` | **new** — multi-phase Dijkstra over linear extensions | sound subset of `P_R` | polynomial per extension |
| `peel-only` | **new** — Dijkstra optimal-subgraph peeling, nothing else | exact when fully peelable | 2 Dijkstras per peeled rule |
| `peel-exact` | peel, then `exact` on the residual | yes | exponential in the residual only |
| `peel-rapex` | peel, then RA\*pex on the residual | ε-covering | — |
| `brute` | path enumeration, validation only | yes | factorial |

---

## 2. `topolex` — multi-phase Dijkstra over linear extensions

**Construction.** Take any linear extension σ of the rulebook's quotient DAG (any
ordering of rules in which nothing appears before something strictly more
important). Run Slutsky et al.'s Algorithm 2 along σ: for rule `σ(1)`, a forward
and a backward Dijkstra identify the edges lying on some `σ(1)`-optimal path and
the graph is shrunk to that optimal subgraph; repeat for `σ(2)` inside the shrunk
graph, and so on. Because each `R_+` is a regular cost monoid, every surviving
`s–t` path is optimal for every rule processed, so the result is the
**lexicographically σ-minimal** solution — at a cost of `2N` Dijkstras.

Doing this for *every* linear extension and taking the union gives `topolex`.

**Soundness lemma.** *If `y <_R x` then `y <_lex,σ x` for every linear extension σ
of the rulebook.*

> Let `j` be the first σ-index at which `x` and `y` differ, and suppose
> `x_j < y_j`. Since `y ≲_R x` and `y_j > x_j`, some strictly more important rule
> `r_i > r_j` has `y_i < x_i`. But `r_i` strictly above `r_j` means `r_i`'s class is
> a strict ancestor of `r_j`'s, so `r_i` precedes `r_j` in σ — and all earlier
> indices are tied. Contradiction. So `y_j < x_j`. ∎

**Consequence.** A lex-σ-minimal solution cannot be strictly rule-dominated, so
**every vector `topolex` returns is genuinely rulebook-optimal**. `topolex`
is a *sound but incomplete* frontier generator: it returns the "lex-supported"
corner points of `P_R` and can miss interior ones, exactly as scalarisation
misses non-supported Pareto points.

This is a different guarantee from RA\*pex's, and complementary:

- RA\*pex: **complete coverage**, elements *need not* be optimal (Definition 6 of
  the paper explicitly allows non-optimal members).
- `topolex`: **every element optimal**, coverage not guaranteed.

**Cost.** At most `M!` extensions for `M` quotient classes, `2N` Dijkstras each —
but extensions are explored as a DFS over the prefix tree, so a shared prefix pays
for its phases once. Cap it with `--max-ext` for an anytime variant.

---

## 3. `peel` — exact Dijkstra reduction of the dominant prefix

**Lemma.** *If the quotient DAG has a unique source class and that class holds
exactly one rule `r`, then every rulebook-optimal solution is `r`-optimal.*

> A unique source of a DAG reaches every node, so `r` is strictly more important
> than every other rule. Let `c*` be the best achievable `r`-cost and let `x` be a
> solution with `r(x) > c*`, `y` one with `r(y) = c*`. For any rule `r_j ≠ r` where
> `y_j > x_j`, the rule `r` itself excuses it (`r > r_j` and `r(y) < r(x)`), so
> `y ≲_R x`. And `x ̸≲_R y`, since `r(x) > r(y)` and nothing sits above `r`. Hence
> `y <_R x` and `x ∉ P_R`. ∎

So that rule can be discharged by **two Dijkstras** instead of a search: shrink to
its optimal subgraph, delete it from the rulebook, and repeat while the condition
still holds. What remains is a much smaller graph whose rulebook is topped by an
antichain — and only that residual needs multi-objective search.

This makes hierarchies of exactly the shape used in the RA\*pex paper's own
experiments (`r4 > r1, r2, r3`; `r1 > r2 > r3 > r4`) collapse almost entirely.
Validated against brute force: `peel-exact` reproduces the exact frontier on every
tiny instance.

**Where it stops.** A source class with two or more *equivalent* rules is a
genuine multi-objective subproblem, not a Dijkstra, so peeling halts there. Two or
more incomparable source classes likewise. Extending the peel to "Pareto-optimal
subgraph of a source antichain" is the obvious next step and is not implemented.

---

## 4. Layout

```
src/rulebook.hpp    quotient DAG, rule-dominance, eps-rule-dominance, extensions
src/graph.hpp       CSR graph, Dijkstra, optimal-subgraph reduction, instance IO
src/solvers.hpp     exact / rapex / topolex / peel / brute
src/main.cpp        CLI, JSON output
bench/gen_instances.py  instance suite (random, grid, layered DAG, tiny)
bench/validate.py       cross-validation against brute force
bench/run_bench.py      the sweep, with the wall-clock cap
bench/analyze.py        scoring and report generation
results/report.md       generated comparison
```

Build and run:

```sh
make
python3 bench/gen_instances.py --out instances --queries 3
python3 bench/validate.py                     # correctness first
python3 bench/run_bench.py --timeout 600      # the 10-minute cap
python3 bench/analyze.py
```

### The 10-minute cap

`--timeout 600` is enforced inside the search loop *and* as a subprocess kill.
`analyze.py` then **drops a case entirely — for every solver — as soon as any
solver exceeded the cap on it**, so every surviving comparison is like-for-like on
identical instances. Difficulty is monotone in graph size within a family, so
`run_bench.py` also stops running a solver on larger sizes once it has been capped
at a smaller one, and records those as capped without burning the wall clock.

### Metrics

Against the exact frontier, per surviving case:

- `runtime` — search time; heuristic preprocessing excluded, as in the RA\*pex paper
- `size` — solutions returned
- `sound` — every returned vector is genuinely rulebook-optimal
- `recall` — fraction of `P_R` returned exactly
- `covered` — fraction of `P_R` weakly rule-dominated by the result
- `eps*` — smallest uniform ε at which the result ε-rule-dominates all of `P_R`

`recall` is the metric `topolex` is expected to lose on and `eps*` the one it
should still do well on; `sound` is the metric RA\*pex is *not* expected to hold.

---

## 5. Correctness

`bench/validate.py` checks, on 12 tiny graphs × 7 rulebooks × 3 queries, against
independent brute-force path enumeration and an independent Python implementation
of the dominance relation:

- `exact` == brute force
- `peel-exact` == brute force (peeling is lossless)
- `topolex` ⊆ brute force (soundness), and how often it is also complete
- RA\*pex variants cover the frontier under rule-dominance

---

## 6. Open threads

1. **Peel past the antichain.** Replace the Dijkstra reduction with a
   Pareto-optimal-subgraph reduction when the source class is an antichain of
   several rules — would extend the exact polynomial prefix much further.
2. **`topolex` as a seed.** Its output is a set of certified-optimal points,
   available in milliseconds. Feeding them to `exact`/RA\*pex as an initial
   solution set should prune hard from the first expansion — not yet tried.
3. **Which points does `topolex` miss?** Characterising the non-lex-supported part
   of `P_R` would say when the cheap sound set is good enough on its own.
4. **Real road networks.** The suite is synthetic (random / grid / layered DAG).
   The DIMACS BAY roadmap used in both papers is not included here.
