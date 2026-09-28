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
bench/ab_seed.py        contention-free A/B for the seeding question
bench/compare_rapex_topolex.py  RA*pex vs multi-level Dijkstra, two-way, no
                        ground truth needed -- see section 7
bench/fetch_dimacs.sh   grab the public DIMACS road networks
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

## 6. Findings

Full sweep: 294 cases (random / grid / layered-DAG graphs x 7 rulebook shapes),
10-minute per-case cap, 13 cases dropped, 281 kept. Generated report in
[`results/report.md`](results/report.md).

The short version: **the two cheap angles work exactly where they are not
needed, and the expensive case stays expensive.**

### TopoLex is sound, fast, and does not cover the frontier

On `flat4`, the hardest shape (mean frontier 949 solutions):

| solver | mean s | max s | size | sound | recall | eps\* |
|---|---|---|---|---|---|---|
| `exact` | 5.51 | 77.40 | 949 | 30/30 | 1.000 | 0.000 |
| `rapex` eps=0.01 | 4.03 | 62.55 | 561 | **23/30** | 0.773 | 0.010 |
| `topolex` | **0.0011** | 0.0075 | 4.0 | **30/30** | **0.158** | 1.552 |

The soundness lemma holds up empirically: `topolex` never returned a
non-optimal solution, in 294/294 benchmark cases and 126/126 brute-force
validation triples. RA\*pex at eps=0.01 returned a non-optimal solution in 7 of
30 `flat4` cases — which is not a bug, Definition 6 of the paper explicitly
permits it, but it is the concrete trade between the two guarantees.

The cost is coverage. TopoLex returns one solution per linear extension — 4 on a
flat 4-rule book, 2 on `diamond4` — so it recovers 16–26% of the frontier and
needs eps ~1.5 to cover the rest. It is a **certified-optimal sample**, not a
frontier approximation, and should not be sold as the latter.

### Peeling is exact and nearly free, and buys nothing

Peeling removes ~97% of edges where it applies:

| rulebook | rules peeled | edges left | frontier size |
|---|---|---|---|
| `chain4` | 4.00 | 2.8% | 1.0 |
| `chain2top` | 2.00 | 2.8% | 1.0 |
| `top1` | 1.00 | 2.9% | 1.0 |
| `diamond4` | **0.00** | 100% | 30.2 |
| `flat4` | **0.00** | 100% | 949.0 |

But look at the last column. A rule that is strictly above everything else
collapses the frontier to (nearly) a single solution, so the rulebooks peeling
can attack are exactly the ones that were already trivial — `exact` solves them
in under 2 ms. And the rulebooks that are actually hard, `diamond4` and the flat
ones, have no unique singleton top at all, so **zero** rules peel and the
residual is the whole problem.

This corrects an earlier reading of the 97% figure as a win. The reduction is
real and the lemma is sound — `peel-exact` reproduces the exact frontier on
every validated case — but on this instance suite it never converts into a
runtime saving, and `peel-exact` is in fact a hair slower than `exact` from the
reduction overhead.

The one place this could still pay off is a rulebook with a genuine hierarchy
*above* a wide antichain — a global top rule that does not collapse the frontier
because several incomparable rules sit below it. `chain2top` was meant to be
that case but its frontier still collapsed to 1. Constructing an instance family
where peeling meets a large residual frontier is the open question.

### Seeding the exact search does not work

Contention-free A/B, single-threaded, both arms back to back with alternating
order and 3 repeats, 30 paired cases: `seed-exact` was faster in **0/30**
(median 0.98x, min 0.78x, max 1.02x). See
[`results/ab_seed.csv`](results/ab_seed.csv).

The mechanism is structural. Seeding could only matter on the flat rulebooks,
whose quotient DAG is an antichain — so TopoLex hands over the N lexicographic
*corner* points, each minimal in one rule and large in the rest. Those dominate
almost nothing. Meanwhile the exact search's queue is already lexicographically
ordered, so it reaches the best of those points in its first expansion anyway.
The warm start pays for N Dijkstra phases and prunes nothing new.

### What this says about the original question

For the hard cases, RA\*pex remains the thing to beat; neither new angle
displaces it. What the new angles add is a **different guarantee**, not a faster
route to the same one:

- need every answer certified optimal, and a sample is enough → `topolex`, four
  orders of magnitude cheaper
- need coverage of the whole frontier → RA\*pex, accepting that some returned
  solutions are not optimal
- rulebook has a singleton global top → peel it, exactly, for free — but expect
  the frontier to be small anyway

---

## 7. Just RA\*pex vs multi-level Dijkstra, on the reference dataset

`bench/compare_rapex_topolex.py` runs only those two and scores them against
each other. No exact ground truth is needed, which is the point: on a road
network the size of BAY the exact frontier is not computable in any reasonable
budget, so the usual reference is unavailable.

What makes the comparison work anyway is that **TopoLex's output is itself a
certificate**. Every solution it returns is provably rulebook-optimal, so:

- an RA\*pex solution strictly rule-dominated by a TopoLex solution is **proven
  non-optimal**, with no frontier required;
- the eps RA\*pex needs to cover TopoLex's certified points is a **live test of
  RA\*pex's Theorem 1** — it must come out at or below the eps it was run with;
- the eps TopoLex needs to cover RA\*pex's set measures the price of TopoLex's
  incompleteness;
- a TopoLex solution strictly dominated by RA\*pex would contradict the
  soundness lemma, and is flagged loudly if it ever happens.

Not reported, deliberately: the plain (eps = 0) coverage fractions between the
two sets. Weak rule-dominance is antisymmetric on distinct vectors, so two
distinct *optimal* points can never weakly dominate one another — those
fractions are structurally zero and just restate `shared`.

### Getting the dataset

The reference repo ships code but no data. Two sources:

```sh
# 1. the public DIMACS road networks (objectives 1 and 2: distance, time)
bench/fetch_dimacs.sh BAY data/dimacs      # or NY for a quick smoke test

# 2. by hand, from the Google Drive folder linked in the reference repo README
#    https://github.com/Infus3d/Rulebook_approximation
#      USA-road-3.BAY.gr     third objective
#      USA-road-4.BAY.gr     fourth objective (BAY only)
#      BAY_instances.txt     queries, one "s,t" per line
#      4_rules_eps_0.01.txt  rulebook
#    -> drop them all in data/dimacs/
```

### Running it

```sh
make

python3 bench/compare_rapex_topolex.py \
  --gr data/dimacs/USA-road-d.BAY.gr,data/dimacs/USA-road-t.BAY.gr,data/dimacs/USA-road-3.BAY.gr,data/dimacs/USA-road-4.BAY.gr \
  --rules data/dimacs/4_rules_eps_0.01.txt \
  --query data/dimacs/BAY_instances.txt \
  --eps 0.01 \
  --timeout 600 \
  --out results/rapex_vs_topolex_BAY.csv
```

`--gr` takes one DIMACS file per objective, comma separated, in rule order — the
same convention as the reference implementation's `-m` flag. The rules file
format is unchanged from that repo, so its own rules files work as they are.
`--eps` applies to RA\*pex only; TopoLex has no approximation parameter.
`--max-ext N` caps TopoLex's linear extensions for an anytime variant (default:
all of them, which is `M!` for `M` quotient classes).

Either solver alone, JSON on stdout:

```sh
./bin/rbsearch --gr <files> --rules <rules> --query <queries> --alg rapex   --eps 0.01
./bin/rbsearch --gr <files> --rules <rules> --query <queries> --alg topolex --eps 0
```

### Shape of the output

From a 4-objective `diamond4` run on a synthetic DIMACS-format grid (1.6k nodes,
6.2k arcs) used to verify the path end to end:

```
       1 -> 1600   rapex   0.338s n=75   topolex 0.001s n=2   shared=0  rapex_bad=0
      40 -> 1561   rapex   0.185s n=72   topolex 0.001s n=2   shared=0  rapex_bad=0

topolex speedup over rapex               : median 210.0x
eps TopoLex needs to cover RA*pex's set  : median 0.549  (max 0.624)
eps RA*pex needs to cover certified pts  : max 0.009896 vs eps=0.01  (within guarantee)
RA*pex solutions PROVEN non-optimal      : 0 across 0/3 queries
TopoLex soundness violations             : 0 (as expected)
```

RA\*pex landing at 0.009896 against a 0.01 budget is Theorem 1 holding with
almost nothing to spare — which is what you would hope to see, and is worth
re-checking on BAY where the frontiers are far larger.

---

## 8. Open threads

1. **Peel past the antichain.** Replace the Dijkstra reduction with a
   Pareto-optimal-subgraph reduction when the source class is an antichain of
   several rules. This is the one change that would let peeling touch
   `diamond4` and the flat books, i.e. the cases that actually cost time.
2. **Instances where peeling meets a large frontier.** Every hierarchical
   rulebook here collapsed the frontier to ~1 solution, which hid whatever
   peeling is worth. Needs a construction where a global top rule admits many
   optimal solutions below it.
3. **More linear extensions, better coverage?** TopoLex currently returns one
   point per extension. Whether interpolating between corner points (or
   perturbing the phase order) recovers interior frontier points cheaply is
   untested.
4. **Characterise what TopoLex misses.** Which part of `P_R` is
   non-lex-supported, and can its size be predicted from the rulebook shape
   before running anything?
5. **Real road networks.** The suite is synthetic. The DIMACS BAY roadmap used
   in both papers is not included here, and RA\*pex is a reimplementation rather
   than the reference binary (no Boost in the build environment), so absolute
   runtimes are not comparable to the published figures — only the relative
   orderings measured here are.
