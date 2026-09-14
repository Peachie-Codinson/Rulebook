// The angles being compared.
//
//   exact     label-setting ground truth for the rulebook-optimal frontier
//   rapex     RA*pex (Muhammetkulyyev et al.), with and without dim. reduction
//   topolex   NEW: multi-phase Dijkstra over linear extensions of the rulebook
//   peel      NEW: Dijkstra optimal-subgraph peeling + residual search
//   brute     path enumeration, only for validating the above on toy graphs
#pragma once

#include <chrono>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>

#include "graph.hpp"
#include "rulebook.hpp"

namespace rbp {

using Clock = std::chrono::steady_clock;

struct Result {
    std::string alg;
    std::vector<Cost> costs;
    double runtime = 0.0;
    bool timed_out = false;
    uint64_t expansions = 0;
    uint64_t generations = 0;
    uint64_t dijkstras = 0;      // phase-based solvers
    uint64_t extensions = 0;     // linear extensions explored
    size_t residual_edges = 0;   // edges left after peeling
    size_t peeled_rules = 0;
    std::string note;
};

struct Budget {
    Clock::time_point deadline;
    bool expired() const { return Clock::now() >= deadline; }
};

inline Budget makeBudget(double seconds) {
    return Budget{Clock::now() +
                  std::chrono::milliseconds(static_cast<long long>(seconds * 1000))};
}

// ---------------------------------------------------------------------------
// Non-dominated cost set under rule-dominance
// ---------------------------------------------------------------------------
//
// Weak rule-dominance is antisymmetric on distinct vectors (if u <~ v and
// v <~ u then u == v: a degradation in a maximal class can never be excused, so
// the two agree there, and induction carries it down the DAG). That makes weak
// dominance safe to prune with -- it never discards a rulebook-optimal vector.
class CostSet {
  public:
    explicit CostSet(const Rulebook &rb) : rb_(&rb) {}

    // Returns false if `c` was dominated by (or equal to) something already in.
    bool insert(const Cost &c) {
        for (const auto &e : items_)
            if (rb_->dominates(e, c)) return false;
        items_.erase(std::remove_if(items_.begin(), items_.end(),
                                    [&](const Cost &e) { return rb_->dominates(c, e); }),
                     items_.end());
        items_.push_back(c);
        return true;
    }

    const std::vector<Cost> &items() const { return items_; }
    size_t size() const { return items_.size(); }

  private:
    const Rulebook *rb_;
    std::vector<Cost> items_;
};

// ---------------------------------------------------------------------------
// Ground truth: exact rulebook-optimal frontier by label setting
// ---------------------------------------------------------------------------
//
// Open is ordered lexicographically in the rulebook's topological rule order.
// That is sound because strict rule-dominance implies strict lexicographic
// dominance under *every* linear extension: at the first differing rule, a
// degradation cannot be excused, since everything strictly more important is
// already tied. Rule-dominance is also translation invariant (adding the same
// edge cost to both sides preserves every componentwise comparison), so
// pruning partial labels never discards an optimal solution.
//
// `seed` pre-loads the goal set with solution costs already known to be
// achievable and rulebook-optimal (TopoLex produces exactly that, in
// milliseconds). Pruning against them is sound: if a seed c weakly
// rule-dominates a partial label's g, then it also weakly rule-dominates
// g plus any suffix, because costs only grow -- a rule at which c is worse than
// g + suffix is a rule at which c is worse than g, and the same more-important
// rule still excuses it.
inline Result exactFrontier(const Graph &g, const Rulebook &rb, uint32_t s,
                            uint32_t t, const Budget &budget,
                            const std::vector<Cost> *seed = nullptr) {
    Result res;
    res.alg = "exact";
    auto t0 = Clock::now();
    const auto &order = rb.orderedRules();

    struct Label {
        Cost g;
        uint32_t v;
        uint64_t id;
    };
    struct Cmp {
        const std::vector<size_t> *order;
        bool operator()(const Label &a, const Label &b) const {
            return lexCompare(a.g, b.g, *order) > 0;  // min-heap
        }
    };
    std::priority_queue<Label, std::vector<Label>, Cmp> open{Cmp{&order}};

    std::vector<std::vector<Cost>> closed(g.n);
    CostSet goal(rb);
    if (seed)
        for (const auto &c : *seed) goal.insert(c);

    Cost zero(g.num_rules, 0);
    open.push({zero, s, 0});
    uint64_t iter = 0;

    while (!open.empty()) {
        if ((++iter & 0x3FF) == 0 && budget.expired()) {
            res.timed_out = true;
            break;
        }
        Label cur = open.top();
        open.pop();

        // Was this label superseded while it sat in the queue?
        bool stale = false;
        for (const auto &e : closed[cur.v])
            if (rb.dominates(e, cur.g)) { stale = true; break; }
        if (stale) continue;
        for (const auto &e : goal.items())
            if (rb.dominates(e, cur.g)) { stale = true; break; }
        if (stale) continue;

        closed[cur.v].push_back(cur.g);
        res.expansions++;

        if (cur.v == t) {
            goal.insert(cur.g);
            continue;
        }

        for (const uint32_t *it = g.outBegin(cur.v); it != g.outEnd(cur.v); ++it) {
            if (!g.active[*it]) continue;
            const Edge &e = g.edges[*it];
            Cost ng = cur.g;
            for (size_t r = 0; r < g.num_rules; ++r) ng[r] += e.cost[r];

            bool dom = false;
            for (const auto &c : closed[e.to])
                if (rb.dominates(c, ng)) { dom = true; break; }
            if (dom) continue;
            for (const auto &c : goal.items())
                if (rb.dominates(c, ng)) { dom = true; break; }
            if (dom) continue;

            res.generations++;
            open.push({std::move(ng), e.to, 0});
        }
    }

    res.costs = goal.items();
    res.runtime = std::chrono::duration<double>(Clock::now() - t0).count();
    return res;
}

// ---------------------------------------------------------------------------
// RA*pex
// ---------------------------------------------------------------------------

struct RPair {
    uint32_t v;
    Cost apex_g, apex_f;
    Cost real_g, real_f;
    bool active = true;
};
using RPairPtr = std::shared_ptr<RPair>;

// Truncated / residual projections used for dimensionality reduction: Tr drops
// the leading rule, alpha additionally drops every rule strictly below it.
struct DimReduce {
    std::vector<char> in_tr;     // per-rule mask
    std::vector<char> in_alpha;

    DimReduce(const Rulebook &rb) {
        size_t N = rb.numRules();
        in_tr.assign(N, 0);
        in_alpha.assign(N, 0);
        const auto &order = rb.orderedRules();
        if (order.empty()) return;
        size_t r1 = order[0];
        size_t c1 = rb.classOf(r1);
        for (size_t i = 1; i < order.size(); ++i) {
            size_t r = order[i];
            in_tr[r] = 1;
            // strictly below r1 iff reachable from r1's class in the quotient DAG
            in_alpha[r] = rb.classOf(r) == c1 ? 1 : 0;
        }
        // Mark everything strictly below r1 as excluded from alpha.
        std::vector<char> below(rb.numClasses(), 0);
        std::function<void(size_t)> mark = [&](size_t c) {
            for (size_t d : rb.classChildren(c)) {
                if (below[d]) continue;
                below[d] = 1;
                mark(d);
            }
        };
        mark(c1);
        for (size_t i = 1; i < order.size(); ++i) {
            size_t r = order[i];
            in_alpha[r] = below[rb.classOf(r)] ? 0 : 1;
        }
    }
};

// Projects a cost vector onto a rule mask by replacing masked-out rules with a
// value that makes them invisible to the dominance test (equal on both sides).
inline void projected(const Cost &a, const Cost &b, const std::vector<char> &mask,
                      Cost &pa, Cost &pb) {
    pa = a;
    pb = b;
    for (size_t r = 0; r < mask.size(); ++r)
        if (!mask[r]) { pa[r] = 0; pb[r] = 0; }
}

inline Result raPex(const Graph &g, const Rulebook &rb, uint32_t s, uint32_t t,
                    const Eps &eps, bool use_dr, const std::vector<Cost> &h,
                    const Budget &budget, uint32_t seed = 1) {
    Result res;
    res.alg = use_dr ? "rapex" : "rapex-nodr";
    auto t0 = Clock::now();
    const auto &order = rb.orderedRules();
    DimReduce dr(rb);
    std::mt19937 rng(seed);

    auto addF = [&](const Cost &gv, uint32_t v) {
        Cost f = gv;
        for (size_t r = 0; r < g.num_rules; ++r) f[r] += h[v][r];
        return f;
    };

    struct Cmp {
        const std::vector<size_t> *order;
        bool operator()(const RPairPtr &a, const RPairPtr &b) const {
            return lexCompare(a->apex_f, b->apex_f, *order) > 0;
        }
    };
    std::priority_queue<RPairPtr, std::vector<RPairPtr>, Cmp> open{Cmp{&order}};
    std::vector<std::list<RPairPtr>> open_at(g.n);

    // Closed sets. Without DR: one list per state. With DR: the G= / G< split.
    std::vector<std::list<RPairPtr>> gcl(g.n), gcl_eq(g.n), gcl_lt(g.n);
    std::vector<RPairPtr> solutions;

    // Merge `b` into `a` if the result stays eps-bounded (Alg. 4).
    auto tryMerge = [&](const RPairPtr &a, const RPairPtr &b) -> bool {
        Cost ng = a->apex_g, nf = a->apex_f;
        for (size_t r = 0; r < g.num_rules; ++r) {
            ng[r] = std::min(ng[r], b->apex_g[r]);
            nf[r] = std::min(nf[r], b->apex_f[r]);
        }
        bool a_ok = rb.dominatesEps(a->real_f, nf, eps);
        bool b_ok = rb.dominatesEps(b->real_f, nf, eps);
        if (!a_ok && !b_ok) return false;
        bool take_b = b_ok && (!a_ok || (rng() & 1));
        a->apex_g = std::move(ng);
        a->apex_f = std::move(nf);
        if (take_b) {
            a->real_g = b->real_g;
            a->real_f = b->real_f;
        }
        return true;
    };

    // Alg. 5 lines 1-6: once the extracted f_1 strictly increases at a state,
    // everything in G= becomes "strictly better in r1" and moves to G<.
    auto transfer = [&](uint32_t v, const RPairPtr &p) {
        size_t r1 = order[0];
        bool any_less = false;
        for (const auto &q : gcl_eq[v])
            if (q->apex_f[r1] < p->apex_f[r1]) { any_less = true; break; }
        if (!any_less) return;
        Cost pa, pb;
        gcl_lt[v].remove_if([&](const RPairPtr &old) {
            for (const auto &q : gcl_eq[v]) {
                projected(q->apex_f, old->apex_f, dr.in_alpha, pa, pb);
                if (rb.dominates(pa, pb)) return true;
            }
            return false;
        });
        for (auto &q : gcl_eq[v]) gcl_lt[v].push_back(q);
        gcl_eq[v].clear();
    };

    auto isDominated = [&](const RPairPtr &p, bool do_transfer) -> bool {
        uint32_t v = p->v;
        Cost pa, pb;
        if (use_dr) {
            if (do_transfer) transfer(v, p);
            for (const auto &q : gcl_eq[v]) {
                projected(q->apex_f, p->apex_f, dr.in_tr, pa, pb);
                if (rb.dominates(pa, pb)) return true;
            }
            for (const auto &q : gcl_lt[v]) {
                projected(q->apex_f, p->apex_f, dr.in_alpha, pa, pb);
                if (rb.dominates(pa, pb)) return true;
            }
        } else {
            for (const auto &q : gcl[v])
                if (rb.dominates(q->apex_f, p->apex_f)) return true;
        }
        // Global check against the solution set (Alg. 2 line 1 / Alg. 5 line 11).
        for (auto &sol : solutions) {
            if (rb.dominatesEps(sol->real_f, p->apex_f, eps)) {
                for (size_t r = 0; r < g.num_rules; ++r) {
                    sol->apex_g[r] = std::min(sol->apex_g[r], p->apex_g[r]);
                    sol->apex_f[r] = std::min(sol->apex_f[r], p->apex_f[r]);
                }
                return true;
            }
        }
        return false;
    };

    auto insertOpen = [&](const RPairPtr &p) {
        for (auto &q : open_at[p->v]) {
            if (!q->active) continue;
            if (tryMerge(p, q)) {
                if (p->apex_g != q->apex_g || p->real_g != q->real_g) {
                    q->active = false;
                    open.push(p);
                    open_at[p->v].push_back(p);
                }
                return;
            }
        }
        open.push(p);
        open_at[p->v].push_back(p);
    };

    auto insertSolution = [&](const RPairPtr &p) {
        for (auto &q : solutions)
            if (tryMerge(q, p)) return;
        solutions.push_back(p);
    };

    auto root = std::make_shared<RPair>();
    root->v = s;
    root->apex_g = root->real_g = Cost(g.num_rules, 0);
    root->apex_f = root->real_f = addF(root->apex_g, s);
    insertOpen(root);

    uint64_t iter = 0;
    while (!open.empty()) {
        if ((++iter & 0x3FF) == 0 && budget.expired()) {
            res.timed_out = true;
            break;
        }
        RPairPtr p = open.top();
        open.pop();
        res.generations++;
        if (!p->active) continue;
        if (isDominated(p, true)) continue;

        if (use_dr) gcl_eq[p->v].push_back(p);
        else gcl[p->v].push_back(p);
        res.expansions++;

        if (p->v == t) {
            insertSolution(p);
            continue;
        }

        for (const uint32_t *it = g.outBegin(p->v); it != g.outEnd(p->v); ++it) {
            if (!g.active[*it]) continue;
            const Edge &e = g.edges[*it];
            auto c = std::make_shared<RPair>();
            c->v = e.to;
            c->apex_g = p->apex_g;
            c->real_g = p->real_g;
            for (size_t r = 0; r < g.num_rules; ++r) {
                c->apex_g[r] += e.cost[r];
                c->real_g[r] += e.cost[r];
            }
            c->apex_f = addF(c->apex_g, e.to);
            c->real_f = addF(c->real_g, e.to);
            if (isDominated(c, false)) continue;
            insertOpen(c);
        }
    }

    for (auto &sol : solutions) res.costs.push_back(sol->real_g);
    res.runtime = std::chrono::duration<double>(Clock::now() - t0).count();
    return res;
}

// ---------------------------------------------------------------------------
// TopoLex: multi-phase Dijkstra over linear extensions
// ---------------------------------------------------------------------------
//
// For one linear extension sigma of the rulebook, running Algorithm 2 of
// Slutsky et al. (iterated optimal-subgraph reduction, one Dijkstra pair per
// rule, top to bottom) yields the lexicographically-sigma-minimal solution.
// Such a solution is always rulebook-optimal, because strict rule-dominance
// implies strict lexicographic dominance under every linear extension. So the
// union over extensions is a guaranteed-sound subset of the frontier -- it can
// miss points, but it never returns a non-optimal one.
//
// Extensions are explored as a DFS over the prefix tree so that a shared prefix
// only pays for its Dijkstra phases once.
inline Result topoLex(Graph &g, const Rulebook &rb, uint32_t s, uint32_t t,
                      long max_extensions, const Budget &budget) {
    Result res;
    res.alg = "topolex";
    auto t0 = Clock::now();

    const auto live = rb.liveRules();
    size_t N = live.size();

    // Rule r is available once every strictly more important rule is placed.
    std::vector<std::vector<size_t>> above(rb.numRules());
    for (size_t r : live) {
        size_t cr = rb.classOf(r);
        for (size_t q : live) {
            if (q == r) continue;
            size_t cq = rb.classOf(q);
            if (cq == cr) continue;
            // q strictly above r iff r is reachable from q in the quotient DAG
            std::vector<char> seen(rb.numClasses(), 0);
            std::function<bool(size_t)> reach = [&](size_t c) {
                if (c == cr) return true;
                if (seen[c]) return false;
                seen[c] = 1;
                for (size_t d : rb.classChildren(c))
                    if (reach(d)) return true;
                return false;
            };
            if (reach(cq)) above[r].push_back(q);
        }
    }

    // Rules already discharged elsewhere (by peeling) are not live in `rb`, but
    // they still occupy a slot in every cost vector. Every surviving s-t path
    // shares the same cost for such a rule -- that is exactly what peeling
    // established -- so its value is read back with one Dijkstra at each leaf.
    // Leaving them at zero would emit cost vectors that no path actually has.
    std::vector<size_t> dead;
    for (size_t r = 0; r < rb.numRules(); ++r)
        if (!rb.isLive(r)) dead.push_back(r);

    std::vector<char> placed(rb.numRules(), 0);
    Cost fixed(rb.numRules(), 0);
    CostSet found(rb);
    bool aborted = false;

    std::function<void(size_t)> dfs = [&](size_t depth) {
        if (aborted) return;
        if (budget.expired()) { aborted = true; res.timed_out = true; return; }
        if (max_extensions > 0 && static_cast<long>(res.extensions) >= max_extensions)
            return;
        if (depth == N) {
            res.extensions++;
            Cost leaf = fixed;
            for (size_t r : dead) {
                leaf[r] = g.dijkstra(s, r, false)[t];
                res.dijkstras++;
            }
            found.insert(leaf);
            return;
        }
        for (size_t r : live) {
            if (placed[r]) continue;
            bool ready = true;
            for (size_t q : above[r])
                if (!placed[q]) { ready = false; break; }
            if (!ready) continue;

            // Apply one phase: keep only edges on a rule-r-optimal path.
            std::vector<char> saved = g.active;
            int64_t c = g.reduceToOptimalSubgraph(s, t, r);
            res.dijkstras += 2;
            if (c >= kInf) {           // target unreachable: nothing to explore
                g.active = saved;
                return;
            }
            fixed[r] = c;
            placed[r] = 1;
            dfs(depth + 1);
            placed[r] = 0;
            g.active = std::move(saved);
            if (aborted) return;
            if (max_extensions > 0 &&
                static_cast<long>(res.extensions) >= max_extensions)
                return;
        }
    };

    dfs(0);
    res.costs = found.items();
    res.runtime = std::chrono::duration<double>(Clock::now() - t0).count();
    return res;
}

// Seeded exact search: spend a few milliseconds on TopoLex to get a handful of
// certified-optimal solutions, then let the exact search start with them
// already in its goal set so it prunes from the very first expansion.
inline Result seededExact(Graph &g, const Rulebook &rb, uint32_t s, uint32_t t,
                          long max_extensions, const Budget &budget) {
    Result res;
    res.alg = "seed-exact";
    auto t0 = Clock::now();

    std::vector<char> saved = g.active;
    Result seed = topoLex(g, rb, s, t, max_extensions, budget);
    g.active = saved;

    Result full = exactFrontier(g, rb, s, t, budget, &seed.costs);
    res.costs = full.costs;
    res.expansions = full.expansions;
    res.generations = full.generations;
    res.dijkstras = seed.dijkstras;
    res.extensions = seed.extensions;
    res.timed_out = seed.timed_out || full.timed_out;
    res.runtime = std::chrono::duration<double>(Clock::now() - t0).count();
    return res;
}

// ---------------------------------------------------------------------------
// Peel: exact Dijkstra reduction of the globally-dominant prefix
// ---------------------------------------------------------------------------
//
// While the remaining quotient DAG has a unique source class holding exactly
// one rule r, that rule is strictly more important than every other remaining
// rule (a unique source of a DAG reaches all of it). Any solution that is
// suboptimal for r is then strictly rule-dominated by an r-optimal one, since r
// excuses every degradation below it. So restricting to the r-optimal subgraph
// is exact, not a heuristic -- and it is two Dijkstras instead of a search.
//
// What is left is an antichain-topped residual problem on a far smaller graph,
// handed to whichever residual solver the caller picks.
enum class Residual { Exact, RApex, None, SeededExact };

inline Result peel(Graph &g, const Rulebook &rb_in, uint32_t s, uint32_t t,
                   Residual residual, const Eps &eps, bool use_dr,
                   const Budget &budget) {
    Result res;
    res.alg = residual == Residual::Exact        ? "peel-exact"
              : residual == Residual::RApex       ? "peel-rapex"
              : residual == Residual::SeededExact ? "peel-seed-exact"
                                                  : "peel-only";
    auto t0 = Clock::now();

    Rulebook rb = rb_in;
    std::vector<char> saved = g.active;

    while (rb.liveRules().size() > 0 && rb.hasUniqueSingletonTop()) {
        size_t r = rb.classRules(rb.sourceClasses()[0])[0];
        int64_t c = g.reduceToOptimalSubgraph(s, t, r);
        res.dijkstras += 2;
        if (c >= kInf) {  // unreachable
            g.active = saved;
            res.runtime = std::chrono::duration<double>(Clock::now() - t0).count();
            return res;
        }
        res.peeled_rules++;
        if (rb.liveRules().size() == 1) {
            rb = rb.without(r);
            break;
        }
        rb = rb.without(r);
    }
    res.residual_edges = g.numActive();

    if (rb.liveRules().empty()) {
        // Fully peeled: every surviving s-t path has the same cost vector.
        auto come = g.dijkstra(s, 0, false);
        Cost c(g.num_rules, 0);
        for (size_t r = 0; r < g.num_rules; ++r) {
            auto d = g.dijkstra(s, r, false);
            c[r] = d[t];
            res.dijkstras++;
        }
        res.costs.push_back(c);
    } else if (residual == Residual::Exact) {
        Result sub = exactFrontier(g, rb, s, t, budget);
        res.costs = sub.costs;
        res.expansions = sub.expansions;
        res.generations = sub.generations;
        res.timed_out = sub.timed_out;
    } else if (residual == Residual::SeededExact) {
        Result sub = seededExact(g, rb, s, t, 0, budget);
        res.costs = sub.costs;
        res.expansions = sub.expansions;
        res.generations = sub.generations;
        res.dijkstras += sub.dijkstras;
        res.extensions = sub.extensions;
        res.timed_out = sub.timed_out;
    } else if (residual == Residual::RApex) {
        auto h = g.idealHeuristic(t);
        Result sub = raPex(g, rb, s, t, eps, use_dr, h, budget);
        res.costs = sub.costs;
        res.expansions = sub.expansions;
        res.generations = sub.generations;
        res.timed_out = sub.timed_out;
    }

    g.active = std::move(saved);
    res.runtime = std::chrono::duration<double>(Clock::now() - t0).count();
    return res;
}

// ---------------------------------------------------------------------------
// Brute force validator (toy graphs only)
// ---------------------------------------------------------------------------
inline Result bruteForce(const Graph &g, const Rulebook &rb, uint32_t s,
                         uint32_t t, size_t max_paths, const Budget &budget) {
    Result res;
    res.alg = "brute";
    auto t0 = Clock::now();
    std::vector<Cost> all;
    std::vector<char> on_path(g.n, 0);
    Cost cur(g.num_rules, 0);

    std::function<void(uint32_t)> rec = [&](uint32_t v) {
        if (res.timed_out || all.size() > max_paths) return;
        if (budget.expired()) { res.timed_out = true; return; }
        if (v == t) {
            all.push_back(cur);
            return;
        }
        on_path[v] = 1;
        for (const uint32_t *it = g.outBegin(v); it != g.outEnd(v); ++it) {
            if (!g.active[*it]) continue;
            const Edge &e = g.edges[*it];
            if (on_path[e.to]) continue;
            for (size_t r = 0; r < g.num_rules; ++r) cur[r] += e.cost[r];
            rec(e.to);
            for (size_t r = 0; r < g.num_rules; ++r) cur[r] -= e.cost[r];
        }
        on_path[v] = 0;
    };
    rec(s);

    res.costs = filterOptimal(rb, all);
    res.note = "paths=" + std::to_string(all.size());
    res.runtime = std::chrono::duration<double>(Clock::now() - t0).count();
    return res;
}

}  // namespace rbp
