// Rulebook algebra: quotient DAG, rule-dominance, eps-rule-dominance,
// topological orders and linear extensions.
//
// Semantics follow Muhammetkulyyev et al., "Approximate Multi-Objective Search
// Under Rulebooks" (Def. 5 and Def. 8) and are cross-checked against the
// reference implementation in Infus3d/Rulebook_approximation
// (src/Utils/RulebookGraph.cpp).
#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace rbp {

using Cost = std::vector<int64_t>;
using Eps = std::vector<double>;

// A rulebook: N rules plus a preorder "is strictly more important than".
// Equivalent rules (mutually reachable) collapse into one quotient class; the
// classes form a DAG that is processed in topological order.
class Rulebook {
  public:
    explicit Rulebook(size_t num_rules = 0)
        : num_rules_(num_rules), adj_(num_rules), live_(num_rules, 1) {}

    size_t numRules() const { return num_rules_; }
    size_t numClasses() const { return classes_.size(); }

    // `from` is strictly more important than `to`.
    void addRelation(size_t from, size_t to) {
        size_t need = std::max(from, to) + 1;
        if (need > num_rules_) {
            num_rules_ = need;
            adj_.resize(num_rules_);
            live_.resize(num_rules_, 1);
        }
        adj_[from].push_back(to);
    }

    // Declares that two rules have equal priority (mutual reachability makes
    // them land in the same SCC, i.e. the same equivalence class).
    void addEquivalence(size_t a, size_t b) {
        addRelation(a, b);
        addRelation(b, a);
    }

    // Builds the quotient DAG (SCC condensation) and its topological order.
    void build() {
        stronglyConnectedComponents();
        topologicalSortClasses();

        ordered_rules_.clear();
        for (size_t c : ordered_classes_)
            for (size_t r : classes_[c])
                ordered_rules_.push_back(r);

        // Descendant closure per class, used to excuse degradations below a
        // class at which one side strictly improves.
        descendants_.assign(numClasses(), std::vector<char>(numClasses(), 0));
        for (size_t c = 0; c < numClasses(); ++c)
            markDescendants(c, descendants_[c]);

        sources_.clear();
        std::vector<int> indeg(numClasses(), 0);
        for (size_t c = 0; c < numClasses(); ++c)
            for (size_t d : class_adj_[c]) indeg[d]++;
        for (size_t c = 0; c < numClasses(); ++c)
            if (indeg[c] == 0) sources_.push_back(c);

        built_ = true;
    }

    bool isBuilt() const { return built_; }

    // Rules in non-increasing priority (a linear extension of the quotient DAG,
    // flattened). This is the order RA*pex uses to sort its priority queue.
    const std::vector<size_t> &orderedRules() const { return ordered_rules_; }
    const std::vector<size_t> &orderedClasses() const { return ordered_classes_; }
    const std::vector<size_t> &classRules(size_t c) const { return classes_[c]; }
    const std::vector<size_t> &classChildren(size_t c) const { return class_adj_[c]; }
    size_t classOf(size_t rule) const { return class_of_[rule]; }
    const std::vector<size_t> &sourceClasses() const { return sources_; }

    // ---------------------------------------------------------------------
    // Dominance
    // ---------------------------------------------------------------------

    // Weak rule-dominance a <~_R b (Def. 8): every rule at which `a` is worse
    // than `b` is excused by a strictly more important rule at which `a` is
    // strictly better.
    bool dominates(const Cost &a, const Cost &b) const {
        return dominatesEps(a, b, nullptr);
    }

    // eps-rule-dominance a <~^eps_R b (Def. 5): a_j <= (1+eps_j) b_j counts as
    // "not worse".
    bool dominatesEps(const Cost &a, const Cost &b, const Eps &eps) const {
        return dominatesEps(a, b, &eps);
    }

    // Strict rule-dominance a <_R b.
    bool strictlyDominates(const Cost &a, const Cost &b) const {
        return dominates(a, b) && !dominates(b, a);
    }

    // Smallest uniform eps >= 0 such that a eps-rule-dominates b, or -1 if no
    // finite eps works. Used to score how well an approximate set covers the
    // exact frontier.
    double minEpsToDominate(const Cost &a, const Cost &b) const {
        // eps only ever helps at rules where a is worse, and the candidate
        // thresholds are the ratios a_j/b_j at those rules.
        std::vector<double> cands{0.0};
        for (size_t j = 0; j < num_rules_; ++j) {
            if (a[j] > b[j]) {
                if (b[j] == 0) continue;  // unreachable with finite eps
                cands.push_back(static_cast<double>(a[j]) / static_cast<double>(b[j]) - 1.0);
            }
        }
        std::sort(cands.begin(), cands.end());
        cands.erase(std::unique(cands.begin(), cands.end()), cands.end());
        for (double c : cands) {
            Eps e(num_rules_, c + 1e-12);
            if (dominatesEps(a, b, e)) return c;
        }
        return -1.0;
    }

    // ---------------------------------------------------------------------
    // Linear extensions
    // ---------------------------------------------------------------------

    // Enumerates linear extensions of the quotient DAG (each an ordering of
    // class ids), stopping after `cap` of them. Cap <= 0 means unbounded.
    std::vector<std::vector<size_t>> linearExtensions(long cap = -1) const {
        std::vector<std::vector<size_t>> out;
        size_t M = numClasses();
        std::vector<int> indeg(M, 0);
        for (size_t c = 0; c < M; ++c)
            for (size_t d : class_adj_[c]) indeg[d]++;
        std::vector<size_t> cur;
        std::vector<char> used(M, 0);
        enumerate(indeg, used, cur, out, cap);
        return out;
    }

    // Flattens an ordering of classes into an ordering of rules.
    std::vector<size_t> flatten(const std::vector<size_t> &class_order) const {
        std::vector<size_t> out;
        out.reserve(num_rules_);
        for (size_t c : class_order)
            for (size_t r : classes_[c]) out.push_back(r);
        return out;
    }

    // True when the remaining rulebook can be peeled one Dijkstra phase at a
    // time: a unique source class that holds exactly one rule is strictly more
    // important than every other remaining rule, so every rulebook-optimal
    // solution must be optimal for it.
    bool hasUniqueSingletonTop() const {
        return sources_.size() == 1 && classes_[sources_[0]].size() == 1;
    }

    // The rulebook with `rule` deleted (relations through it are transitively
    // reconnected).  Used by the peeling solver.
    Rulebook without(size_t rule) const {
        Rulebook rb(0);
        rb.num_rules_ = num_rules_;  // rule ids stay stable so cost indices do too
        rb.adj_.assign(num_rules_, {});
        rb.live_ = live_;
        rb.live_[rule] = 0;
        for (size_t r = 0; r < num_rules_; ++r) {
            if (!rb.live_[r]) continue;
            for (size_t s : adj_[r])
                if (rb.live_[s]) rb.adj_[r].push_back(s);
        }
        // Transitive reconnection through the removed rule, so priorities it
        // used to mediate are preserved.
        for (size_t p = 0; p < num_rules_; ++p) {
            if (!rb.live_[p]) continue;
            if (std::find(adj_[p].begin(), adj_[p].end(), rule) == adj_[p].end())
                continue;
            for (size_t s : adj_[rule])
                if (rb.live_[s] && s != p) rb.adj_[p].push_back(s);
        }
        rb.dropped_ = dropped_;
        rb.dropped_.push_back(rule);
        rb.build();
        return rb;
    }

    // Rules removed by `without`, in removal order.
    const std::vector<size_t> &droppedRules() const { return dropped_; }

    // Rules still present in the quotient DAG.
    std::vector<size_t> liveRules() const {
        std::vector<size_t> out;
        for (size_t r = 0; r < num_rules_; ++r)
            if (live_[r]) out.push_back(r);
        return out;
    }

    bool isLive(size_t rule) const { return live_[rule] != 0; }

  private:
    bool dominatesEps(const Cost &a, const Cost &b, const Eps *eps) const {
        assert(built_ && "Rulebook::build() must be called first");
        std::vector<char> excused(numClasses(), 0);
        for (size_t c : ordered_classes_) {
            if (excused[c]) continue;
            bool strictly_better = false;
            for (size_t r : classes_[c]) {
                double lim = eps ? (1.0 + (*eps)[r]) * static_cast<double>(b[r])
                                 : static_cast<double>(b[r]);
                double av = static_cast<double>(a[r]);
                if (av > lim + kTol) return false;   // worse, and nothing excuses it
                if (av < lim - kTol) strictly_better = true;
            }
            if (strictly_better) {
                // A strict improvement here excuses every strictly less
                // important rule.
                const auto &d = descendants_[c];
                for (size_t k = 0; k < numClasses(); ++k)
                    if (d[k]) excused[k] = 1;
            }
        }
        return true;
    }

    void enumerate(std::vector<int> &indeg, std::vector<char> &used,
                   std::vector<size_t> &cur,
                   std::vector<std::vector<size_t>> &out, long cap) const {
        if (cap > 0 && static_cast<long>(out.size()) >= cap) return;
        if (cur.size() == numClasses()) {
            out.push_back(cur);
            return;
        }
        for (size_t c = 0; c < numClasses(); ++c) {
            if (used[c] || indeg[c] != 0) continue;
            used[c] = 1;
            for (size_t d : class_adj_[c]) indeg[d]--;
            cur.push_back(c);
            enumerate(indeg, used, cur, out, cap);
            cur.pop_back();
            for (size_t d : class_adj_[c]) indeg[d]++;
            used[c] = 0;
            if (cap > 0 && static_cast<long>(out.size()) >= cap) return;
        }
    }

    void dfsOrder(size_t v, const std::vector<std::vector<size_t>> &g,
                  std::vector<char> &vis, std::vector<size_t> &out) const {
        // Iterative to keep deep rule graphs off the C stack.
        std::vector<std::pair<size_t, size_t>> stk{{v, 0}};
        vis[v] = 1;
        while (!stk.empty()) {
            auto &[u, i] = stk.back();
            if (i < g[u].size()) {
                size_t w = g[u][i++];
                if (!vis[w]) {
                    vis[w] = 1;
                    stk.push_back({w, 0});
                }
            } else {
                out.push_back(u);
                stk.pop_back();
            }
        }
    }

    void stronglyConnectedComponents() {
        size_t n = num_rules_;
        classes_.clear();
        class_adj_.clear();
        class_of_.assign(n, 0);

        std::vector<size_t> order;
        std::vector<char> vis(n, 0);
        for (size_t i = 0; i < n; ++i)
            if (!live_[i]) vis[i] = 1;  // dead rules never enter a class
        for (size_t i = 0; i < n; ++i)
            if (!vis[i]) dfsOrder(i, adj_, vis, order);

        std::vector<std::vector<size_t>> rev(n);
        for (size_t v = 0; v < n; ++v)
            for (size_t u : adj_[v]) rev[u].push_back(v);

        vis.assign(n, 0);
        for (size_t i = 0; i < n; ++i)
            if (!live_[i]) vis[i] = 1;
        std::reverse(order.begin(), order.end());
        for (size_t v : order) {
            if (vis[v]) continue;
            std::vector<size_t> comp;
            dfsOrder(v, rev, vis, comp);
            size_t id = classes_.size();
            for (size_t u : comp) class_of_[u] = id;
            std::sort(comp.begin(), comp.end());
            classes_.push_back(comp);
        }

        class_adj_.assign(classes_.size(), {});
        for (size_t v = 0; v < n; ++v) {
            if (!live_[v]) continue;
            for (size_t u : adj_[v])
                if (live_[u] && class_of_[v] != class_of_[u])
                    class_adj_[class_of_[v]].push_back(class_of_[u]);
        }
        for (auto &e : class_adj_) {
            std::sort(e.begin(), e.end());
            e.erase(std::unique(e.begin(), e.end()), e.end());
        }
    }

    void topologicalSortClasses() {
        size_t M = classes_.size();
        std::vector<char> vis(M, 0);
        ordered_classes_.clear();
        for (size_t i = 0; i < M; ++i)
            if (!vis[i]) dfsOrder(i, class_adj_, vis, ordered_classes_);
        std::reverse(ordered_classes_.begin(), ordered_classes_.end());
    }

    void markDescendants(size_t c, std::vector<char> &mark) const {
        if (mark[c]) return;
        mark[c] = 1;
        for (size_t d : class_adj_[c]) markDescendants(d, mark);
    }

    static constexpr double kTol = 1e-9;

    size_t num_rules_;
    std::vector<std::vector<size_t>> adj_;  // strict-priority edges
    std::vector<std::vector<size_t>> classes_;
    std::vector<std::vector<size_t>> class_adj_;
    std::vector<size_t> class_of_;
    std::vector<size_t> ordered_classes_;
    std::vector<size_t> ordered_rules_;
    std::vector<std::vector<char>> descendants_;
    std::vector<size_t> sources_;
    std::vector<size_t> dropped_;
    std::vector<char> live_;
    bool built_ = false;
};

// Lexicographic comparison of cost vectors under a rule permutation.
inline int lexCompare(const Cost &a, const Cost &b, const std::vector<size_t> &order) {
    for (size_t r : order) {
        if (a[r] != b[r]) return a[r] < b[r] ? -1 : 1;
    }
    return 0;
}

// Keeps only the rule-dominance-minimal (i.e. rulebook-optimal) vectors, after
// de-duplication.
inline std::vector<Cost> filterOptimal(const Rulebook &rb, std::vector<Cost> v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    std::vector<Cost> out;
    for (size_t i = 0; i < v.size(); ++i) {
        bool dominated = false;
        for (size_t j = 0; j < v.size() && !dominated; ++j) {
            if (i == j) continue;
            if (rb.strictlyDominates(v[j], v[i])) dominated = true;
        }
        if (!dominated) out.push_back(v[i]);
    }
    return out;
}

}  // namespace rbp
