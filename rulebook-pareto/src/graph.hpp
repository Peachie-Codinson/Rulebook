// Graph representation, instance IO, and the single-rule Dijkstra primitives
// the phase-based solvers are built out of.
#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "rulebook.hpp"

namespace rbp {

constexpr int64_t kInf = std::numeric_limits<int64_t>::max() / 4;

struct Edge {
    uint32_t from;
    uint32_t to;
    Cost cost;
};

// Forward/backward CSR over an *active subset* of the edge list. The phase
// solvers never copy the graph; they just flip edges off in `active`.
class Graph {
  public:
    uint32_t n = 0;             // number of vertices (ids 0..n-1)
    size_t num_rules = 0;
    std::vector<Edge> edges;

    void build() {
        out_start.assign(n + 1, 0);
        in_start.assign(n + 1, 0);
        for (const auto &e : edges) {
            out_start[e.from + 1]++;
            in_start[e.to + 1]++;
        }
        for (uint32_t v = 0; v < n; ++v) {
            out_start[v + 1] += out_start[v];
            in_start[v + 1] += in_start[v];
        }
        out_idx.resize(edges.size());
        in_idx.resize(edges.size());
        auto ocur = out_start, icur = in_start;
        for (size_t i = 0; i < edges.size(); ++i) {
            out_idx[ocur[edges[i].from]++] = static_cast<uint32_t>(i);
            in_idx[icur[edges[i].to]++] = static_cast<uint32_t>(i);
        }
        active.assign(edges.size(), 1);
    }

    // Edge-index ranges. Callers must skip inactive edges themselves.
    const uint32_t *outBegin(uint32_t v) const { return out_idx.data() + out_start[v]; }
    const uint32_t *outEnd(uint32_t v) const { return out_idx.data() + out_start[v + 1]; }
    const uint32_t *inBegin(uint32_t v) const { return in_idx.data() + in_start[v]; }
    const uint32_t *inEnd(uint32_t v) const { return in_idx.data() + in_start[v + 1]; }

    std::vector<char> active;  // per-edge on/off, used by subgraph reduction

    size_t numActive() const {
        size_t k = 0;
        for (char c : active) k += (c != 0);
        return k;
    }

    // Single-rule Dijkstra over the active subgraph.
    // `reverse` walks in-edges instead of out-edges (cost-to-go).
    std::vector<int64_t> dijkstra(uint32_t src, size_t rule, bool reverse) const {
        std::vector<int64_t> dist(n, kInf);
        using QE = std::pair<int64_t, uint32_t>;
        std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pq;
        dist[src] = 0;
        pq.push({0, src});
        while (!pq.empty()) {
            auto [d, v] = pq.top();
            pq.pop();
            if (d > dist[v]) continue;
            const uint32_t *b = reverse ? inBegin(v) : outBegin(v);
            const uint32_t *e = reverse ? inEnd(v) : outEnd(v);
            for (const uint32_t *it = b; it != e; ++it) {
                if (!active[*it]) continue;
                const Edge &ed = edges[*it];
                uint32_t u = reverse ? ed.from : ed.to;
                int64_t nd = d + ed.cost[rule];
                if (nd < dist[u]) {
                    dist[u] = nd;
                    pq.push({nd, u});
                }
            }
        }
        return dist;
    }

    // Componentwise ideal-point heuristic: h[v][r] = cheapest cost from v to
    // `target` under rule r alone. Consistent, which RA*pex requires.
    std::vector<Cost> idealHeuristic(uint32_t target) const {
        std::vector<Cost> h(n, Cost(num_rules, 0));
        for (size_t r = 0; r < num_rules; ++r) {
            auto d = dijkstra(target, r, /*reverse=*/true);
            for (uint32_t v = 0; v < n; ++v) h[v][r] = d[v];
        }
        return h;
    }

    // Algorithm 1 of Slutsky et al.: restrict `active` to the edges that lie on
    // some rule-`rule`-optimal s-t path. Returns that optimal cost, or kInf if
    // the target is unreachable (in which case `active` is left alone).
    int64_t reduceToOptimalSubgraph(uint32_t s, uint32_t t, size_t rule) {
        auto come = dijkstra(s, rule, false);
        auto go = dijkstra(t, rule, true);
        int64_t c = come[t];
        if (c >= kInf) return kInf;
        for (size_t i = 0; i < edges.size(); ++i) {
            if (!active[i]) continue;
            const Edge &e = edges[i];
            if (come[e.from] >= kInf || go[e.to] >= kInf ||
                come[e.from] + e.cost[rule] + go[e.to] != c) {
                active[i] = 0;
            }
        }
        return c;
    }

  private:
    std::vector<size_t> out_start, in_start;
    std::vector<uint32_t> out_idx, in_idx;
};

// ---------------------------------------------------------------------------
// Instance IO
// ---------------------------------------------------------------------------

// Text instance format (one self-contained file, unlike the split DIMACS .gr
// layout the reference implementation uses):
//
//   n m N            vertices, edges, rules
//   <m lines>        u v c_0 c_1 ... c_{N-1}      (0-indexed vertices)
//
// The rulebook and the queries live in separate files so one graph can be
// reused across rulebooks.
inline bool loadGraph(const std::string &path, Graph &g) {
    std::ifstream f(path);
    if (!f) {
        std::cerr << "cannot open graph file " << path << "\n";
        return false;
    }
    size_t m;
    f >> g.n >> m >> g.num_rules;
    g.edges.resize(m);
    for (size_t i = 0; i < m; ++i) {
        Edge &e = g.edges[i];
        f >> e.from >> e.to;
        e.cost.resize(g.num_rules);
        for (size_t r = 0; r < g.num_rules; ++r) f >> e.cost[r];
    }
    g.build();
    return true;
}

// Rulebook file, matching the reference repo's layout:
//
//   N                     number of rules
//   eps_0 ... eps_{N-1}
//   K                     number of equivalence classes
//   <K lines>             size followed by that many rule ids
//   P                     number of priority relations
//   <P lines>             "a b"  meaning rule a is strictly more important
//
inline bool loadRulebook(const std::string &path, Rulebook &rb, Eps &eps) {
    std::ifstream f(path);
    if (!f) {
        std::cerr << "cannot open rules file " << path << "\n";
        return false;
    }
    size_t N;
    f >> N;
    eps.resize(N);
    for (size_t i = 0; i < N; ++i) f >> eps[i];
    rb = Rulebook(N);

    size_t K;
    f >> K;
    for (size_t i = 0; i < K; ++i) {
        size_t cnt;
        f >> cnt;
        std::vector<size_t> cls(cnt);
        for (size_t j = 0; j < cnt; ++j) f >> cls[j];
        for (size_t j = 1; j < cnt; ++j) rb.addEquivalence(cls[0], cls[j]);
    }
    size_t P;
    f >> P;
    for (size_t i = 0; i < P; ++i) {
        size_t a, b;
        f >> a >> b;
        rb.addRelation(a, b);
    }
    rb.build();
    return true;
}

inline bool loadQueries(const std::string &path,
                        std::vector<std::pair<uint32_t, uint32_t>> &q) {
    std::ifstream f(path);
    if (!f) {
        std::cerr << "cannot open query file " << path << "\n";
        return false;
    }
    uint32_t s, t;
    while (f >> s >> t) q.push_back({s, t});
    return true;
}

}  // namespace rbp
