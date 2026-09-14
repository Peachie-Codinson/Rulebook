// CLI driver. One binary, one algorithm per invocation, JSON on stdout so the
// Python harness can enforce the wall-clock cap from the outside and still get
// structured results back.
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include "graph.hpp"
#include "rulebook.hpp"
#include "solvers.hpp"

using namespace rbp;

namespace {

void usage() {
    std::cerr <<
        "usage: rbsearch --graph G --rules R --query Q --alg A [options]\n"
        "  --alg      exact | rapex | rapex-nodr | topolex | seed-exact |\n"
        "             peel-exact | peel-rapex | peel-seed-exact | peel-only |\n"
        "             brute\n"
        "  --timeout  seconds per query (default 600)\n"
        "  --eps      override every eps in the rules file with this scalar\n"
        "  --max-ext  cap on linear extensions for topolex (default 0 = all)\n"
        "  --seed     RNG seed for RA*pex merge tie-breaks (default 1)\n";
}

std::string jsonCosts(const std::vector<Cost> &cs) {
    std::string out = "[";
    for (size_t i = 0; i < cs.size(); ++i) {
        out += "[";
        for (size_t j = 0; j < cs[i].size(); ++j) {
            out += std::to_string(cs[i][j]);
            if (j + 1 < cs[i].size()) out += ",";
        }
        out += "]";
        if (i + 1 < cs.size()) out += ",";
    }
    return out + "]";
}

std::string arg(int argc, char **argv, const std::string &key,
                const std::string &def) {
    for (int i = 1; i + 1 < argc; ++i)
        if (key == argv[i]) return argv[i + 1];
    return def;
}

}  // namespace

int main(int argc, char **argv) {
    std::string gfile = arg(argc, argv, "--graph", "");
    std::string rfile = arg(argc, argv, "--rules", "");
    std::string qfile = arg(argc, argv, "--query", "");
    std::string alg = arg(argc, argv, "--alg", "");
    double timeout = std::stod(arg(argc, argv, "--timeout", "600"));
    std::string eps_override = arg(argc, argv, "--eps", "");
    long max_ext = std::stol(arg(argc, argv, "--max-ext", "0"));
    uint32_t seed = static_cast<uint32_t>(std::stoul(arg(argc, argv, "--seed", "1")));

    if (gfile.empty() || rfile.empty() || qfile.empty() || alg.empty()) {
        usage();
        return 2;
    }

    Graph g;
    Rulebook rb;
    Eps eps;
    std::vector<std::pair<uint32_t, uint32_t>> queries;
    if (!loadGraph(gfile, g) || !loadRulebook(rfile, rb, eps) ||
        !loadQueries(qfile, queries))
        return 2;

    if (!eps_override.empty())
        std::fill(eps.begin(), eps.end(), std::stod(eps_override));

    if (rb.numRules() != g.num_rules) {
        std::cerr << "rule count mismatch: graph has " << g.num_rules
                  << ", rulebook has " << rb.numRules() << "\n";
        return 2;
    }

    // Heuristic preprocessing is excluded from the reported runtime, matching
    // the convention in the RA*pex paper.
    std::vector<Cost> h;
    bool needs_h = (alg == "rapex" || alg == "rapex-nodr");

    std::cout << "{\"alg\":\"" << alg << "\",\"graph\":\"" << gfile
              << "\",\"rules\":\"" << rfile << "\",\"n\":" << g.n
              << ",\"m\":" << g.edges.size() << ",\"N\":" << g.num_rules
              << ",\"classes\":" << rb.numClasses() << ",\"eps\":" << eps[0]
              << ",\"results\":[";

    bool first = true;
    for (auto [s, t] : queries) {
        if (needs_h) h = g.idealHeuristic(t);
        auto budget = makeBudget(timeout);
        std::fill(g.active.begin(), g.active.end(), 1);

        Result r;
        if (alg == "exact") {
            r = exactFrontier(g, rb, s, t, budget);
        } else if (alg == "rapex") {
            r = raPex(g, rb, s, t, eps, true, h, budget, seed);
        } else if (alg == "rapex-nodr") {
            r = raPex(g, rb, s, t, eps, false, h, budget, seed);
        } else if (alg == "topolex") {
            r = topoLex(g, rb, s, t, max_ext, budget);
        } else if (alg == "seed-exact") {
            r = seededExact(g, rb, s, t, max_ext, budget);
        } else if (alg == "peel-seed-exact") {
            r = peel(g, rb, s, t, Residual::SeededExact, eps, true, budget);
        } else if (alg == "peel-exact") {
            r = peel(g, rb, s, t, Residual::Exact, eps, true, budget);
        } else if (alg == "peel-rapex") {
            r = peel(g, rb, s, t, Residual::RApex, eps, true, budget);
        } else if (alg == "peel-only") {
            r = peel(g, rb, s, t, Residual::None, eps, true, budget);
        } else if (alg == "brute") {
            r = bruteForce(g, rb, s, t, 2000000, budget);
        } else {
            usage();
            return 2;
        }

        if (!first) std::cout << ",";
        first = false;
        std::cout << std::setprecision(9) << std::fixed
                  << "{\"s\":" << s << ",\"t\":" << t
                  << ",\"runtime\":" << r.runtime
                  << ",\"timed_out\":" << (r.timed_out ? "true" : "false")
                  << ",\"num_solutions\":" << r.costs.size()
                  << ",\"expansions\":" << r.expansions
                  << ",\"generations\":" << r.generations
                  << ",\"dijkstras\":" << r.dijkstras
                  << ",\"extensions\":" << r.extensions
                  << ",\"residual_edges\":" << r.residual_edges
                  << ",\"peeled_rules\":" << r.peeled_rules
                  << ",\"costs\":" << jsonCosts(r.costs) << "}";
        std::cout.flush();
    }
    std::cout << "]}" << std::endl;
    return 0;
}
