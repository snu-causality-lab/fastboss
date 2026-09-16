#include "fastboss.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace {

struct Args {
    int p = 100;
    int n = 1000;
    int degree = 5;
    unsigned int seed = 20260602U;
    double alpha = 2.0;
    int sweeps = 0;
    int use_gst = 1;
    int random_start = 1;
    double min_gain = 1e-6;
    int score_only = 0;
};

void usage(const char* argv0) {
    std::fprintf(stderr,
                 "usage: %s [--p N] [--n N] [--degree D] [--seed S] "
                 "[--alpha A] [--sweeps K] [--use-gst 0|1] "
                 "[--random-start 0|1] [--min-gain X] [--score-only 0|1]\n",
                 argv0);
}

int parse_int(const char* value, const char* name) {
    char* end = nullptr;
    const long out = std::strtol(value, &end, 10);
    if (!end || *end != '\0') {
        std::fprintf(stderr, "invalid integer for %s: %s\n", name, value);
        std::exit(2);
    }
    return static_cast<int>(out);
}

double parse_double(const char* value, const char* name) {
    char* end = nullptr;
    const double out = std::strtod(value, &end);
    if (!end || *end != '\0') {
        std::fprintf(stderr, "invalid double for %s: %s\n", name, value);
        std::exit(2);
    }
    return out;
}

Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (key == "--help" || key == "-h") {
            usage(argv[0]);
            std::exit(0);
        }
        if (i + 1 >= argc) {
            usage(argv[0]);
            std::exit(2);
        }
        const char* value = argv[++i];
        if (key == "--p") {
            args.p = parse_int(value, key.c_str());
        } else if (key == "--n") {
            args.n = parse_int(value, key.c_str());
        } else if (key == "--degree") {
            args.degree = parse_int(value, key.c_str());
        } else if (key == "--seed") {
            args.seed = static_cast<unsigned int>(parse_int(value, key.c_str()));
        } else if (key == "--alpha") {
            args.alpha = parse_double(value, key.c_str());
        } else if (key == "--sweeps") {
            args.sweeps = parse_int(value, key.c_str());
        } else if (key == "--use-gst") {
            args.use_gst = parse_int(value, key.c_str());
        } else if (key == "--random-start") {
            args.random_start = parse_int(value, key.c_str());
        } else if (key == "--min-gain") {
            args.min_gain = parse_double(value, key.c_str());
        } else if (key == "--score-only") {
            args.score_only = parse_int(value, key.c_str());
        } else {
            std::fprintf(stderr, "unknown option: %s\n", key.c_str());
            usage(argv[0]);
            std::exit(2);
        }
    }
    if (args.p < 2 || args.n < 2 || args.degree < 0 || args.sweeps < 0 ||
        args.alpha <= 0.0 || args.min_gain < 0.0) {
        usage(argv[0]);
        std::exit(2);
    }
    return args;
}

std::vector<double> simulate_data(const Args& args, int* true_edges) {
    std::mt19937 rng(args.seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::uniform_real_distribution<double> magnitude(0.3, 1.0);

    std::vector<double> weights(static_cast<std::size_t>(args.p) * args.p, 0.0);
    *true_edges = 0;
    for (int child = 1; child < args.p; ++child) {
        const double edge_prob = std::min(1.0, static_cast<double>(args.degree) / child);
        std::vector<int> parents;
        for (int parent = 0; parent < child; ++parent) {
            if (uniform(rng) < edge_prob) {
                parents.push_back(parent);
            }
        }
        const double scale = parents.empty() ? 1.0 : 1.0 / std::sqrt(static_cast<double>(parents.size()));
        for (int parent : parents) {
            const double sign = uniform(rng) < 0.5 ? -1.0 : 1.0;
            weights[static_cast<std::size_t>(parent) * args.p + child] = sign * magnitude(rng) * scale;
            ++(*true_edges);
        }
    }

    std::vector<double> x(static_cast<std::size_t>(args.n) * args.p, 0.0);
    for (int j = 0; j < args.p; ++j) {
        for (int i = 0; i < args.n; ++i) {
            x[static_cast<std::size_t>(i) + static_cast<std::size_t>(args.n) * j] = normal(rng);
        }
        for (int parent = 0; parent < j; ++parent) {
            const double beta = weights[static_cast<std::size_t>(parent) * args.p + j];
            if (beta == 0.0) {
                continue;
            }
            for (int i = 0; i < args.n; ++i) {
                x[static_cast<std::size_t>(i) + static_cast<std::size_t>(args.n) * j] +=
                        beta * x[static_cast<std::size_t>(i) + static_cast<std::size_t>(args.n) * parent];
            }
        }
    }
    return x;
}

void print_stat(const char* key, double value) {
    std::printf("%s,%.17g\n", key, value);
}

} // namespace

int main(int argc, char** argv) {
    const Args args = parse_args(argc, argv);
    int true_edges = 0;
    std::vector<double> x = simulate_data(args, &true_edges);

    fastboss_options options = fastboss_default_options();
    options.alpha = args.alpha;
    options.max_sweeps = args.sweeps;
    options.use_gst = args.use_gst;
    options.random_start = args.random_start;
    options.seed = args.seed;
    options.min_gain = args.min_gain;

    fastboss_result* result = nullptr;
    const fastboss_status status = args.score_only
            ? fastboss_score_order_column_major(x.data(), args.n, args.p, &options, nullptr, &result)
            : fastboss_fit_column_major(x.data(), args.n, args.p, &options, nullptr, &result);
    if (status != FASTBOSS_OK) {
        std::fprintf(stderr, "fastboss failed: %s\n", fastboss_last_error());
        return 1;
    }

    const fastboss_statistics* stats = fastboss_result_statistics(result);
    const int* adjacency = fastboss_result_adjacency(result);
    int learned_edges = 0;
    for (int i = 0; i < args.p * args.p; ++i) {
        learned_edges += adjacency[i] != 0;
    }

    print_stat("p", args.p);
    print_stat("n", args.n);
    print_stat("degree_target", args.degree);
    print_stat("seed", args.seed);
    print_stat("alpha", args.alpha);
    print_stat("sweeps_requested", args.sweeps);
    print_stat("use_gst", args.use_gst);
    print_stat("true_edges", true_edges);
    print_stat("learned_edges", learned_edges);
    print_stat("score", fastboss_result_score(result));
    print_stat("runtime_s", stats->runtime_s);
    print_stat("sweeps", stats->sweeps);
    print_stat("accepted_moves", stats->accepted_moves);
    print_stat("local_score_calls", stats->local_score_calls);
    print_stat("gst_trace_calls", stats->gst_trace_calls);
    print_stat("gst_cached_nodes_allocated", stats->gst_cached_nodes_allocated);
    print_stat("gst_cached_distinct_node_visits", stats->gst_cached_distinct_node_visits);
    print_stat("gst_cached_child_materializations", stats->gst_cached_child_materializations);
    print_stat("gst_cached_branch_entries", stats->gst_cached_branch_entries);
    print_stat("gst_cached_branch_considerations", stats->gst_cached_branch_considerations);
    print_stat("gst_cached_branch_prefix_hits", stats->gst_cached_branch_prefix_hits);
    print_stat("gst_cached_grow_candidate_checks", stats->gst_cached_grow_candidate_checks);
    print_stat("gst_nocache_trace_calls", stats->gst_nocache_trace_calls);

    fastboss_free_result(result);
    return 0;
}
