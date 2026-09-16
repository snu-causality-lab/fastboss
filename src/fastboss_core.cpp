#include "fastboss.h"
#include "tie_break.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

struct fastboss_result {
    int p = 0;
    double score = 0.0;
    std::vector<int> order;
    std::vector<int> initial_order;
    std::vector<int> adjacency; /* row-major: parent * p + child */
    fastboss_statistics stats{};
};

namespace {

using Clock = std::chrono::high_resolution_clock;
constexpr double kPi = 3.141592653589793238462643383279502884;
thread_local std::string g_last_error;

void set_error(const char* message) {
    g_last_error = message ? message : "unknown fastboss error";
}

double elapsed_s(const Clock::time_point& start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

struct FlatSet {
    std::vector<int> values;

    bool contains(int x) const {
        return std::binary_search(values.begin(), values.end(), x);
    }

    void insert(int x) {
        auto it = std::lower_bound(values.begin(), values.end(), x);
        if (it == values.end() || *it != x) {
            values.insert(it, x);
        }
    }

    void erase(int x) {
        auto it = std::lower_bound(values.begin(), values.end(), x);
        if (it != values.end() && *it == x) {
            values.erase(it);
        }
    }

    void clear() {
        values.clear();
    }
};

bool validate_permutation(const int* order, int p) {
    if (!order) {
        return true;
    }
    std::vector<int> seen(p, 0);
    for (int i = 0; i < p; ++i) {
        if (order[i] < 0 || order[i] >= p || seen[order[i]]) {
            return false;
        }
        seen[order[i]] = 1;
    }
    return true;
}

bool validate_dimensions(int n, int p) {
    if (n < 2 || p < 2) {
        return false;
    }
    const auto max_size = std::numeric_limits<std::size_t>::max();
    const auto n_size = static_cast<std::size_t>(n);
    const auto p_size = static_cast<std::size_t>(p);
    const auto int_max = static_cast<std::size_t>(std::numeric_limits<int>::max());
    return p_size <= max_size / p_size &&
           n_size <= max_size / p_size &&
           p_size <= int_max / p_size &&
           n_size <= int_max / p_size;
}

bool load_options(const fastboss_options* options, fastboss_options& out) {
    out = fastboss_default_options();
    if (!options) {
        return true;
    }

    const std::size_t header_size = offsetof(fastboss_options, alpha);
    const std::size_t provided_size = options->struct_size;
    if (provided_size < header_size ||
        provided_size < sizeof(fastboss_options) ||
        options->abi_version != FASTBOSS_ABI_VERSION) {
        set_error("fastboss_options has incompatible struct_size or abi_version");
        return false;
    }

    std::memcpy(&out, options, sizeof(fastboss_options));
    return true;
}

bool validate_data(const double* data, int n, int p) {
    const std::size_t count = static_cast<std::size_t>(n) * static_cast<std::size_t>(p);
    for (std::size_t i = 0; i < count; ++i) {
        if (!std::isfinite(data[i])) {
            return false;
        }
    }
    return true;
}

void erase_value(std::vector<int>& values, int x) {
    auto it = std::find(values.begin(), values.end(), x);
    if (it != values.end()) {
        values.erase(it);
    }
}

void insert_at_slot(std::vector<int>& order, int v, int slot) {
    erase_value(order, v);
    slot = std::max(0, std::min(slot, static_cast<int>(order.size())));
    order.insert(order.begin() + slot, v);
}

int find_pos(const std::vector<int>& order, int v) {
    auto it = std::find(order.begin(), order.end(), v);
    return it == order.end() ? -1 : static_cast<int>(std::distance(order.begin(), it));
}

int score_index_to_slot(int score_index, int old_pos) {
    return score_index - static_cast<int>(score_index > old_pos);
}

std::vector<int> identity_order(int p) {
    std::vector<int> out(p);
    std::iota(out.begin(), out.end(), 0);
    return out;
}

bool cholesky_solve(std::vector<double> a, std::vector<double>& b, int n) {
    for (int j = 0; j < n; ++j) {
        double diag = a[j * n + j];
        for (int k = 0; k < j; ++k) {
            const double l_jk = a[j * n + k];
            diag -= l_jk * l_jk;
        }
        if (diag <= 1e-12 || !std::isfinite(diag)) {
            return false;
        }
        a[j * n + j] = std::sqrt(diag);
        for (int i = j + 1; i < n; ++i) {
            double value = a[i * n + j];
            for (int k = 0; k < j; ++k) {
                value -= a[i * n + k] * a[j * n + k];
            }
            a[i * n + j] = value / a[j * n + j];
        }
    }

    for (int i = 0; i < n; ++i) {
        double value = b[i];
        for (int k = 0; k < i; ++k) {
            value -= a[i * n + k] * b[k];
        }
        b[i] = value / a[i * n + i];
    }
    for (int i = n - 1; i >= 0; --i) {
        double value = b[i];
        for (int k = i + 1; k < n; ++k) {
            value -= a[k * n + i] * b[k];
        }
        b[i] = value / a[i * n + i];
    }
    return true;
}

class GaussianBICScorer {
public:
    GaussianBICScorer(const double* x, int n, int p, double alpha)
        : n(n), p(p), alpha(alpha), cross(static_cast<std::size_t>(p) * p, 0.0) {
        std::vector<double> means(p, 0.0);
        for (int j = 0; j < p; ++j) {
            double sum = 0.0;
            for (int i = 0; i < n; ++i) {
                sum += x[static_cast<std::size_t>(i) +
                         static_cast<std::size_t>(n) * static_cast<std::size_t>(j)];
            }
            means[j] = sum / static_cast<double>(n);
        }
        for (int a = 0; a < p; ++a) {
            for (int b = 0; b <= a; ++b) {
                double value = 0.0;
                for (int i = 0; i < n; ++i) {
                    const auto ia = static_cast<std::size_t>(i) +
                            static_cast<std::size_t>(n) * static_cast<std::size_t>(a);
                    const auto ib = static_cast<std::size_t>(i) +
                            static_cast<std::size_t>(n) * static_cast<std::size_t>(b);
                    value += (x[ia] - means[a]) * (x[ib] - means[b]);
                }
                cross[static_cast<std::size_t>(a) * p + b] = value;
                cross[static_cast<std::size_t>(b) * p + a] = value;
            }
        }
    }

    double local_score(int target, const FlatSet& parents) {
        local_score_calls += 1.0;
        const int m = static_cast<int>(parents.values.size());
        const double syy = std::max(cross[static_cast<std::size_t>(target) * p + target], 1e-12);
        double sse = syy;
        if (m > 0) {
            std::vector<double> a(static_cast<std::size_t>(m) * m, 0.0);
            std::vector<double> b(m, 0.0);
            for (int i = 0; i < m; ++i) {
                const int pi = parents.values[i];
                b[i] = cross[static_cast<std::size_t>(pi) * p + target];
                for (int j = 0; j <= i; ++j) {
                    const int pj = parents.values[j];
                    a[static_cast<std::size_t>(i) * m + j] =
                            cross[static_cast<std::size_t>(pi) * p + pj];
                    a[static_cast<std::size_t>(j) * m + i] =
                            a[static_cast<std::size_t>(i) * m + j];
                }
                a[static_cast<std::size_t>(i) * m + i] += 1e-10;
            }
            std::vector<double> beta = b;
            if (cholesky_solve(a, beta, m)) {
                double explained = 0.0;
                for (int i = 0; i < m; ++i) {
                    explained += b[i] * beta[i];
                }
                sse = syy - explained;
            }
        }
        sse = std::max(sse, 1e-12);
        const double variance = std::max(sse / static_cast<double>(n), 1e-12);
        const double loglik = -0.5 * static_cast<double>(n) *
                (std::log(2.0 * kPi) + 1.0 + std::log(variance));
        const double parameters = static_cast<double>(m + 1);
        return loglik - 0.5 * alpha * parameters * std::log(static_cast<double>(n));
    }

    int n;
    int p;
    double alpha;
    std::vector<double> cross;
    double local_score_calls = 0.0;
};

class BossGST;

class BossGSTNode {
public:
    BossGSTNode(BossGST* tree, double score);
    double trace(const std::vector<char>& prefix_mask, std::vector<int>& available, FlatSet& parents);

private:
    struct Branch {
        int add;
        double grow_score;
        std::unique_ptr<BossGSTNode> child;
        Branch(int add, double grow_score) : add(add), grow_score(grow_score) {}
    };

    BossGST* tree;
    double grow_score;
    double shrink_score;
    bool ever_visited = false;
    bool branches_built = false;
    bool shrink_built = false;
    std::vector<Branch> branches;
    std::vector<int> remove_sequence;

    void grow(const std::vector<int>& available, FlatSet& parents);
    void shrink(FlatSet& parents);
};

class BossGST {
public:
    BossGST(int vertex, GaussianBICScorer* scorer, int p, bool use_cache)
        : vertex(vertex), scorer(scorer), use_cache(use_cache) {
        available_template.reserve(p - 1);
        for (int i = 0; i < p; ++i) {
            if (i != vertex) {
                available_template.push_back(i);
            }
        }
        FlatSet empty;
        root_score = scorer->local_score(vertex, empty);
        if (use_cache) {
            stats.gst_cached_nodes_allocated += 1.0;
            root = std::make_unique<BossGSTNode>(this, root_score);
        }
    }

    double trace(const std::vector<char>& prefix_mask, FlatSet& parents) {
        stats.gst_trace_calls += 1.0;
        parents.clear();
        available.assign(available_template.begin(), available_template.end());
        if (use_cache) {
            return root->trace(prefix_mask, available, parents);
        }
        return trace_without_cache(root_score, prefix_mask, available, parents);
    }

    double greedy_shrink(double current_score, FlatSet& parents) {
        while (true) {
            int best_remove = -1;
            double best_score = current_score;
            std::vector<int> parent_values = parents.values;
            for (int candidate : parent_values) {
                parents.erase(candidate);
                const double score = scorer->local_score(vertex, parents);
                parents.insert(candidate);
                if (score > current_score &&
                    (best_remove < 0 || fastboss::detail::better_indexed_score(
                            score, candidate, best_score, best_remove))) {
                    best_score = score;
                    best_remove = candidate;
                }
            }
            if (best_remove < 0) {
                break;
            }
            current_score = best_score;
            parents.erase(best_remove);
        }
        return current_score;
    }

    double trace_without_cache(double current_score,
                               const std::vector<char>& prefix_mask,
                               std::vector<int>& available_values,
                               FlatSet& parents) {
        stats.gst_nocache_trace_calls += 1.0;
        std::vector<std::pair<int, double>> branches;
        branches.reserve(available_values.size());
        for (int candidate : available_values) {
            if (parents.contains(candidate)) {
                continue;
            }
            const double before = current_score;
            parents.insert(candidate);
            const double score = scorer->local_score(vertex, parents);
            parents.erase(candidate);
            if (score > before) {
                branches.emplace_back(candidate, score);
            }
        }
        std::sort(branches.begin(), branches.end(),
                  [](const auto& lhs, const auto& rhs) {
                      return fastboss::detail::better_indexed_score(
                              lhs.second, lhs.first, rhs.second, rhs.first);
                  });
        for (const auto& branch : branches) {
            erase_value(available_values, branch.first);
            if (prefix_mask[branch.first]) {
                parents.insert(branch.first);
                return trace_without_cache(branch.second, prefix_mask, available_values, parents);
            }
        }
        return greedy_shrink(current_score, parents);
    }

    int vertex;
    GaussianBICScorer* scorer;
    bool use_cache;
    double root_score;
    std::vector<int> available_template;
    std::vector<int> available;
    std::unique_ptr<BossGSTNode> root;
    fastboss_statistics stats{};
};

BossGSTNode::BossGSTNode(BossGST* tree, double score)
    : tree(tree), grow_score(score), shrink_score(score) {}

void BossGSTNode::grow(const std::vector<int>& available, FlatSet& parents) {
    for (int candidate : available) {
        if (parents.contains(candidate)) {
            continue;
        }
        tree->stats.gst_cached_grow_candidate_checks += 1.0;
        parents.insert(candidate);
        const double score = tree->scorer->local_score(tree->vertex, parents);
        parents.erase(candidate);
        if (score > grow_score) {
            tree->stats.gst_cached_branch_entries += 1.0;
            branches.emplace_back(candidate, score);
        }
    }
    std::sort(branches.begin(), branches.end(),
              [](const Branch& lhs, const Branch& rhs) {
                  return fastboss::detail::better_indexed_score(
                          lhs.grow_score, lhs.add, rhs.grow_score, rhs.add);
              });
    branches_built = true;
}

void BossGSTNode::shrink(FlatSet& parents) {
    remove_sequence.clear();
    shrink_score = grow_score;
    while (true) {
        int best_remove = -1;
        double best_score = shrink_score;
        std::vector<int> parent_values = parents.values;
        for (int candidate : parent_values) {
            parents.erase(candidate);
            const double score = tree->scorer->local_score(tree->vertex, parents);
            parents.insert(candidate);
            if (score > shrink_score &&
                (best_remove < 0 || fastboss::detail::better_indexed_score(
                        score, candidate, best_score, best_remove))) {
                best_score = score;
                best_remove = candidate;
            }
        }
        if (best_remove < 0) {
            break;
        }
        shrink_score = best_score;
        remove_sequence.push_back(best_remove);
        parents.erase(best_remove);
    }
    shrink_built = true;
}

double BossGSTNode::trace(const std::vector<char>& prefix_mask,
                          std::vector<int>& available,
                          FlatSet& parents) {
    if (!ever_visited) {
        ever_visited = true;
        tree->stats.gst_cached_distinct_node_visits += 1.0;
    }
    if (!branches_built) {
        grow(available, parents);
    }
    for (auto& branch : branches) {
        tree->stats.gst_cached_branch_considerations += 1.0;
        erase_value(available, branch.add);
        if (prefix_mask[branch.add]) {
            tree->stats.gst_cached_branch_prefix_hits += 1.0;
            if (!branch.child) {
                tree->stats.gst_cached_nodes_allocated += 1.0;
                tree->stats.gst_cached_child_materializations += 1.0;
                branch.child = std::make_unique<BossGSTNode>(tree, branch.grow_score);
            }
            parents.insert(branch.add);
            return branch.child->trace(prefix_mask, available, parents);
        }
    }
    if (!shrink_built) {
        shrink(parents);
    } else {
        for (int remove : remove_sequence) {
            parents.erase(remove);
        }
    }
    return shrink_score;
}

struct Fit {
    std::vector<int> order;
    std::vector<int> initial_order;
    std::vector<FlatSet> parents;
    fastboss_statistics stats{};
};

class BossSearch {
public:
    BossSearch(GaussianBICScorer* scorer, int p, bool use_gst)
        : scorer(scorer), p(p) {
        for (int v = 0; v < p; ++v) {
            gsts.push_back(std::make_unique<BossGST>(v, scorer, p, use_gst));
        }
    }

    double trace_order(const std::vector<int>& order, std::vector<FlatSet>& parents) {
        parents.assign(p, FlatSet{});
        std::vector<char> prefix_mask(p, 0);
        FlatSet pa;
        double total = 0.0;
        for (int node : order) {
            const double score = gsts[node]->trace(prefix_mask, pa);
            parents[node] = pa;
            total += score;
            prefix_mask[node] = 1;
        }
        return total;
    }

    std::pair<double, int> propose_move(const std::vector<int>& order, int v, double min_gain) {
        const int old_pos = find_pos(order, v);
        if (old_pos < 0) {
            return {0.0, old_pos};
        }
        std::vector<double> scores(p + 1, -std::numeric_limits<double>::infinity());
        std::vector<char> prefix_mask(p, 0);
        FlatSet pa;
        double score_prefix = 0.0;
        for (int j = 0; j < p; ++j) {
            const int w = order[j];
            scores[j] = gsts[v]->trace(prefix_mask, pa) + score_prefix;
            if (v != w) {
                score_prefix += gsts[w]->trace(prefix_mask, pa);
                prefix_mask[w] = 1;
            }
        }
        scores[p] = gsts[v]->trace(prefix_mask, pa) + score_prefix;

        int best = p;
        prefix_mask[v] = 1;
        double suffix_score = 0.0;
        for (int j = p - 1; j >= 0; --j) {
            const int w = order[j];
            if (v != w) {
                prefix_mask[w] = 0;
                suffix_score += gsts[w]->trace(prefix_mask, pa);
            }
            scores[j] += suffix_score;
            const int candidate_slot = score_index_to_slot(j, old_pos);
            const int best_slot = score_index_to_slot(best, old_pos);
            if (fastboss::detail::better_relocation(
                    scores[j], candidate_slot,
                    scores[best], best_slot, old_pos)) {
                best = j;
            }
        }
        const double gain = scores[best] - scores[old_pos];
        const int slot = gain > min_gain ? score_index_to_slot(best, old_pos) : old_pos;
        return {gain, slot};
    }

    Fit run(std::vector<int> order, int max_sweeps, int seed, double min_gain) {
        const auto start = Clock::now();
        Fit fit;
        fit.initial_order = order;
        std::vector<int> variables = identity_order(p);
        std::mt19937 rng(static_cast<std::uint32_t>(seed));
        while (true) {
            bool improved = false;
            std::shuffle(variables.begin(), variables.end(), rng);
            for (int v : variables) {
                const auto proposal = propose_move(order, v, min_gain);
                if (proposal.first > min_gain) {
                    insert_at_slot(order, v, proposal.second);
                    improved = true;
                    fit.stats.accepted_moves += 1;
                }
            }
            fit.stats.sweeps += 1;
            if (!improved || (max_sweeps > 0 && fit.stats.sweeps >= max_sweeps)) {
                break;
            }
        }
        fit.order = order;
        fit.stats.score = trace_order(order, fit.parents);
        fit.stats.runtime_s = elapsed_s(start);
        collect_stats(fit);
        return fit;
    }

    Fit score(std::vector<int> order) {
        const auto start = Clock::now();
        Fit fit;
        fit.initial_order = order;
        fit.order = order;
        fit.stats.score = trace_order(order, fit.parents);
        fit.stats.runtime_s = elapsed_s(start);
        collect_stats(fit);
        return fit;
    }

private:
    void collect_stats(Fit& fit) {
        fit.stats.local_score_calls = scorer->local_score_calls;
        for (const auto& gst : gsts) {
            const auto& s = gst->stats;
            fit.stats.gst_trace_calls += s.gst_trace_calls;
            fit.stats.gst_cached_nodes_allocated += s.gst_cached_nodes_allocated;
            fit.stats.gst_cached_distinct_node_visits += s.gst_cached_distinct_node_visits;
            fit.stats.gst_cached_child_materializations += s.gst_cached_child_materializations;
            fit.stats.gst_cached_branch_entries += s.gst_cached_branch_entries;
            fit.stats.gst_cached_branch_considerations += s.gst_cached_branch_considerations;
            fit.stats.gst_cached_branch_prefix_hits += s.gst_cached_branch_prefix_hits;
            fit.stats.gst_cached_grow_candidate_checks += s.gst_cached_grow_candidate_checks;
            fit.stats.gst_nocache_trace_calls += s.gst_nocache_trace_calls;
        }
    }

    GaussianBICScorer* scorer;
    int p;
    std::vector<std::unique_ptr<BossGST>> gsts;
};

std::vector<int> make_order(int p, const int* initial_order, bool random_start, unsigned seed) {
    std::vector<int> order = identity_order(p);
    if (initial_order) {
        std::copy(initial_order, initial_order + p, order.begin());
    } else if (random_start) {
        std::mt19937 rng(seed);
        std::shuffle(order.begin(), order.end(), rng);
    }
    return order;
}

fastboss_result* result_from_fit(Fit&& fit, int p) {
    auto result = std::make_unique<fastboss_result>();
    result->p = p;
    result->score = fit.stats.score;
    result->order = std::move(fit.order);
    result->initial_order = std::move(fit.initial_order);
    result->adjacency.assign(static_cast<std::size_t>(p) * p, 0);
    for (int child = 0; child < p; ++child) {
        for (int parent : fit.parents[child].values) {
            result->adjacency[static_cast<std::size_t>(parent) * p + child] = 1;
        }
    }
    result->stats = fit.stats;
    result->stats.struct_size = sizeof(fastboss_statistics);
    result->stats.abi_version = FASTBOSS_ABI_VERSION;
    return result.release();
}

fastboss_status run_impl(const double* data,
                         int n,
                         int p,
                         const fastboss_options* options,
                         const int* initial_order,
                         fastboss_result** out,
                         bool score_only) {
    if (!out) {
        set_error("out must not be null");
        return FASTBOSS_ERROR_INVALID_ARGUMENT;
    }
    *out = nullptr;
    if (!data || !validate_dimensions(n, p)) {
        set_error("data must be a non-null column-major matrix with n >= 2 and p >= 2");
        return FASTBOSS_ERROR_INVALID_ARGUMENT;
    }
    if (!validate_data(data, n, p)) {
        set_error("data must contain only finite values");
        return FASTBOSS_ERROR_INVALID_ARGUMENT;
    }
    if (!validate_permutation(initial_order, p)) {
        set_error("initial_order must be a zero-based permutation");
        return FASTBOSS_ERROR_INVALID_ARGUMENT;
    }
    fastboss_options opt;
    if (!load_options(options, opt)) {
        return FASTBOSS_ERROR_INVALID_ARGUMENT;
    }
    if (!std::isfinite(opt.alpha) || !std::isfinite(opt.min_gain) ||
        opt.alpha <= 0.0 || opt.max_sweeps < 0 || opt.min_gain < 0.0) {
        set_error("invalid fastboss options");
        return FASTBOSS_ERROR_INVALID_ARGUMENT;
    }
    try {
        GaussianBICScorer scorer(data, n, p, opt.alpha);
        BossSearch search(&scorer, p, opt.use_gst != 0);
        std::vector<int> order = make_order(p, initial_order, opt.random_start != 0, opt.seed);
        Fit fit = score_only
                ? search.score(order)
                : search.run(order, opt.max_sweeps, static_cast<int>(opt.seed), opt.min_gain);
        *out = result_from_fit(std::move(fit), p);
        g_last_error.clear();
        return FASTBOSS_OK;
    } catch (const std::bad_alloc&) {
        set_error("allocation failed");
        return FASTBOSS_ERROR_ALLOCATION;
    } catch (const std::exception& e) {
        set_error(e.what());
        return FASTBOSS_ERROR_INTERNAL;
    }
}

} // namespace

extern "C" fastboss_options fastboss_default_options(void) {
    fastboss_options options;
    options.struct_size = sizeof(fastboss_options);
    options.abi_version = FASTBOSS_ABI_VERSION;
    options.alpha = 2.0;
    options.max_sweeps = 0;
    options.use_gst = 1;
    options.random_start = 0;
    options.seed = 20260525U;
    options.min_gain = 1e-6;
    std::fill(std::begin(options.reserved), std::end(options.reserved), 0);
    return options;
}

extern "C" fastboss_status fastboss_fit_column_major(
        const double* data,
        int n,
        int p,
        const fastboss_options* options,
        const int* initial_order,
        fastboss_result** out) {
    return run_impl(data, n, p, options, initial_order, out, false);
}

extern "C" fastboss_status fastboss_score_order_column_major(
        const double* data,
        int n,
        int p,
        const fastboss_options* options,
        const int* order,
        fastboss_result** out) {
    if (out) {
        *out = nullptr;
    }
    if (!order) {
        set_error("order must not be null");
        return FASTBOSS_ERROR_INVALID_ARGUMENT;
    }
    return run_impl(data, n, p, options, order, out, true);
}

extern "C" void fastboss_free_result(fastboss_result* result) {
    delete result;
}

extern "C" const char* fastboss_last_error(void) {
    return g_last_error.c_str();
}

extern "C" double fastboss_result_score(const fastboss_result* result) {
    return result ? result->score : std::numeric_limits<double>::quiet_NaN();
}

extern "C" int fastboss_result_num_variables(const fastboss_result* result) {
    return result ? result->p : 0;
}

extern "C" const int* fastboss_result_order(const fastboss_result* result) {
    return result ? result->order.data() : nullptr;
}

extern "C" const int* fastboss_result_initial_order(const fastboss_result* result) {
    return result ? result->initial_order.data() : nullptr;
}

extern "C" const int* fastboss_result_adjacency(const fastboss_result* result) {
    return result ? result->adjacency.data() : nullptr;
}

extern "C" const fastboss_statistics* fastboss_result_statistics(const fastboss_result* result) {
    return result ? &result->stats : nullptr;
}
