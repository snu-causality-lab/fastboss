#ifndef FASTBOSS_H
#define FASTBOSS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) && defined(FASTBOSS_STATIC)
#  define FASTBOSS_API
#elif defined(_WIN32) && defined(FASTBOSS_BUILD_SHARED)
#  define FASTBOSS_API __declspec(dllexport)
#elif defined(_WIN32)
#  define FASTBOSS_API __declspec(dllimport)
#else
#  define FASTBOSS_API __attribute__((visibility("default")))
#endif

typedef enum {
    FASTBOSS_OK = 0,
    FASTBOSS_ERROR_INVALID_ARGUMENT = 1,
    FASTBOSS_ERROR_ALLOCATION = 2,
    FASTBOSS_ERROR_INTERNAL = 3
} fastboss_status;

#define FASTBOSS_ABI_VERSION 1

typedef struct {
    size_t struct_size;      /* Must be sizeof(fastboss_options). */
    int abi_version;         /* Must be FASTBOSS_ABI_VERSION. */
    double alpha;
    int max_sweeps;          /* 0 means run to convergence. */
    int use_gst;             /* nonzero enables lazy GST. */
    int random_start;        /* ignored when initial_order is non-null. */
    unsigned int seed;
    double min_gain;
    int reserved[8];
} fastboss_options;

typedef struct {
    size_t struct_size;
    int abi_version;
    double runtime_s;
    int sweeps;
    int accepted_moves;
    double score;
    double local_score_calls;
    double gst_trace_calls;
    double gst_cached_nodes_allocated;
    double gst_cached_distinct_node_visits;
    double gst_cached_child_materializations;
    double gst_cached_branch_entries;
    double gst_cached_branch_considerations;
    double gst_cached_branch_prefix_hits;
    double gst_cached_grow_candidate_checks;
    double gst_nocache_trace_calls;
    double reserved[8];
} fastboss_statistics;

typedef struct fastboss_result fastboss_result;

FASTBOSS_API fastboss_options fastboss_default_options(void);

/*
 * Run BOSS on a column-major n by p matrix, matching R/Fortran layout.
 * `initial_order` may be null. If supplied, it must contain a zero-based
 * permutation of 0..p-1.
 *
 * On success, `*out` owns the returned result. The caller must release it with
 * `fastboss_free_result()`. Result accessors return pointers owned by that
 * result and valid only until it is freed.
 */
FASTBOSS_API fastboss_status fastboss_fit_column_major(
        const double* data,
        int n,
        int p,
        const fastboss_options* options,
        const int* initial_order,
        fastboss_result** out);

/*
 * Score a supplied zero-based order without relocation.
 */
FASTBOSS_API fastboss_status fastboss_score_order_column_major(
        const double* data,
        int n,
        int p,
        const fastboss_options* options,
        const int* order,
        fastboss_result** out);

FASTBOSS_API void fastboss_free_result(fastboss_result* result);

/*
 * Returns a thread-local diagnostic string for the most recent error on the
 * calling thread. The pointer is invalidated by the next fastboss call on that
 * same thread.
 */
FASTBOSS_API const char* fastboss_last_error(void);

FASTBOSS_API double fastboss_result_score(const fastboss_result* result);
FASTBOSS_API int fastboss_result_num_variables(const fastboss_result* result);
FASTBOSS_API const int* fastboss_result_order(const fastboss_result* result);
FASTBOSS_API const int* fastboss_result_initial_order(const fastboss_result* result);
FASTBOSS_API const int* fastboss_result_adjacency(const fastboss_result* result);
FASTBOSS_API const fastboss_statistics* fastboss_result_statistics(const fastboss_result* result);

#ifdef __cplusplus
}
#endif

#endif
