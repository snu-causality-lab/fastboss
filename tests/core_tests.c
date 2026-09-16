#include "fastboss.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void fill_data(double* x, int n, int p) {
    (void) p;
    for (int i = 0; i < n; ++i) {
        const double t = (double) i / (double) n;
        const double x0 = sin(6.283185307179586 * t);
        const double x1 = cos(6.283185307179586 * t);
        x[i + n * 0] = x0;
        x[i + n * 1] = 0.8 * x0 + 0.1 * sin(19.0 * t);
        x[i + n * 2] = -0.7 * x0 + 0.2 * x1;
        x[i + n * 3] = 0.6 * x[i + n * 1] + 0.5 * x[i + n * 2] + 0.05 * cos(17.0 * t);
        x[i + n * 4] = 0.9 * x[i + n * 3] + 0.1 * sin(23.0 * t);
    }
}

static int check_tied_parent_choice(void) {
    const int n = 128;
    const int p = 66;
    const int target = p - 1;
    double* x = (double*) calloc((size_t) n * p, sizeof(double));
    int* order = (int*) calloc((size_t) p, sizeof(int));
    if (!x || !order) {
        free(x);
        free(order);
        return 0;
    }

    for (int j = 0; j < p - 1; ++j) {
        for (int i = 0; i < n; ++i) {
            const double t = (double) i / (double) n;
            x[i + n * j] = sin(6.283185307179586 * t);
        }
    }
    for (int i = 0; i < n; ++i) {
        const double t = (double) i / (double) n;
        x[i + n * target] =
                0.8 * sin(6.283185307179586 * t) + 0.1 * cos(17.0 * t);
    }
    for (int j = 0; j < p; ++j) {
        order[j] = j;
    }

    fastboss_options gst = fastboss_default_options();
    gst.use_gst = 1;
    fastboss_options nocache = gst;
    nocache.use_gst = 0;
    fastboss_result* cached = NULL;
    fastboss_result* uncached = NULL;
    const fastboss_status cached_status =
            fastboss_score_order_column_major(x, n, p, &gst, order, &cached);
    const fastboss_status uncached_status =
            fastboss_score_order_column_major(x, n, p, &nocache, order, &uncached);
    int ok = cached_status == FASTBOSS_OK && uncached_status == FASTBOSS_OK;
    if (ok) {
        const int* cached_adjacency = fastboss_result_adjacency(cached);
        const int* uncached_adjacency = fastboss_result_adjacency(uncached);
        int cached_parent = -1;
        int uncached_parent = -1;
        int cached_count = 0;
        int uncached_count = 0;
        for (int parent = 0; parent < target; ++parent) {
            if (cached_adjacency[parent * p + target]) {
                cached_parent = parent;
                cached_count += 1;
            }
            if (uncached_adjacency[parent * p + target]) {
                uncached_parent = parent;
                uncached_count += 1;
            }
        }
        ok = cached_count == 1 && uncached_count == 1 &&
             cached_parent == 0 && uncached_parent == 0;
    }

    fastboss_free_result(cached);
    fastboss_free_result(uncached);
    free(order);
    free(x);
    return ok;
}

static int same_int_array(const int* a, const int* b, int n) {
    for (int i = 0; i < n; ++i) {
        if (a[i] != b[i]) {
            return 0;
        }
    }
    return 1;
}

int main(void) {
    if (!check_tied_parent_choice()) {
        fprintf(stderr, "equal-score GST branches did not choose the smallest parent ID\n");
        return 1;
    }

    const int n = 128;
    const int p = 5;
    double* x = (double*) calloc((size_t) n * p, sizeof(double));
    if (!x) {
        return 2;
    }
    fill_data(x, n, p);

    int initial_order[] = {0, 1, 2, 3, 4};
    fastboss_options gst = fastboss_default_options();
    gst.max_sweeps = 3;
    gst.use_gst = 1;
    gst.seed = 99U;
    fastboss_options nocache = gst;
    nocache.use_gst = 0;

    fastboss_result* a = NULL;
    fastboss_result* b = NULL;
    fastboss_status sa = fastboss_fit_column_major(x, n, p, &gst, initial_order, &a);
    fastboss_status sb = fastboss_fit_column_major(x, n, p, &nocache, initial_order, &b);
    if (sa != FASTBOSS_OK || sb != FASTBOSS_OK) {
        fprintf(stderr, "fit failed: %s\n", fastboss_last_error());
        fastboss_free_result(a);
        fastboss_free_result(b);
        free(x);
        return 1;
    }

    const double score_diff = fabs(fastboss_result_score(a) - fastboss_result_score(b));
    if (score_diff > 1e-7) {
        fprintf(stderr, "GST/no-GST score mismatch: %.12g\n", score_diff);
        return 1;
    }
    if (!same_int_array(fastboss_result_order(a), fastboss_result_order(b), p)) {
        fprintf(stderr, "GST/no-GST final order mismatch\n");
        return 1;
    }
    if (!same_int_array(fastboss_result_adjacency(a), fastboss_result_adjacency(b), p * p)) {
        fprintf(stderr, "GST/no-GST adjacency mismatch\n");
        return 1;
    }

    fastboss_result* scored = NULL;
    fastboss_status ss = fastboss_score_order_column_major(x, n, p, &gst, initial_order, &scored);
    if (ss != FASTBOSS_OK || scored == NULL) {
        fprintf(stderr, "valid score_order failed: %s\n", fastboss_last_error());
        return 1;
    }
    if (!same_int_array(fastboss_result_order(scored), initial_order, p) ||
        !isfinite(fastboss_result_score(scored))) {
        fprintf(stderr, "valid score_order returned invalid result\n");
        return 1;
    }
    fastboss_free_result(scored);

    fastboss_result* stale = (fastboss_result*) 0x1;
    fastboss_status bad = fastboss_score_order_column_major(x, n, p, &gst, NULL, &stale);
    if (bad == FASTBOSS_OK || stale != NULL) {
        fprintf(stderr, "invalid order did not clear out pointer\n");
        return 1;
    }

    double x_bad[10];
    for (int i = 0; i < 10; ++i) {
        x_bad[i] = 0.0;
    }
    x_bad[3] = NAN;
    fastboss_result* bad_result = NULL;
    bad = fastboss_fit_column_major(x_bad, 5, 2, &gst, NULL, &bad_result);
    if (bad == FASTBOSS_OK || bad_result != NULL) {
        fprintf(stderr, "non-finite data was accepted\n");
        return 1;
    }

    fastboss_options invalid = gst;
    invalid.alpha = 0.0;
    bad_result = (fastboss_result*) 0x1;
    bad = fastboss_fit_column_major(x, n, p, &invalid, NULL, &bad_result);
    if (bad == FASTBOSS_OK || bad_result != NULL) {
        fprintf(stderr, "invalid alpha was accepted\n");
        return 1;
    }

    invalid = gst;
    invalid.struct_size = 1;
    bad_result = (fastboss_result*) 0x1;
    bad = fastboss_fit_column_major(x, n, p, &invalid, NULL, &bad_result);
    if (bad == FASTBOSS_OK || bad_result != NULL) {
        fprintf(stderr, "invalid options struct_size was accepted\n");
        return 1;
    }

    invalid = gst;
    invalid.abi_version = FASTBOSS_ABI_VERSION + 1;
    bad_result = (fastboss_result*) 0x1;
    bad = fastboss_fit_column_major(x, n, p, &invalid, NULL, &bad_result);
    if (bad == FASTBOSS_OK || bad_result != NULL) {
        fprintf(stderr, "invalid options abi_version was accepted\n");
        return 1;
    }

    bad = fastboss_fit_column_major(x, n, p, &gst, NULL, NULL);
    if (bad == FASTBOSS_OK) {
        fprintf(stderr, "null out pointer was accepted\n");
        return 1;
    }

    fastboss_free_result(a);
    fastboss_free_result(b);
    free(x);
    printf("fastboss core tests passed\n");
    return 0;
}
