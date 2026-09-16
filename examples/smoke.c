#include "fastboss.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double normal01(unsigned int* state) {
    *state = 1664525U * (*state) + 1013904223U;
    double u1 = ((*state >> 8) + 1.0) / 16777217.0;
    *state = 1664525U * (*state) + 1013904223U;
    double u2 = ((*state >> 8) + 1.0) / 16777217.0;
    return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

int main(void) {
    const int n = 200;
    const int p = 5;
    double* x = (double*) calloc((size_t) n * p, sizeof(double));
    if (!x) {
        return 2;
    }
    unsigned int rng = 7U;
    for (int i = 0; i < n; ++i) {
        double e0 = normal01(&rng);
        double e1 = normal01(&rng);
        double e2 = normal01(&rng);
        double e3 = normal01(&rng);
        double e4 = normal01(&rng);
        x[i + n * 0] = e0;
        x[i + n * 1] = 0.8 * e0 + e1;
        x[i + n * 2] = -0.7 * e0 + e2;
        x[i + n * 3] = 0.6 * x[i + n * 1] + 0.5 * x[i + n * 2] + e3;
        x[i + n * 4] = 0.9 * x[i + n * 3] + e4;
    }

    fastboss_options options = fastboss_default_options();
    options.max_sweeps = 2;
    options.random_start = 1;
    options.seed = 42U;

    fastboss_result* result = NULL;
    fastboss_status status = fastboss_fit_column_major(x, n, p, &options, NULL, &result);
    if (status != FASTBOSS_OK) {
        fprintf(stderr, "fastboss failed: %s\n", fastboss_last_error());
        free(x);
        return 1;
    }

    const fastboss_statistics* stats = fastboss_result_statistics(result);
    printf("score %.6f\n", fastboss_result_score(result));
    printf("p %d sweeps %d accepted %d local_score_calls %.0f\n",
           fastboss_result_num_variables(result),
           stats->sweeps,
           stats->accepted_moves,
           stats->local_score_calls);
    printf("order");
    const int* order = fastboss_result_order(result);
    for (int j = 0; j < fastboss_result_num_variables(result); ++j) {
        printf(" %d", order[j]);
    }
    printf("\n");

    fastboss_free_result(result);
    free(x);
    return 0;
}
