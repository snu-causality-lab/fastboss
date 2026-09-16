# Algorithm implemented here

FastBOSS evaluates an order of `p` variables by assigning each variable a parent set from its predecessors and summing local Gaussian scores. It then relocates variables to improve that total. Both the order search and the parent search are greedy; the result is not guaranteed to be a globally optimal DAG.

## Local score

For target `j` and parent set `P`, the kernel centers each column and fits a linear regression of `j` on `P`. With `n` observations, its local score is

```text
S(j, P) = -n/2 * [log(2*pi) + 1 + log(variance)]
          - alpha/2 * (|P| + 1) * log(n)
variance = max(SSE/n, 1e-12)
```

The `+1` counts residual variance. The kernel adds `1e-10` to the diagonal of the parent cross-product matrix before its Cholesky solve. It also floors the centered target sum of squares and residual sum of squares at `1e-12`; if the solve fails, it uses the target sum of squares. These fixed constants are implementation choices, so scores need not equal those from another Gaussian BIC implementation. `alpha` defaults to `2.0` and must be positive. Continuous, finite, reasonably scaled inputs are expected.

## Parent trace and Lazy GST

For each target, a grow-shrink trace starts with no parents. At each grow step, it scores adding each remaining candidate variable. It sorts strictly improving candidates by score, with smaller variable ID first on ties, and follows the first candidate that precedes the target in the current order. Candidates visited before that choice are removed from the remaining list for that trace. When no eligible improving grow branch remains, it greedily removes parents while the local score strictly improves. Equal-score shrink choices remove the smaller ID first. The final parents are all predecessors, so the resulting adjacency matrix is acyclic.

With `use_gst=True`, the tree for each target stores improving grow branches once. It allocates a child only if a later trace follows that branch; it caches the shrink sequence only at a reached terminal node. `use_gst=False` recomputes the same trace without a tree. For the same data, initial order, and seed, the two modes should return the same order, score, and adjacency; their score-call counts and runtime can differ. The tree is a parent-trace cache, not a cache of Gaussian local scores.

## Order search

The default initial order is `0, 1, ..., p-1`. `random_start=True` shuffles it with `seed`; an explicit `order` overrides that setting. Each sweep shuffles the variable visit order with the same seeded generator. For each visited variable, the search evaluates its possible insertion slots and accepts a relocation only if total score improves by more than `min_gain`. Equal-score relocation choices prefer the shortest move, then the smaller insertion slot. `sweeps=0` continues until a full sweep accepts no move; a positive value caps the number of sweeps. A fixed seed makes repeated runs deterministic for the same build and input.

`score_order` only traces a supplied order; it does not relocate variables. The final `adjacency[parent, child]` value is `1` exactly when the trace selected that parent for that child. The score is the sum of the final local scores.

## Boundaries

The implementation is serial and uses a linear Gaussian observational score. It does not model missing values, categorical variables, interventions, latent confounding, or a Markov equivalence class. Neither the tree cache nor the greedy search promises a particular runtime or graph-recovery rate. The result is a candidate DAG to validate against the assumptions and data of the application.
