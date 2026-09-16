# FastBOSS

FastBOSS searches for a directed acyclic graph from continuous, observational data. It implements a serial Best Order Score Search (BOSS) with a linear Gaussian score. Lazy grow-shrink trees (Lazy GST) reuse parent-search traces across order evaluations. The same search can run without that cache by setting `use_gst=False`.

The repository contains a C++17 library with a versioned C interface and a small Python wrapper. Start with Python if you want to analyze a NumPy matrix; use the C interface if you are integrating the kernel into another program. [Algorithm and numerical details](docs/algorithm.md) explain what the implementation actually computes.

## Run the Python example

From the repository root on macOS or Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install numpy
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
PYTHONPATH=wrappers/python python examples/python_quickstart.py
```

You need CMake 3.16 or newer, a C++17 compiler, Python 3.9 or newer, and NumPy. The native build creates `build/libfastboss.dylib` on macOS or `build/libfastboss.so` on Linux. The Python wrapper finds it there automatically. If you build elsewhere, set `FASTBOSS_LIBRARY` to the full path of the shared library.

The example generates five variables, fits a graph, prints its score and order, and lists edges as `parent -> child`. Its output is an estimate from one simulated dataset, not a guarantee of recovering every generating edge.

## Use your own data

Pass a two-dimensional array with **rows = observations** and **columns = variables**. Data must be finite and real-valued; the kernel requires at least two rows and two columns. There is no missing-value imputation or categorical-data model. If column scales differ sharply, consider standardizing nonconstant columns because the numerical ridge and variance floors are fixed absolute values.

```python
import numpy as np
from pyfastboss import fit, score_order

X = np.loadtxt("observations.csv", delimiter=",")  # shape: (n, p), no header
result = fit(X, alpha=2.0, seed=7)
print(result.order)                         # zero-based variable IDs
print(result.adjacency)                     # adjacency[parent, child] == 1
print(result.score)
print(result.statistics["local_score_calls"])

fixed = score_order(X, result.order)        # score this order; do not relocate
assert np.isclose(fixed.score, result.score)
```

Run it with `PYTHONPATH=wrappers/python python your_script.py` after building the shared library. If your CSV has headers, missing values, or nonnumeric columns, prepare a numeric matrix first. [Python wrapper notes](wrappers/python/README.md) describe library discovery and result ownership.

| `fit` argument | Default | Meaning |
| --- | --- | --- |
| `alpha` | `2.0` | Positive multiplier on the Gaussian BIC penalty. |
| `sweeps` | `0` | Maximum relocation sweeps; `0` runs until a sweep accepts no move. |
| `order` | `None` | Optional zero-based initial permutation of `0, ..., p-1`. |
| `random_start` | `False` | Shuffle the initial order; ignored if `order` is supplied. |
| `seed` | `20260525` | Seed for initial shuffling and variable visit order; integer in `[0, 2**32 - 1]`. |
| `use_gst` | `True` | Cache lazy grow-shrink traces. `False` recomputes them. |
| `min_gain` | `1e-6` | A relocation must improve the total score by more than this amount. |

`score_order(X, order, alpha=2.0, use_gst=True)` evaluates one supplied order without relocation. `FastBossResult` contains `score`, `order`, `initial_order`, a `p × p` integer `adjacency` matrix, and a `statistics` dictionary. Scores are comparable only when the data and score settings are the same.

Start with the defaults: `sweeps=0` searches until no move is accepted, and larger `alpha` values penalize additional parents more strongly. Set `sweeps` to a positive cap for a quicker exploratory run; change `min_gain` only when you intend to change the move-acceptance threshold. Higher scores are preferred for the same matrix and `alpha`; score magnitude alone does not measure graph quality. `local_score_calls` counts scorer evaluations for profiling, not statistical confidence. Observational edges are candidate relationships, not validated causal effects.

## C interface and development

The public interface is [`include/fastboss.h`](include/fastboss.h). Input is a column-major `double` matrix: observation `i` of variable `j` is `data[i + n*j]`. Adjacency output is row-major: `adjacency[parent*p + child]`. Initialize options with `fastboss_default_options()`, check the returned status and `fastboss_last_error()`, and release every successful result with `fastboss_free_result()`. [`examples/smoke.c`](examples/smoke.c) shows the full ownership pattern. A static library can be built with `-DFASTBOSS_BUILD_SHARED=OFF`.

Run the native tests with the CMake and CTest commands above. To run the optional Python tests after building:

```bash
PYTHONPATH=wrappers/python python -m unittest discover -s tests -p test_python.py -v
```

The [test guide](tests/README.md) describes the score oracle and cache checks. Source lives in [`src/fastboss_core.cpp`](src/fastboss_core.cpp); [`AGENTS.md`](AGENTS.md) gives coding agents a repository map and verification rules, with the same instructions in `CLAUDE.md`.

This repository implements the serial Gaussian search only. It does not include asynchronous or parallel BOSS, an observational equivalence-class conversion, interventional scoring, or the experimental pipelines from which this core was extracted. Numerical conventions and limits are documented in [the algorithm note](docs/algorithm.md).

## License

MIT; see [LICENSE](LICENSE).
