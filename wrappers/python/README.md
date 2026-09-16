# Python wrapper

`pyfastboss` is a thin NumPy/`ctypes` interface to the C ABI. It is imported directly from this checkout; there is no separate package installation step. Follow the repository [quick start](../../README.md#run-the-python-example) to build the native shared library and run the example.

From the repository root:

```bash
PYTHONPATH=wrappers/python python examples/python_quickstart.py
```

The wrapper searches ancestor `build/` directories for `libfastboss.dylib` on macOS or `libfastboss.so` on Linux. When building to another directory, set `FASTBOSS_LIBRARY=/absolute/path/to/libfastboss.dylib` (or `.so`). If loading fails, the exception lists the paths searched.

`fit(X, ...)` returns a `FastBossResult` with a score, zero-based order, initial order, adjacency matrix, and statistics. Input rows are observations and columns are variables. The wrapper copies returned arrays and statistics before freeing the native result, so they remain valid after the call. Invalid native calls raise `FastBossError` with the C library's diagnostic. The full option list and data requirements are in the [main README](../../README.md#use-your-own-data).

`score_order(X, order, alpha=2.0, use_gst=True)` evaluates a given order and returns the same result fields; it does not search for a better order.
