# FastBOSS repository guide for coding agents

Read `README.md` for the user workflow, `docs/algorithm.md` for the implemented score and search, and `include/fastboss.h` for the public C ABI before changing code. This is a standalone serial Gaussian BOSS implementation. Do not infer features from older async or experimental BOSS repositories.

## Map

- `src/fastboss_core.cpp`: Gaussian scorer, grow-shrink traces, Lazy GST, order search, and C ABI implementation.
- `src/tie_break.h`: deterministic tie rules.
- `wrappers/python/pyfastboss/core.py`: NumPy validation, native-library loading, and result copying.
- `examples/`: runnable C and Python usage; `tests/`: native and Python checks.

## Change rules

- Preserve the versioned C ABI and result ownership contract. Update both the C header and Python `ctypes` layout if any public struct changes.
- For a C/C++ behavior change, add a failing regression first. Run its focused test and the full native CTest suite. For score or graph behavior, also compare a small case against the independent NumPy least-squares oracle in `tests/test_python.py`; cached-versus-uncached agreement alone is not an independent correctness check.
- Verify a release build with warnings treated as errors and an ASan/UBSan build when the host toolchain supports them. Keep generated builds, datasets, and large experiment outputs out of Git.
- For a Python wrapper change, run `PYTHONPATH=wrappers/python python -m unittest discover -s tests -p test_python.py -v` after building the shared library. `FASTBOSS_LIBRARY` can point to a build outside the checkout.
- Keep `README.md`, `docs/algorithm.md`, examples, and tests aligned with actual behavior. Do not claim speedups or graph recovery without a reproducible comparison.
- `AGENTS.md` and `CLAUDE.md` are twins: make the same instruction change in both files.
