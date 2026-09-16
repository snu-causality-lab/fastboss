# Tests

`core_tests.c` exercises the public C ABI, compares cached and uncached GST runs, verifies ownership and invalid-input behavior, and checks deterministic parent selection under score ties.

`determinism_tests.cpp` fixes the internal tie-breaking contract for GST choices and relocation slots.

`test_python.py` checks the wrapper's input validation, compares cached and uncached runs across seeded random matrices, and independently computes small-order Gaussian scores and parent sets with NumPy least squares. The independent score check does not use the C++ scorer.

Run all tests with:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
PYTHONPATH=wrappers/python python3 -m unittest discover -s tests -p test_python.py -v
```
