"""Checks for the optional Python interface against a built native library."""

import unittest
import math

import numpy as np

from pyfastboss import fit, score_order


def independent_order_score(data, order, alpha=2.0):
    """Small NumPy reference: least-squares BIC plus uncached greedy tracing."""
    n, p = data.shape
    centered = data - data.mean(axis=0)

    def local_score(target, parents):
        y = centered[:, target]
        if parents:
            design = centered[:, sorted(parents)]
            coefficients = np.linalg.lstsq(design, y, rcond=None)[0]
            residual = y - design @ coefficients
        else:
            residual = y
        variance = max(float(residual @ residual) / n, 1e-12)
        loglik = -0.5 * n * (math.log(2 * math.pi) + 1 + math.log(variance))
        return loglik - 0.5 * alpha * (len(parents) + 1) * math.log(n)

    adjacency = np.zeros((p, p), dtype=np.int32)
    total = 0.0
    prefix = set()
    for target in order:
        parents = set()
        available = set(range(p)) - {target}
        current = local_score(target, parents)
        while True:
            improvements = []
            for candidate in available:
                trial = local_score(target, parents | {candidate})
                if trial > current:
                    improvements.append((trial, candidate))
            improvements.sort(key=lambda item: (-item[0], item[1]))
            chosen = next((item for item in improvements if item[1] in prefix), None)
            if chosen is None:
                break
            # All higher-ranked branches were explored before this choice.
            available -= {candidate for _, candidate in improvements[: improvements.index(chosen) + 1]}
            current, parent = chosen
            parents.add(parent)
        while parents:
            improvements = []
            for candidate in parents:
                trial = local_score(target, parents - {candidate})
                if trial > current:
                    improvements.append((trial, candidate))
            if not improvements:
                break
            current, removed = min(improvements, key=lambda item: (-item[0], item[1]))
            parents.remove(removed)
        for parent in parents:
            adjacency[parent, target] = 1
        total += current
        prefix.add(target)
    return total, adjacency


class PythonInterfaceTests(unittest.TestCase):
    def test_order_score_matches_independent_least_squares_oracle(self):
        for seed in range(5):
            with self.subTest(seed=seed):
                data = np.random.default_rng(seed).normal(size=(80, 4))
                order = [2, 0, 3, 1]
                expected_score, expected_adjacency = independent_order_score(data, order)
                actual = score_order(data, order)
                self.assertAlmostEqual(actual.score, expected_score, places=6)
                np.testing.assert_array_equal(actual.adjacency, expected_adjacency)

    def test_cached_and_uncached_runs_agree(self):
        for p in (3, 5, 8):
            for seed in range(5):
                with self.subTest(p=p, seed=seed):
                    data = np.random.default_rng(seed).normal(size=(64, p))
                    cached = fit(data, sweeps=2, seed=seed, use_gst=True)
                    uncached = fit(data, sweeps=2, seed=seed, use_gst=False)
                    self.assertAlmostEqual(cached.score, uncached.score, places=7)
                    np.testing.assert_array_equal(cached.order, uncached.order)
                    np.testing.assert_array_equal(cached.adjacency, uncached.adjacency)
                    scored = score_order(data, cached.order)
                    self.assertAlmostEqual(cached.score, scored.score, places=7)

    def test_complex_data_is_rejected(self):
        data = np.ones((8, 3), dtype=np.complex128)
        with self.assertRaisesRegex(ValueError, "real"):
            fit(data, sweeps=1)

    def test_complex_order_is_rejected(self):
        data = np.random.default_rng(1).normal(size=(8, 3))
        with self.assertRaisesRegex(ValueError, "real"):
            score_order(data, np.array([0, 1 + 2j, 2]))

    def test_seed_must_fit_native_unsigned_integer(self):
        data = np.random.default_rng(1).normal(size=(8, 3))
        for seed in (-1, 2**32):
            with self.subTest(seed=seed):
                with self.assertRaisesRegex(ValueError, "seed"):
                    fit(data, sweeps=1, seed=seed)


if __name__ == "__main__":
    unittest.main()
