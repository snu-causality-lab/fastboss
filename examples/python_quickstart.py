"""Runnable FastBOSS example; run from the repository root after building."""

import numpy as np

from pyfastboss import fit, score_order


def main():
    rng = np.random.default_rng(7)
    noise = rng.normal(size=(300, 5))
    data = np.empty_like(noise)
    data[:, 0] = noise[:, 0]
    data[:, 1] = 0.8 * data[:, 0] + noise[:, 1]
    data[:, 2] = -0.7 * data[:, 0] + noise[:, 2]
    data[:, 3] = 0.6 * data[:, 1] + 0.5 * data[:, 2] + noise[:, 3]
    data[:, 4] = 0.9 * data[:, 3] + noise[:, 4]

    result = fit(data, seed=7)
    print(f"score: {result.score:.3f}")
    print(f"order: {result.order.tolist()}")
    print("estimated edges:")
    for parent, child in np.argwhere(result.adjacency):
        print(f"  {parent} -> {child}")
    print(f"local score calls: {result.statistics['local_score_calls']:.0f}")

    fixed = score_order(data, result.order)
    print(f"final-order score check: {np.isclose(fixed.score, result.score)}")


if __name__ == "__main__":
    main()
