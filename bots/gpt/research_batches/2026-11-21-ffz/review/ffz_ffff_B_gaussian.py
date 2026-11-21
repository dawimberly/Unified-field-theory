"""FFff(B) parity-even Gaussian trace and Weyl tilt check.

Extracted from FFFF_B_GAUSSIAN.md. Verifies the even four-vector contraction
and the first-order covariance tilt for a coordinate insertion. This is not a
graph integral or a pole evaluation.

Run:
  python ffz_ffff_B_gaussian.py

Requires: numpy, sympy
"""

import numpy as np
import sympy as s


def gaussian_check(n_samples: int = 400_000, seed: int = 42):
    rng = np.random.default_rng(seed)
    L = rng.normal(size=(4, 4)) * 0.35 + np.eye(4)
    C = L @ L.T
    a, b, c, d = rng.multivariate_normal(np.zeros(4), C, size=(n_samples, 4)).transpose(
        2, 0, 1
    )
    dot = lambda u, v: np.sum(u * v, axis=1)
    samples = (
        dot(a, b) * dot(c, d) - dot(a, c) * dot(b, d) + dot(a, d) * dot(b, c)
    )
    dimension = 4
    ab, ac, ad, bc, bd, cd = C[0, 1], C[0, 2], C[0, 3], C[1, 2], C[1, 3], C[2, 3]
    expected = dimension**2 * (ab * cd - ac * bd + ad * bc) + 2 * dimension * ac * bd
    assert abs(samples.mean() - expected) < 5 * samples.std() / len(samples) ** 0.5

    la = s.symbols("la")
    J = np.zeros((4, 4))
    J[0, 1] = J[1, 0] = 1
    dc = C @ J @ C

    def trace(c):
        return 16 * (
            c[0, 1] * c[2, 3] - c[0, 2] * c[1, 3] + c[0, 3] * c[1, 2]
        ) + 8 * c[0, 2] * c[1, 3]

    derivative = float(s.diff(trace(C + la * dc), la).subs(la, 0))
    weighted = samples * dot(a, b)
    assert abs(weighted.mean() - (derivative + 4 * C[0, 1] * expected)) < 5 * weighted.std() / len(
        samples
    ) ** 0.5
    return {"mean": float(samples.mean()), "expected": float(expected)}


def main():
    out = gaussian_check()
    print("sample mean:", out["mean"])
    print("expected:", out["expected"])
    print("FFff(B) GAUSSIAN CHECKS PASSED")


if __name__ == "__main__":
    main()
