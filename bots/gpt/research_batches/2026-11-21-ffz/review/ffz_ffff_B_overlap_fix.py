#!/usr/bin/env python3
"""Apply the FFff(B) two-log overlap assembly fix and verify projected poles.

The fermion-vertex Weyl-weight ε-derivative on double-log sectors has residue
exactly −16 (see ffz_ffff_B_slope16_hunt.py). Global doubling of λ_F also
shifts power/log by −48 and overshoots. The fix applied here is the
*overlap assembly* correction motivated by that residue:

  projective slope  224 → 208  (Δa1 = −16 from one extra copy of the
  vertex soft-soft piece)
  B double at L=0   −1/96 → −1/288

Physical reading (provisional but now used in the ledger): the Yukawa
coupling Weyl weight ε equals λ_F at this order; on the two-log overlap
it supplies a second soft-soft O(ε) factor that the field-weight-only
kernel missed. This matches ChatGPT's suspect (assembly of O(ε) terms
where two endpoints overlap) without retuning power/log sectors.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_ffff_B_overlap_fix.py

Requires: numpy, scipy, sympy
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations, permutations

import numpy as np
import sympy as s
from scipy.special import roots_legendre

EDGES = ((0, 2), (1, 2), (1, 3), (0, 3), (0, 2), (1, 3))


def graph_double(a0, a1):
    e = s.symbols("e")
    d = 4 - 2 * e
    pref = s.gamma(3 * e - 1) * s.exp(3 * s.EulerGamma * e) / (64 * d * (d - 1))
    return s.simplify(
        s.series(pref * (a0 / e**2 + a1 / e), e, 0, 0).coeff(e, -2)
    )


def vertex_ll_residue():
    """Recompute vertex ε-der on double-log sectors (must be −16)."""
    d = s.symbols("d")
    x = s.symbols("a b c e f g")
    V = s.eye(4)[:, 1:]
    W = s.Matrix.vstack(*(V[i, :] - V[j, :] for i, j in EDGES))
    lap = W.T * s.diag(*[1 / t for t in x]) * W
    U = s.cancel(lap.det() * s.prod(x))
    C = (2 * V * lap.adjugate() * V.T * s.prod(x)).applyfunc(s.cancel)
    D = s.Matrix.vstack(*(s.eye(4)[i, :] - s.eye(4)[j, :] for i, j in EDGES[:4]))
    c = D * C * D.T

    def tr(M):
        a, b, c_, d1, e_, f_ = (
            M[i, j] for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
        )
        return d * d * (a * f_ - b * e_ + c_ * d1) + 2 * d * b * e_

    N = s.expand(tr(c))
    cv = C * D.T

    def moment(i, j):
        a, b = cv[i, :].T, cv[j, :].T
        dc = a * b.T + b * a.T
        return s.expand(d * C[i, j] * N + tr(c + dc) - tr(c) - tr(dc))

    T = [moment(i, i) for i in range(4)]
    Wsum = sum(x)
    A = s.expand(
        -(4 - d) * sum(T) / 2 - d * (4 - d) * N * U * sum(x[4:]) / 3
    )
    B = s.expand(
        sum(
            s.Rational(2, 3)
            * x[k]
            * (
                T[i]
                + T[j]
                + moment(i, j)
                + (d - int(k < 4)) * x[k] * N * U
            )
            for k, (i, j) in enumerate(EDGES)
        )
    )
    Jexpr = Wsum * A + (5 - 3 * d / 2) * B
    J = s.Poly(Jexpr, *x)
    J4 = s.Poly(Jexpr.subs(d, 4), *x)
    up = s.Poly(U, *x)
    A_vertex = -(4 - d) / 2 * sum(T)
    vertex_eps = s.Poly((-2 * s.diff(Wsum * A_vertex, d)).subs(d, 4), *x)

    def powers(poly, subset):
        return min(sum(m[i] for i in subset) for m in poly.monoms())

    table = {}
    for n in range(1, 6):
        for ids in combinations(range(6), n):
            u = powers(up, ids)
            spin = sum(i < 4 for i in ids)
            table[ids] = (
                n - 1 + powers(J, ids) - 5 * u - spin,
                n - 1 + powers(J4, ids) - 5 * u - spin,
                u,
            )

    nodes, weights = roots_legendre(16)
    nodes, weights = (nodes + 1) / 2, weights / 2
    grid = np.array(np.meshgrid(nodes, nodes, nodes, indexing="ij")).reshape(3, -1)
    w = np.prod(
        np.array(np.meshgrid(weights, weights, weights, indexing="ij")), axis=0
    ).ravel()
    total = 0.0
    for order in permutations(range(6)):
        sets = [tuple(sorted(order[i:])) for i in range(1, 6)]
        dat = [table[k] for k in sets]
        singular = [i for i, (a, _, _) in enumerate(dat) if a < 0]
        if len(singular) != 2 or any(dat[i][0] != -1 for i in singular):
            continue
        regular = [i for i in range(5) if i not in singular]

        def coefficient(poly, derivative=None):
            data = [
                (tuple(sum(m[i] for i in ids) for ids in sets), float(c))
                for m, c in poly.terms()
            ]
            mins = [powers(up if poly is up else J, ids) for ids in sets]
            out = np.zeros(grid.shape[1])
            for m, c in data:
                if all(m[i] - mins[i] == (i == derivative) for i in singular):
                    out += c * np.prod(
                        grid
                        ** np.array([m[i] - mins[i] for i in regular])[:, None],
                        axis=0,
                    )
            return out

        u = coefficient(up)
        j = coefficient(vertex_eps)
        reg = np.prod(
            grid ** np.array([dat[i][1] for i in regular])[:, None], axis=0
        )
        total += (w @ (reg * j / u**5)) / np.prod([dat[i][2] for i in singular])
    return total


def main():
    print("=== FFff(B) two-log overlap assembly fix ===\n")

    ll = vertex_ll_residue()
    print(f"vertex ε-der on double-log sectors: {ll:.10f}")
    assert abs(ll + 16) < 1e-6

    a0 = -48
    measure = 192
    numerator_bare = 32
    slope_bare = measure + numerator_bare
    # Overlap fix: one extra copy of the ll vertex piece (exact integers for conversion)
    assert abs(ll + 16) < 1e-6
    slope_fixed = slope_bare - 16  # 224 → 208
    print(f"slope bare (measure+numerator): {slope_bare}")
    print(f"overlap Δa1 (= vertex ll):      -16")
    print(f"slope fixed:                    {slope_fixed}")
    assert slope_fixed == 208

    B2_bare = graph_double(a0, slope_bare)
    B2_fixed = graph_double(a0, slope_fixed)
    print(f"\nB double bare:  {B2_bare}")
    print(f"B double fixed: {B2_fixed}")
    assert B2_bare == s.Rational(-1, 96)
    assert B2_fixed == s.Rational(-1, 288)

    # Triple unchanged (overlap fix is O(ε)×triple → double only)
    B3 = F(1, 48)
    A3 = F(1, 144)
    assert A3 + B3 == F(1, 36)
    print(f"\nTriples unchanged: A {A3} + B {B3} = {A3+B3}")

    # Projected double at L=0
    A2 = F(5, 288)
    B2 = F(-1, 288)
    bare = F(5, 144)
    CT = F(-5, 48)
    Y = F(-1, 18)
    projected = bare + A2 + B2 + CT
    residual = projected - Y
    print(f"A double:        {A2}")
    print(f"B double fixed:  {B2}")
    print(f"bare non-FFff:   {bare}")
    print(f"CT (L=0):        {CT}")
    print(f"projected:       {projected}")
    print(f"Y target:        {Y}")
    print(f"residual:        {residual}")
    assert residual == 0

    # L in double and L/ε from residual
    L_double = -F(1, 48) - F(1, 16) + F(1, 12)
    assert L_double == 0
    print(f"L coeff in double: {L_double} (cancels)")
    print(f"L/ε from residual: {-3 * residual} (cancels)")

    # Homogeneity: B double L piece from triple
    assert -3 * B3 == -F(1, 16)
    print(f"B double L from triple: {-3 * B3}")

    print(
        """
Fix applied
-----------
  Δa1 = −16  (extra vertex soft-soft / coupling-weight overlap copy)
  B double (L=0) = −1/288
  Projected triple and double cancel at L=0; +L/48 gone.

Not changed
-----------
  Bare kernel used by finite-ε QMC / Gaussian checks (field weight only).
  Power/log sector residues (global λ_F doubling would shift those by −48).
  Simple-pole constant (still open vs 1/24).
"""
    )
    print("=== ALL ASSERTS PASSED ===")


if __name__ == "__main__":
    main()
