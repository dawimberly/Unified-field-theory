#!/usr/bin/env python3
"""Projected Newton simple-pole ledger at a_X=0, a_Y=-1 (L=0).

Assembles known L=0 simple constants, converts FFff(B) projective
(a0, a1, c0) → graph simple, and evaluates the gap vs target 1/24.

Overlap-fixed slope a1=208 is used. The field-weight regulated finite
projective coefficient is taken as c0≈192 (FFFF_B_REGULATED.md). The
extra vertex-ll copy that supplied Δa1=−16 also contributes an exact
Δc0=64/3 from the W^{-3ε} and U^ε expansions on double-log sectors.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_ffff_B_simple_ledger.py

Requires: numpy, scipy, sympy, mpmath
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations, permutations

import mpmath as mp
import numpy as np
import sympy as s
from scipy.special import roots_legendre

EDGES = ((0, 2), (1, 2), (1, 3), (0, 3), (0, 2), (1, 3))


def graph_simple(a0, a1, c0):
    e = s.symbols("e")
    d = 4 - 2 * e
    pref = s.gamma(3 * e - 1) * s.exp(3 * s.EulerGamma * e) / (64 * d * (d - 1))
    ser = s.series(pref * (a0 / e**2 + a1 / e + c0), e, 0, 1)
    return s.simplify(ser.coeff(e, -1))


def sunset_C0():
    mp.mp.dps = 40
    return (
        -mp.mpf(15) / 2
        + 2
        * mp.sqrt(3)
        * mp.im(mp.polylog(2, mp.exp(1j * mp.pi / 3)))
    )


def ct_simple_L0():
    """Combined one+two-loop curvature counterterm constant at L=0."""
    C0 = sunset_C0()
    return (C0 + 6) / 18 + (1 - mp.zeta(2)) / 24


def known_pieces():
    """L=0 simple constants in units h²z/K (packet + review QMC)."""
    return {
        "scalar_QQ_QCc_CCcc": mp.mpf("-0.05121964552"),
        "QFf": mp.mpf("-0.07466151326053"),
        "CcFf": mp.mpf("-0.02480"),  # review QMC; checkpoint quoted −0.0247
        "FFff_A": 4 * mp.mpf("0.03724051498984436"),  # ×4 spin normalization
        "CT": ct_simple_L0(),
    }


def vertex_ll_WU_finite(n_quad=16):
    """Finite projective Δc0 of one extra vertex-ll copy from W^{-3ε} and U^ε.

    Returns (residue, delta_c0). Residue must be −16; delta_c0 is 64/3.
    """
    d = s.symbols("d")
    x = s.symbols("a b c e f g")
    V = s.eye(4)[:, 1:]
    W = s.Matrix.vstack(*(V[i, :] - V[j, :] for i, j in EDGES))
    lap = W.T * s.diag(*[1 / t for t in x]) * W
    U = s.cancel(lap.det() * s.prod(x))
    C = (2 * V * lap.adjugate() * V.T * s.prod(x)).applyfunc(s.cancel)
    D = s.Matrix.vstack(
        *(s.eye(4)[i, :] - s.eye(4)[j, :] for i, j in EDGES[:4])
    )
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
    A_vertex = -(4 - d) / 2 * sum(T)
    A = s.expand(
        A_vertex - d * (4 - d) * N * U * sum(x[4:]) / 3
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
    vertex_eps = s.Poly((-2 * s.diff(Wsum * A_vertex, d)).subs(d, 4), *x)
    Upoly = s.Poly(s.expand(U), *x)
    Wpoly = s.Poly(Wsum, *x)

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

    nodes, weights = roots_legendre(n_quad)
    nodes, weights = (nodes + 1) / 2, weights / 2
    grid = np.array(np.meshgrid(nodes, nodes, nodes, indexing="ij")).reshape(
        3, -1
    )
    w = np.prod(
        np.array(np.meshgrid(weights, weights, weights, indexing="ij")), axis=0
    ).ravel()

    residue = 0.0
    delta = 0.0
    for order in permutations(range(6)):
        sets = [tuple(sorted(order[i:])) for i in range(1, 6)]
        dat = [table[k] for k in sets]
        singular = [i for i, (a, _, _) in enumerate(dat) if a < 0]
        if len(singular) != 2 or any(dat[i][0] != -1 for i in singular):
            continue
        regular = [i for i in range(5) if i not in singular]

        def coefficient(poly):
            data = [
                (tuple(sum(m[i] for i in ids) for ids in sets), float(c_))
                for m, c_ in poly.terms()
            ]
            mins = [powers(up if poly is up else J, ids) for ids in sets]
            out = np.zeros(grid.shape[1])
            for m, c_ in data:
                if all(m[i] - mins[i] == 0 for i in singular):
                    out += c_ * np.prod(
                        grid
                        ** np.array(
                            [m[i] - mins[i] for i in regular]
                        )[:, None],
                        axis=0,
                    )
            return out

        def soft_value(poly):
            data = [
                (tuple(sum(m[i] for i in ids) for ids in sets), float(c_))
                for m, c_ in poly.terms()
            ]
            mins = [powers(poly, ids) for ids in sets]
            out = np.zeros(grid.shape[1])
            for m, c_ in data:
                if all(m[i] - mins[i] == 0 for i in singular):
                    out += c_ * np.prod(
                        grid
                        ** np.array(
                            [m[i] - mins[i] for i in regular]
                        )[:, None],
                        axis=0,
                    )
            return out

        u = coefficient(up)
        j = coefficient(vertex_eps)
        reg = np.prod(
            grid ** np.array([dat[i][1] for i in regular])[:, None], axis=0
        )
        dens = (reg * j / u**5) / np.prod([dat[i][2] for i in singular])
        W0 = np.clip(soft_value(Wpoly), 1e-30, None)
        u_safe = np.clip(u, 1e-30, None)
        residue += w @ dens
        delta += w @ (dens * (-3 * np.log(W0) + np.log(u_safe)))
    return residue, delta


def main():
    print("=== FFff(B) projected simple-pole ledger (L=0) ===\n")

    pieces = known_pieces()
    subtotal = sum(pieces.values())
    target = mp.mpf(1) / 24
    B_needed = target - subtotal

    print("Known L=0 simple pieces:")
    for k, v in pieces.items():
        print(f"  {k:22s}  {float(v):+.10f}")
    print(f"  {'subtotal':22s}  {float(subtotal):+.10f}")
    print(f"  {'target 1/24':22s}  {float(target):+.10f}")
    print(f"  {'B simple needed':22s}  {float(B_needed):+.10f}")

    # Exact conversion identities
    a0, a1 = -48, 208
    s0 = graph_simple(a0, a1, 0)
    ds = graph_simple(a0, a1, 1) - s0
    assert ds == s.Rational(-1, 2304)
    print(f"\nConversion: simple = ({s0}) + c0*({ds})")

    # Vertex-ll finite from W,U
    ll, d_c0 = vertex_ll_WU_finite()
    print(f"\nvertex-ll residue:     {ll:.10f}  (expect −16)")
    print(f"vertex-ll Δc0 (W+U):   {d_c0:.10f}  (expect 64/3)")
    assert abs(ll + 16) < 1e-6
    assert abs(d_c0 - float(F(64, 3))) < 1e-6
    d_c0 = F(64, 3)

    c0_bare = 192  # field-weight regulated quadratic (FFFF_B_REGULATED)
    scenarios = {
        "c0=192 (field-weight QMC only)": c0_bare,
        "c0=192+64/3 (add vertex-ll W/U)": c0_bare + d_c0,
        "c0=192+2*(64/3) (empirical 2×; not derived)": c0_bare + 2 * d_c0,
    }

    print("\nB simple under c0 scenarios (a0=-48, a1=208):")
    print(f"  {'scenario':48s}  {'c0':>10}  {'B simple':>14}  {'gap':>14}")
    for name, c0 in scenarios.items():
        B = graph_simple(a0, a1, c0)
        gap = float(B) - float(B_needed)
        print(
            f"  {name:48s}  {float(c0):10.4f}  {float(B):+14.10f}  {gap:+14.10f}"
        )

    # Homogeneity logs (not in L=0 constant)
    assert -3 * F(1, 48) == -F(1, 16)
    assert -3 * F(-1, 288) == F(1, 96)
    assert F(9, 2) * F(1, 48) == F(3, 32)
    print("\nHomogeneity in bare B simple (for L≠0 checks):")
    print("  L from double −1/288:  +L/96")
    print("  L² from triple 1/48:   +3L²/32")

    print(
        """
Reading
-------
  Higher poles cancel after the overlap fix. The simple test is open.
  Adding the derived vertex-ll W/U finite 64/3 halves the gap but does
  not close it. Soft-linear and other O(ε)×ll finite pieces are not yet
  controlled. FFff(A)'s simple constant is a numerical extrapolation and
  remains the other large uncertainty. The empirical 2×(64/3) scenario
  lands inside CcFf/A noise — coincidence until derived.
"""
    )
    print("=== LEDGER COMPLETE ===")


if __name__ == "__main__":
    main()
