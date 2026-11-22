#!/usr/bin/env python3
"""Audit soft-linear Δc0 of the extra FFff(B) vertex-ll copy.

Computes, on double-log sectors only:

  1. Leading soft-soft residue (must be −16)
  2. W^{-3ε}+U^ε finite of that density (must be 64/3)
  3. Soft-linear finite with the U^{-5} product rule
     (dj − 5 j du/u), which is quadrature-stable at ≈ −105.304
     with U-rule shift exactly +30 relative to the raw linear

Physical scope of the overlap assembly fix
------------------------------------------
The Δa1=−16 correction is a *two-log overlap* assembly term: both soft
axes logarithmic. Soft-linear is a *single*-soft Taylor channel (one log
endpoint × regular). Global doubling of λ_F is already excluded because
it shifts power/log by −48. By the same region logic, soft-linear is
outside the overlap fix and must not be added to c0.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_ffff_B_soft_linear_audit.py

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
    pref = s.gamma(3 * e - 1) * s.exp(3 * s.EulerGamma * e) / (
        64 * d * (d - 1)
    )
    return s.simplify(
        s.series(pref * (a0 / e**2 + a1 / e + c0), e, 0, 1).coeff(e, -1)
    )


def build_and_integrate(n_quad=20):
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
            M[i, j]
            for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
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
    A = s.expand(A_vertex - d * (4 - d) * N * U * sum(x[4:]) / 3)
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
        np.array(np.meshgrid(weights, weights, weights, indexing="ij")),
        axis=0,
    ).ravel()

    residue = wu = lin = lin_raw = u_corr = 0.0
    for order in permutations(range(6)):
        sets = [tuple(sorted(order[i:])) for i in range(1, 6)]
        dat = [table[k] for k in sets]
        singular = [i for i, (a, _, _) in enumerate(dat) if a < 0]
        if len(singular) != 2 or any(dat[i][0] != -1 for i in singular):
            continue
        regular = [i for i in range(5) if i not in singular]
        s0, s1 = singular

        def coeff(poly, deriv=None):
            data = [
                (tuple(sum(m[i] for i in ids) for ids in sets), float(c_))
                for m, c_ in poly.terms()
            ]
            mins = [powers(up if poly is up else J, ids) for ids in sets]
            out = np.zeros(grid.shape[1])
            for m, c_ in data:
                if all(m[i] - mins[i] == (i == deriv) for i in singular):
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

        u = coeff(up)
        j0 = coeff(vertex_eps)
        du0, du1 = coeff(up, s0), coeff(up, s1)
        jv0, jv1 = coeff(vertex_eps, s0), coeff(vertex_eps, s1)
        j_s0 = jv0 - 5 * j0 * du0 / u
        j_s1 = jv1 - 5 * j0 * du1 / u
        u0, u1 = dat[s0][2], dat[s1][2]
        reg = np.prod(
            grid ** np.array([dat[i][1] for i in regular])[:, None], axis=0
        )
        dens0 = (reg * j0 / u**5) / (u0 * u1)
        W0 = np.clip(soft_value(Wpoly), 1e-30, None)
        u_safe = np.clip(u, 1e-30, None)
        residue += w @ dens0
        wu += w @ (dens0 * (-3 * np.log(W0) + np.log(u_safe)))
        lin += w @ ((reg / u**5) * (j_s0 / u1 + j_s1 / u0))
        lin_raw += w @ ((reg / u**5) * (jv0 / u1 + jv1 / u0))
        u_corr += w @ (
            (reg / u**5)
            * ((-5 * j0 * du0 / u) / u1 + (-5 * j0 * du1 / u) / u0)
        )
    return dict(
        residue=residue,
        wu=wu,
        lin=lin,
        lin_raw=lin_raw,
        u_corr=u_corr,
    )


def B_needed():
    mp.mp.dps = 40
    C0 = (
        -mp.mpf(15) / 2
        + 2
        * mp.sqrt(3)
        * mp.im(mp.polylog(2, mp.exp(1j * mp.pi / 3)))
    )
    ct = (C0 + 6) / 18 + (1 - mp.zeta(2)) / 24
    sub = (
        mp.mpf("-0.05121964552")
        + mp.mpf("-0.07466151326053")
        + mp.mpf("-0.02480")
        + 4 * mp.mpf("0.03724051498984436")
        + ct
    )
    return mp.mpf(1) / 24 - sub


def main():
    print("=== FFff(B) soft-linear audit (vertex-ll extra copy) ===\n")
    r = build_and_integrate(20)
    print(f"leading residue:     {r['residue']:.10f}  (expect −16)")
    print(f"W/U finite Δc0:      {r['wu']:.10f}  (expect 64/3)")
    print(f"soft-linear (U rule):{r['lin']:.10f}")
    print(f"soft-linear raw:     {r['lin_raw']:.10f}")
    print(f"U-rule shift:        {r['u_corr']:.10f}  (expect +30)")

    assert abs(r["residue"] + 16) < 1e-6
    assert abs(r["wu"] - float(F(64, 3))) < 1e-6
    assert abs(r["u_corr"] - 30) < 1e-6

    need = float(B_needed())
    scenarios = {
        "overlap only: 192+64/3": 192 + F(64, 3),
        "wrongly + soft-linear": 192 + F(64, 3) + r["lin"],
        "empirical 2×(64/3)": 192 + 2 * F(64, 3),
    }
    print(f"\nB needed at L=0: {need:+.10f}")
    print(f"{'scenario':32s}  {'B simple':>14}  {'gap':>14}")
    for name, c0 in scenarios.items():
        B = float(graph_simple(-48, 208, c0))
        print(f"  {name:30s}  {B:+14.10f}  {B-need:+14.10f}")

    # A projective c0 that would close with overlap-only B
    B_overlap = float(graph_simple(-48, 208, 192 + F(64, 3)))
    mp.mp.dps = 40
    C0 = (
        -mp.mpf(15) / 2
        + 2
        * mp.sqrt(3)
        * mp.im(mp.polylog(2, mp.exp(1j * mp.pi / 3)))
    )
    ct = (C0 + 6) / 18 + (1 - mp.zeta(2)) / 24
    base = (
        mp.mpf("-0.05121964552")
        + mp.mpf("-0.07466151326053")
        + mp.mpf("-0.02480")
        + ct
    )
    A_raw_needed = (mp.mpf(1) / 24 - base - B_overlap) / 4
    c0_A_needed = 1152 * (
        A_raw_needed - mp.mpf(79) / 6912 - mp.pi**2 / 768
    )
    print(
        f"\nTo close with overlap-only B, FFff(A) projective c0 must be "
        f"{float(c0_A_needed):.4f} (reported ~14.93; Δ≈{float(c0_A_needed-14.93):+.2f})."
    )

    print(
        """
Verdict
-------
  Soft-linear is real and stable, but it is a single-soft channel.
  The overlap assembly fix only corrects two-log soft-soft. Do NOT add
  soft-linear to c0. Remaining gap after 192+64/3 is ~+0.009 — same size
  as a ~2.6 shift in A's projective finite coefficient. Next: re-extract
  FFff(A) simple / projective c0.
"""
    )
    print("=== AUDIT COMPLETE ===")


if __name__ == "__main__":
    main()
