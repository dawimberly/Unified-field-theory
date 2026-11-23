#!/usr/bin/env python3
"""Locate the concrete −16 projective-slope term in FFff(B).

Finds that the fermion-vertex Weyl-weight ε-derivative, restricted to
double-log (two overlapping logarithmic endpoints) sectors, has residue
exactly −16 — the unique shift that takes slope 224→208 and B double
−1/96→−1/288, cancelling the projected −1/144 residual and the uncanceled
+L/48.

Hypothesis (not yet a derivation): the Yukawa *coupling* Weyl weight ε
(packet: cubic/quartic/Yukawa weights 1+ε, 2ε, ε) is omitted beside the
field weight λ_F=(4−d)/2=ε, so the O(ε)×soft-soft vertex piece is
undercounted by a factor of two on the two-log overlap.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_ffff_B_slope16_hunt.py

Requires: numpy, scipy, sympy
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations, permutations

import numpy as np
import sympy as s
from scipy.special import roots_legendre

EDGES = ((0, 2), (1, 2), (1, 3), (0, 3), (0, 2), (1, 3))


def build():
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
    A_vertex = -(4 - d) / 2 * sum(T)
    A_scalar = -d * (4 - d) * N * U * sum(x[4:]) / 3
    A = s.expand(A_vertex + A_scalar)
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
    return dict(
        d=d,
        x=x,
        T=T,
        Wsum=Wsum,
        A_vertex=A_vertex,
        A_scalar=A_scalar,
        B=B,
        J=J,
        J4=J4,
        up=up,
        Jexpr=Jexpr,
    )


def powers(poly, subset):
    return min(sum(m[i] for i in subset) for m in poly.monoms())


def sector_table(J, J4, up):
    table = {}
    for n in range(1, 6):
        for ids in combinations(range(6), n):
            u = powers(up, ids)
            spin = sum(i < 4 for i in ids)
            E0 = n - 1 + powers(J, ids) - 5 * u - spin
            E4 = n - 1 + powers(J4, ids) - 5 * u - spin
            table[ids] = (E0, E4, u)
    return table


def residue(n, numerator, table, up, J, sector="all"):
    nodes, weights = roots_legendre(n)
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
        if len(singular) != 2:
            continue
        is_ll = all(dat[i][0] == -1 for i in singular)
        is_pl = any(dat[i][0] == -2 for i in singular) and any(
            dat[i][0] == -1 for i in singular
        )
        if sector == "ll" and not is_ll:
            continue
        if sector == "pl" and not is_pl:
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
        j = coefficient(numerator)
        power = next((i for i in singular if dat[i][0] == -2), None)
        if power is not None:
            j = coefficient(numerator, power) - 5 * j * coefficient(up, power) / u
        reg = np.prod(
            grid ** np.array([dat[i][1] for i in regular])[:, None], axis=0
        )
        total += (w @ (reg * j / u**5)) / np.prod([dat[i][2] for i in singular])
    return total


def graph_double(a0, a1):
    e = s.symbols("e")
    d = 4 - 2 * e
    pref = s.gamma(3 * e - 1) * s.exp(3 * s.EulerGamma * e) / (64 * d * (d - 1))
    return s.simplify(
        s.series(pref * (a0 / e**2 + a1 / e), e, 0, 0).coeff(e, -2)
    )


def main():
    print("=== FFff(B) −16 slope hunt ===\n")
    P = build()
    d, x, Wsum, J, J4, up, Jexpr, B = (
        P[k] for k in ("d", "x", "Wsum", "J", "J4", "up", "Jexpr", "B")
    )
    table = sector_table(J, J4, up)

    Jeps = s.Poly((-2 * s.diff(Jexpr, d)).subs(d, 4), *x)
    vertex_eps = s.Poly(
        (-2 * s.diff(Wsum * P["A_vertex"], d)).subs(d, 4), *x
    )
    scalar_eps = s.Poly(
        (-2 * s.diff(Wsum * P["A_scalar"], d)).subs(d, 4), *x
    )
    explicit = s.Poly((3 * B).subs(d, 4), *x)
    Bd = s.Poly((2 * s.diff(B, d)).subs(d, 4), *x)

    print("Numerator ε-der residues by sector:")
    splits = {}
    for name, poly in (
        ("full", Jeps),
        ("vertex", vertex_eps),
        ("scalar", scalar_eps),
        ("explicit (k-1)", explicit),
        ("B_d", Bd),
    ):
        ll = residue(16, poly, table, up, J, "ll")
        pl = residue(16, poly, table, up, J, "pl")
        splits[name] = (ll, pl, ll + pl)
        print(f"  {name:16s}  ll={ll:12.6f}  pl={pl:12.6f}  sum={ll+pl:12.6f}")

    assert abs(splits["full"][2] - 32) < 1e-6
    assert abs(splits["vertex"][0] + 16) < 1e-6  # ll vertex == -16
    assert abs(splits["scalar"][0]) < 1e-6  # scalar absent on ll
    print("\nSMOKING GUN: vertex ε-der on double-log sectors =", splits["vertex"][0])
    print("             required Δa1 for budget double −1/288 = -16")

    # Doubling the ll vertex piece: full numerator 32 + (-16) = 16
    num_fixed = 32 + splits["vertex"][0]  # add one more copy of ll vertex
    assert abs(num_fixed - 16) < 1e-6
    slope_fixed = 192 + 16
    print(f"\nIf ll vertex is undercounted by ×2: numerator {num_fixed:.0f}, slope {slope_fixed}")
    print("  double(−48, 224) =", graph_double(-48, 224))
    print("  double(−48, 208) =", graph_double(-48, 208))
    assert graph_double(-48, 208) == s.Rational(-1, 288)

    # Projected ledger
    A2, B2, bare, CT, Y = F(5, 288), F(-1, 288), F(5, 144), F(-5, 48), F(-1, 18)
    residual = bare + A2 + B2 + CT - Y
    print(f"\nProjected residual with B=−1/288: {residual}")
    assert residual == 0
    print("L/ε from residual:", -3 * residual, "(cancels)")

    print(
        """
Hypothesis (derived match; physical justification provisional):
  Packet assigns Yukawa *coupling* Weyl weight ε, equal at this order to
  λ_F=(4−d)/2. The kernel includes the field-weight insertion −λ_F Σ M_ii
  but not a parallel coupling-weight insertion. On the two-log overlap the
  vertex O(ε)×soft-soft residue is −16; a second equal copy supplies the
  missing −16 of projective slope.

Status: first concrete term that accounts for the full double-pole gap and
the uncanceled +L/48. Still needs an independent derivation that the
coupling-weight insertion has the same soft-soft residue as λ_F (not merely
the same engineering weight ε).
"""
    )
    print("=== ALL ASSERTS PASSED ===")


if __name__ == "__main__":
    main()
