#!/usr/bin/env python3
"""FFff(B) Wick multiplicity (PSE framing) and ε×triple → double contamination.

Answers two follow-ups from the parallel ChatGPT audit:

1. Allowed contractions: Yukawa vertex = one scalar + two Weyl fermions, so a
   connected four-scalar response needs four vertices; pairing the four scalar
   legs into two internal scalars yields L = 6−4+1 = 3. Reconstruct the signed
   Wick weight used in the pole conversion.

2. Why the double can fail while the triple passes: an O(ε) relative error on the
   triple density shifts only the double. Quantify the shift that would make the
   B double equal the budget target −1/288, and split the numerator ε-derivative
   residue (+32) into explicit (k−1) vs ∂A/∂d, ∂B/∂d pieces.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_ffff_B_contraction_contamination.py

Requires: numpy, scipy, sympy
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations, permutations

import numpy as np
import sympy as s
from scipy.special import roots_legendre


EDGES = ((0, 2), (1, 2), (1, 3), (0, 3), (0, 2), (1, 3))


def topology_and_wick():
    print("=== Topology (Tong / PSE framing → Weyl Yukawa) ===")
    print("Vertex content: one scalar + two Weyl fermions (no φ²ψ vertex).")
    print("Four external/internal scalar legs ⇒ four Yukawa vertices.")
    print("Pairing those into two internal scalars + four fermion lines:")
    print("  V=4, I=6, L=I−V+1=3  (FFff(B) three-loop graph)")
    assert 6 - 4 + 1 == 3

    # Fermion wirings: match 4 outs to 4 inns.
    n4 = n22 = 0
    for perm in permutations(range(4)):
        seen = [False] * 4
        lengths = []
        for start in range(4):
            if seen[start]:
                continue
            i, L = start, 0
            while not seen[i]:
                seen[i] = True
                i = perm[i]
                L += 1
            lengths.append(L)
        if lengths == [4]:
            n4 += 1
        elif sorted(lengths) == [2, 2]:
            n22 += 1
    print(f"Fermion pairings 4!={24}: single 4-cycles={n4}, double 2-cycles={n22}")
    assert n4 == 6 and n22 == 3
    print("FFff(B) keeps the single 4-cycle family; FFff(A) is the 2+2 family.")

    print("\n=== Signed Wick weight ===")
    # Packet: 32 matchings, vertex halves, factorial symmetry, cumulant/trace signs.
    # Runnable pole convention uses overall +1 after z|y|^4.
    wick = 32
    half4 = F(1, 2) ** 4
    one_sym = wick * half4 / 2
    two_sym = wick * half4 / 4
    print(f"32 × (1/2)^4 / 2!  = {one_sym}   ← recovers runnable +1")
    print(f"32 × (1/2)^4 / 2!² = {two_sym}   ← packet '2!²' wording")
    assert one_sym == 1
    assert two_sym == F(1, 2)
    print(
        "FLAG: packet text says 2!² (would give +1/2); the pole scripts use +1."
        " A global 1/2 would break the triple sum 1/36, so the runnable +1 is the"
        " one consistent with the passing triple check. Treat 2!² as a wording"
        " ambiguity, not as evidence to rescale the double alone."
    )


def contamination_arithmetic():
    print("\n=== ε×triple contamination of the double ===")
    e = s.symbols("e")
    d = 4 - 2 * e
    pref = s.gamma(3 * e - 1) * s.exp(3 * s.EulerGamma * e) / (64 * d * (d - 1))
    a1 = s.symbols("a1")
    I = -48 / e**2 + a1 / e
    dbl = s.simplify(s.series(pref * I, e, 0, 0).coeff(e, -2))
    # dbl = 25/288 - a1/2304
    a1_needed = s.solve(s.Eq(dbl, s.Rational(-1, 288)), a1)[0]
    print(f"projective I = -48/ε² + a1/ε + …")
    print(f"graph double = {dbl}")
    print(f"a1 for double -1/96 (current): 224 → {s.simplify(dbl.subs(a1, 224))}")
    print(f"a1 for double -1/288 (budget): {a1_needed}")
    print(f"Δa1 required: {a1_needed - 224}")
    assert a1_needed == 208
    # Relative O(ε) on the triple density: (1+δ ε)(-48/ε²)= -48/ε² -48δ/ε
    delta = s.Rational(1, 3)
    assert -48 * delta == -16
    print(f"(1+δ ε)×triple density with δ=1/3 shifts a1 by {-48 * delta} — exactly the need.")
    print(
        f"Equivalently: numerator ε-der residue +32 instead of +16 "
        f"(measure 192 + 16 = 208). Gap in graph units: "
        f"{F(-1, 96) - F(-1, 288)} = -1/144 = -δ × B_triple."
    )
    assert F(-1, 96) - F(-1, 288) == F(-1, 144)
    assert F(-1, 144) == -delta * F(1, 48)


def split_numerator_slope():
    print("\n=== Split of numerator ε-derivative residue (+32) ===")
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
    km1 = 5 - 3 * d / 2
    Jexpr = sum(x) * A + km1 * B
    J = s.Poly(Jexpr, *x)
    J4 = s.Poly(Jexpr.subs(d, 4), *x)
    up = s.Poly(U, *x)

    def powers(poly, subset):
        return min(sum(m[i] for i in subset) for m in poly.monoms())

    table = {}
    for n in range(1, 6):
        for ids in combinations(range(6), n):
            u = powers(up, ids)
            spin = sum(i < 4 for i in ids)
            E0 = n - 1 + powers(J, ids) - 5 * u - spin
            E4 = n - 1 + powers(J4, ids) - 5 * u - spin
            table[ids] = (E0, E4, u)

    def residue(n, numerator):
        nodes, weights = roots_legendre(n)
        nodes, weights = (nodes + 1) / 2, weights / 2
        grid = np.array(np.meshgrid(nodes, nodes, nodes, indexing="ij")).reshape(
            3, -1
        )
        w = np.prod(
            np.array(np.meshgrid(weights, weights, weights, indexing="ij")),
            axis=0,
        ).ravel()
        total = 0.0
        for order in permutations(range(6)):
            sets = [tuple(sorted(order[i:])) for i in range(1, 6)]
            dat = [table[k] for k in sets]
            singular = [i for i, (a, _, _) in enumerate(dat) if a < 0]
            if len(singular) != 2:
                continue
            regular = [i for i in range(5) if i not in singular]

            def coefficient(poly, derivative=None):
                data = [
                    (tuple(sum(m[i] for i in ids) for ids in sets), float(c))
                    for m, c in poly.terms()
                ]
                mins = [powers(poly, ids) for ids in sets]
                out = np.zeros(grid.shape[1])
                for m, c in data:
                    if all(m[i] - mins[i] == (i == derivative) for i in singular):
                        out += c * np.prod(
                            grid
                            ** np.array([m[i] - mins[i] for i in regular])[
                                :, None
                            ],
                            axis=0,
                        )
                return out

            u = coefficient(up)
            j = coefficient(numerator)
            power = next((i for i in singular if dat[i][0] == -2), None)
            if power is not None:
                j = (
                    coefficient(numerator, power)
                    - 5 * j * coefficient(up, power) / u
                )
            reg = np.prod(
                grid ** np.array([dat[i][1] for i in regular])[:, None], axis=0
            )
            total += (w @ (reg * j / u**5)) / np.prod(
                [dat[i][2] for i in singular]
            )
        return total

    # ∂J/∂ε = -2 ∂J/∂d; at d=4 the explicit ∂(k−1)/∂ε piece is +3 B.
    full = s.Poly((-2 * s.diff(Jexpr, d)).subs(d, 4), *x)
    explicit = s.Poly((3 * B).subs(d, 4), *x)
    B4 = s.Poly(B.subs(d, 4), *x)
    r_full = residue(16, full)
    r_expl = residue(16, explicit)
    r_B = residue(16, B4)
    print(f"full numerator ε-der residue:     {r_full:.8f}  (expect 32)")
    print(f"explicit (k−1) piece (+3 B):      {r_expl:.8f}  (expect 144)")
    print(f"residue of B at d=4:              {r_B:.8f}  (expect 48 = −J4 residue)")
    assert abs(r_full - 32) < 1e-6
    assert abs(r_expl - 144) < 1e-6
    assert abs(r_B - 48) < 1e-6
    print(
        "Dropping the explicit (k−1) piece is NOT the δ=1/3 fix: it removes 144,"
        " not 16, and drives the slope far from 208."
    )
    print(
        "The budget's 'need +16 not +32' from FFFF_A_ANALYTIC_DOUBLE.md is the"
        " cancellation requirement, not an independent evaluation of the residue."
    )
    print(
        "Radial check: numeric A+(k−1)B with Γ(k−1) equals Γ(k−1)A+Γ(k)B, so the"
        " explicit (k−1) in the projective numerator is required once, not twice."
    )


def why_double_fails():
    print("\n=== Why this review says the double is failing ===")
    print("Not a vague suspicion — concrete ledger:")
    print("  • Triple passes: A 1/144 + B 1/48 = 1/36.")
    print("  • Independent B double (residue/regulated) = −1/96, not −1/288.")
    print("  • Residual vs Y at L=0: −1/144, first bad term = that B constant.")
    print("  • e^{−3εL} maps the residual into uncanceled +L/48 in the simple pole.")
    print("  • Scalar-mass (+1/18) and E-split (−1/54) pieces passed as internals.")
    print("  • ChatGPT's 59/648 is a different assembly (1/8−11/324); same qualitative")
    print("    failure mode (cutoff logs cancel, L/ε does not).")
    print(
        "Leading mechanical suspect consistent with 'ε×triple leaves double':"
        " a δ=1/3 relative O(ε) on the triple density, or anything else that removes"
        " 16 from the projective slope (224→208). Not yet identified in the kernel."
    )


def main():
    topology_and_wick()
    contamination_arithmetic()
    split_numerator_slope()
    why_double_fails()
    print("\n=== ALL ASSERTS PASSED ===")


if __name__ == "__main__":
    main()
