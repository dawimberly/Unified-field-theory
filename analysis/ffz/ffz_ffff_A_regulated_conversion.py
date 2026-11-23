"""Exact conversion of the FFff(A) projective Laurent series to graph poles.

Takes the numerical projective inference
    I = 2/ε² − (10/3)/ε + c₀ + O(ε)
from FFFF_A_REGULATED.md and multiplies by
    −Γ(3ε−1) exp(3γ_E ε − 3ε L) / [32 d(d−1)]
at d = 4−2ε. With γ_E identified with EulerGamma this yields the claimed
graph poles in units h²z/K:
    1/(576 ε³) + (5/1152 − L/192)/ε² + (c₀/1152 + 79/6912 + π²/768 − 5L/384 + L²/128)/ε
and the implied FFff(B) higher-pole targets against the combined bare budget.

Run:
  python ffz_ffff_A_regulated_conversion.py

Requires: sympy
"""

import sympy as s


def conversion_check():
    e, L = s.symbols("eps L")
    gE = s.EulerGamma
    d = 4 - 2 * e
    # Γ(3ε−1) = Γ(1+3ε)/((3ε−1) 3ε)
    G = s.gamma(1 + 3 * e) / ((3 * e - 1) * 3 * e)
    pref = -G * s.exp(3 * gE * e - 3 * e * L) / (32 * d * (d - 1))
    ser = s.series(pref, e, 0, 2)
    coeffs = {n: s.simplify(ser.coeff(e, n)) for n in (-1, 0, 1)}

    assert coeffs[-1] == s.Rational(1, 1152)

    c0 = s.symbols("c0")
    # I = 2/e^2 - (10/3)/e + c0
    triple = s.simplify(coeffs[-1] * 2)
    double = s.simplify(coeffs[-1] * (-s.Rational(10, 3)) + coeffs[0] * 2)
    simple = s.simplify(
        coeffs[-1] * c0 + coeffs[0] * (-s.Rational(10, 3)) + coeffs[1] * 2
    )

    assert triple == s.Rational(1, 576)
    assert s.simplify(double - (s.Rational(5, 1152) - L / 192)) == 0

    poly = s.Poly(s.expand(simple), L)
    assert poly.coeff_monomial(L**2) == s.Rational(1, 128)
    assert s.simplify(poly.coeff_monomial(L) - (-s.Rational(5, 384))) == 0
    assert s.simplify(simple.subs(L, 0) - (c0 / 1152 + s.Rational(79, 6912) + s.pi**2 / 768)) == 0

    # Combined bare FFff budget minus this A inference => B targets
    b_triple = s.Rational(1, 36) - s.Rational(1, 576)
    b_double = s.Rational(1, 72) - s.Rational(5, 1152) - L / 12 + L / 192
    assert b_triple == s.Rational(5, 192)
    assert s.simplify(b_double - (s.Rational(11, 1152) - s.Rational(5, 64) * L)) == 0

    # Analytic two-tail leading coefficient (Dirichlet form in FFFF_A_REGULATED).
    # Rewrite Beta through Gamma so the double pole is visible to series.
    k = 3 * e
    d_sym = 4 - 2 * e

    def D(i, j):
        a = e - 1 + i
        b = 1 - 2 * e + j
        c = 4 - e
        return s.gamma(a) * s.gamma(b) * s.gamma(c - a - b) / s.gamma(c)

    def B(a, b):
        return s.gamma(a) * s.gamma(b) / s.gamma(a + b)

    one_tail = (
        (k - 1)
        * 2
        * (d_sym - 1)
        / 3
        * d_sym
        / 4
        * B(2 - e, 2 - e)
        * B(2 * e, 2 * e)
        * (
            d_sym * (D(0, 0) + D(1, 0) + D(0, 1) + D(1, 1))
            + 2 * D(1, 1)
        )
    )
    two_tails = s.series(2 * one_tail, e, 0, 0)
    assert s.simplify(two_tails.coeff(e, -2) - 2) == 0

    return {
        "triple": triple,
        "double": double,
        "simple": simple,
        "B_triple_target": b_triple,
        "B_double_target": s.simplify(b_double),
    }


def main():
    out = conversion_check()
    print("FFff(A) triple:", out["triple"])
    print("FFff(A) double:", out["double"])
    print("FFff(A) simple:", out["simple"])
    print("implied FFff(B) triple target:", out["B_triple_target"])
    print("implied FFff(B) double target:", out["B_double_target"])
    print("FFFF(A) REGULATED CONVERSION CHECKS PASSED")


if __name__ == "__main__":
    main()
