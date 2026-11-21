"""FFff(A) soft-corner poles for the FFZ three-loop Newton sum.

Exact SymPy averages of the short-pair Taylor coefficients of the curvature
density H, and the resulting triple-pole coefficient. The double pole still
needs the once-soft faces (and any long-pair tails); those are not claimed here.

Run:
  python ffz_ffff_A_poles.py

Requires: sympy
"""

from fractions import Fraction as F
import sympy as sp


def taylor_averages():
    """Angular averages of H = a00 + bq*q + br*r + bqr*q*r + ... at d=4.

    H = (q N1)(r N2) S^{-2} B multiplies the projective measure
    q^{ε-2} r^{ε-2} dq dr dg dh dv on the gauge v+w=1.
    """
    g, h, v = sp.symbols("g h v")
    # Coefficients in effective lengths Q=q g(1-g), R=r h(1-h), from the
    # independent Wick expansion at d=4 (see review notes).
    a00 = sp.Rational(32, 3)
    a10 = 0  # coeff of Q
    a01 = (
        32
        * v
        * (3 * h**2 * v - 3 * h**2 - 3 * h * v + 3 * h + 2 * v - 2)
        / (h * (h - 1))
    )
    poly = (
        -81 * h**2 * v**2
        + 54 * h**2 * v
        - h**2
        + 81 * h * v**2
        - 54 * h * v
        + h
        - 54 * v**2
        + 36 * v
    )
    a11 = 16 * poly / (3 * h * (h - 1))

    # Convert to physical q,r: ∂_r H = a01 * h(1-h), ∂_q∂_r H = a11*g(1-g)*h(1-h).
    br = sp.simplify(a01 * h * (1 - h))
    # h(1-h)/(h(h-1)) = -1
    assert sp.simplify(h * (1 - h) / (h * (h - 1)) + 1) == 0
    br = sp.simplify(-32 * v * (3 * h**2 * v - 3 * h**2 - 3 * h * v + 3 * h + 2 * v - 2))
    bqr = sp.simplify(a11 * g * (1 - g) * h * (1 - h))
    # a11*h(1-h) = -16/3 * poly
    bqr = sp.Rational(-16, 3) * poly * g * (1 - g)

    avg_br = sp.integrate(sp.integrate(br, (h, 0, 1)), (v, 0, 1))
    avg_bqr = sp.integrate(
        sp.integrate(sp.integrate(bqr, (g, 0, 1)), (h, 0, 1)), (v, 0, 1)
    )
    assert avg_br == 8
    assert avg_bqr == sp.Rational(-4, 27)
    return {
        "a00": a00,
        "bq": 0,
        "br": avg_br,
        "bqr": avg_bqr,
    }


def poles_from_corner(avg):
    """Projective soft-corner poles times Wick and curvature normalization.

    I = a00/(ε-1)^2 + bq/(ε(ε-1)) + br/(ε(ε-1)) + bqr/ε^2 + (faces...)
    The constant a00 produces no pole. The bilinear bqr produces 1/ε^2.
    Gamma(3ε)~1/(3ε), weight -1/8, 1/[4d(d-1)]=1/48 at d=4 => -1/384.
    """
    e = sp.symbols("e")
    I_corner = (
        avg["a00"] / (e - 1) ** 2
        + avg["bq"] / (e * (e - 1))
        + avg["br"] / (e * (e - 1))
        + avg["bqr"] / e**2
    )
    # Leading Laurent of the corner alone (faces omitted).
    Iser = sp.series(I_corner, e, 0, 1).removeO()
    pre = sp.Rational(-1, 384)
    graph = sp.series(pre * (1 / (3 * e)) * I_corner, e, 0, 0)
    triple = pre * sp.Rational(1, 3) * avg["bqr"]  # 1/(3ε) * bqr/ε^2
    # Linear corner piece alone (incomplete for the double pole without faces):
    double_corner_only = pre * sp.Rational(1, 3) * (-(avg["bq"] + avg["br"]))
    assert triple == sp.Rational(1, 7776)
    assert double_corner_only == sp.Rational(1, 144)
    return {
        "I_corner_series": Iser,
        "triple_pole": triple,
        "double_from_corner_linear_only": double_corner_only,
        "note": "double pole incomplete without once-soft faces",
    }


def assembly_snapshot():
    """Known higher poles versus the Y-matching targets. Not a finish-line claim."""
    target_triple = F(1, 36)
    target_double = F(-1, 18)
    mass_ct_triple = F(1, 12)
    mass_ct_double = F(-1, 16)  # constant piece; also -L/12
    fffa_triple = F(1, 7776)
    # Verified doubles among bare graphs so far (faces for FFff(A) open):
    qq = F(1, 72)  # packet, not re-derived in this review
    qcc = F(1, 72)  # independently verified
    ccff = F(1, 144)  # independently verified
    # FFff(A) double not yet complete; corner-only piece would be +1/144
    bare_double_verified = qq + qcc + ccff
    return {
        "target_triple": target_triple,
        "target_double": target_double,
        "mass_ct_triple": mass_ct_triple,
        "mass_ct_double_const": mass_ct_double,
        "fffa_triple": fffa_triple,
        "bare_double_verified_without_ffff": bare_double_verified,
        "triple_gap_after_mass_ct_and_fffa": target_triple
        - mass_ct_triple
        - fffa_triple,
    }


def main():
    avg = taylor_averages()
    print("FFff(A) Taylor averages:", avg)
    poles = poles_from_corner(avg)
    print("triple pole:", poles["triple_pole"])
    print("corner-only double (incomplete):", poles["double_from_corner_linear_only"])
    print("I_corner:", poles["I_corner_series"])
    snap = assembly_snapshot()
    print("assembly snapshot:", snap)
    print("ALL FFff(A) CORNER CHECKS PASSED")


if __name__ == "__main__":
    main()
