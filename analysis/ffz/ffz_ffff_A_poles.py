"""FFff(A) soft-corner check (corrected).

Supersedes the earlier bilinear claim -4/27. The connected-cycle Wick formula
and finite-d numerator are those of FFFF_A_CORNER_CORRECTION.md.

At d=4 the angular averages of H are (-32/3, 0, -8, 0) for
(constant, q, r, qr). There is no short-short triple pole. The unit-cutoff
corner Laurent series is 10/(3ε)-169/3, giving graph double pole +5/(1728ε²)
from this corner alone. Once-soft faces and FFff(B) remain open.

Run:
  python ffz_ffff_A_poles.py

Requires: sympy
"""

import sympy as s


def corner_check():
    q, r, f, g, v, d = s.symbols("q r f g v d")
    P, Q = q * f, r * g
    S = 1 + P + Q
    t = [0, P, P + v, P + v + Q]
    C = s.Matrix(4, 4, lambda i, j: 2 * t[min(i, j)] * (S - t[max(i, j)]) / S)
    a = C[:, 1]
    b = C[:, 3] - C[:, 2]
    c = 2 * P * (S - P) / S
    e = 2 * Q * (S - Q) / S
    j = -2 * P * Q / S
    N = d * (d * (S - P) * (S - Q) + 2 * P * Q) / (4 * S * S)
    T = (
        d * C * N
        + (d * (a * a.T * e + b * b.T * c) + 2 * j * (a * b.T + b * a.T))
        * d
        / (8 * P * Q)
    ).applyfunc(s.cancel)
    A = -(4 - d) * s.trace(T) / 2 - d * (4 - d) * N / 3
    B = s.Rational(2, 3) * (
        q * T[1, 1]
        + r * (T[2, 2] + T[3, 3] + T[2, 3])
        + v * (T[1, 1] + T[2, 2] + T[1, 2])
        + (1 - v) * T[3, 3]
        + N
        * (
            (d - 1) * (q * q * (1 - 2 * f) + r * r * (1 - 2 * g))
            + d * (v * v + (1 - v) ** 2)
        )
    )
    k = 6 - 3 * d / 2
    W = 1 + q + r
    H = S ** (-d / 2) * W ** (1 - k) * (A + (k - 1) * B / W)

    def at(expr):
        return s.factor(expr.subs({q: 0, r: 0}, simultaneous=True))

    def avg(expr):
        x, y = s.symbols("x y")
        return s.factor(
            s.integrate(
                s.integrate(
                    s.integrate(
                        s.expand(expr.subs({f: x * (1 - x), g: y * (1 - y)})),
                        (x, 0, 1),
                    ),
                    (y, 0, 1),
                ),
                (v, 0, 1),
            )
        )

    b_avgs = [
        avg(at(h))
        for h in (
            B.subs(d, 4) / S**2,
            s.diff(B.subs(d, 4) / S**2, q),
            s.diff(B.subs(d, 4) / S**2, r),
            s.diff(B.subs(d, 4) / S**2, q, r),
        )
    ]
    assert b_avgs == [s.Rational(32, 3), 0, 8, 0]

    vals = [avg(at(h)) for h in (H, s.diff(H, q), s.diff(H, r), s.diff(H, q, r))]
    assert [h.subs(d, 4) for h in vals] == [-s.Rational(32, 3), 0, -8, 0]

    ep = s.symbols("eps")
    corner = (
        vals[0] / (ep - 1) ** 2
        + (vals[1] + vals[2]) / (ep * (ep - 1))
        + vals[3] / ep**2
    )
    laurent = s.series(corner.subs(d, 4 - 2 * ep), ep, 0, 1)
    assert s.simplify(laurent.coeff(ep, -1) - s.Rational(10, 3)) == 0

    # Gamma(3ε-1)~-1/(3ε), Wick -1/8, 1/[4d(d-1)]=1/48
    pre = s.Rational(-1, 8) / 48
    graph_double = pre * (-s.Rational(1, 3)) * s.Rational(10, 3)
    assert graph_double == s.Rational(5, 1728)

    return {
        "B_averages": b_avgs,
        "H_averages_d4": [h.subs(d, 4) for h in vals],
        "corner_laurent": laurent,
        "corner_double_pole": graph_double,
        "triple_pole_from_corner": 0,
    }


def main():
    out = corner_check()
    print("B averages (const,q,r,qr):", out["B_averages"])
    print("H averages at d=4:", out["H_averages_d4"])
    print("corner projective Laurent:", out["corner_laurent"])
    print("corner graph double pole:", out["corner_double_pole"])
    print("corner triple pole:", out["triple_pole_from_corner"])
    print("CORRECTED FFff(A) CORNER CHECKS PASSED")


if __name__ == "__main__":
    main()
