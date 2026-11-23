"""Independent review of the FFZ QCc, CcFf, and FFff(A) poles.

Re-derives the short-pair overlaps and the long-pair endpoints at d=4-2*epsilon
before the Laurent step. Does not import the checkpoint scripts.

Run from this directory:
  python ffz_independent_pole_review.py

Requires Python 3, NumPy, SciPy, SymPy, and mpmath.
The CcFf simple-pole block is a quasi-Monte Carlo estimate of an integrand
built in this file. Set OPENBLAS_NUM_THREADS=1 if BLAS oversubscription appears.
"""

from math import exp
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.special import beta, gamma
from scipy.stats import qmc

mp.mp.dps = 40


def section(title):
    print("\n==", title)


# ---------------------------------------------------------------------------
# QCc
# ---------------------------------------------------------------------------

def qcc_corner_and_overlap():
    """Short-pair overlap of the doubled QCc graph, before setting epsilon to 0."""
    e, d = sp.symbols("e d")
    # Homogeneous numerator on the gauge with the single edge fixed at length 1.
    # At p=q=0 every non-root covariance vanishes, W=1, and only the curvature-mass
    # term and the single-edge piece of B survive.
    A = -d * (4 - d) / 3
    B = 2 * d / 3
    # k+1 = 3*epsilon for five scalar edges.
    N = sp.simplify(A + 3 * e * B)
    d_sub = 4 - 2 * e
    N_eps = sp.factor(N.subs(d, d_sub))
    assert sp.simplify(N.subs(d, d_sub) - 4 * d_sub * e / 3) == 0
    # The superseded checkpoint used 2*d*epsilon/3. The difference of the two
    # constant corners, integrated on p,q < 1 against p^{e-1} q^{e-1}, is exact.
    c_new = 4 * d_sub * e / 3
    c_old = 2 * d_sub * e / 3
    diff = sp.expand(c_new / e**2 - c_old / e**2)
    assert sp.simplify(diff - (sp.Rational(8, 3) / e - sp.Rational(4, 3))) == 0
    # Previous projective residue 16/3 plus the added 8/3.
    residue = sp.Rational(16, 3) + sp.Rational(8, 3)
    assert residue == 8
    # Gamma(3e)~1/(3e), Wick +1/4, and 1/[4 d (d-1)] combine into 1/192 at d=4.
    double = residue / (3 * 4 * 48)
    assert double == sp.Rational(1, 72)
    return {
        "corner": N_eps,
        "overlap_difference": diff,
        "projective_residue": residue,
        "double_pole": double,
    }


def qcc_faces_and_tail():
    """Long-pair tail and the two face integrals that fix the simple pole."""
    # Geometric Dirichlet integral, W^{-3e} omitted. The beta-function form is
    # re-derived in the review text. Pole and finite part only; O(e) does not
    # enter the simple pole.
    def g0_fun(z):
        zz = mp.mpf(z)
        dd = 4 - 2 * zz
        pref_m = -dd * (6 - dd)
        br = (
            2 * mp.gamma(1 + zz) * mp.gamma(zz) * mp.gamma(2 - 3 * zz)
            + 2 * mp.gamma(1 + zz) ** 2 * mp.gamma(1 - 3 * zz)
        ) / mp.gamma(3 - zz)
        return pref_m * mp.beta(1 - zz, 1 - zz) ** 2 * br

    # e*G0 -> -8, and G0 + 8/e -> -32.
    z = mp.mpf("1e-10")
    pole = z * g0_fun(z)
    finite = g0_fun(z) + 8 / z
    assert abs(pole + 8) < mp.mpf("1e-8")
    assert abs(finite + 32) < mp.mpf("1e-6")

    def tail_fun(z):
        zz = mp.mpf(z)
        dd = 4 - 2 * zz
        return (
            (8 * dd * zz / 3)
            * mp.beta(1 - 2 * zz, 1 + zz)
            * mp.beta(2 * zz, 2 * zz)
            * mp.beta(zz, 1)
            * mp.beta(1 - zz, 1 - zz)
        )

    tail_pole = z * tail_fun(z)
    tail_finite = tail_fun(z) - mp.mpf(32) / (3 * z)
    assert abs(tail_pole - mp.mpf(32) / 3) < mp.mpf("1e-7")
    assert abs(tail_finite - mp.mpf(80) / 3) < mp.mpf("1e-5")

    # Evanescent once-subtracted face, divided by d*epsilon. Closed form checked
    # against the rational integrand, still at finite r, before d->4.
    q, r = sp.symbols("q r", positive=True)
    W = 1 + q
    S = 1 + q * r
    C = 2 * q * r / S
    bracket = -2 * W / 3 + 2 * (W * C + 1 + q**2 * (1 - 2 * r)) / W - sp.Rational(4, 3) * q
    Y = sp.together(bracket / S**2)
    integrand_short = sp.together((Y - sp.Rational(4, 3)) / q)
    I_face = sp.simplify(
        sp.integrate(integrand_short, (q, 0, 1)) + sp.integrate(sp.together(Y / q), (q, 1, sp.oo))
    )
    claimed = sp.Rational(2, 3) - 4 * (r + 2) * sp.log(r) / (3 * (r - 1))
    assert sp.simplify(I_face - claimed) == 0

    def geom(rr):
        return -mp.log(rr) / (1 - rr) ** 2 - 1 / (1 - rr)

    def evan(rr):
        return mp.mpf(2) / 3 - 4 * (rr + 2) * mp.log(rr) / (3 * (rr - 1))

    g1 = 24 * mp.quad(lambda x: geom(x * (1 - x)), [0, mp.mpf("0.5"), 1])
    v = 8 * mp.quad(lambda x: evan(x * (1 - x)), [0, mp.mpf("0.5"), 1])
    # Old overlap finite piece -4/3, corrected piece -8/3.
    I0_old = -32 + g1 + v + mp.mpf(80) / 3 - mp.mpf(4) / 3
    I0_new = -32 + g1 + v + mp.mpf(80) / 3 - mp.mpf(8) / 3
    # 1/[d(d-1)] = (1/12)(1 + (7/6) e + ...); Gamma(3e) e^{3 gamma e} = 1/(3e) + O(e).
    def simple(I0, residue):
        return (I0 / 3 + residue * mp.mpf(7) / 18) / 192

    c_new = simple(I0_new, 8)
    # z^{-3e} = 1 - 3 e L + ... hits 1/(72 e^2) and produces -L/(24 e).
    return {
        "I0_old": I0_old,
        "I0_new": I0_new,
        "simple_pole_L0": c_new,
        "log_coefficient": sp.Rational(-1, 24),
    }


# ---------------------------------------------------------------------------
# CcFf
# ---------------------------------------------------------------------------

def ccff_corner():
    """Laplacian derivation of the p=0 Taylor coefficients, d kept symbolic."""
    d, q, g, v = sp.symbols("d q g v")
    r = g * (1 - g)
    Q = q * r
    # Rooted triangle: scalar edges v and 1-v, fermion effective edge Q.
    cA, cB, cF = 1 / (1 - v), 1 / v, 1 / Q
    L = sp.Matrix([[cA + cF, -cF], [-cF, cB + cF]])
    C = sp.simplify(2 * L.inv())
    cAA, cAB, cBB = C[0, 0], C[0, 1], C[1, 1]
    cdel = sp.factor(cAA + cBB - 2 * cAB)
    vA, vB = sp.factor(cAA - cAB), sp.factor(cAB - cBB)
    Maa = sp.factor(d * cAA + 2 * vA**2 / cdel)
    Mbb = sp.factor(d * cBB + 2 * vB**2 / cdel)
    Mab = sp.factor(d * cAB + 2 * vA * vB / cdel)
    # q N = d c / (4 Q) = d / (2 S), S=1+Q. Then one more S^{-d/2}.
    lam = (4 - d) / 2
    Apoly = -lam * (Maa + Mbb) - d * (4 - d) / 3
    B = (sp.Integer(2) / 3) * (
        (1 - v) * (Maa + d * (1 - v))
        + v * (Mbb + d * v)
        + q * (Maa + Mbb + Mab)
        + (d - 1) * q**2 * (1 - 2 * r)
    )
    k = 6 - sp.Rational(3, 2) * d
    S, W = 1 + Q, 1 + q
    raw = (d / 2) * (Apoly + k * B / W)
    H = raw * S ** (-1 - d / 2) * W ** (-k)
    H0 = sp.factor(sp.limit(sp.simplify(H), q, 0))
    H1 = sp.factor(sp.limit(sp.diff(sp.simplify(H), q), q, 0))

    def average(expr):
        return sp.factor(sp.integrate(sp.integrate(sp.expand(expr), (v, 0, 1)), (g, 0, 1)))

    a0, a1 = average(H0), average(H1)
    claimed0 = -d**2 * (d - 4) / 6
    claimed1 = -d * (d - 4) * (17 * d**2 - 72 * d + 4) / 72
    assert sp.factor(a0 - claimed0) == 0
    assert sp.simplify(a1 - claimed1) == 0
    ee = sp.symbols("ee")
    h1 = sp.series(a1.subs(d, 4 - 2 * ee), ee, 0, 2).removeO()
    # Leading piece -4*epsilon/3. Integrated with p^{e-1} q^{e-1} this is -4/(3 epsilon).
    assert sp.simplify(h1 + sp.Rational(4, 3) * ee) == 0 or sp.degree(sp.expand(h1 + sp.Rational(4, 3) * ee)) >= 2
    leading = sp.series(a1.subs(d, 4 - 2 * ee), ee, 0, 2).coeff(ee, 1)
    assert leading == sp.Rational(-4, 3)
    return {"H0": a0, "H1": a1, "H1_linear": leading}


def ccff_tail_and_geometric():
    e = sp.symbols("e", positive=True)
    d = 4 - 2 * e
    G = sp.gamma
    # Regularized form of (2 d^2 e / 3) * four beta values. The Q integral is the
    # continuation of Beta(e-1, 4-2e), taken before e->0.
    tail_reg = (
        -2 * d * d / 3
        * G(1 + 2 * e) ** 2 / G(1 + 4 * e)
        * G(1 + e) * G(4 - 2 * e) / ((1 - e) * G(3 - e))
        * G(2 - e) ** 2 / G(4 - 2 * e)
        / ((1 - 2 * e) * (2 - 2 * e))
    )
    # This regularized expression equals e * (model integral) near e=0? The
    # checkpoint's shifted product equals the model times e, up to the sign
    # convention checked by the value -8/3.
    assert sp.simplify(tail_reg.subs(e, 0)) == sp.Rational(-8, 3)
    deriv = sp.simplify(sp.diff(tail_reg, e).subs(e, 0))
    assert deriv == sp.Rational(-20, 3)

    P, Q = sp.symbols("P Q", positive=True)
    S = 1 + P + Q
    M = 8 * P * (1 + Q) / S + 4 * P**2 * Q / (S * (1 + P))
    H = -2 * (1 + P) * M / S**3
    assert sp.simplify(H.subs(Q, 0) + 16 * P / (1 + P) ** 3) == 0
    linear = sp.factor(sp.diff(H, Q).subs(Q, 0) / 6)
    assert sp.integrate(sp.cancel(linear / P), (P, 0, sp.oo)) == 2

    # Residue: geometric +2, overlap -4/3, tail -8/3.
    residue = 2 - sp.Rational(4, 3) - sp.Rational(8, 3)
    assert residue == -2
    # Weight -1/2, Gamma(3e)~1/(3e), 1/[4 d (d-1)] = 1/48 at d=4.
    double = residue * (sp.Rational(-1, 2)) / (3 * 48)
    assert double == sp.Rational(1, 144)
    # z^{-3e} on 1/(144 e^2) produces -L/(48 e).
    return {"residue": residue, "double_pole": double, "log_coefficient": sp.Rational(-1, 48)}


def H_ccff(P, Q, f, g, v, e):
    """CcFf projective density. f,g are split-fraction products a(1-a).

    The long-scalar-pair model is removed here and restored analytically by the caller.
    """
    d = 4 - 2 * e
    k = 3 * e
    p = P / np.maximum(f, 1e-300)
    q = Q / np.maximum(g, 1e-300)
    S = 1 + P + Q
    W = 1 + p + q
    ts = [np.zeros_like(P), P, P + (1 - v), P + (1 - v) + Q]
    ls = [np.zeros_like(P), P, P + (1 - v), -v]

    def C(i, j):
        ti, tj = ts[i], ts[j]
        mn = np.minimum(ti, tj)
        mx = np.maximum(ti, tj)
        return 2 * mn * (S - mx) / S

    coef = 4 * Q / (S * (1 + P))

    def M(i, j):
        return d * C(i, j) + coef * ls[i] * ls[j]

    lam_c = (6 - d) / 2
    lam_f = (4 - d) / 2
    A = -lam_c * M(1, 1) - lam_f * (M(2, 2) + M(3, 3)) - d * (4 - d) / 3 * (1 + p)
    B = (2 / 3) * (
        p * M(1, 1)
        + d * p * p * (1 - 2 * f)
        + (1 - v) * (M(1, 1) + M(2, 2) + M(1, 2) + d * (1 - v))
        + v * (M(3, 3) + d * v)
        + q * (M(2, 2) + M(3, 3) + M(2, 3))
        + (d - 1) * q * q * (1 - 2 * g)
    )
    qN = d * (1 + P) / (2 * S)
    bare = qN * S ** (-0.5 * d) * W ** (-k) * (A + k * B / W)
    model = qN * S ** (-0.5 * d) * (4 * d * e / 3) * p ** (1 - 3 * e)
    return bare - model


def ccff_simple_pole(power=17, seeds=(11, 22)):
    """QMC of the independently coded CcFf integrand.

    Returns epsilon*I (projective) and the graph simple-pole estimator at L=0.
    """

    def run(e, seed):
        n = 2**power
        x, y, a, b, v = qmc.Sobol(5, scramble=True, seed=seed).random(n).T
        P, Q = x / (1 - x), y / (1 - y)
        f, g = a * (1 - a), b * (1 - b)
        short_p, short_q = P < 1, Q < 1
        wp = P ** (e - 1) / (1 - x) ** 2 * f ** (-e)
        wq = Q ** (e - 2) / (1 - y) ** 2 * g ** (1 - e)
        step = 1e-7

        def faces(Pv):
            C = H_ccff(Pv, np.zeros_like(Pv), f, g, v, e)
            D = (H_ccff(Pv, np.full_like(Pv, step), f, g, v, e) - C) / step
            return C, D

        C, D = faces(P)
        C0, D0 = faces(np.zeros_like(P))
        F = H_ccff(P, Q, f, g, v, e) - np.where(short_q, C + Q * D, 0.0)
        F0 = H_ccff(np.zeros_like(P), Q, f, g, v, e) - np.where(short_q, C0 + Q * D0, 0.0)
        mixed = np.mean(wp * wq * (F - np.where(short_p, F0, 0.0)))
        face_q = np.mean(wq * F0)
        c = np.mean(wp * g ** (1 - e) * (C - np.where(short_p, C0, 0.0)))
        dlt = np.mean(wp * g ** (1 - e) * (D - np.where(short_p, D0, 0.0)))
        c0 = np.mean(g ** (1 - e) * C0)
        d0 = np.mean(g ** (1 - e) * D0)
        fac = beta(1 - e, 1 - e) / e
        dd = 4 - 2 * e
        tail = (2 * dd * dd * e / 3) * beta(2 * e, 2 * e) * beta(2 - e, 2 - e) * beta(e - 1, 4 - 2 * e) * beta(1 - 2 * e, 2)
        I = mixed + fac * face_q + (c + fac * c0) / (e - 1) + (dlt + fac * d0) / e + tail
        graph = -I * gamma(3 * e) * exp(3 * 0.5772156649015329 * e) / (8 * dd * (dd - 1))
        return e * I, e * (graph - 1 / (144 * e * e))

    rows = []
    for e in (0.01, 0.005, 0.0025):
        vals = [run(e, seed) for seed in seeds]
        rows.append(np.mean(vals, axis=0))
        print("  epsilon", e, "mean of streams", rows[-1])
    y = np.array(rows)
    # Remove O(e) and O(e^2): samples at 4h, 2h, h with h=0.0025.
    extrap = (8 * y[2] - 6 * y[1] + y[0]) / 3
    return {"samples": rows, "extrapolated": extrap}


# ---------------------------------------------------------------------------
# FFff(A)
# ---------------------------------------------------------------------------

def ffff_A_structure():
    """Curvature insertion for the class-A graph: two parallel fermion pairs.

    Edges, with vertex 0 fixed at the origin:
      fermion pair 0-2 of total q and effective length Q,
      fermion pair 1-3 of total r and effective length R,
      scalar 0-3 of length v, scalar 1-2 of length w.
    Class A is the multiplicity-8 matching. Its flat Wick weight is -1/8.
    """
    q, r, g, h, v, w = sp.symbols("q r g h v w", positive=True)
    Q, R = q * g * (1 - g), r * h * (1 - h)
    fa, fb, fc, fe = q * g, q * (1 - g), r * h, r * (1 - h)
    idx = {1: 0, 2: 1, 3: 2}
    L = sp.zeros(3)
    for i, j, length in (
        (0, 2, fa), (0, 2, fb), (1, 3, fc), (1, 3, fe), (0, 3, v), (1, 2, w)
    ):
        cond = 1 / length
        if i and j:
            a, b = idx[i], idx[j]
            L[a, a] += cond
            L[b, b] += cond
            L[a, b] -= cond
            L[b, a] -= cond
        elif i == 0:
            L[idx[j], idx[j]] += cond
        else:
            L[idx[i], idx[i]] += cond
    U = sp.factor(fa * fb * fc * fe * v * w * L.det())
    assert sp.factor(U - q * r * (v + w + Q + R)) == 0

    # Each spin weight scales as rho^{-1}, so the pair product scales as rho^{-2}.
    # With six edges the radial degrees are alpha_A = 5 - 3d/2 and alpha_B = 6 - 3d/2.
    # At d=4-2e these are 3e-1 and 3e, and Gamma(3e)=(3e-1) Gamma(3e-1), so the
    # gauge-fixed combination is [A + (3e-1) B] Gamma(3e-1) z^{1-3e}.
    e, d = sp.symbols("e d")
    alpha_A = sp.simplify((5 - sp.Rational(3, 2) * (4 - 2 * e)))
    alpha_B = sp.simplify((6 - sp.Rational(3, 2) * (4 - 2 * e)))
    assert alpha_A == 3 * e - 1
    assert alpha_B == 3 * e

    # Corner: both effective lengths vanish and v+w=1. Each qN -> d/2.
    # The surviving mass profile is the two scalar edges between the merged vertices.
    # C = 2 v w, M = d C on the non-root supervertex, fermion edges drop out.
    vv = sp.symbols("vv")
    Cslow = 2 * vv * (1 - vv)
    Mslow = d * Cslow
    Bcorner = sp.simplify(
        (sp.Integer(2) / 3)
        * (vv * (Mslow + d * vv) + (1 - vv) * (Mslow + d * (1 - vv)))
    )
    assert sp.simplify(Bcorner - 2 * d / 3) == 0
    # Constant term of the projective measure is q^{-2+e} r^{-2+e}. Its integral
    # on the unit square is 1/(e-1)^2, finite at e=0. It does not make a triple
    # or double projective pole. Those sit in the linear and bilinear Taylor terms.
    return {
        "U": "q*r*(v+w+Q+R)",
        "radial": "Gamma(3*e-1) * z**(1-3*e) * [A + (3*e-1)*B]",
        "weight": sp.Rational(-1, 8),
        "corner_B": sp.Rational(2, 3) * d,
        "corner_spin": (d / 2) ** 2,
    }


def main():
    section("QCc overlap, before d->4")
    qcc = qcc_corner_and_overlap()
    for key, value in qcc.items():
        print(f"  {key}: {value}")

    section("QCc faces, tail, simple pole")
    faces = qcc_faces_and_tail()
    print("  I0 corrected", faces["I0_new"])
    print("  simple pole at L=0", faces["simple_pole_L0"])
    print("  log coefficient", faces["log_coefficient"])
    assert abs(faces["simple_pole_L0"] + mp.mpf("0.01649742330")) < mp.mpf("1e-11")

    section("CcFf short-pair corner")
    corner = ccff_corner()
    print("  averaged H0", corner["H0"])
    print("  averaged H1", corner["H1"])
    print("  H1 coefficient of epsilon", corner["H1_linear"])

    section("CcFf long pair, geometric face, double pole")
    tail = ccff_tail_and_geometric()
    for key, value in tail.items():
        print(f"  {key}: {value}")

    section("CcFf simple pole, this file's QMC")
    num = ccff_simple_pole()
    proj, simple = num["extrapolated"]
    print("  extrapolated epsilon*I", proj)
    print("  extrapolated simple pole", simple)
    assert abs(proj + 2) < 0.01
    assert abs(simple + 0.0247) < 0.0015

    section("FFff(A) insertion")
    info = ffff_A_structure()
    for key, value in info.items():
        print(f"  {key}: {value}")
    print("ALL INDEPENDENT CHECKS PASSED")


if __name__ == "__main__":
    main()
