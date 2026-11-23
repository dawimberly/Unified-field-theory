#!/usr/bin/env python3
"""Re-extract FFff(A) projective c0 / simple pole with controlled error.

Uses the regulated integrand of FFFF_A_REGULATED.md (long tails analytic,
short-face Taylor subtractions). Does NOT tune to the 1/24 Newton target.

Fit model on y(ε)=ε(I − 2/ε²):
  y(ε) = −10/3 + c0·ε + b·ε² + …
Reports c0 ± error from cut/seed/ε systematics, then the graph simple
constant (pre- and post- ×4 spin normalization).

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_ffff_A_simple_reextract.py
  OPENBLAS_NUM_THREADS=1 python3 ffz_ffff_A_simple_reextract.py --quick

Requires: numpy, scipy, sympy, mpmath
"""

from __future__ import annotations

import argparse
from itertools import product

import mpmath as mp
import numpy as np
from scipy.special import beta, gamma
from scipy.stats import qmc


def H(P, Q, f, g, v, e):
    d = 4 - 2 * e
    k = 3 * e
    q, r = P / f, Q / g
    S, W = 1 + P + Q, 1 + q + r
    t = [0, P, P + v, P + v + Q]

    def C(i, j):
        return 2 * t[min(i, j)] * (S - t[max(i, j)]) / S

    N = d * (d * (1 + Q) * (1 + P) + 2 * P * Q) / (4 * S * S)
    a = [0, S - P, S - t[2], S - t[3]]
    b = [0, -P, -t[2], S - t[3]]

    def T(i, j):
        return d * C(i, j) * N + d / S**3 * (
            d
            * (
                P * (1 + P) * a[i] * a[j]
                + Q * (1 + Q) * b[i] * b[j]
            )
            - 2 * P * Q * (a[i] * b[j] + b[i] * a[j])
        )

    A = -(4 - d) * (T(1, 1) + T(2, 2) + T(3, 3)) / 2 - d * (4 - d) * N / 3
    B = (
        2
        / 3
        * (
            q * T(1, 1)
            + r * (T(2, 2) + T(3, 3) + T(2, 3))
            + v * (T(1, 1) + T(2, 2) + T(1, 2))
            + (1 - v) * T(3, 3)
            + N
            * (
                (d - 1) * (q * q * (1 - 2 * f) + r * r * (1 - 2 * g))
                + d * (v * v + (1 - v) ** 2)
            )
        )
    )
    return S ** (-d / 2) * (
        W ** (1 - k) * A
        + (k - 1)
        * (
            W ** (-k) * B
            - 2 * (d - 1) * N * (q ** (2 - k) + r ** (2 - k)) / 3
        )
    )


def tail(e):
    d = 4 - 2 * e
    k = 3 * e

    def D(i, j):
        a = e - 1 + i
        b = 1 - 2 * e + j
        c = 4 - e
        return gamma(a) * gamma(b) * gamma(c - a - b) / gamma(c)

    return (
        (k - 1)
        * 2
        * (d - 1)
        / 3
        * d
        / 4
        * beta(2 - e, 2 - e)
        * beta(2 * e, 2 * e)
        * (
            d * (D(0, 0) + D(1, 0) + D(0, 1) + D(1, 1))
            + 2 * D(1, 1)
        )
    )


def run_once(e, seed, power, cut=1.0):
    z = qmc.Sobol(5, scramble=True, seed=seed)
    out = []
    h = 1e-5

    def dp(P, Q, f, g, v):
        return H(P + 1e-27j, Q, f, g, v, e).imag / 1e-27

    def dq(P, Q, f, g, v):
        return H(P, Q + 1e-27j, f, g, v, e).imag / 1e-27

    nblocks = max(1, 2 ** (power - 14))
    for _ in range(nblocks):
        x, y, a, b, v = z.random(2**14).T
        P, Q = x / (1 - x), y / (1 - y)
        f, g = a * (1 - a), b * (1 - b)
        X, Y = P < cut, Q < cut
        wp = P ** (e - 2) / (1 - x) ** 2
        wq = Q ** (e - 2) / (1 - y) ** 2
        c = H(0, 0, f, g, v, e)
        cp = dp(0, 0, f, g, v)
        cq = dq(0, 0, f, g, v)

        def regular(p, q):
            d = 4 - 2 * e
            k = 3 * e
            S = 1 + p + q
            N = d * (d * (1 + p) * (1 + q) + 2 * p * q) / (4 * S * S)
            return H(p, q, f, g, v, e) + (k - 1) * 2 * (
                d - 1
            ) / 3 * N * S ** (-d / 2) * (
                (p / f) ** (2 - k) + (q / g) ** (2 - k)
            )

        step = h * g
        dp0 = lambda q: regular(1e-25j, q).imag / 1e-25
        mixed = (-3 * dp0(0) + 4 * dp0(step) - dp0(2 * step)) / (2 * step)
        U = H(P, 0, f, g, v, e)
        V = H(0, Q, f, g, v, e)
        D = dp(0, Q, f, g, v)
        E = dq(P, 0, f, g, v)
        M = (
            H(P, Q, f, g, v, e)
            - X * (V + P * D)
            - Y * (U + Q * E)
            + (X & Y) * (c + P * cp + Q * cq + P * Q * mixed)
        )
        a0 = cut ** (e - 1) / (e - 1)
        a1 = cut**e / e
        out.append(
            np.mean(
                (f * g) ** (1 - e)
                * (
                    wp * wq * M
                    + wq * (V - Y * (c + Q * cq)) * a0
                    + wq * (D - Y * (cp + Q * mixed)) * a1
                    + wp * (U - X * (c + P * cp)) * a0
                    + wp * (E - X * (cq + P * mixed)) * a1
                    + c * a0 * a0
                    + (cp + cq) * a0 * a1
                    + mixed * a1 * a1
                )
            )
        )
    return np.mean(out) + 2 * tail(e)


def fit_c0(eps, ys):
    """Fit y = -10/3 + c0 ε + b ε² to points (ε, y)."""
    eps = np.asarray(eps, float)
    ys = np.asarray(ys, float)
    A = np.column_stack([np.ones_like(eps), eps, eps**2])
    # Fix intercept near -10/3 by using free fit then reporting both
    coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
    intercept, c0, b = coef
    # Also constrained intercept = -10/3
    A2 = np.column_stack([eps, eps**2])
    y2 = ys - (-10.0 / 3.0)
    coef2, *_ = np.linalg.lstsq(A2, y2, rcond=None)
    c0_c, b_c = coef2
    return {
        "intercept_free": intercept,
        "c0_free": c0,
        "b_free": b,
        "c0_constrained": c0_c,
        "b_constrained": b_c,
    }


def graph_simple_from_c0(c0):
    """Pre-×4 graph simple at L=0 from projective c0 (FFFF_A conversion)."""
    return float(c0) / 1152 + 79 / 6912 + np.pi**2 / 768


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--quick",
        action="store_true",
        help="fewer samples (smoke test)",
    )
    args = ap.parse_args()
    power = 15 if args.quick else 17
    eps_list = [0.0075, 0.00375, 0.001875] if args.quick else [
        0.01,
        0.0075,
        0.005,
        0.00375,
        0.0025,
        0.001875,
    ]
    cuts = [0.5, 1.0, 2.0]
    seeds = [12, 41] if args.quick else [12, 41, 77]

    print("=== FFff(A) simple re-extraction (no 1/24 tuning) ===\n")
    print(f"power={power}, eps={eps_list}, cuts={cuts}, seeds={seeds}\n")

    # Collect y(ε)=ε(I-2/ε²) for each (cut, seed)
    records = []  # (cut, seed, eps, y, e2I)
    for cut, seed in product(cuts, seeds):
        ys = []
        for e in eps_list:
            I = run_once(e, seed, power, cut)
            y = e * (I - 2 / (e * e))
            e2I = e * e * I
            records.append((cut, seed, e, y, e2I))
            ys.append(y)
            print(
                f"  cut={cut:g} seed={seed} eps={e:.5f}  "
                f"ε²I={e2I:.6f}  y={y:.6f}",
                flush=True,
            )
        fit = fit_c0(eps_list, ys)
        print(
            f"    -> c0_constrained={fit['c0_constrained']:.6f}  "
            f"c0_free={fit['c0_free']:.6f}  "
            f"intercept_free={fit['intercept_free']:.6f}",
            flush=True,
        )

    # Aggregate c0 estimates: one constrained fit per (cut, seed)
    c0_vals = []
    for cut, seed in product(cuts, seeds):
        ys = [r[3] for r in records if r[0] == cut and r[1] == seed]
        es = [r[2] for r in records if r[0] == cut and r[1] == seed]
        fit = fit_c0(es, ys)
        c0_vals.append(fit["c0_constrained"])

    c0_vals = np.array(c0_vals, float)
    c0_mean = float(np.mean(c0_vals))
    c0_std = float(np.std(c0_vals, ddof=1)) if len(c0_vals) > 1 else 0.0
    # Half-range as a conservative systematics envelope
    c0_half = 0.5 * (float(np.max(c0_vals)) - float(np.min(c0_vals)))
    c0_err = max(c0_std, c0_half)

    # Also fit using only the three smallest ε (closer to continuum)
    c0_small = []
    small = sorted(eps_list)[:3]
    for cut, seed in product(cuts, seeds):
        # Keep (ε, y) paired; do not reorder ys against a differently sorted ε list.
        pairs = sorted(
            [
                (r[2], r[3])
                for r in records
                if r[0] == cut and r[1] == seed and r[2] in small
            ]
        )
        es = [p[0] for p in pairs]
        ys = [p[1] for p in pairs]
        fit = fit_c0(es, ys)
        c0_small.append(fit["c0_constrained"])
    c0_small = np.array(c0_small, float)

    print("\n--- Aggregate ---")
    print(f"c0 (all ε, constrained −10/3): {c0_mean:.6f} ± {c0_err:.6f}")
    print(
        f"  spread over cut×seed: min={c0_vals.min():.6f} "
        f"max={c0_vals.max():.6f}"
    )
    print(
        f"c0 (three smallest ε): {c0_small.mean():.6f} ± "
        f"{max(c0_small.std(ddof=1), 0.5*(c0_small.max()-c0_small.min())):.6f}"
    )
    print(f"packet prior (FFFF_A_REGULATED): ~14.93")

    s_raw = graph_simple_from_c0(c0_mean)
    s_lo = graph_simple_from_c0(c0_mean - c0_err)
    s_hi = graph_simple_from_c0(c0_mean + c0_err)
    print(f"\nA simple L=0 pre-×4:  {s_raw:.8f}  [{s_lo:.8f}, {s_hi:.8f}]")
    print(
        f"A simple L=0 post-×4: {4*s_raw:.8f}  "
        f"[{4*s_lo:.8f}, {4*s_hi:.8f}]"
    )

    # Gap diagnostic vs overlap-only B (informational; not a fit target)
    mp.mp.dps = 30
    C0 = (
        -mp.mpf(15) / 2
        + 2
        * mp.sqrt(3)
        * mp.im(mp.polylog(2, mp.exp(1j * mp.pi / 3)))
    )
    ct = (C0 + 6) / 18 + (1 - mp.zeta(2)) / 24
    import sympy as sy

    e = sy.symbols("e")
    d = 4 - 2 * e
    pref = sy.gamma(3 * e - 1) * sy.exp(3 * sy.EulerGamma * e) / (
        64 * d * (d - 1)
    )
    B = float(
        sy.simplify(
            sy.series(
                pref * (-48 / e**2 + 208 / e + 192 + sy.Rational(64, 3)),
                e,
                0,
                1,
            ).coeff(e, -1)
        )
    )
    base = (
        mp.mpf("-0.05121964552")
        + mp.mpf("-0.07466151326053")
        + mp.mpf("-0.02480")
        + ct
    )
    projected = float(base) + 4 * s_raw + B
    target = 1.0 / 24.0
    print(f"\nInformational (overlap-only B, this A):")
    print(f"  projected simple = {projected:.8f}")
    print(f"  target 1/24      = {target:.8f}")
    print(f"  gap              = {projected - target:+.8f}")
    print(
        f"  c0 that would null the gap ≈ 12.33 "
        f"(distance from this c0: {c0_mean - 12.33:+.2f})"
    )
    print(
        "\nNOTE: gap diagnostic is not used to adjust c0. "
        "Reported c0 is from the regulated integrand alone."
    )
    print("=== RE-EXTRACTION COMPLETE ===")


if __name__ == "__main__":
    main()
