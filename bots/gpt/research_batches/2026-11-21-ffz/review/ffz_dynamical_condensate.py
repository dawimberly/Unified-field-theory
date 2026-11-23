#!/usr/bin/env python3
"""Dynamical condensate / BCS–NJL gap equation for <X> in the WZ model.

New-idea rung (stacked on the ξv² contact analysis): v should come from a
gap equation, not be chosen by hand. Compare any predicted v²/z to the
gap-absorbing target 3.703e-4 h² from FFFF_CONDENSATE_R.md.

Truncations and SUSY-breaking assumptions are stated explicitly.
No tuning to 1/24.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_dynamical_condensate.py

Requires: numpy, sympy, mpmath
"""

from __future__ import annotations

from fractions import Fraction as F

import mpmath as mp
import numpy as np
import sympy as s

mp.mp.dps = 40

# Gap-absorbing target from tree ξv² contact (PR #16)
TARGET_V2_OVER_Z_PREF = 3.703108e-4  # v²/z = PREF × h²
GAP_SIMPLE = 0.00974619
K_NUM = float(16 * mp.pi**2)  # ≈ 157.914


# ---------------------------------------------------------------------------
# (1) Gap-equation setup
# ---------------------------------------------------------------------------

def model_and_truncation():
    """WZ + truncations for a dynamical <X>.

    Superpotential (packet): W = m X Y + y X³/6.
    F-term potential: V_F = |m Y + (y/2) X²|² + |m X|².
    Quartic after eliminating auxiliaries: |y|²|X|⁴/4 (+ mixed terms).

    Truncations used below:
      (T1) Holomorphic / SUSY F-term stationarity (no soft terms).
      (T2) Soft scalar mass m_s²|X|² added by hand; minimize tree V_F+V_soft
           along Y=0 (or Y integrated at tree). No CW loops in the
           minimization (tree soft-broken gap).
      (T3) NJL/BCS four-fermion truncation: integrate out a heavy messenger
           or freeze the radial mode to get G (ψψ)(ψ̄ψ̄) with cutoff Λ;
           gap equation for fermion bilinear Δ, then v ~ Δ/y.
           ASSUMED: cutoff regularization; large-N / mean-field.

    SUSY breaking:
      T1: none (exact SUSY).
      T2: soft mass m_s² — ASSUMED external.
      T3: cutoff Λ breaks SUSY in the UV completion — ASSUMED.
    """
    return {
        "W": "m X Y + y X³/6",
        "h": "|y|²/K",
        "z": "|m|²",
        "truncations": {
            "T1": "SUSY F-term stationarity (exact)",
            "T2": "tree V_F + soft m_s²|X|²; Y=0 slice",
            "T3": "mean-field NJL/BCS with cutoff Λ; G ~ y²/z",
        },
        "labels": {
            "T1": "derived",
            "T2": "derived (min); assumed (soft m_s² present)",
            "T3": "assumed (NJL truncation + cutoff)",
        },
    }


def T1_SUSY_F_term_gap():
    """Exact SUSY: F_X = F_Y = 0.

    F_Y = m X = 0  ⇒  X = 0 (m≠0).
    Non-renormalization: W not corrected ⇒ F-term equations uncorrected
    ⇒ no dynamical shift of the SUSY vacuum by loops.
    Only solution: v = 0.
    """
    m, y, X, Y = s.symbols("m y X Y")
    F_X = m * Y + y * X**2 / 2
    F_Y = m * X
    assert F_X.subs({X: 0, Y: 0}) == 0
    # Unique solution when m≠0
    sol_Y_from_FY = s.solve(F_Y, X)
    assert sol_Y_from_FY == [0]
    return {
        "nontrivial_solution": False,
        "v": 0,
        "v2_over_z": 0.0,
        "reason": (
            "F_Y=mX=0 forces X=0; non-renormalization of W forbids loop-generated "
            "shift of the F-flat locus. SUSY forbids a dynamical condensate."
        ),
        "critical_coupling": None,
        "label": "derived",
    }


def T2_soft_tree_gap(r_soft):
    """Tree soft-broken potential on Y=0:
        V = (|m|² + m_s²)|X|² + (|y|²/4)|X|⁴
    with r = m_s²/z, z=|m|², |y|² = h K.

    Stationarity: |m|²+m_s² + (|y|²/2)|X|² = 0
      ⇒ v² ≡ |X|² = -2(|m|²+m_s²)/|y|² = -2 z(1+r)/(h K)
      ⇒ v²/z = -2(1+r)/(h K)

    Nontrivial only for r < -1 (tachyonic soft mass overcoming |m|²).
    Soft ratio r is an INPUT — not predicted by the Newton test.
    """
    # Symbolic
    h, r = s.symbols("h r", real=True)
    v2_over_z = s.simplify(-2 * (1 + r) / (h * (16 * s.pi**2)))
    # Numeric for given r_soft
    if r_soft >= -1:
        return {
            "r": r_soft,
            "nontrivial_solution": False,
            "v2_over_z": 0.0,
            "reason": f"r=m_s²/z={r_soft} ≥ -1 ⇒ X=0 is the minimum",
            "formula": "v²/z = -2(1+r)/(h K)   [only for r<-1]",
            "label": "derived (min); assumed (soft r input)",
        }
    # Return symbolic coefficient: v²/z = C(r)/h with C = -2(1+r)/K
    C = -2 * (1 + r_soft) / K_NUM
    return {
        "r": r_soft,
        "nontrivial_solution": True,
        "v2_over_z_form": f"{C:.6e} / h",
        "C_over_h": C,  # v²/z = C/h
        "formula": "v²/z = -2(1+r)/(h K)",
        "soft_input_required": True,
        "label": "derived (min); assumed (soft r input)",
    }


def T3_NJL_gap_equation():
    """Mean-field NJL/BCS with hard cutoff Λ (4d).

    Four-fermion coupling from Yukawa at scale z: G = κ y² / z = κ h K / z
    with κ = O(1) (ASSUMED κ=1).

    Standard 4d NJL gap for Dirac mass Δ (Weyl×2 ~ Dirac):
        1 = (G/(4π²)) [Λ² - Δ² ln((Λ²+Δ²)/Δ²)]
    Critical: G_c = 4π²/Λ²  ⇒  nontrivial iff G > G_c
        i.e. κ h K / z > 4π²/Λ²  ⇒  (Λ/√z)² > 4π²/(κ h K)

    Near criticality Δ ≪ Λ:
        Δ²/Λ² ≈ (G-G_c)/G * (something) — expand:
        1 = (G Λ²)/(4π²) [1 - (Δ²/Λ²) ln(Λ²/Δ²) + …]
    Let g = G Λ²/(4π²); nontrivial for g>1:
        1/g = 1 - x ln(1/x) roughly with x=Δ²/Λ²  (leading)
    More precisely from 1 = g (1 - x ln(1+1/x)) with x=Δ²/Λ².

    Then v ~ Δ/|y|  (Yukawa converts bilinear to scalar vev; ASSUMED)
      ⇒ v²/z = (Δ/|y|)² / z = (Δ²/z) / (h K)

    The cutoff ratio λ = Λ/√z remains a FREE UV parameter.
    """
    # Critical coupling
    # G = h K / z (κ=1), G_c = 4π²/Λ²
    # g = G Λ²/(4π²) = h K λ² / (4π²) with λ=Λ/√z
    # Nontrivial iff g > 1 iff λ² > 4π²/(h K)

    def g_of(h, lam):
        return h * K_NUM * lam**2 / (4 * np.pi**2)

    def solve_x(g, n_grid=5000):
        """Solve 1 = g (1 - x ln(1+1/x)) for x=Δ²/Λ² ∈ (0,1)."""
        if g <= 1:
            return 0.0
        # f(x) = g(1 - x ln(1+1/x)) - 1
        xs = np.linspace(1e-12, 0.999, n_grid)
        vals = g * (1 - xs * np.log(1 + 1 / xs)) - 1
        # find root (vals starts positive at tiny x if g>1, crosses to neg)
        for i in range(len(xs) - 1):
            if vals[i] == 0:
                return float(xs[i])
            if vals[i] * vals[i + 1] < 0:
                # linear interpolate
                t = vals[i] / (vals[i] - vals[i + 1])
                return float(xs[i] + t * (xs[i + 1] - xs[i]))
        return float(xs[np.argmin(np.abs(vals))])

    # Scan: for given h, minimal λ for nontrivial, and sample solutions
    samples = []
    for h in (0.1, 1.0):
        lam_crit = np.sqrt(4 * np.pi**2 / (h * K_NUM))
        for lam in (lam_crit * 1.01, lam_crit * 1.5, lam_crit * 3.0, 10.0):
            g = g_of(h, lam)
            x = solve_x(g)
            # Δ² = x Λ² = x λ² z  ⇒  Δ²/z = x λ²
            # v²/z = (Δ²/z)/(h K) = x λ² /(h K)
            v2z = x * lam**2 / (h * K_NUM) if x > 0 else 0.0
            samples.append(
                {
                    "h": h,
                    "lam": lam,
                    "lam_crit": lam_crit,
                    "g": g,
                    "x_Delta2_over_Lam2": x,
                    "v2_over_z": v2z,
                }
            )
    return {
        "G": "h K / z  (κ=1 assumed)",
        "G_c": "4π²/Λ²",
        "gap_eq": "1 = (G/(4π²))[Λ² - Δ² ln((Λ²+Δ²)/Δ²)]",
        "v_from_Delta": "v = Δ/|y|  (assumed)",
        "nontrivial_requires": "Λ/√z > 2π/√(h K)  and soft/cutoff SUSY breaking",
        "samples": samples,
        "free_parameter": "Λ/√z (and κ)",
        "label": "assumed (NJL+cutoff+κ=1+v=Δ/y)",
    }


# ---------------------------------------------------------------------------
# (2)–(3) Existence, prediction vs target
# ---------------------------------------------------------------------------

def compare_to_target(v2_over_z, h):
    """predicted / target with target = PREF × h²."""
    target = TARGET_V2_OVER_Z_PREF * h**2
    if target == 0:
        return float("nan")
    if v2_over_z == 0:
        return 0.0
    return v2_over_z / target


def soft_r_that_hits_target(h):
    """Invert T2: v²/z = -2(1+r)/(h K) = PREF h²
    ⇒ 1+r = -PREF h² * h K / 2 = -PREF h³ K / 2
    ⇒ r = -1 - PREF h³ K / 2
    """
    r = -1 - TARGET_V2_OVER_Z_PREF * (h**3) * K_NUM / 2
    return r


def main():
    print("=== Dynamical condensate / BCS–NJL gap equation ===\n")
    print("No tuning to 1/24. Target v²/z = 3.703e-4 × h² (from ξv² contact).\n")

    info = model_and_truncation()
    print("-- (1) Setup --")
    print(f"  W = {info['W']};  h={info['h']};  z={info['z']}")
    for k, v in info["truncations"].items():
        print(f"  {k}: {v}  [{info['labels'][k]}]")

    # T1
    print("\n-- (2a) T1: exact SUSY F-term --")
    t1 = T1_SUSY_F_term_gap()
    print(f"  Nontrivial solution? {t1['nontrivial_solution']}")
    print(f"  v²/z = {t1['v2_over_z']}")
    print(f"  {t1['reason']}")
    print(f"  [{t1['label']}]")
    print(f"  Ratio predicted/target at h=1: {compare_to_target(0.0, 1.0):.3f}")

    # T2
    print("\n-- (2b) T2: soft tree gap --")
    print("  formula: v²/z = -2(1+r)/(h K),  r=m_s²/z, nontrivial iff r<-1")
    for r in (0.0, -0.5, -1.0, -1.5, -2.0):
        t2 = T2_soft_tree_gap(r)
        if not t2["nontrivial_solution"]:
            print(f"  r={r:+.1f}: v=0  ({t2['reason']})")
        else:
            for h in (0.1, 1.0):
                v2z = t2["C_over_h"] / h
                ratio = compare_to_target(v2z, h)
                print(
                    f"  r={r:+.1f}, h={h}: v²/z={v2z:.6e}, "
                    f"target={TARGET_V2_OVER_Z_PREF*h**2:.6e}, "
                    f"ratio={ratio:.4e}"
                )
    print("  Soft r required to HIT target:")
    for h in (0.1, 1.0):
        r_star = soft_r_that_hits_target(h)
        print(f"    h={h}: r_* = m_s²/z = {r_star:.6e}  (tuned input)")
    print("  [derived min; assumed soft r — FREE PARAMETER remains]")

    # T3
    print("\n-- (2c) T3: NJL/BCS with cutoff --")
    t3 = T3_NJL_gap_equation()
    print(f"  G = {t3['G']}")
    print(f"  Gap eq: {t3['gap_eq']}")
    print(f"  v from Δ: {t3['v_from_Delta']}")
    print(f"  Nontrivial requires: {t3['nontrivial_requires']}")
    print(f"  Free parameter: {t3['free_parameter']}")
    print(f"  {'h':>4}  {'Λ/√z':>10}  {'g':>8}  {'v²/z':>12}  {'target':>12}  {'ratio':>10}")
    for srow in t3["samples"]:
        h = srow["h"]
        tgt = TARGET_V2_OVER_Z_PREF * h**2
        ratio = compare_to_target(srow["v2_over_z"], h)
        print(
            f"  {h:4.1f}  {srow['lam']:10.4f}  {srow['g']:8.3f}  "
            f"{srow['v2_over_z']:12.4e}  {tgt:12.4e}  {ratio:10.3e}"
        )
    print(f"  [{t3['label']}]")

    # Summary ratios at "natural" points
    print("\n-- (3) Predicted / target summary --")
    print("  T1 (SUSY):           ratio = 0  (v=0; SUSY forbids condensate)")
    print("  T2 (soft):           ratio fixed only after choosing r=m_s²/z — not predicted")
    print("  T3 (NJL):            ratio fixed only after choosing Λ/√z — not predicted")
    # Example: pick T2 r=-2, h=1
    t2_ex = T2_soft_tree_gap(-2.0)
    v2z_ex = t2_ex["C_over_h"] / 1.0
    print(
        f"  Example T2 r=-2, h=1: v²/z={v2z_ex:.6e}, "
        f"ratio={compare_to_target(v2z_ex, 1.0):.4e}"
    )
    # Example T3: first sample with h=1 and nontrivial
    for srow in t3["samples"]:
        if srow["h"] == 1.0 and srow["v2_over_z"] > 0:
            print(
                f"  Example T3 h=1, Λ/√z={srow['lam']:.3f}: "
                f"v²/z={srow['v2_over_z']:.6e}, "
                f"ratio={compare_to_target(srow['v2_over_z'], 1.0):.4e}"
            )
            break

    print(
        f"""
Verdict
-------
(1) Gap equations set up: T1 SUSY F-term; T2 soft tree; T3 NJL/BCS+cutoff.
    Truncations and SUSY-breaking assumptions stated above.

(2) Nontrivial v:
    • Exact SUSY (T1): FORBIDDEN (non-renormalization / F_Y=mX=0). No critical
      Yukawa that generates <X>≠0.
    • Soft tree (T2): exists only for tachyonic soft mass r=m_s²/z < -1.
      Soft ratio is an INPUT.
    • NJL (T3): exists for Λ/√z above critical ~ 2π/√(hK). Cutoff is an INPUT.

(3) Predicted v²/z vs target 3.703e-4 h²:
    • T1: predicted=0 ⇒ ratio=0.
    • T2/T3: no unique prediction; ratio depends on free soft/cutoff parameter.
      Hitting the target means tuning that parameter (e.g. r_*={soft_r_that_hits_target(1.0):.4e} at h=1).

(4) Free parameters remaining: soft m_s²/z (T2) and/or Λ/√z, κ (T3).
    Dynamical condensate does NOT remove the tuning — it relocates it into
    the SUSY-breaking / UV input. This rung does not close the Newton gap
    without tuning.
"""
    )
    assert t1["v2_over_z"] == 0.0
    assert not t1["nontrivial_solution"]
    assert T2_soft_tree_gap(0.0)["nontrivial_solution"] is False
    assert T2_soft_tree_gap(-2.0)["nontrivial_solution"] is True
    assert soft_r_that_hits_target(1.0) < -1
    print("=== DYNAMICAL CONDENSATE RUNG COMPLETE ===")


if __name__ == "__main__":
    main()
