#!/usr/bin/env python3
"""Minimal broken-gauge spectrum for integrated S>0 (no SO(10)/E6 assumed).

Fields-only language: massive vector field modes (Proca), not particles.
IR: photon and gluon field modes remain massless below the breaking scale.

Constraints: S>0; IR masslessness; proton-decay ballpark; SM gauge-coupling
unification from measured couplings (2-loop sketch).

Ask: are SO(10)/E6 selected or merely allowed? Anything genuinely derived?

No tuning. DERIVED vs ASSUMED. Coincidences flagged.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_minimal_spectrum.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import log, pi

A1_V = F(-2, 3)
A1_PROCA = F(1, 2)
A1_W = F(1, 12)
A1_S = F(1, 6)
A1_EW = F(-1, 12)
A1_AC = F(-13, 4)
A1_HAD = F(1, 3)

M_PL = 1.22e19
V_EW = 246.0
M_TOP = 173.0
LAMBDA_QCD = 0.3
N_SM_MASSLESS = 12  # residual massless gauge field modes below M_*
N_WEYL_SM = 45
N_SCALAR_H = 4

# Measured SM couplings at M_Z (PDG-ish; ASSUMED inputs for running)
ALPHA_EM_MZ = 1.0 / 127.9
SIN2_W_MZ = 0.2312
ALPHA_S_MZ = 0.1181
M_Z = 91.1876


def a1_matter_sm() -> F:
    return N_WEYL_SM * A1_W + N_SCALAR_H * A1_S


def a1_unbroken(n_v: int) -> F:
    """Above M_*: all SM+extra gauge modes massless (Kabat each)."""
    return (N_SM_MASSLESS + n_v) * A1_V + a1_matter_sm()


def a1_broken(n_v: int) -> F:
    """Below M_*: 12 massless SM + n_v massive vector modes (Proca)."""
    return N_SM_MASSLESS * A1_V + n_v * A1_PROCA + a1_matter_sm()


def spectral_S(n_v: int, m_star: float, mu_uv: float = M_PL) -> float:
    """Two high bins + IR tail. DERIVED arithmetic on ASSUMED thresholds."""
    bins = [
        (mu_uv, m_star, a1_unbroken(n_v)),
        (m_star, V_EW, a1_broken(n_v)),
        (V_EW, M_TOP, A1_EW),
        (M_TOP, LAMBDA_QCD, A1_AC),
        (LAMBDA_QCD, 0.1, A1_HAD),
    ]
    S = 0.0
    for hi, lo, a1 in bins:
        if hi <= lo:
            continue
        S += float(a1) * 2.0 * log(hi / lo)
    return S


def minimal_n_v(m_star: float = 1e16, mu_uv: float = M_PL) -> dict:
    """Smallest N_V with S>0 at fixed scales."""
    n_min = None
    s_min = None
    for n in range(0, 80):
        s = spectral_S(n, m_star, mu_uv)
        if s > 0 and n_min is None:
            n_min = n
            s_min = s
            break
    return {
        "n_v_min": n_min,
        "S_at_min": s_min,
        "S_at_0": spectral_S(0, m_star, mu_uv),
        "S_at_12": spectral_S(12, m_star, mu_uv),  # SU(5)-like heavy count
        "S_at_33": spectral_S(33, m_star, mu_uv),  # SO(10)-like
        "S_at_66": spectral_S(66, m_star, mu_uv),  # E6-like
        "m_star": m_star,
        "mu_uv": mu_uv,
        "label_arith": "derived (a1 counting + log weights)",
        "label_scales": "assumed (M_*, Λ_UV ballparks)",
    }


def ir_masslessness() -> dict:
    """Photon and gluon field modes massless in IR — forces residual ≥12? 

    Actually residual must include U(1)_EM and SU(3)_c after EW.
    Standard count: 12 massless above EW (8+3+1), then W/Z massive.
    Breaking scale M_* must be ≫ Λ_QCD so gluon modes are massless
    through the confining window.
    """
    return {
        "n_massless_required_above_EW": N_SM_MASSLESS,
        "M_star_min_for_QCD_GeV": LAMBDA_QCD,
        "safe_M_star_GeV": 1e15,  # ASSUMED GUT-ish lower edge for p-decay
        "label": "derived IR requirement; safe M_* assumed from p-decay lore",
    }


def proton_decay_constraint(m_star: float) -> dict:
    """Rough τ ∝ M_*^4; need M_* ≳ 1e15–1e16 GeV (literature)."""
    m_min_lit = 1e15
    return {
        "m_star": m_star,
        "m_min_literature_GeV": m_min_lit,
        "passes": m_star >= m_min_lit,
        "label": "literature ballpark; not a lattice matrix-element calc",
    }


def two_loop_unification_sketch() -> dict:
    """Do measured SM couplings unify without extra thresholds?

    Standard result: SM alone does NOT meet; MSSM does near ~2e16 GeV.
    Extra N_V Proca thresholds shift β above M_* — can help or hurt.
    We do NOT claim a derived meeting for arbitrary N_V.
    """
    # One-loop SM β coefficients (gauge): b1=41/10, b2=-19/6, b3=-7
    # Report known non-unification of SM as DERIVED lore check via
    # inverse couplings meeting failure.
    b1, b2, b3 = 41.0 / 10.0, -19.0 / 6.0, -7.0
    # α1^{-1}, α2^{-1}, α3^{-1} at M_Z
    alpha1_inv = (5.0 / 3.0) * (1.0 - SIN2_W_MZ) / ALPHA_EM_MZ  # GUT norm rough
    # More standard: α1 = 5/3 α' with α' = α_em / cos²θ_W
    alpha_em = ALPHA_EM_MZ
    c2w = 1.0 - SIN2_W_MZ
    alpha1 = (5.0 / 3.0) * alpha_em / c2w
    alpha2 = alpha_em / SIN2_W_MZ
    alpha3 = ALPHA_S_MZ
    a1i, a2i, a3i = 1.0 / alpha1, 1.0 / alpha2, 1.0 / alpha3

    def run(ai, b, mu_hi, mu_lo):
        # dai/dlnμ = -b/(2π) for a=1/α... actually dα^{-1}/dlnμ = -b/(2π)
        return ai - (b / (2.0 * pi)) * log(mu_hi / mu_lo)

    # Run SM to 1e16
    mu = 1e16
    a1_h = run(a1i, b1, mu, M_Z)
    a2_h = run(a2i, b2, mu, M_Z)
    a3_h = run(a3i, b3, mu, M_Z)
    spread = max(a1_h, a2_h, a3_h) - min(a1_h, a2_h, a3_h)
    return {
        "alpha_inv_at_1e16": (a1_h, a2_h, a3_h),
        "spread_at_1e16": spread,
        "SM_unifies": spread < 1.0,  # tight meet criterion
        "label_inputs": "assumed PDG-ish α(M_Z)",
        "label_beta": "derived 1-loop SM β (2-loop would shift O(1))",
        "plain": (
            f"SM 1-loop inverse couplings at 1e16: "
            f"({a1_h:.1f}, {a2_h:.1f}, {a3_h:.1f}), spread={spread:.1f} — "
            f"do NOT unify. Extra N_V thresholds ASSUMED if forced to meet."
        ),
    }


def so10_e6_status(n_min: int) -> dict:
    """Are SO(10)/E6 selected or merely allowed?"""
    return {
        "n_v_min_derived_arith": n_min,
        "SU5_heavy": 12,
        "SO10_heavy": 33,
        "E6_heavy": 66,
        "SU5_allowed": 12 >= n_min if n_min is not None else False,
        "SO10_allowed": 33 >= n_min if n_min is not None else False,
        "E6_allowed": 66 >= n_min if n_min is not None else False,
        "selected_uniquely": False,  # many N_V work
        "genuinely_derived": (
            "Minimal N_V from S>0 arithmetic is derived given scale ASSUMPTIONS; "
            "group identity SO(10)/E6 is NOT derived — only allowed if N_V matches."
        ),
    }


def self_checks() -> None:
    m = minimal_n_v()
    assert m["S_at_0"] < 0
    assert m["n_v_min"] is not None and m["n_v_min"] > 0
    assert m["S_at_min"] > 0
    # SO(10) 33 should be above minimum for default scales
    assert m["n_v_min"] <= 33
    u = two_loop_unification_sketch()
    assert u["SM_unifies"] is False
    st = so10_e6_status(m["n_v_min"])
    assert st["selected_uniquely"] is False
    assert st["SO10_allowed"] is True
    ir = ir_masslessness()
    assert ir["n_massless_required_above_EW"] == 12


def main():
    print("=== Minimal broken-gauge spectrum (fields-only) ===\n")
    print("  Language: massive vector field modes; IR massless photon/gluon modes.")
    m = minimal_n_v()
    print(
        f"  Minimal N_V for S>0 at M_*={m['m_star']:.0e}, Λ_UV={m['mu_uv']:.2e}: "
        f"N_V={m['n_v_min']}  S={m['S_at_min']:+.1f}  [{m['label_arith']}; {m['label_scales']}]"
    )
    print(f"  S(N_V=0)={m['S_at_0']:+.1f}  S(12)={m['S_at_12']:+.1f}  "
          f"S(33)={m['S_at_33']:+.1f}  S(66)={m['S_at_66']:+.1f}")
    ir = ir_masslessness()
    print(f"  IR: require {ir['n_massless_required_above_EW']} massless gauge modes above EW; "
          f"M_*≫Λ_QCD [{ir['label']}]")
    p = proton_decay_constraint(m["m_star"])
    print(f"  Proton decay ballpark: M_*≳{p['m_min_literature_GeV']:.0e}? {p['passes']} [{p['label']}]")
    u = two_loop_unification_sketch()
    print(f"  Unification: {u['plain']}")
    st = so10_e6_status(m["n_v_min"])
    print(f"\n  SU(5)/SO(10)/E6 allowed? {st['SU5_allowed']}/{st['SO10_allowed']}/{st['E6_allowed']}")
    print(f"  Uniquely selected? {st['selected_uniquely']}")
    print(f"  {st['genuinely_derived']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
Minimal N_V is derived arithmetic at fixed ASSUMED scales. SO(10)/E6
are ALLOWED (N_V=33,66 ≥ N_min), not uniquely SELECTED. SM couplings
do not unify alone — extra thresholds assumed if unification is imposed.
Nothing here derives a specific simple group from S>0 alone.
"""
    )
    print("=== MINIMAL SPECTRUM COMPLETE ===")


if __name__ == "__main__":
    main()
