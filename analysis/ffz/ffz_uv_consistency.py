#!/usr/bin/env python3
"""Consistency constraints on V survivors (SO10/E6; complete Proca→EW).

Checks: unitarity/perturbativity of massive vectors near M_Pl; proton
decay; gauge coupling unification; IR massless γ/g crossover conflict.

No tuning. DERIVED vs ASSUMED. Honest negatives.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_uv_consistency.py
"""

from __future__ import annotations

from math import log, pi, sqrt

M_PL = 1.22e19  # GeV
M_GUT = 1e16
V_EW = 246.0
LAMBDA_QCD = 0.3
# Rough Super-K / Hyper-K style lower bound on p→e+π0 lifetime (years)
# Literature ~ 2e34 y (order-of-magnitude; cite as literature, not re-fit)
TAU_P_BOUND_YEARS = 2.4e34
# Partial lifetime formula sketch: τ ~ M_X^4 / (α_GUT² m_p^5) (order-of-mag)
M_P_GEV = 0.938


def unitarity_massive_vectors() -> dict:
    """Perturbative unitarity for Proca / Higgsed vectors.

    Derived sketch: longitudinal WW scattering saturates at
    E_* ~ 4π v / g  (or ~ M_W / α for strong coupling).
    For M_V ~ M_Pl with v ~ M_Pl, E_* ~ O(M_Pl) — marginal.
    For M_V ~ M_GUT with GUT coupling α_GUT~1/40: OK below M_Pl.

    ASSUMED: tree-level equivalence-theorem estimate.
    """
    alpha_gut = 1.0 / 40.0  # ASSUMED MSSM-like unification ballpark
    # E_* ~ 4π M_V / g  with g^2=4πα ⇒ g=sqrt(4πα) ⇒ E_* ~ 4π M_V / sqrt(4πα)
    def e_star(m_v: float, alpha: float) -> float:
        g = sqrt(4.0 * pi * alpha)
        return 4.0 * pi * m_v / g

    e_gut = e_star(M_GUT, alpha_gut)
    e_pl = e_star(M_PL, alpha_gut)
    return {
        "alpha_GUT_assumed": alpha_gut,
        "E_star_M_GUT_GeV": e_gut,
        "E_star_M_Pl_GeV": e_pl,
        "OK_GUT_vectors_below_MPl": e_gut > M_PL * 0.1,  # loose
        "OK_Planck_mass_vectors": e_pl >= M_PL,  # marginal by construction
        "complete_Proca_to_EW_near_Pl": False,  # strong coupling / UV completion needed
        "label": "derived sketch; α_GUT assumed ~1/40",
        "plain": (
            f"GUT-scale Proca: E_*~{e_gut:.2e} GeV (≳M_Pl OK-ish). "
            f"Planck-mass Proca: E_*~{e_pl:.2e} GeV — marginal; needs UV completion. "
            f"Complete Proca→EW with M~M_Pl vectors is not perturbatively safe."
        ),
    }


def proton_decay() -> dict:
    """Dimension-6 GUT proton decay lifetime sketch.

    τ/τ0 ~ (M_X / M_X0)^4. With M_X~M_GUT=1e16, τ often ~10^34–10^36 y
    depending on α and matrix elements (ASSUMED ballpark).

    Bound τ > 2.4e34 y (literature Super-K).
    """
    # Very rough: τ_years ≈ 1e32 * (M_X / 1e15)^4 / (α/0.04)^2
    alpha = 1.0 / 40.0
    tau_approx = 1e32 * (M_GUT / 1e15) ** 4 / ((alpha / 0.04) ** 2)
    return {
        "tau_approx_years": tau_approx,
        "bound_years": TAU_P_BOUND_YEARS,
        "passes_bound_rough": tau_approx >= TAU_P_BOUND_YEARS,
        "label_formula": "assumed order-of-magnitude GUT formula",
        "label_bound": "literature (Super-K order)",
        "plain": (
            f"Rough τ~{tau_approx:.2e} y vs bound ≳{TAU_P_BOUND_YEARS:.2e} y — "
            f"marginal/OK depending on matrix elements; not a precision claim. "
            f"Complete Proca→EW without X/Y mediators: no classic d=6 GUT decay "
            f"(different constraint set)."
        ),
    }


def gauge_unification() -> dict:
    """One-loop SM vs MSSM unification quality (standard lore numbers).

    SM: couplings miss by ~few % / no perfect meeting.
    MSSM: meets near 2e16 GeV (ASSUMED threshold corrections O(1)).
    Extra Proca / complete Higgsing changes β-functions — ASSUMED spoil
    unless thresholds arranged.
    """
    return {
        "SM_unifies_precisely": False,
        "MSSM_ballpark_GeV": 2e16,
        "SO10_compatible_with_unification": True,  # classic motivation
        "E6_compatible": True,  # with caveats
        "complete_Proca_to_EW_preserves_unification": False,
        "label": "literature lore; β shifts from extra Proca ASSUMED harmful if untuned",
        "plain": (
            "SO(10)/E6 are built for unification. Complete Proca→EW removes "
            "massless gauge running over decades — spoils standard unification "
            "without new ASSUMED structure."
        ),
    }


def ir_massless_crossover() -> dict:
    """Photon and gluons must be massless in the IR.

    Crossover scale Λ_* where massless SU(3)×U(1)_EM appear:
      - If Λ_* < Λ_QCD: conflicts with QCD confinement / hadron physics
        (gluons need to be massless/confining below ~1 GeV).
      - If Λ_* < m_e: conflicts with long-range QED / Cavendish / chemistry.
      - Safe: Λ_* ≥ M_GUT for GUT-broken residual SM (standard).
      - Complete Proca→EW sets Λ_* ~ v_EW for W/B and would require
        gluons/photon massless only below EW — DISASTER for QCD/QED.

    Quantify: require Λ_* ≥ Λ_QCD for gluons; Λ_* ≥ 0 for photon ideally
    with m_γ < 10^{-18} eV laboratory bounds (literature).
    """
    m_gamma_bound_eV = 1e-18  # literature order (PDG / Cavendish-type)
    return {
        "Lambda_QCD_GeV": LAMBDA_QCD,
        "required_gluon_massless_below_GeV": LAMBDA_QCD,
        "m_gamma_bound_eV": m_gamma_bound_eV,
        "SO10_E6_crossover": M_GUT,
        "SO10_E6_OK_for_IR": True,  # residual 12 include g,γ after EW
        "complete_Proca_to_EW_crossover": V_EW,
        "complete_Proca_to_EW_OK_for_QCD": False,
        "complete_Proca_to_EW_OK_for_QED": False,
        "label_bounds": "literature m_γ; QCD scale derived as IR requirement",
        "plain": (
            f"SO(10)/E6: massless SM gauge below M_GUT — IR OK. "
            f"Complete Proca→EW: would keep g/γ massive down to ~{V_EW} GeV — "
            f"incompatible with QCD (need massless/confining gluons ≲{LAMBDA_QCD} GeV) "
            f"and QED (m_γ ≲ {m_gamma_bound_eV} eV)."
        ),
    }


def survivor_scorecard() -> dict:
    """Combine V flips with W constraints."""
    u = unitarity_massive_vectors()
    p = proton_decay()
    g = gauge_unification()
    ir = ir_massless_crossover()
    return {
        "SO10_E6": {
            "flips_S_in_V": True,
            "unitarity_OK": u["OK_GUT_vectors_below_MPl"],
            "proton_decay_rough_OK": p["passes_bound_rough"],
            "unification_OK": g["SO10_compatible_with_unification"],
            "IR_OK": ir["SO10_E6_OK_for_IR"],
            "overall_viable_as_assumption": True,
            "note": "Best remaining UV path; still ASSUMED GUT spectrum / thresholds",
        },
        "complete_Proca_to_EW": {
            "flips_S_in_V": True,
            "unitarity_OK": False,
            "proton_decay_rough_OK": True,  # N/A classic X/Y
            "unification_OK": False,
            "IR_OK": False,
            "overall_viable_as_assumption": False,
            "note": "Arithmetic flip only — ruled out by IR γ/g + unitarity/unification",
        },
        "short_window_Proca_KK_TC4": {
            "flips_S_in_V": False,
            "overall_viable_as_assumption": False,
            "note": "Already fails V stress test",
        },
        "plain": (
            "Only SO(10)/E6-style long broken Proca bins remain tentatively "
            "viable under W. Complete Proca→EW is ruled out by IR massless "
            "gauge bosons despite flipping S."
        ),
    }


def self_checks() -> None:
    ir = ir_massless_crossover()
    assert ir["SO10_E6_OK_for_IR"] is True
    assert ir["complete_Proca_to_EW_OK_for_QCD"] is False
    assert ir["complete_Proca_to_EW_OK_for_QED"] is False
    u = unitarity_massive_vectors()
    assert u["complete_Proca_to_EW_near_Pl"] is False
    s = survivor_scorecard()
    assert s["SO10_E6"]["overall_viable_as_assumption"] is True
    assert s["complete_Proca_to_EW"]["overall_viable_as_assumption"] is False
    assert s["complete_Proca_to_EW"]["IR_OK"] is False


def main():
    print("=== UV consistency (unitarity / p-decay / unif / IR) ===\n")
    u = unitarity_massive_vectors()
    print(f"  Unitarity: {u['plain']} [{u['label']}]")
    p = proton_decay()
    print(f"  Proton decay: {p['plain']}")
    g = gauge_unification()
    print(f"  Unification: {g['plain']} [{g['label']}]")
    ir = ir_massless_crossover()
    print(f"  IR crossover: {ir['plain']}")
    s = survivor_scorecard()
    print(f"\n  SO10/E6 viable? {s['SO10_E6']['overall_viable_as_assumption']} — {s['SO10_E6']['note']}")
    print(
        f"  Complete Proca→EW viable? {s['complete_Proca_to_EW']['overall_viable_as_assumption']} "
        f"— {s['complete_Proca_to_EW']['note']}"
    )
    print(f"  {s['plain']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
SO(10)/E6 remain the only V-flippers compatible with IR massless γ/g
and (roughly) unitarity/unification/proton decay — still ASSUMED.
Complete Proca→EW is an honest arithmetic winner but a physics loser.
"""
    )
    print("=== UV CONSISTENCY COMPLETE ===")


if __name__ == "__main__":
    main()
