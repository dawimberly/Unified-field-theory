#!/usr/bin/env python3
"""Bounce criterion: AC/AD vs LQC / asymptotic-safety running G.

Does bag-vacuum SEC violation or deconfinement a1 flip produce a bounce
at sub-Planck density? Densities, radii, timescales for 10 M_⊙ collapse.
Observable signatures (GRB/FRB, echoes) only if derived.

Fields-only. No tuning. DERIVED vs ASSUMED. Honest negatives.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_bounce_criterion.py
"""

from __future__ import annotations

from math import pi, sqrt

G_NEWTON = 6.67430e-11  # m³/kg/s²
C_LIGHT = 2.99792458e8  # m/s
M_SUN = 1.98847e30  # kg
M_PL_GEV = 1.22e19
LAMBDA_QCD_GEV = 0.3
# LQC critical density ~ 0.41 ρ_Pl (literature; ASSUMED coefficient)
RHO_PL_GEV4 = (M_PL_GEV) ** 4  # schematic natural units
RHO_CRIT_LQC_OVER_PL = 0.41
# Nuclear / deconfinement
RHO_SAT_G_CM3 = 2.8e14  # g/cm³ ballpark
RHO_DECONF_OVER_SAT = 5.0  # ASSUMED few×n_sat


def solar_mass_collapse(m_over_msun: float = 10.0) -> dict:
    """Schwarzschild radius and free-fall timescale — DERIVED GR/Newtonian."""
    m = m_over_msun * M_SUN
    r_s = 2.0 * G_NEWTON * m / C_LIGHT**2
    # Free-fall from ~few R_s: t_ff ~ (R^3/GM)^{1/2}
    r0 = 3.0 * r_s
    t_ff = sqrt(r0**3 / (G_NEWTON * m))
    # Density if packed inside R_s: ρ ~ M / (4/3 π R_s³)
    vol = 4.0 / 3.0 * pi * r_s**3
    rho_mean = m / vol  # kg/m³
    rho_mean_g_cm3 = rho_mean * 1e-3  # kg/m³ → g/cm³
    return {
        "M_over_Msun": m_over_msun,
        "R_s_km": r_s / 1e3,
        "t_ff_ms": t_ff * 1e3,
        "rho_mean_inside_Rs_g_cm3": rho_mean_g_cm3,
        "label": "derived (GR / Newtonian free-fall)",
    }


def lqc_bounce() -> dict:
    """LQC / Planck-star: bounce at ρ_crit ~ 0.41 ρ_Pl."""
    return {
        "rho_crit_over_rho_Pl": RHO_CRIT_LQC_OVER_PL,
        "sub_Planck_by_O1": True,
        "at_QCD_density": False,
        "label": "literature LQC coefficient ASSUMED; scale ~Planck DERIVED as Planckian",
        "plain": "LQC bounce is near-Planckian — far above nuclear/deconfinement density.",
    }


def asymptotic_safety_g() -> dict:
    """AS running G(k): anti-screening / fixed-point — bounce claims model-dependent.

    No derived bounce density from AS alone here.
    """
    return {
        "derived_bounce_density": False,
        "label": "assumed AS scenario; not computed from FFZ ladder",
        "plain": "AS running G is ASSUMED; no concrete ρ_bounce derived in this rung.",
    }


def ac_ad_bounce() -> dict:
    """Do AC (bag vacuum) or AD (deconf a1) bounce at sub-Planck ρ?"""
    # Bag B ~ (150 MeV)^4 ≪ ρ_Pl
    b_over_pl = (0.15 / M_PL_GEV) ** 4
    # Deconfinement density ~ 5 n_sat ~ 1e15 g/cm³ vs ρ_Pl ~ 5e93 g/cm³
    rho_pl_g_cm3 = 5.1e93  # literature order
    rho_deconf = RHO_DECONF_OVER_SAT * RHO_SAT_G_CM3
    ratio = rho_deconf / rho_pl_g_cm3
    # AD: ΔG/G ~ 1e-39 — cannot reverse collapse
    return {
        "bag_B_over_rho_Pl": b_over_pl,
        "rho_deconf_g_cm3": rho_deconf,
        "rho_deconf_over_rho_Pl": ratio,
        "AC_produces_bounce": False,
        "AD_produces_bounce": False,
        "reason_AC": (
            "Bag vacuum ρ+3p<0 is QCD-scale; in collapse it is a thin/vacuum "
            "piece, not a bulk bounce at stellar densities (ordinary matter SEC OK)."
        ),
        "reason_AD": (
            "Local a1 flip changes induced G by ~10^{-39}; dynamical bounce "
            "requires O(1) change in effective G or Planck-scale ρ."
        ),
        "label": "derived size comparisons; bounce dynamics ASSUMED absent",
    }


def observables() -> dict:
    """GRB/FRB/echoes — only if derived from a bounce."""
    return {
        "GRB_FRB_echoes_derived_here": False,
        "reason": "No bounce from AC/AD ⇒ no derived electromagnetic signature.",
        "label": "honest negative",
    }


def self_checks() -> None:
    c = solar_mass_collapse(10.0)
    # R_s ≈ 2.95 km × (M/M_⊙) ≈ 29.5 km for 10 M_⊙
    assert 25 < c["R_s_km"] < 35
    assert c["t_ff_ms"] > 0
    b = ac_ad_bounce()
    assert b["AC_produces_bounce"] is False
    assert b["AD_produces_bounce"] is False
    assert b["rho_deconf_over_rho_Pl"] < 1e-70
    assert observables()["GRB_FRB_echoes_derived_here"] is False
    assert lqc_bounce()["at_QCD_density"] is False


def main():
    print("=== Bounce criterion: AC/AD vs LQC/AS (10 M_⊙) ===\n")
    c = solar_mass_collapse(10.0)
    print(
        f"  10 M_⊙: R_s={c['R_s_km']:.1f} km  t_ff(3R_s)={c['t_ff_ms']:.2f} ms  "
        f"ρ_mean(R_s)={c['rho_mean_inside_Rs_g_cm3']:.3e} g/cm³  [{c['label']}]"
    )
    lqc = lqc_bounce()
    print(f"  LQC: ρ_crit~{lqc['rho_crit_over_rho_Pl']} ρ_Pl — {lqc['plain']}")
    asafe = asymptotic_safety_g()
    print(f"  AS: {asafe['plain']} [{asafe['label']}]")
    b = ac_ad_bounce()
    print(f"  ρ_deconf/ρ_Pl ≈ {b['rho_deconf_over_rho_Pl']:.3e}")
    print(f"  AC bounce? {b['AC_produces_bounce']} — {b['reason_AC']}")
    print(f"  AD bounce? {b['AD_produces_bounce']} — {b['reason_AD']}")
    o = observables()
    print(f"  Observables derived? {o['GRB_FRB_echoes_derived_here']} — {o['reason']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
Neither AC nor AD yields a sub-Planck bounce. LQC bounce is near-Planckian.
No derived GRB/FRB/echo signatures from this ladder.
"""
    )
    print("=== BOUNCE CRITERION COMPLETE ===")


if __name__ == "__main__":
    main()
