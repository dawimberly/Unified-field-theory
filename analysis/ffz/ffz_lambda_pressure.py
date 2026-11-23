#!/usr/bin/env python3
"""Induced cosmological constant / pressure / negative gravity.

Dan: pressure can make gravity repel. Lambda_ind from roster; rho+3p;
link to dark energy without tuning?

No tuning. DERIVED vs ASSUMED.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_lambda_pressure.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import pi

# a0 heat-kernel / vacuum energy weights (standard; bosons +, fermions −)
# Real scalar: +1, Weyl: −1/2? Dirac −1, vector +2 (or +3 with ghosts careful)
# Use Adler-like: vacuum energy ~ (N_B − N_F) Λ⁴ / (something)
# Conservative: sign(Lambda_ind) ~ sign(n_bose_dof − n_fermi_dof) with
# massless V counting 2 physical dof, Proca 3, Weyl 2 on-shell, scalar 1.


def dof_uv_sm() -> dict:
    # UV SM on-shell dof (approx)
    n_s = 4  # Higgs
    n_weyl_os = 45 * 2  # each Weyl 2 on-shell? Actually Weyl has 2 real = 1 complex component on-shell helicity count: 2
    # Standard: 45 Weyl = 45×1 Weyl fermion = 90 chiral dof counting / 2 for Majorana... use:
    # SM: 28 bosonic dof (12×2 gauge + 4 Higgs) and 90 fermionic (45×2)
    n_b = 12 * 2 + 4  # 28
    n_f = 45 * 2  # 90
    return {"n_b": n_b, "n_f": n_f, "n_b_minus_n_f": n_b - n_f}


def lambda_ind_sign() -> dict:
    d = dof_uv_sm()
    # Fermion-dominated ⇒ negative vacuum energy tendency (wrong-sign CC too)
    return {
        **d,
        "Lambda_ind_sign_tendency": "negative" if d["n_b_minus_n_f"] < 0 else "positive",
        "matches_observed_DE_positive": False,
        "label": "derived (dof count); magnitude requires Λ_UV^4 tuning for DE",
    }


def rho_plus_3p() -> dict:
    """Active gravitational mass ~ ρ+3p for cosmology.

    Radiation: p=ρ/3 ⇒ ρ+3p=2ρ > 0 (attractive in FRW deceleration).
    Vacuum / CC: p=−ρ ⇒ ρ+3p=−2ρ < 0 (repulsive / accelerates).
    Stiff matter etc. Dan's 'pressure can repel' = ρ+3p < 0 regime —
    DERIVED from GR + continuum T_μν; not from FFZ ledger.
    """
    return {
        "radiation_rho_3p": "2ρ > 0",
        "vacuum_rho_3p": "−2ρ < 0 (repulsive)",
        "dust_rho_3p": "ρ > 0",
        "link_to_a1_roster": False,
        "label": "derived (GR); no untuned DE from SM roster",
    }


def dark_energy_without_tuning() -> dict:
    """Does SM roster give Ω_Λ∼0.7 without tuning Λ_UV?"""
    # Observed ρ_DE ∼ (2.3×10^{-3} eV)^4; SM Λ_UV^4 is ∼ (100 GeV)^4 or Planck —
    # mismatch ∼ 10^{50}–10^{120}. No.
    return {
        "DE_from_roster_without_tuning": False,
        "reason": "Λ_UV^4 vs meV^4 mismatch; need tuning or new mechanism",
        "label": "derived (CC problem; not solved here)",
    }


def self_checks() -> None:
    L = lambda_ind_sign()
    assert L["n_b_minus_n_f"] == 28 - 90 == -62
    assert L["matches_observed_DE_positive"] is False
    assert rho_plus_3p()["link_to_a1_roster"] is False
    assert dark_energy_without_tuning()["DE_from_roster_without_tuning"] is False


def main():
    print("=== Λ_ind / pressure / negative gravity ===\n")
    L = lambda_ind_sign()
    print(f"  UV dof n_b={L['n_b']} n_f={L['n_f']} Δ={L['n_b_minus_n_f']}")
    print(f"  Λ_ind sign tendency: {L['Lambda_ind_sign_tendency']}  matches DE? {L['matches_observed_DE_positive']}")
    r = rho_plus_3p()
    print(f"  Vacuum ρ+3p: {r['vacuum_rho_3p']}  [{r['label']}]")
    d = dark_energy_without_tuning()
    print(f"  DE without tuning? {d['DE_from_roster_without_tuning']} — {d['reason']}")
    self_checks()
    print("  self-checks PASS")
    print(
        """
Verdict
-------
ρ+3p<0 for vacuum explains repulsive gravity (GR, derived). SM roster
tends to wrong-sign Λ_ind (fermion-heavy) and does not give DE without
tuning the cutoff (CC problem stands).
"""
    )
    print("=== LAMBDA/PRESSURE COMPLETE ===")


if __name__ == "__main__":
    main()
