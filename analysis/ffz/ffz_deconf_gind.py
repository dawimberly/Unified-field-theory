#!/usr/bin/env python3
"""Induced-gravity shift across deconfinement in collapse.

Using #36/#40: hadronic a1=1/3±0.29 (gluon Kabats absent); at
deconfinement free-gluon Kabats return (8×−2/3). Compute local shift
in induced 1/G. Sign/magnitude vs M_Pl²? Likely ~ Λ_QCD²/M_Pl² — quantify.

Fields-only. No tuning. DERIVED vs ASSUMED. Honest negatives.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_deconf_gind.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import pi

A1_HAD_CENTRAL = F(1, 3)
A1_HAD_LOW = F(1, 6)
A1_HAD_HIGH = F(3, 4)
A1_V = F(-2, 3)
# Restoring 8 massless gluon modes adds 8×(−2/3)=−16/3 to a1.
DELTA_A1_RESTORE_8 = 8 * A1_V  # −16/3

LAMBDA_QCD_GEV = 0.3
M_PL_GEV = 1.22e19
N_SAT_GEV_FM3 = 0.15  # rough GeV/fm³ order for n_sat energy density
# Deconfinement crossover ~ few × n_sat (ASSUMED lore)


def a1_across_deconfinement() -> dict:
    a1_had = A1_HAD_CENTRAL
    a1_dec = a1_had + DELTA_A1_RESTORE_8  # 1/3 − 16/3 = −15/3 = −5
    assert DELTA_A1_RESTORE_8 == F(-16, 3)
    assert a1_dec == F(-5)
    return {
        "a1_hadronic": a1_had,
        "a1_deconfined": a1_dec,
        "delta_a1": DELTA_A1_RESTORE_8,
        "sign_flip_local_a1": (a1_had > 0) and (a1_dec < 0),
        "label": "derived from #36/#40 Kabat restoration arithmetic",
        "plain": (
            f"Hadronic a1={a1_had}>0; deconfined a1={a1_dec}<0 "
            f"(Δ={DELTA_A1_RESTORE_8}). Local a1 sign flips across deconfinement."
        ),
    }


def induced_g_shift() -> dict:
    """Sakharov sketch: 1/(16πG_ind) ~ a1 Λ²/(12π) with Λ~Λ_QCD for local IR.

    Absolute 1/G_ind from QCD-scale cutoff is tiny vs M_Pl²:
      (Λ_QCD / M_Pl)² ~ (0.3/1.22e19)² ~ 6e-40.
    Local fractional change of Planck G is negligible unless the UV
    Sakharov piece is somehow dominated by the IR a1 (it is not — #41).
    """
    ratio = (LAMBDA_QCD_GEV / M_PL_GEV) ** 2
    # Relative shift of IR-only piece: Δ(1/G)_IR / (1/G)_Pl ~ |Δa1/a1_eff| * ratio
    # With |Δa1|=16/3≈5.33 and a1_UV~O(1), still O(ratio)
    delta_frac = abs(float(DELTA_A1_RESTORE_8)) * ratio  # crude
    return {
        "Lambda_QCD_over_MPl_sq": ratio,
        "delta_a1": float(DELTA_A1_RESTORE_8),
        "frac_shift_vs_Planck_G": delta_frac,
        "local_a1_sign_flips": True,
        "G_ind_Planck_sign_flips": False,  # UV Kabats dominate integrated (#41)
        "honest_size": (
            f"(Λ_QCD/M_Pl)² ≈ {ratio:.3e}; fractional IR-induced shift "
            f"≲ {delta_frac:.3e} — utterly negligible vs Planck G."
        ),
        "label_ratio": "derived",
        "label_dominance": "derived from #41 integrated S",
    }


def deconfinement_density() -> dict:
    """Crossover at few × n_sat (ASSUMED lore / lattice QCD)."""
    return {
        "n_over_n_sat_ballpark": (2.0, 10.0),
        "T_c_MeV_ballpark": (150.0, 170.0),
        "label": "assumed lattice/heavy-ion lore",
    }


def self_checks() -> None:
    a = a1_across_deconfinement()
    assert a["sign_flip_local_a1"] is True
    assert a["a1_deconfined"] == F(-5)
    g = induced_g_shift()
    assert g["G_ind_Planck_sign_flips"] is False
    assert g["Lambda_QCD_over_MPl_sq"] < 1e-35
    assert g["frac_shift_vs_Planck_G"] < 1e-35


def main():
    print("=== Induced G across deconfinement (collapse) ===\n")
    a = a1_across_deconfinement()
    print(f"  {a['plain']}  [{a['label']}]")
    g = induced_g_shift()
    print(f"  {g['honest_size']}")
    print(f"  Local a1 sign flips? {g['local_a1_sign_flips']}")
    print(f"  Planck-scale G_ind flips? {g['G_ind_Planck_sign_flips']}  [{g['label_dominance']}]")
    d = deconfinement_density()
    print(f"  Crossover ~ {d['n_over_n_sat_ballpark']} n_sat, "
          f"T~{d['T_c_MeV_ballpark']} MeV  [{d['label']}]")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
Local a1 flips − at deconfinement (Kabat restoration). Effect on
physical G is ~10^{-39} relative to M_Pl — negligible. No bounce from
this alone.
"""
    )
    print("=== DECONF G_IND COMPLETE ===")


if __name__ == "__main__":
    main()
