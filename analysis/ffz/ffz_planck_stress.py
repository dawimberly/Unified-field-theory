#!/usr/bin/env python3
"""Planck→GUT stress test (#46 style) for U survivors.

For each surviving UV-mass mechanism: spectral integral S with
Λ_UV∼M_Pl. Report sign, magnitude, ASSUMED inputs.

No tuning. DERIVED vs ASSUMED. Coincidences flagged.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_planck_stress.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import log

A1_V = F(-2, 3)
A1_PROCA = F(1, 2)
A1_W = F(1, 12)
A1_S = F(1, 6)
A1_EW = F(-1, 12)
A1_AC = F(-13, 4)
A1_HAD = F(1, 3)
A1_SM_UV = F(-43, 12)  # 4 scalars ξ=0 + 45 Weyl + 12 massless

M_PL = 1.22e19
M_GUT = 1e16
V_EW = 246.0
M_TOP = 173.0
LAMBDA_QCD = 0.3


def a1_sm_matter() -> F:
    """Fermions + Higgs scalars (ξ=0), no gauge."""
    return 45 * A1_W + 4 * A1_S  # 15/4 + 2/3 = 53/12


def a1_gauge(n_massless: int, n_proca: int) -> F:
    return n_massless * A1_V + n_proca * A1_PROCA


def spectral_S(bins: list[dict]) -> dict:
    S = 0.0
    rows = []
    for b in bins:
        if b["mu_hi"] <= b["mu_lo"]:
            continue
        w = 2.0 * log(b["mu_hi"] / b["mu_lo"])
        c = float(b["a1"]) * w
        S += c
        rows.append({**b, "weight": w, "contrib": c})
    return {"rows": rows, "S": S, "sign_positive": S > 0}


def low_energy_tail() -> list[dict]:
    return [
        {"name": "EW_to_top", "mu_hi": V_EW, "mu_lo": M_TOP, "a1": A1_EW},
        {"name": "top_to_QCD", "mu_hi": M_TOP, "mu_lo": LAMBDA_QCD, "a1": A1_AC},
        {"name": "hadronic", "mu_hi": LAMBDA_QCD, "mu_lo": 0.1, "a1": A1_HAD},
    ]


def baseline_SM() -> dict:
    """#41-style: SM UV from M_Pl→v (all SM massless gauge in one fat bin)."""
    bins = [
        {
            "name": "SM_UV_Planck_to_EW",
            "mu_hi": M_PL,
            "mu_lo": V_EW,
            "a1": A1_SM_UV,
        }
    ] + low_energy_tail()
    r = spectral_S(bins)
    r["mechanism"] = "baseline_SM"
    r["assumed"] = ["free-field a1 thresholds", "Λ_UV=M_Pl"]
    return r


def stress_complete_higgsing_short() -> dict:
    """All 12 Proca only in short Planck→GUT; SM massless restored below.

    Honest #46-style stress: short high-scale window. Likely FAILS because
    the long M_GUT→EW bin still has SM a1=−43/12.
    """
    a1_uv = a1_gauge(0, 12) + a1_sm_matter()  # 6 + 53/12 = 125/12
    bins = [
        {"name": "all_Proca_Planck_GUT", "mu_hi": M_PL, "mu_lo": M_GUT, "a1": a1_uv},
        {"name": "SM_GUT_to_EW", "mu_hi": M_GUT, "mu_lo": V_EW, "a1": A1_SM_UV},
    ] + low_energy_tail()
    r = spectral_S(bins)
    r["mechanism"] = "complete_Higgsing_short_window"
    r["a1_Planck_GUT"] = a1_uv
    r["assumed"] = [
        "full 12 Proca only above M_GUT",
        "SM massless restored below M_GUT",
        "matter content SM-like",
    ]
    return r


def stress_complete_higgsing_to_EW() -> dict:
    """All 12 Proca from M_Pl all the way to EW (#45 scenario).

    ASSUMED: no massless SM gauge bosons until EW (conflicts with QCD/QED
    above Λ_QCD — flagged for W).
    """
    a1_uv = a1_gauge(0, 12) + a1_sm_matter()
    bins = [
        {"name": "all_Proca_Planck_to_EW", "mu_hi": M_PL, "mu_lo": V_EW, "a1": a1_uv},
    ] + low_energy_tail()
    r = spectral_S(bins)
    r["mechanism"] = "complete_Higgsing_to_EW"
    r["assumed"] = [
        "12 Proca from M_Pl→EW",
        "SM massless only below EW (IR conflict — see W)",
    ]
    return r


def stress_so10_heavy() -> dict:
    """Unbroken SO(10) M_Pl→M_GUT; broken→SM (33 Proca+12 massless) M_GUT→v.

    Gauge+fermion only (conservative, #46); 48 Weyl.
    """
    a1_unb = a1_gauge(45, 0) + 48 * A1_W  # −30 + 4 = −26
    a1_brk = a1_gauge(12, 33) + 48 * A1_W  # 17/2 + 4 = 25/2
    bins = [
        {"name": "SO10_unbroken", "mu_hi": M_PL, "mu_lo": M_GUT, "a1": a1_unb},
        {"name": "SO10_broken_SM", "mu_hi": M_GUT, "mu_lo": V_EW, "a1": a1_brk},
    ] + low_energy_tail()
    r = spectral_S(bins)
    r["mechanism"] = "SO10_heavy_Proca"
    r["a1_unbroken"] = a1_unb
    r["a1_broken"] = a1_brk
    r["assumed"] = ["SO(10) content", "breaking at M_GUT", "no extra Higgs scalars"]
    return r


def stress_e6_heavy() -> dict:
    """E6: dim=78; unbroken M_Pl→M_GUT; 66 Proca+12 massless below.

    Fermions: ASSUMED 3×27 = 81 Weyl (one 27 per gen).
    """
    a1_unb = a1_gauge(78, 0) + 81 * A1_W  # −52 + 27/4 = −181/4
    a1_brk = a1_gauge(12, 66) + 81 * A1_W  # −8 + 33 + 27/4 = 121/4
    bins = [
        {"name": "E6_unbroken", "mu_hi": M_PL, "mu_lo": M_GUT, "a1": a1_unb},
        {"name": "E6_broken_SM", "mu_hi": M_GUT, "mu_lo": V_EW, "a1": a1_brk},
    ] + low_energy_tail()
    r = spectral_S(bins)
    r["mechanism"] = "E6_heavy_Proca"
    r["assumed"] = ["E6 dim=78", "3×27 Weyl", "break at M_GUT", "no Higgs scalars"]
    return r


def stress_kk_extras() -> dict:
    """SM zero modes massless everywhere; +12 KK Proca in Planck→GUT only.

    ASSUMED: one KK level active only above M_GUT (truncated tower).
    """
    a1_uv = A1_SM_UV + 12 * A1_PROCA  # SM + 12 extras
    a1_mid = A1_SM_UV
    bins = [
        {"name": "SM_plus_KK_Planck_GUT", "mu_hi": M_PL, "mu_lo": M_GUT, "a1": a1_uv},
        {"name": "SM_GUT_to_EW", "mu_hi": M_GUT, "mu_lo": V_EW, "a1": a1_mid},
    ] + low_energy_tail()
    r = spectral_S(bins)
    r["mechanism"] = "string_KK_level1"
    r["assumed"] = ["exactly 1 KK level", "only above M_GUT", "zero modes massless"]
    return r


def stress_technicolor_n4() -> dict:
    """Minimal n=4 dynamical Proca above M_GUT; SM otherwise.

    ASSUMED: 4 of 12 SM-like vectors massed only in Planck→GUT bin
    (Δ vs SM = 4×7/6 = 14/3 in that bin only).
    """
    # In UV bin: 8 massless + 4 Proca + matter
    a1_uv = a1_gauge(8, 4) + a1_sm_matter()
    a1_mid = A1_SM_UV
    bins = [
        {"name": "TC4_Planck_GUT", "mu_hi": M_PL, "mu_lo": M_GUT, "a1": a1_uv},
        {"name": "SM_GUT_to_EW", "mu_hi": M_GUT, "mu_lo": V_EW, "a1": a1_mid},
    ] + low_energy_tail()
    r = spectral_S(bins)
    r["mechanism"] = "technicolor_n4"
    r["assumed"] = ["exactly 4 dynamical Proca above M_GUT", "restored SM below"]
    return r


def stress_technicolor_n12_short() -> dict:
    """TC-12 ≡ complete Higgsing in short Planck→GUT window only."""
    return {**stress_complete_higgsing_short(), "mechanism": "technicolor_n12_short"}


def all_stress_tests() -> dict:
    tests = [
        baseline_SM(),
        stress_complete_higgsing_short(),
        stress_complete_higgsing_to_EW(),
        stress_so10_heavy(),
        stress_e6_heavy(),
        stress_kk_extras(),
        stress_technicolor_n4(),
        stress_technicolor_n12_short(),
    ]
    flips = [t for t in tests if t["sign_positive"] and t["mechanism"] != "baseline_SM"]
    return {
        "tests": tests,
        "n_flip": len(flips),
        "flip_names": [t["mechanism"] for t in flips],
        "coincidence_flag": (
            "M_GUT=1e16 GeV and Λ_UV=M_Pl are standard ballparks — structural"
        ),
        "plain": (
            "Short Planck→GUT-only Proca (complete/TC-12/KK/TC-4) does NOT "
            "flip S — long SM bin still wins. SO(10)/E6 flip because the "
            "long M_GUT→EW broken bin has many heavy Proca. "
            "Complete Higgsing M_Pl→EW (#45) flips but conflicts with IR "
            "massless γ/g (W)."
        ),
    }


def self_checks() -> None:
    b = baseline_SM()
    assert b["sign_positive"] is False
    assert stress_complete_higgsing_short()["sign_positive"] is False
    assert stress_complete_higgsing_to_EW()["sign_positive"] is True
    assert stress_so10_heavy()["sign_positive"] is True
    assert stress_e6_heavy()["sign_positive"] is True
    assert stress_kk_extras()["sign_positive"] is False
    assert stress_technicolor_n4()["sign_positive"] is False
    a = all_stress_tests()
    assert "complete_Higgsing_to_EW" in a["flip_names"]
    assert "SO10_heavy_Proca" in a["flip_names"]
    assert "string_KK_level1" not in a["flip_names"]
    assert "complete_Higgsing_short_window" not in a["flip_names"]


def main():
    print("=== Planck→GUT stress test for U survivors ===\n")
    a = all_stress_tests()
    for t in a["tests"]:
        print(
            f"  {t['mechanism']:28s}  S={t['S']:+8.1f}  "
            f"sign+={t['sign_positive']}  assumed={t.get('assumed')}"
        )
    print(f"\n  Flip ({a['n_flip']}): {a['flip_names']}")
    print(f"  coincidence: {a['coincidence_flag']}")
    print(f"  {a['plain']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
Short Planck→GUT Proca windows fail (SM bin dominates). SO(10)/E6
survive via long broken-to-SM Proca bin. Complete Higgsing M_Pl→EW
flips arithmetically but needs IR consistency (W).
"""
    )
    print("=== PLANCK STRESS COMPLETE ===")


if __name__ == "__main__":
    main()
