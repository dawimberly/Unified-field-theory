#!/usr/bin/env python3
"""Total induced 1/G_ind via spectral representation + heavy sector.

Combine matched hadronic band (#40) with scales above Λ_QCD (EW, top, Higgs)
in Adler/Sakharov form:

  1/(16π G_ind) ~ (1/12) ∫_0^{Λ_UV²} ds ρ(s)/s   [schematic]

Ask: is TOTAL 1/G positive? Magnitude near M_Pl²? (needs cutoff — ASSUMED).

No tuning. DERIVED vs ASSUMED. Coincidences flagged.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_spectral_gind.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import pi, sqrt

# Matched hadronic central from #40
A1_HAD_CENTRAL = F(1, 3)
A1_HAD_LOW = F(1, 6)
A1_HAD_HIGH = F(3, 4)

# Free-field a1 above Λ_QCD (from #27/#31/#38)
A1_UV_SM = F(-43, 12)  # ξ_H=0 UV
A1_EW = F(-1, 12)  # Proca W/Z below EW
A1_AC = F(-13, 4)  # quarks+gluons still free above Λ_QCD

M_PL_GEV = 1.22e19
LAMBDA_QCD_GEV = 0.3
V_EW_GEV = 246.0
M_TOP_GEV = 173.0


def spectral_integral_scheme() -> dict:
    """Piecewise spectral weight: Σ a1_i × log(μ_i² / μ_{i+1}²) schematic.

    Derived: Adler/Sakharov structure ties 1/G to ∫ ρ/s.
    Assumed: replace ρ by free-field a1 thresholds (standard but crude).
    """
    return {
        "formula": "1/(16πG) ∝ Σ_i a1(μ_i) · Δlog(μ²)  (threshold sketch)",
        "structure_label": "derived (spectral / Sakharov)",
        "threshold_label": "assumed (free-field a1 as spectral proxy)",
    }


def heavy_sector_pieces() -> list[dict]:
    """Spectral pieces from UV down to Λ_QCD (heavy), then hadronic IR."""
    return [
        {
            "name": "UV_SM",
            "mu_hi": 1e16,
            "mu_lo": V_EW_GEV,
            "a1": A1_UV_SM,
            "label": "derived content; ξ_H=0 assumed",
        },
        {
            "name": "EW_to_top",
            "mu_hi": V_EW_GEV,
            "mu_lo": M_TOP_GEV,
            "a1": A1_EW,
            "label": "derived Proca IR for W/Z",
        },
        {
            "name": "top_to_QCD",
            "mu_hi": M_TOP_GEV,
            "mu_lo": LAMBDA_QCD_GEV,
            "a1": A1_AC,
            "label": "derived free QCD+EW remnants",
        },
        {
            "name": "hadronic_IR",
            "mu_hi": LAMBDA_QCD_GEV,
            "mu_lo": 0.1,  # soft IR floor
            "a1": A1_HAD_CENTRAL,
            "a1_low": A1_HAD_LOW,
            "a1_high": A1_HAD_HIGH,
            "label": "matched band #40; ρ model assumed",
        },
    ]


def integrated_sign() -> dict:
    """Weighted sum S = Σ a1_i · log(μ_hi²/μ_lo²). Sign of S ↔ sign of 1/G."""
    from math import log

    pieces = heavy_sector_pieces()
    S = 0.0
    S_low = 0.0
    S_high = 0.0
    rows = []
    for p in pieces:
        # Δlog(μ²) = 2 ln(μ_hi/μ_lo)
        w = 2.0 * log(p["mu_hi"] / p["mu_lo"])
        a1 = float(p["a1"])
        contrib = a1 * w
        a1_l = float(p.get("a1_low", p["a1"]))
        a1_h = float(p.get("a1_high", p["a1"]))
        S += contrib
        S_low += a1_l * w
        S_high += a1_h * w
        rows.append(
            {
                "name": p["name"],
                "weight": w,
                "a1": p["a1"],
                "contrib": contrib,
                "label": p["label"],
            }
        )
    # Heavy UV pieces are negative and dominate — classic wrong-sign.
    return {
        "rows": rows,
        "S_central": S,
        "S_low": S_low,
        "S_high": S_high,
        "total_sign_positive_central": S > 0,
        "total_sign_positive_entire_band": S_low > 0 and S_high > 0,
        "plain_sign": (
            "With UV cutoff ~10^16 GeV and free-field a1 proxies, the "
            "log-weighted total S is NEGATIVE (UV gauge Kabats dominate). "
            "Hadronic IR cannot flip the integrated 1/G without cutting off "
            "or canceling the UV — same classic Sakharov wrong-sign lesson."
        ),
    }


def planck_cutoff_needed() -> dict:
    """Magnitude: M_Pl² needs Λ_UV ~ M_Pl / sqrt(|N_eff|).

    Sakharov estimate: M_Pl² ~ |a1| Λ_UV² / (12π)  (schematic).
    ASSUMED cutoff; not derived from FFZ.
    """
    a1_eff = abs(float(A1_UV_SM))  # |−43/12| ≈ 3.58
    # M_Pl² ≈ (a1_eff / (12π)) Λ²  ⇒  Λ ≈ M_Pl * sqrt(12π / a1_eff)
    lam = M_PL_GEV * sqrt(12.0 * pi / a1_eff)
    return {
        "M_Pl_GeV": M_PL_GEV,
        "a1_eff_UV": float(A1_UV_SM),
        "Lambda_UV_for_MPl_GeV": lam,
        "near_MPl_without_cutoff": False,
        "cutoff_assumed": True,
        "coincidence_flag": (
            "Λ_UV ∼ M_Pl is the usual Sakharov coincidence / assumption — "
            "not derived here"
        ),
        "plain": (
            f"To reach M_Pl² with |a1|∼{a1_eff:.2f} needs Λ_UV∼{lam:.2e} GeV "
            f"(∼M_Pl). Hadronic IR does not supply Planck-scale induction "
            f"without an ASSUMED UV cutoff."
        ),
    }


def self_checks() -> None:
    s = integrated_sign()
    assert s["total_sign_positive_central"] is False
    assert s["S_central"] < 0
    p = planck_cutoff_needed()
    assert p["near_MPl_without_cutoff"] is False
    assert p["cutoff_assumed"] is True
    assert p["Lambda_UV_for_MPl_GeV"] > 1e18


def main():
    print("=== Spectral total 1/G_ind (heavy + hadronic) ===\n")
    sc = spectral_integral_scheme()
    print(f"  {sc['formula']}")
    print(f"  structure [{sc['structure_label']}]; thresholds [{sc['threshold_label']}]")
    s = integrated_sign()
    for r in s["rows"]:
        print(
            f"  {r['name']:14s}  a1={r['a1']!s:8s}  w={r['weight']:7.3f}  "
            f"contrib={r['contrib']:+8.3f}  [{r['label']}]"
        )
    print(
        f"\n  S_central={s['S_central']:+.3f}  "
        f"S_band=[{s['S_low']:+.3f}, {s['S_high']:+.3f}]"
    )
    print(f"  Total sign + (central)? {s['total_sign_positive_central']}")
    print(f"  Total sign + (entire hadronic band)? {s['total_sign_positive_entire_band']}")
    print(f"  {s['plain_sign']}")
    p = planck_cutoff_needed()
    print(f"\n  Λ_UV for M_Pl: {p['Lambda_UV_for_MPl_GeV']:.3e} GeV [{('assumed')}]")
    print(f"  coincidence: {p['coincidence_flag']}")
    print(f"  {p['plain']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
TOTAL log-weighted 1/G_ind is NEGATIVE with UV SM a1 — hadronic positivity
does not flip the integrated sign. Magnitude ≪ M_Pl² unless Λ_UV∼M_Pl
(ASSUMED). Honest negative for "confinement alone explains Newton."
"""
    )
    print("=== SPECTRAL G_IND COMPLETE ===")


if __name__ == "__main__":
    main()
