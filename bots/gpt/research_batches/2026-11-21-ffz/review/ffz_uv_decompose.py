#!/usr/bin/env python3
"""Why integrated S~-265 (#41) is negative while IR a1 is positive.

Decompose by scale; identify dominant negative UV sector; quantify
minimal physically motivated UV content (not tuning) that would flip
the integral — compare Visser / Adler ⟨TT⟩ / asymptotic safety / ξ.

No tuning. DERIVED vs ASSUMED. Coincidences flagged.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_uv_decompose.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import log

A1_UV_SM = F(-43, 12)
A1_EW = F(-1, 12)
A1_AC = F(-13, 4)
A1_HAD = F(1, 3)
A1_V = F(-2, 3)
A1_PROCA = F(1, 2)
A1_W = F(1, 12)
A1_S = F(1, 6)

LAMBDA_QCD = 0.3
V_EW = 246.0
M_TOP = 173.0
MU_UV = 1e16


def pieces() -> list[dict]:
    return [
        {"name": "UV_SM", "mu_hi": MU_UV, "mu_lo": V_EW, "a1": A1_UV_SM},
        {"name": "EW_to_top", "mu_hi": V_EW, "mu_lo": M_TOP, "a1": A1_EW},
        {"name": "top_to_QCD", "mu_hi": M_TOP, "mu_lo": LAMBDA_QCD, "a1": A1_AC},
        {"name": "hadronic_IR", "mu_hi": LAMBDA_QCD, "mu_lo": 0.1, "a1": A1_HAD},
    ]


def decompose() -> dict:
    rows = []
    S = 0.0
    for p in pieces():
        w = 2.0 * log(p["mu_hi"] / p["mu_lo"])
        contrib = float(p["a1"]) * w
        S += contrib
        rows.append(
            {
                "name": p["name"],
                "a1": p["a1"],
                "weight": w,
                "contrib": contrib,
                "frac_of_abs_S": None,  # filled below
            }
        )
    abs_sum = sum(abs(r["contrib"]) for r in rows)
    for r in rows:
        r["frac_of_abs_S"] = abs(r["contrib"]) / abs_sum if abs_sum else 0.0
    # Dominant negative
    neg = [r for r in rows if r["contrib"] < 0]
    dominant = min(neg, key=lambda r: r["contrib"])
    # UV gauge Kabat decomposition of SM a1
    n_gauge = 12
    a1_gauge = n_gauge * A1_V  # −8
    a1_ferm = 45 * A1_W  # +45/12 = +15/4
    a1_scalar = 4 * A1_S  # +2/3 at ξ=0
    assert a1_gauge + a1_ferm + a1_scalar == A1_UV_SM
    return {
        "rows": rows,
        "S": S,
        "dominant_negative": dominant["name"],
        "dominant_contrib": dominant["contrib"],
        "uv_gauge_a1": a1_gauge,
        "uv_fermion_a1": a1_ferm,
        "uv_scalar_a1": a1_scalar,
        "plain": (
            f"S≈{S:.1f}. Dominant negative: {dominant['name']} "
            f"({dominant['contrib']:.1f}). Within UV SM, gauge Kabats "
            f"a1={a1_gauge} dominate over fermions (+{a1_ferm}) and "
            f"scalars (+{a1_scalar})."
        ),
    }


def minimal_flip() -> dict:
    """What minimal UV change flips S>0? Not tuning — physically motivated options.

    Options compared (literature):
      A. Remove / Higgs all 12 massless Kabats → Proca (+1/2 each): Δa1 = 12×(1/2−(−2/3))=14
      B. Conformal ξ=1/6 on all scalars (Visser / conformal): kills scalar +; doesn't flip alone
      C. Adler ⟨TT⟩ sum-rule positivity on spectral density (assumed UV completion)
      D. Asymptotic safety: UV a1 effectively screened (assumed)
      E. Extra positive species (large N_f or scalars at ξ=0)

    Compute Δa1 needed in the UV_SM bin alone to make S≥0.
    """
    d = decompose()
    S = d["S"]
    uv = next(r for r in d["rows"] if r["name"] == "UV_SM")
    w_uv = uv["weight"]
    # Need S + δa1 * w_uv ≥ 0  ⇒  δa1 ≥ −S / w_uv
    delta_a1_needed = -S / w_uv
    a1_uv_target = float(A1_UV_SM) + delta_a1_needed

    # Option A: all 12 vectors → Proca
    delta_proca_12 = 12 * (A1_PROCA - A1_V)  # 12*(1/2+2/3)=12*(7/6)=14
    assert delta_proca_12 == 14
    S_if_proca12 = S + float(delta_proca_12) * w_uv

    # Option B: ξ_H=1/6 only on 4 Higgs → Δa1 = −4/6 = −2/3 (more negative!)
    delta_xi = 4 * (0 - A1_S)  # conformal vs minimal: each loses 1/6
    assert delta_xi == F(-2, 3)
    S_if_xi = S + float(delta_xi) * w_uv

    # Option E: how many extra minimal scalars needed in UV bin
    # δa1 = n_s / 6 ≥ delta_a1_needed ⇒ n_s ≥ 6 * delta_a1_needed
    n_scalars_needed = 6.0 * delta_a1_needed
    n_weyl_needed = 12.0 * delta_a1_needed  # each Weyl +1/12

    return {
        "delta_a1_UV_needed": delta_a1_needed,
        "a1_UV_target": a1_uv_target,
        "option_A_all_Proca_12": {
            "delta_a1": float(delta_proca_12),
            "S_after": S_if_proca12,
            "flips": S_if_proca12 > 0,
            "label": "derived arithmetic; physical motivation = complete Higgsing (ASSUMED spectrum)",
        },
        "option_B_conformal_Higgs": {
            "delta_a1": float(delta_xi),
            "S_after": S_if_xi,
            "flips": S_if_xi > 0,
            "label": "derived; makes UV more negative — does NOT flip",
        },
        "option_C_Adler_sum_rule": {
            "flips": "unknown",
            "label": "assumed UV completion with positive ⟨TT⟩ spectral weight",
            "refs": "Adler induced gravity; Visser gr-qc/0204062",
        },
        "option_D_asymptotic_safety": {
            "flips": "unknown",
            "label": "assumed UV screening of a1; not computed here",
        },
        "option_E_extra_species": {
            "n_minimal_scalars": n_scalars_needed,
            "n_weyl": n_weyl_needed,
            "label": "derived count; adding O(10²) scalars is ASSUMED BSM — flag as not SM",
        },
        "coincidence_flag": (
            "δa1_needed ≈ |a1_UV| scale — structural from log(M_UV/v) weight, not miraculous"
        ),
        "plain": (
            f"Need Δa1_UV ≳ {delta_a1_needed:.2f} (target a1_UV≳{a1_uv_target:.2f}). "
            f"Higgsing all 12 vectors to Proca (Δ=+14) DOES flip S "
            f"(→{S_if_proca12:.1f}). Conformal Higgs alone does not. "
            f"Extra species need ~{n_scalars_needed:.0f} minimal scalars "
            f"or ~{n_weyl_needed:.0f} Weyl — not SM content."
        ),
    }


def self_checks() -> None:
    d = decompose()
    assert d["S"] < 0
    assert d["dominant_negative"] == "UV_SM"
    assert d["uv_gauge_a1"] == F(-8)
    m = minimal_flip()
    assert m["option_A_all_Proca_12"]["flips"] is True
    assert m["option_B_conformal_Higgs"]["flips"] is False
    assert m["delta_a1_UV_needed"] > 0


def main():
    print("=== UV decompose S (#41) + minimal flip ===\n")
    d = decompose()
    for r in d["rows"]:
        print(
            f"  {r['name']:14s}  a1={r['a1']!s:8s}  w={r['weight']:7.3f}  "
            f"contrib={r['contrib']:+8.3f}  (|S| share {100*r['frac_of_abs_S']:.1f}%)"
        )
    print(f"\n  S={d['S']:+.3f}")
    print(f"  Dominant negative: {d['dominant_negative']} ({d['dominant_contrib']:+.3f})")
    print(
        f"  UV SM split: gauge={d['uv_gauge_a1']}  ferm={d['uv_fermion_a1']}  "
        f"scalar={d['uv_scalar_a1']}"
    )
    print(f"  {d['plain']}")
    m = minimal_flip()
    print(f"\n  Δa1_UV needed to flip: {m['delta_a1_UV_needed']:.3f}")
    print(
        f"  A all-Proca-12: Δ={m['option_A_all_Proca_12']['delta_a1']:+.1f}  "
        f"S→{m['option_A_all_Proca_12']['S_after']:+.1f}  "
        f"flip={m['option_A_all_Proca_12']['flips']}  [{m['option_A_all_Proca_12']['label']}]"
    )
    print(
        f"  B conformal Higgs: Δ={m['option_B_conformal_Higgs']['delta_a1']:+.3f}  "
        f"flip={m['option_B_conformal_Higgs']['flips']}  [{m['option_B_conformal_Higgs']['label']}]"
    )
    print(f"  C Adler ⟨TT⟩: [{m['option_C_Adler_sum_rule']['label']}]")
    print(f"  D asymp. safety: [{m['option_D_asymptotic_safety']['label']}]")
    print(
        f"  E extra species: ~{m['option_E_extra_species']['n_minimal_scalars']:.0f} scalars "
        f"or ~{m['option_E_extra_species']['n_weyl']:.0f} Weyl  [{m['option_E_extra_species']['label']}]"
    )
    print(f"  coincidence: {m['coincidence_flag']}")
    print(f"  {m['plain']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
Negative S is dominated by UV SM gauge Kabats (12×−2/3). IR positivity
is a rounding error on the log-weighted integral. Minimal flip that is
arithmetically clean: Higgs all 12 vectors (Proca) — physical content
ASSUMED. Conformal ξ alone fails. Extra SM-like species need O(10²).
"""
    )
    print("=== UV DECOMPOSE COMPLETE ===")


if __name__ == "__main__":
    main()
