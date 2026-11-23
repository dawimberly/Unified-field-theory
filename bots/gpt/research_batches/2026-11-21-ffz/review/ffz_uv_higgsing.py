#!/usr/bin/env python3
"""Survey UV mechanisms that mass ALL gauge bosons above the GUT scale.

For each: does it remove UV gauge Kabats with Δa1_UV ≥ +4.2 (#45)?
Explicit counting. No tuning. DERIVED vs ASSUMED.

Mechanisms: complete Higgsing; SO(10)/E6 full chains; Stueckelberg;
string/KK towers; asymptotic safety / gauge condensation; technicolor.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_uv_higgsing.py
"""

from __future__ import annotations

from fractions import Fraction as F

A1_V = F(-2, 3)
A1_PROCA = F(1, 2)
DELTA_PER_VECTOR = A1_PROCA - A1_V  # 1/2 + 2/3 = 7/6
assert DELTA_PER_VECTOR == F(7, 6)

# From #45: need Δa1_UV ≳ +4.232 in the UV_SM log bin to flip S
DELTA_A1_NEEDED = 4.232
N_SM_GAUGE = 12


def delta_from_proca(n_massless_to_proca: int) -> F:
    """Δa1 when n massless vectors (−2/3) become Proca (+1/2)."""
    return n_massless_to_proca * DELTA_PER_VECTOR


def mechanism_complete_higgsing() -> dict:
    """All SM (or GUT) gauge bosons massive above Λ_* — no residual massless.

    Arithmetic: 12 SM → Proca ⇒ Δa1 = 14 ≥ 4.2.
    ASSUMED: such a phase exists and the UV spectral bin sees it.
    IR conflict (γ,g massless) deferred to W.
    """
    d = delta_from_proca(N_SM_GAUGE)
    return {
        "name": "complete_Higgsing_all_12",
        "n_vectors_massed": N_SM_GAUGE,
        "delta_a1": d,
        "delta_a1_float": float(d),
        "meets_4_2": float(d) >= DELTA_A1_NEEDED,
        "survives_U": True,
        "label_arithmetic": "derived",
        "label_existence": "assumed (phase / spectrum not in SM)",
        "plain": "Δ=+14 from 12×(7/6); clears +4.2 arithmetically.",
    }


def mechanism_so10_e6_full() -> dict:
    """SO(10)/E6 broken fully to SM: heavy GUT bosons Proca; 12 SM massless remain.

    Relative to unbroken GUT (all massless Kabats), heavy Proca help.
    Relative to SM UV bin (already 12 massless): GUT extras being Proca
    does NOT remove the 12 SM Kabats in the SM-like bin.
    Δa1 vs SM UV (12 massless): 0 for the SM 12; extras are additional
    Proca (+1/2) not present in SM — net a1 change depends on comparison
    baseline.

    Vs SM baseline: adding N_heavy Proca adds +N_heavy/2 (ASSUMED content).
    SO(10): 33 heavy → +33/2 = +16.5; E6 dim=78, broken to SM ⇒ 66 heavy
    → +33. Cleared vs +4.2 IF those states sit in the UV bin.
    Does NOT by itself mass the photon/gluons.
    """
    n_so10_heavy = 45 - 12  # 33
    n_e6_heavy = 78 - 12  # 66
    d_so10 = n_so10_heavy * A1_PROCA  # added Proca vs SM
    d_e6 = n_e6_heavy * A1_PROCA
    # Also: if baseline were unbroken SO(10), Δ from Higgsing heavies:
    d_so10_vs_unbroken = delta_from_proca(n_so10_heavy)
    return {
        "name": "SO10_E6_full_to_SM",
        "n_so10_heavy_Proca": n_so10_heavy,
        "n_e6_heavy_Proca": n_e6_heavy,
        "delta_a1_vs_SM_SO10": d_so10,
        "delta_a1_vs_SM_E6": d_e6,
        "delta_a1_SO10_vs_unbroken": d_so10_vs_unbroken,
        "meets_4_2_vs_SM": float(d_so10) >= DELTA_A1_NEEDED,
        "removes_SM_Kabats": False,
        "survives_U": True,  # adds enough + vs SM via heavy Proca
        "label_arithmetic": "derived (group dims)",
        "label_higgs": "assumed breaking chain / VEVs",
        "plain": (
            "Heavy GUT Proca add +16.5 (SO10) / +33 (E6) vs SM — clears +4.2 "
            "as extra species, but leaves 12 SM Kabats intact."
        ),
    }


def mechanism_stueckelberg() -> dict:
    """Stueckelberg: Abelian vectors can eat axions → massive without Higgs.

    Non-Abelian Stueckelberg is constrained / typically incomplete.
    SM has 1 Abelian (hypercharge) massless above EW; Stueckelberg can
    mass U(1) factors. Max honest SM Abelian count above EW: 1 (B).
    Δa1 ≤ 7/6 ≈ 1.17 < 4.2 — does NOT clear alone.
    Multiple U(1) string factors ASSUMED could add more.
    """
    n_abelian_sm = 1
    d_sm = delta_from_proca(n_abelian_sm)
    n_needed = int(__import__("math").ceil(DELTA_A1_NEEDED / float(DELTA_PER_VECTOR)))
    return {
        "name": "Stueckelberg",
        "n_abelian_SM": n_abelian_sm,
        "delta_a1_SM_U1": d_sm,
        "n_U1_needed_for_4_2": n_needed,
        "meets_4_2_SM_only": float(d_sm) >= DELTA_A1_NEEDED,
        "survives_U": False,  # SM-only insufficient
        "survives_if_many_U1": True,  # ASSUMED string U(1)^n
        "label_SM": "derived (one U(1)_Y)",
        "label_many_U1": "assumed (string/hidden U(1) stack)",
        "plain": (
            f"SM Stueckelberg Δ=+7/6 < 4.2. Need ≥{n_needed} U(1) factors "
            f"massed — ASSUMED beyond SM."
        ),
    }


def mechanism_string_kk() -> dict:
    """String / KK massive vector towers.

    Zero modes of unbroken SM gauge group remain massless (Kabats stay)
    unless the compactification / Wilson lines break them.
    Massive KK levels: each massive vector ~ Proca +1/2.
    Infinite tower: ASSUMED regulated by string scale; species-like +.

    Explicit finite count ASSUMED: first KK level for each of 12 SM
    adjoint directions → 12 Proca extras ⇒ Δ=+6 vs SM zero-mode-only.
    Meets +4.2 as extra species; does not remove zero-mode Kabats.
    Full Wilson-line breaking of all zero modes ≡ complete Higgsing.
    """
    n_kk_level1 = 12  # ASSUMED one massive KK vector per SM gauge dof
    d_kk = n_kk_level1 * A1_PROCA
    d_wilson_full = delta_from_proca(N_SM_GAUGE)  # if zero modes also massed
    return {
        "name": "string_KK_tower",
        "delta_a1_KK_level1_extras": d_kk,
        "delta_a1_if_zero_modes_also_massed": d_wilson_full,
        "meets_4_2_extras": float(d_kk) >= DELTA_A1_NEEDED,
        "removes_zero_mode_Kabats": False,
        "survives_U": True,  # extras clear +4.2; full Wilson = complete Higgsing
        "label_extras": "assumed (1 KK level / truncation)",
        "label_wilson": "assumed (≡ complete breaking)",
        "plain": (
            "KK extras +6 (ASSUMED 1 level) clear +4.2 as species; "
            "zero-mode Kabats remain unless Wilson lines fully break them."
        ),
    }


def mechanism_asymp_safety_condensation() -> dict:
    """Asymptotic safety or gauge-boson condensation.

    No explicit free-field a1 shift derived here. Could screen Kabats
    (ASSUMED). Cannot claim Δa1 ≥ 4.2 with counting.
    """
    return {
        "name": "asymp_safety_or_condensation",
        "delta_a1": None,
        "meets_4_2": "unknown",
        "survives_U": False,  # no explicit counting → does not survive as derived
        "label": "assumed screening; not an explicit a1 count",
        "plain": "No derived Δa1; cannot claim ≥4.2 without a model calculation.",
    }


def mechanism_technicolor() -> dict:
    """Technicolor-like dynamical breaking: condensate eats gauge bosons.

    For N_TC vectors fully dynamical-massed: same Proca Δ=+7/6 each.
    To clear +4.2 need n ≥ 4 (since 3×7/6=3.5 < 4.2; 4×7/6=14/3≈4.67).
    Full SM 12 → Δ=+14. Existence of UV technicolor ASSUMED.
    """
    n_min = 4
    d_min = delta_from_proca(n_min)
    d_12 = delta_from_proca(12)
    assert float(d_min) >= DELTA_A1_NEEDED
    return {
        "name": "technicolor_dynamical",
        "n_min_for_4_2": n_min,
        "delta_a1_min": d_min,
        "delta_a1_all_12": d_12,
        "meets_4_2_min": True,
        "survives_U": True,
        "label_arithmetic": "derived",
        "label_existence": "assumed UV strong sector",
        "plain": f"n≥{n_min} dynamical Proca clears +4.2; n=12 → Δ=+14.",
    }


def survey() -> dict:
    mechs = [
        mechanism_complete_higgsing(),
        mechanism_so10_e6_full(),
        mechanism_stueckelberg(),
        mechanism_string_kk(),
        mechanism_asymp_safety_condensation(),
        mechanism_technicolor(),
    ]
    survivors = [m for m in mechs if m.get("survives_U")]
    return {
        "delta_a1_needed": DELTA_A1_NEEDED,
        "delta_per_vector": DELTA_PER_VECTOR,
        "mechanisms": mechs,
        "survivors": [m["name"] for m in survivors],
        "n_survivors": len(survivors),
        "plain": (
            "Survivors with explicit Δa1≥4.2 path: complete Higgsing, "
            "SO10/E6 heavy Proca extras, string/KK extras (or Wilson), "
            "technicolor n≥4. Stueckelberg SM-only fails. "
            "Asymp. safety/condensation: no derived count."
        ),
    }


def self_checks() -> None:
    assert DELTA_PER_VECTOR == F(7, 6)
    assert float(delta_from_proca(12)) == 14.0
    assert mechanism_complete_higgsing()["meets_4_2"] is True
    assert mechanism_stueckelberg()["meets_4_2_SM_only"] is False
    assert mechanism_stueckelberg()["survives_U"] is False
    assert mechanism_asymp_safety_condensation()["survives_U"] is False
    assert mechanism_technicolor()["n_min_for_4_2"] == 4
    s = survey()
    assert "complete_Higgsing_all_12" in s["survivors"]
    assert "Stueckelberg" not in s["survivors"]
    assert s["n_survivors"] == 4


def main():
    print("=== UV full-mass mechanisms vs Δa1 ≥ +4.2 ===\n")
    print(f"  Need Δa1_UV ≥ {DELTA_A1_NEEDED}  (per vector Δ={DELTA_PER_VECTOR})")
    s = survey()
    for m in s["mechanisms"]:
        print(f"\n  [{m['name']}]")
        print(f"    {m['plain']}")
        print(f"    survives_U? {m.get('survives_U')}")
    print(f"\n  Survivors ({s['n_survivors']}): {s['survivors']}")
    print(f"  {s['plain']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
Four mechanisms clear +4.2 with explicit counting (existence ASSUMED).
Stueckelberg (SM U(1) only) and asymp. safety/condensation do not
survive as derived Δa1 claims. IR massless γ/g deferred to W.
"""
    )
    print("=== UV HIGGSING SURVEY COMPLETE ===")


if __name__ == "__main__":
    main()
