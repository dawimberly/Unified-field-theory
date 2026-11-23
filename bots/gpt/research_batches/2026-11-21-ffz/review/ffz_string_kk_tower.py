#!/usr/bin/env python3
"""String / KK tower: do massive modes flip Sakharov a1 sign?

Proca +1/2 each; massive spin-2; species bound. Quote known results.

No tuning. DERIVED vs ASSUMED.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_string_kk_tower.py
"""

from __future__ import annotations

from fractions import Fraction as F

A1_PROCA = F(1, 2)
A1_SM_UV = F(-43, 12)


def proca_tower_flip() -> dict:
    """How many Proca modes needed to cancel SM UV −43/12?

    Each Proca adds +1/2. Need N such that −43/12 + N/2 > 0 ⇒ N/2 > 43/12
    ⇒ N > 43/6 ≈ 7.17 ⇒ N ≥ 8 extra Proca (beyond SM) — or convert
    enough massless vectors.
    Converting all 12 SM vectors to Proca: a1 = +125/12 (#29 fantasy).
    A KK tower of N_KK Proca: Δ = N_KK/2.
    """
    N_min = 8
    assert -A1_SM_UV < N_min * A1_PROCA
    assert -A1_SM_UV > (N_min - 1) * A1_PROCA
    return {
        "N_Proca_min_to_flip_UV": N_min,
        "Delta_per_Proca": A1_PROCA,
        "all_12_Proca_a1": F(125, 12),
        "label_arith": "derived",
        "label_tower_exists": "assumed (model-dependent KK/string spectrum)",
    }


def species_bound_note() -> dict:
    """Dvali species bound: M_Pl² ≳ N_sp Λ² / (something) — relates N to weak gravity.

    Induced gravity with many species can give M_Pl² ∼ N Λ² /(16π²).
    Sign still requires the weighted a1 sum positive; large N of *scalars*
    or Proca helps, large N of massless vectors hurts.
    Known: string theory often has correct-sign discussions in AdS/stringy
    induced gravity literature, but no model-independent flip from SM alone.
    """
    return {
        "species_bound": "M_Pl² ≳ N Λ² / 16π² (schematic)",
        "guarantees_positive_a1": False,
        "string_induced_G_model_independent_flip": False,
        "refs": "Dvali species; Sakharov; string-inspired induced gravity reviews",
        "label": "literature schematic; not a derived SM fix",
    }


def massive_spin2() -> dict:
    """Massive graviton / KK graviton: contributes positively to some counts
    (5 dofs) but also ghosts/vDVZ issues — not a free lunch for a1 sign.
    """
    return {
        "dofs": 5,
        "a1_naive": F(5, 6),  # 5×1/6 schematic ASSUMED
        "clean_SM_fix": False,
        "label": "assumed schematic; vDVZ/ghost caveats",
    }


def self_checks() -> None:
    p = proca_tower_flip()
    assert p["N_Proca_min_to_flip_UV"] == 8
    assert species_bound_note()["guarantees_positive_a1"] is False


def main():
    print("=== String/KK tower vs a1 sign ===\n")
    p = proca_tower_flip()
    print(f"  Need ≥{p['N_Proca_min_to_flip_UV']} Proca to flip UV −43/12 [{p['label_arith']}]")
    print(f"  All-12-Proca a1={p['all_12_Proca_a1']} (not SM)")
    s = species_bound_note()
    print(f"  Species bound guarantees +a1? {s['guarantees_positive_a1']} [{s['label']}]")
    m = massive_spin2()
    print(f"  Massive spin-2 naive a1∼{m['a1_naive']} clean fix? {m['clean_SM_fix']}")
    self_checks()
    print("  self-checks PASS")
    print(
        """
Verdict
-------
A KK/Proca tower can flip a1 in arithmetic (≥8 Proca), but existence of
that tower is model-assumed. Species bound relates N to M_Pl, does not
guarantee positive a1. No derived SM-only string fix.
"""
    )
    print("=== STRING/KK COMPLETE ===")


if __name__ == "__main__":
    main()
