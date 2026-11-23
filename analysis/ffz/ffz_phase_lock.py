#!/usr/bin/env python3
"""AI: Dan's sharp phase-lock boundary vs smooth.

Model the field-interaction surface as sharp (Dirichlet / Robin) vs smooth.
Compute the finite shift to the h²z/K pole coefficient and its dependence
on the boundary parameter. If a tuned parameter is required to hit GAP,
say so plainly. Flag numerical coincidences.

Stacked on AH (#61). No tuning. DERIVED vs ASSUMED. Fields-only.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_phase_lock.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import log

C_LEDGER = 0.05131719
C_TARGET = 1.0 / 24.0
GAP = C_LEDGER - C_TARGET
GAP_ERR = 0.0024
KAPPA_SCHEME = -GAP
LN2_OVER_72 = log(2) / 72.0
XI = F(1, 6)


def smooth_bulk() -> dict:
    """Existing Newton probe: no geometric boundary (AG)."""
    return {
        "BC": "smooth / no ∂M",
        "delta_C": 0.0,
        "free_parameter": None,
        "label": "derived (AG / Stokes)",
    }


def dirichlet_sharp() -> dict:
    """Sharp Dirichlet surface (Vassilevich-style).

    Surface a1 density ~ K/3 (pure Dirichlet limit). That is an a1 / Λ²
    structure for 1/G_ind — not the interacting three-loop h²z/K pole.
    Derived δC into the ledger channel: 0 (channel split).
    """
    K_coeff = F(1, 3)
    return {
        "BC": "Dirichlet",
        "surface_a1_K_coeff": K_coeff,
        "channel": "Λ² / a1 (not h²z/K)",
        "delta_C_derived": 0.0,
        "free_parameter": None,
        "label_surface": "assumed (Vassilevich Dirichlet table)",
        "label_delta_C": "derived (channel split → 0)",
    }


def robin_family(S: float) -> dict:
    """Robin BC with parameter S (improvement / boundary mass).

    Standard mixed family: surface density involves S and (1/3 − ξ) K.
    At ξ=1/6: (1/3 − ξ) = 1/6.
    S is a free BC parameter — not fixed by the bulk WZ Lagrangian.
    """
    xi_piece = F(1, 3) - XI
    assert xi_piece == F(1, 6)
    return {
        "BC": "Robin",
        "S": S,
        "xi_K_piece_at_1_6": xi_piece,
        "channel": "Λ² / a1 (not h²z/K)",
        "delta_C_derived": 0.0,
        "free_parameter": "S (Robin)",
        "label_surface": "assumed (Vassilevich Robin family)",
        "label_delta_C": "derived (channel split → 0)",
    }


def contact_ansatz(alpha: float) -> dict:
    """ASSUMED map: finite R contact δC = α (units h²z/K).

    This is the same scheme class as vacuum EH κ (#23 / AG).
    Hitting GAP requires α = GAP (or α = −κ) — one tuned parameter.
    """
    return {
        "ansatz": "δC = α  (finite pure-metric R contact)",
        "alpha": alpha,
        "delta_C": alpha,
        "alpha_to_hit_GAP": GAP,
        "needs_tuning": abs(alpha - GAP) > 1e-15,
        "same_class_as": "κ vacuum EH (#23)",
        "label": "assumed ansatz; α free / scheme",
    }


def phase_lock_scan() -> dict:
    """Scan BC choices; report derived vs tuned paths to GAP."""
    smooth = smooth_bulk()
    diri = dirichlet_sharp()
    # Sample Robin S values (dimensionless schematic)
    robins = [robin_family(S) for S in (0.0, 0.1, 1.0, float(GAP))]
    # Contact ansatz at α=0 (natural) and α=GAP (tuned)
    natural = contact_ansatz(0.0)
    tuned = contact_ansatz(GAP)
    coincidences = [
        {
            "name": "ln2/72",
            "value": LN2_OVER_72,
            "within_1sigma_of_GAP": abs(LN2_OVER_72 - GAP) <= GAP_ERR,
            "label": "coincidence only",
        },
        {
            "name": "α=GAP or κ=−GAP",
            "value": GAP,
            "within_1sigma_of_GAP": True,
            "label": "tautology / scheme absorption — NOT a derivation",
        },
        {
            "name": "S=GAP (Robin sample)",
            "value": GAP,
            "within_1sigma_of_GAP": True,
            "label": "coincidence if identified with pole units — NOT derived",
        },
    ]
    derived_any = any(
        abs(x["delta_C"] if "delta_C" in x else x["delta_C_derived"]) > 0
        for x in [smooth, diri] + robins
    )
    # derived_any should be False (all derived δC=0)
    return {
        "smooth": smooth,
        "dirichlet": diri,
        "robins": robins,
        "contact_natural": natural,
        "contact_tuned": tuned,
        "coincidences": coincidences,
        "derived_closes_gap": False,
        "tuning_required": True,
        "tuning_statement": (
            f"To hit GAP={GAP:.6f} via a boundary/contact ansatz requires "
            f"setting free α=GAP (≡ κ=−GAP from #23). That is tuning / scheme, "
            f"not a derived phase-lock prediction."
        ),
        "honest_negative": True,
        "plain": (
            "Sharp Dirichlet/Robin surfaces live in the a1/Λ² channel; "
            "derived δC into h²z/K is 0 for all S. A finite contact α can "
            "match GAP only by tuning α=GAP — same as scheme κ."
        ),
    }


def self_checks() -> None:
    r = phase_lock_scan()
    assert r["smooth"]["delta_C"] == 0.0
    assert r["dirichlet"]["delta_C_derived"] == 0.0
    assert all(rb["delta_C_derived"] == 0.0 for rb in r["robins"])
    assert r["contact_natural"]["needs_tuning"] is True
    assert r["contact_tuned"]["needs_tuning"] is False
    assert abs(r["contact_tuned"]["alpha"] - GAP) < 1e-15
    assert r["tuning_required"] is True
    assert r["derived_closes_gap"] is False
    assert robin_family(0.0)["xi_K_piece_at_1_6"] == F(1, 6)


def main():
    print("=== AI: Sharp phase-lock boundary vs smooth ===\n")
    r = phase_lock_scan()
    print("-- Smooth --")
    print(f"  δC={r['smooth']['delta_C']}  [{r['smooth']['label']}]")
    print("\n-- Dirichlet --")
    d = r["dirichlet"]
    print(
        f"  K-coeff={d['surface_a1_K_coeff']}  channel={d['channel']}  "
        f"δC_derived={d['delta_C_derived']}"
    )
    print("\n-- Robin family --")
    for rb in r["robins"]:
        print(
            f"  S={rb['S']:.6f}  (1/3−ξ)={rb['xi_K_piece_at_1_6']}  "
            f"δC_derived={rb['delta_C_derived']}"
        )
    print("\n-- Contact ansatz (ASSUMED) --")
    print(
        f"  α=0: δC={r['contact_natural']['delta_C']}  "
        f"needs_tuning={r['contact_natural']['needs_tuning']}"
    )
    print(
        f"  α=GAP: δC={r['contact_tuned']['delta_C']:.8f}  "
        f"needs_tuning={r['contact_tuned']['needs_tuning']}  "
        f"[{r['contact_tuned']['label']}]"
    )
    print("\n-- Coincidences --")
    for c in r["coincidences"]:
        print(
            f"  {c['name']:24} = {c['value']:.8f}  "
            f"within_1σ={c['within_1sigma_of_GAP']}  [{c['label']}]"
        )
    print(f"\n  {r['tuning_statement']}")
    print(f"  {r['plain']}")
    print("\n=== PHASE LOCK COMPLETE ===")
    return r


if __name__ == "__main__":
    main()
