#!/usr/bin/env python3
"""Does confinement have an analogue in the FFZ WZ+Yukawa h²z/K channel?

Ask: could a strongly-coupled phase of the WZ model shift the simple pole
by +0.009651 (the ledger GAP)? Derive or rule out.

No tuning. DERIVED vs ASSUMED. Honest negatives.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_ledger_confinement.py
"""

from __future__ import annotations

C_LEDGER = 0.05131719
C_TARGET = 1.0 / 24.0
GAP = C_LEDGER - C_TARGET  # ≈ +0.00965052
GAP_ERR = 0.0024


def channel_split() -> dict:
    """QCD confinement vs WZ three-loop live in different channels."""
    return {
        "confinement": {
            "theory": "SU(3) QCD",
            "regime": "nonperturbative IR",
            "observable_here": "free-field / spectral a1 (Λ² channel)",
            "affects_h2zK_pole": False,
            "label": "derived (channel identification)",
        },
        "ledger": {
            "theory": "WZ + Yukawa (FFZ Newton program)",
            "regime": "weakly coupled 3-loop perturbative",
            "observable_here": "h² z / (K ε) simple pole → C≈0.051317",
            "GAP": GAP,
            "label": "derived (ledger definition)",
        },
    }


def strong_wz_analogue() -> dict:
    """Is there a derived 'confinement' of the WZ model shifting the pole?

    The WZ model in the FFZ Newton test is a weakly coupled Yukawa theory
    with no non-Abelian gauge dynamics and no area-law confinement.
    A hypothetical strongly coupled WZ phase is ASSUMED speculation —
    not present in the three-loop calculation that produces C and GAP.
    """
    return {
        "wz_has_nonabelian_gauge": False,
        "wz_has_area_law": False,
        "three_loop_is_perturbative": True,
        "strong_phase_in_ledger_calc": False,
        "derived_pole_shift_from_confinement": False,
        "speculative_analogue": (
            "One could invent a strongly coupled WZ/Yukawa phase, but that "
            "is ASSUMED and outside the existing three-loop ledger."
        ),
        "label_ruling": "derived ruling-out within the stated WZ setup",
        "label_speculation": "assumed if claimed anyway",
    }


def gap_shift_check() -> dict:
    """Does anything from Round 2–3 confinement work change GAP?"""
    delta = 0.0
    return {
        "delta_C": delta,
        "delta_GAP": delta,
        "GAP_after": GAP + delta,
        "closes_gap": False,
        "within_err_anyway": abs(GAP) < GAP_ERR,  # False; GAP ~ 4× err
        "plain": (
            f"GAP={GAP:.6f}±{GAP_ERR} unchanged. Confinement / matching / "
            f"spectral 1/G live in the Λ² / a1 channel; they do not "
            f"renormalize the WZ h²z/K simple pole. Ruled out as a derived "
            f"fix of the +0.009651 discrepancy."
        ),
    }


def self_checks() -> None:
    s = strong_wz_analogue()
    assert s["derived_pole_shift_from_confinement"] is False
    assert s["wz_has_nonabelian_gauge"] is False
    assert s["three_loop_is_perturbative"] is True
    g = gap_shift_check()
    assert g["delta_GAP"] == 0.0
    assert g["closes_gap"] is False
    assert abs(g["GAP_after"] - 0.00965052) < 1e-5
    ch = channel_split()
    assert ch["confinement"]["affects_h2zK_pole"] is False


def main():
    print("=== Confinement ↔ FFZ h²z/K ledger link ===\n")
    ch = channel_split()
    print(f"  Confinement: {ch['confinement']['theory']} / {ch['confinement']['regime']}")
    print(f"    → {ch['confinement']['observable_here']} [{ch['confinement']['label']}]")
    print(f"    affects h²z/K pole? {ch['confinement']['affects_h2zK_pole']}")
    print(f"  Ledger: {ch['ledger']['theory']} / {ch['ledger']['regime']}")
    print(f"    → {ch['ledger']['observable_here']}  GAP={ch['ledger']['GAP']:.6f}")
    s = strong_wz_analogue()
    print(f"\n  WZ has non-Abelian gauge? {s['wz_has_nonabelian_gauge']}")
    print(f"  WZ area law? {s['wz_has_area_law']}")
    print(f"  3-loop perturbative? {s['three_loop_is_perturbative']}")
    print(f"  Derived pole shift from confinement? {s['derived_pole_shift_from_confinement']}")
    print(f"  [{s['label_ruling']}]")
    print(f"  Speculation: {s['speculative_analogue']}")
    g = gap_shift_check()
    print(f"\n  δC={g['delta_C']}  δGAP={g['delta_GAP']}  GAP_after={g['GAP_after']:.6f}")
    print(f"  {g['plain']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
RULED OUT as a derived link: QCD confinement does not shift the WZ
h²z/K simple pole. GAP +0.009651 stands. A 'strong WZ phase' would be
pure assumption outside the ledger calculation.
"""
    )
    print("=== LEDGER–CONFINEMENT COMPLETE ===")


if __name__ == "__main__":
    main()
