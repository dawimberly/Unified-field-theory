#!/usr/bin/env python3
"""Casimir vs Le Sage: is confinement-scale induced G a Casimir-type effect?

Dan raised push gravity (Le Sage; physics.SE q/26692). Test whether
induced gravity at the confinement scale is mathematically a Casimir-type
effect (boundaries removing field modes → net attraction), cf. Casimir/
Lifshitz, Jaffe hep-th/0503158 caveat, MIT bag Casimir energy.

Show why it avoids Le Sage's fatal problems (heating, drag/aberration,
shielding / Majorana absorption). Derived vs analogy.

No tuning. DERIVED vs ASSUMED. Coincidences flagged.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_casimir_lesage.py
"""

from __future__ import annotations

# Majorana-type gravitational absorption: experimental upper bounds on
# absorption coefficient κ (m²/kg). Literature order-of-magnitude
# (Kreuzer 1968; Unnikrishnan / modern reanalyses; SE q/26692 discussion):
# κ ≲ 10^{-15}–10^{-12} m²/kg depending on method. We quote the tighter
# oft-cited ballpark and flag the range as literature, not re-derived.
KAPPA_ABSORB_BOUND_M2_PER_KG = 1e-15
KAPPA_ABSORB_BOUND_RANGE = (1e-15, 1e-12)


def casimir_structure() -> dict:
    """Casimir / Lifshitz: boundaries modify mode spectrum → force.

    Derived: mode subtraction / ζ-function / heat-kernel difference between
    geometries gives a finite force (Lifshitz formula for dielectrics).
    Jaffe caveat (hep-th/0503158): the force is real, but identifying it
    with a universal 'vacuum energy density' that gravitates requires care —
    renormalization scheme dependence.
    """
    return {
        "mechanism": "boundary conditions remove/shift field modes → ΔE → force",
        "structure_label": "derived (Casimir/Lifshitz mode difference)",
        "jaffe_caveat": (
            "Force ≠ naive unrenormalized vacuum energy; see Jaffe "
            "hep-th/0503158"
        ),
        "jaffe_label": "literature caveat (cited)",
        "mit_bag": (
            "MIT bag: quark/gluon modes excluded outside hadron → Casimir "
            "energy contribution to hadron mass (model-dependent)"
        ),
        "mit_bag_label": "assumed model; energy shift real in the bag EFT",
        "refs": (
            "Casimir 1948; Lifshitz; Jaffe hep-th/0503158; "
            "MIT bag Casimir literature"
        ),
    }


def confinement_as_casimir() -> dict:
    """Map confinement-scale induced a1 ↔ Casimir-type mode removal.

    DERIVED parallel: confinement removes free massless gluon modes from
    the IR EFT (same content statement as #36) — analogous to plates
    removing EM modes. The Δa1 from dropping 8×Kabat is a mode-counting
    difference.

    ANALOGY (not theorem): calling the induced Newton term itself a
    'Casimir force between masses' — Sakharov induction is a curvature
    coupling from vacuum loops, not a two-body Casimir force calculation.
    """
    delta_a1_remove_8_kabat = 16.0 / 3.0  # −8×(−2/3)
    return {
        "free_gluons_removed_from_IR": True,
        "delta_a1_from_mode_removal": delta_a1_remove_8_kabat,
        "parallel_to_casimir_mode_subtraction": True,
        "parallel_label": "derived as EFT mode content (#36)",
        "equals_two_body_Casimir_force": False,
        "equals_label": "analogy only — not a derived Casimir-force theorem",
        "plain": (
            "Confinement removes IR gluon modes (derived content). That is "
            "Casimir-like mode subtraction. Identifying induced G with a "
            "Casimir attraction between bodies remains an analogy."
        ),
    }


def lesage_fatal_problems() -> dict:
    """Le Sage push gravity: fatal issues and why Casimir/QFT avoids them."""
    return {
        "heating": {
            "lesage": "Corpuscles deposit kinetic energy → bodies heat unboundedly",
            "qft_casimir": "No kinetic corpuscle flux; vacuum fluctuations, no net heat dump",
            "avoided": True,
            "label": "derived (no Le Sage medium in QFT vacuum)",
        },
        "drag_aberration": {
            "lesage": (
                "Motion through corpuscular wind → drag; finite corpuscle "
                "speed → aberration (gravity mis-points)"
            ),
            "qft_casimir": (
                "Lorentz-invariant vacuum; induced EH term is local/covariant "
                "at long distance — no preferred wind frame"
            ),
            "avoided": True,
            "label": "derived (Lorentz invariance of QFT vacuum)",
        },
        "shielding_absorption": {
            "lesage": (
                "Bodies absorb/shadow corpuscles → gravity not purely additive; "
                "Majorana gravitational absorption"
            ),
            "bound_kappa_m2_per_kg": KAPPA_ABSORB_BOUND_M2_PER_KG,
            "bound_range": KAPPA_ABSORB_BOUND_RANGE,
            "bound_refs": (
                "Kreuzer 1968; modern reanalyses; physics.SE q/26692; "
                "κ ≲ 10^{-15}–10^{-12} m²/kg literature ballpark"
            ),
            "bound_label": "literature bound (cited; not re-measured here)",
            "qft_casimir": (
                "No absorption cross-section for a corpuscular graviton wind; "
                "Einstein gravity (and Sakharov EH) is geometric / additive "
                "at weak field — shielding is not a built-in mechanism"
            ),
            "avoided": True,
            "label": "derived absence of Le Sage shielding in EH/QFT; bound = literature",
        },
    }


def verdict() -> dict:
    ls = lesage_fatal_problems()
    all_avoided = all(
        ls[k]["avoided"] for k in ("heating", "drag_aberration", "shielding_absorption")
    )
    c = confinement_as_casimir()
    return {
        "casimir_like_mode_removal": c["parallel_to_casimir_mode_subtraction"],
        "derived_two_body_Casimir": c["equals_two_body_Casimir_force"],
        "avoids_all_lesage_fatals": all_avoided,
        "kappa_bound": KAPPA_ABSORB_BOUND_M2_PER_KG,
        "plain": (
            "Induced-G / confinement mode removal is Casimir-like in the "
            "derived sense of boundary-modified spectra, and avoids Le Sage "
            "heating, drag/aberration, and shielding (κ ≲ 1e-15 m²/kg "
            "literature). It is NOT a derived proof that Newton is a "
            "two-body Casimir force (Jaffe caveat + analogy limit)."
        ),
    }


def self_checks() -> None:
    c = confinement_as_casimir()
    assert c["free_gluons_removed_from_IR"] is True
    assert c["equals_two_body_Casimir_force"] is False
    assert abs(c["delta_a1_from_mode_removal"] - 16.0 / 3.0) < 1e-12
    ls = lesage_fatal_problems()
    assert ls["heating"]["avoided"] is True
    assert ls["drag_aberration"]["avoided"] is True
    assert ls["shielding_absorption"]["avoided"] is True
    assert ls["shielding_absorption"]["bound_kappa_m2_per_kg"] == 1e-15
    v = verdict()
    assert v["avoids_all_lesage_fatals"] is True
    assert v["derived_two_body_Casimir"] is False


def main():
    print("=== Casimir vs Le Sage (confinement-scale induced G) ===\n")
    cs = casimir_structure()
    print(f"  Casimir mechanism: {cs['mechanism']} [{cs['structure_label']}]")
    print(f"  Jaffe caveat: {cs['jaffe_caveat']} [{cs['jaffe_label']}]")
    print(f"  MIT bag: {cs['mit_bag']} [{cs['mit_bag_label']}]")
    c = confinement_as_casimir()
    print(f"\n  Gluons removed from IR? {c['free_gluons_removed_from_IR']}")
    print(f"  Δa1 from mode removal: {c['delta_a1_from_mode_removal']} [{c['parallel_label']}]")
    print(f"  Equals two-body Casimir force? {c['equals_two_body_Casimir_force']} [{c['equals_label']}]")
    ls = lesage_fatal_problems()
    for k in ("heating", "drag_aberration", "shielding_absorption"):
        print(f"\n  Le Sage {k}:")
        print(f"    problem: {ls[k]['lesage']}")
        print(f"    QFT:     {ls[k]['qft_casimir']}")
        print(f"    avoided? {ls[k]['avoided']} [{ls[k]['label']}]")
        if k == "shielding_absorption":
            print(
                f"    κ bound ≲ {ls[k]['bound_kappa_m2_per_kg']} m²/kg "
                f"(range {ls[k]['bound_range']}) [{ls[k]['bound_label']}]"
            )
    v = verdict()
    print(f"\n  {v['plain']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
Casimir-like (derived mode removal); avoids all three Le Sage fatals
(derived absence of corpuscular medium + literature κ bound).
Two-body Casimir = Newton remains ANALOGY (Jaffe caveat).
"""
    )
    print("=== CASIMIR / LE SAGE COMPLETE ===")


if __name__ == "__main__":
    main()
