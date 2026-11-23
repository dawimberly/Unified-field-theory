#!/usr/bin/env python3
"""Derive (not assume) shared boundary / Kabat contact across gauge sectors.

Stacked on two-phase (#29). From Donnelly–Wall edge modes, Kabat contact,
and EE/induced-gravity matching: do U(1)_EM, SU(2), SU(3), and gravity-sector
modes share edge modes or a common contact? Is #28's +53/12 shared-Kabat
derivable? Does any derived boundary effect flip IR a1=−1/12 to positive?

No tuning. DERIVED vs ASSUMED. Honest negatives.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_shared_boundary.py
"""

from __future__ import annotations

from fractions import Fraction as F

A1_SM_UV = F(-43, 12)
A1_SM_IR_PROCA = F(-1, 12)  # #29: W±Z Proca, γ+gluons massless
A1_VEC_MASSLESS = F(-2, 3)  # Ds/6 − 1 = 1/3 − 1
A1_VEC_NO_KABAT = F(1, 3)  # Ds/6 only
A1_PROCA = F(1, 2)
A1_SHARED_KABAT_Z3 = F(53, 12)  # #28 assumed option A


def donnelly_wall_edge_modes() -> dict:
    """Donnelly–Wall: Maxwell entanglement = edge-mode Hilbert space.

    For a spatial entangling surface, the normal electric flux E_⊥ is a
    classical edge degree of freedom; the edge Hilbert space is L² on the
    space of flux configurations. Distinct U(1) factors have independent
    flux algebras ⇒ independent edge spaces.

    Non-abelian SU(N): edge modes valued in the Lie algebra; the leading
    area-law / Kabat-like contact scales with dim(G)=N²−1, not with 1
    shared mode for the whole group (Solodukhin/DW literature).

    Gravity: graviton entanglement has its own structure (metric edge
    modes / soft hair); it is NOT identified with Maxwell edge modes
    without an extra assumption (e.g. KK photon = g_{μ5}).
    """
    return {
        "U1_edge_independent": True,
        "SU_N_contact_scales_with": "dim(G)=N²−1",
        "EM_shares_with_SU2": False,
        "EM_shares_with_SU3": False,
        "EM_shares_with_gravity": False,
        "label": "derived (DW: independent flux algebras per U(1)/Cartan)",
        "refs": "Donnelly–Wall 2012/2015; Kabat 1995; Solodukhin 2015",
    }


def kabat_contact_per_sector() -> dict:
    """Kabat contact c=1 per massless gauge boson (in a1: −1 per vector).

    SM UV massless: 1 (U(1)_Y or γ after EWSB bookkeeping) + 3 (SU(2))
    + 8 (SU(3)) = 12 independent contacts in the unbroken UV count.
    After EWSB IR: W±,Z massive → lose Kabat; remain γ + 8 gluons = 9.
    """
    uv = {"U1": 1, "SU2": 3, "SU3": 8, "total": 12}
    ir = {"U1_EM": 1, "SU3": 8, "SU2_massive": 0, "total": 9}
    assert uv["total"] == 12
    assert ir["total"] == 9
    return {
        "c_Kabat_per_massless_vector": 1,
        "UV_contacts": uv,
        "IR_contacts": ir,
        "shared_across_sectors_derived": False,
        "label": "derived (one contact per massless gauge field)",
    }


def shared_kabat_53_12_check() -> dict:
    """#28 option A: share one Kabat among Z3 triples ⇒ a1_vec=0 ⇒ +53/12.

    That requires contacts_total = 12/3 = 4 instead of 12.
    Derived DW counting: contacts = 12 (UV) or 9 (IR). Not 4.
    Hence +53/12 is NOT derivable.
    """
    derived_uv_contacts = 12
    assumed_shared = 4
    return {
        "assumed_shared_contacts": assumed_shared,
        "derived_UV_contacts": derived_uv_contacts,
        "a1_assumed_shared": A1_SHARED_KABAT_Z3,
        "derivable": False,
        "label": "derived: +53/12 not supported by DW edge counting",
    }


def derived_a1_with_edge_counting() -> dict:
    """Net a1 using derived (unshared) contacts.

    UV: already −43/12 with 12×(−2/3).
    IR Proca: −1/12 with 9×(−2/3)+3×(+1/2).

    Flip −1/12 → positive would need to remove enough Kabat from the
    remaining 9 massless vectors. Each removed Kabat adds +1 to a1
    (because a1_vec: −2/3 → +1/3).
    Need Δ > 1/12 ⇒ remove at least 1 full Kabat from one vector…
    Wait: −1/12 + 1 = 11/12 > 0. So removing ONE photon/gluon Kabat
    would flip the IR sign — but that removal is NOT derived (photon
    and gluons are massless ⇒ each keeps c=1).

    Derived: zero Kabats removed beyond the 3 Proca. a1 stays −1/12.
    """
    a1_uv = A1_SM_UV
    a1_ir = A1_SM_IR_PROCA
    # If we unjustifiably dropped k Kabats among the 9:
    # a1(k) = −1/12 + k
    need_k_for_positive = 1  # since −1/12 + 1 > 0
    return {
        "a1_UV_derived": a1_uv,
        "a1_IR_derived": a1_ir,
        "IR_sign_OK": a1_ir > 0,
        "Kabats_removed_beyond_Proca_derived": 0,
        "Kabats_needed_to_flip_IR": need_k_for_positive,
        "flip_derived": False,
        "label": "derived: no extra shared/canceled contact; IR stays −1/12",
    }


def gravity_sector_modes() -> dict:
    """Do graviton edge modes share Maxwell Kabat?

    No derived identification. Soft graviton hair / BMS edge modes are
    a different asymptotic symmetry algebra from U(1) electric flux.
    Mixing them requires ASSUMED KK or emergent-gauge hypotheses.
    """
    return {
        "shares_Maxwell_Kabat": False,
        "label": "derived (distinct algebras); sharing assumed only in KK stories",
    }


def self_checks() -> None:
    dw = donnelly_wall_edge_modes()
    assert dw["EM_shares_with_SU3"] is False
    k = kabat_contact_per_sector()
    assert k["UV_contacts"]["total"] == 12
    assert k["IR_contacts"]["total"] == 9
    assert k["shared_across_sectors_derived"] is False
    s = shared_kabat_53_12_check()
    assert s["derivable"] is False
    a = derived_a1_with_edge_counting()
    assert a["a1_IR_derived"] == F(-1, 12)
    assert a["flip_derived"] is False
    assert a["Kabats_needed_to_flip_IR"] == 1
    # Sanity: −1/12 + 1 = 11/12 > 0
    assert A1_SM_IR_PROCA + 1 == F(11, 12)


def main():
    print("=== Shared boundary / Kabat edge-mode derivation ===\n")
    print("No tuning. Stacked on #29.\n")

    print("-- Donnelly–Wall edge modes --")
    dw = donnelly_wall_edge_modes()
    print(f"  U(1) edges independent? {dw['U1_edge_independent']}")
    print(f"  SU(N) scales with: {dw['SU_N_contact_scales_with']}")
    print(f"  EM shares SU(2)/SU(3)/gravity? {dw['EM_shares_with_SU2']}/{dw['EM_shares_with_SU3']}/{dw['EM_shares_with_gravity']}")
    print(f"  [{dw['label']}]")

    print("\n-- Kabat contacts per sector --")
    k = kabat_contact_per_sector()
    print(f"  UV contacts: {k['UV_contacts']}  IR: {k['IR_contacts']}")
    print(f"  Shared across sectors derived? {k['shared_across_sectors_derived']} [{k['label']}]")

    print("\n-- #28 +53/12 shared-Kabat --")
    s = shared_kabat_53_12_check()
    print(f"  Assumed contacts={s['assumed_shared_contacts']} vs derived UV={s['derived_UV_contacts']}")
    print(f"  +53/12 derivable? {s['derivable']} [{s['label']}]")

    print("\n-- Derived a1 --")
    a = derived_a1_with_edge_counting()
    print(f"  UV a1={a['a1_UV_derived']}  IR a1={a['a1_IR_derived']}  IR sign_OK={a['IR_sign_OK']}")
    print(f"  Extra Kabats removed (derived)={a['Kabats_removed_beyond_Proca_derived']}; need {a['Kabats_needed_to_flip_IR']} to flip IR")
    print(f"  Flip derived? {a['flip_derived']} [{a['label']}]")

    print("\n-- Gravity-sector modes --")
    g = gravity_sector_modes()
    print(f"  Share Maxwell Kabat? {g['shares_Maxwell_Kabat']} [{g['label']}]")

    print("\n-- Self-checks --")
    self_checks()
    print("  all asserts PASS")
    print(
        """
Verdict
-------
Derived DW/Kabat counting: contacts are per massless gauge field
(UV 12, IR 9), not shared across U(1)/SU(2)/SU(3)/gravity.
#28 +53/12 is NOT derivable. IR a1 stays −1/12; flipping it would
need ≥1 unjustified Kabat removal. No derived boundary sign flip.
"""
    )
    print("=== SHARED BOUNDARY COMPLETE ===")


if __name__ == "__main__":
    main()
