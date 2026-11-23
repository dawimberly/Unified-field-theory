#!/usr/bin/env python3
"""Surface / boundary / contact terms in the induced-R simple pole.

Dan's original idea: gravity emergent from mode phase/frequency/amplitude at
field-interaction surfaces. New-idea rung stacked on the self-consistent loop
(#19). Question: is the residual gap +0.00965±0.0024 a surface/boundary or
contact contribution invisible to the bulk Wick/Schwinger ledger?

Tasks:
  (1) Sakharov spectral / zeta form of the induced R coefficient; split bulk
      vs boundary (GHY-like, □R, total-derivative, vertex δ-contacts) at 3
      loops for WZ+Yukawa, ξ=1/6.
  (2) Check whether ξ=1/6, □R, or coincident-vertex contacts produce a finite
      simple-pole piece the bulk ledger dropped; compute it.
  (3) Compare to GAP; report ratio; derived vs scheme/BC-dependent.
  (4) Optional: Jacobson / entanglement-entropy area-law cross-check for the
      same field content (contact-term ambiguity for nonminimal scalars).

No tuning to 1/24. Derived vs assumed labeled.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_surface_term.py

Requires: sympy, mpmath
"""

from __future__ import annotations

from fractions import Fraction as F
from math import log, pi

import mpmath as mp
import sympy as s

mp.mp.dps = 40

# Standing R* ledger (PR #18)
C_LEDGER = 0.05131719
C_TARGET = 1.0 / 24.0
GAP = C_LEDGER - C_TARGET  # +0.00965052
GAP_ERR = 0.0024
K = 16.0 * pi**2
XI = F(1, 6)


# ---------------------------------------------------------------------------
# (1) Sakharov spectral / zeta form; bulk vs boundary dictionary
# ---------------------------------------------------------------------------

def sakharov_spectral_sum():
    """Induced EH from mode frequency shifts under curvature (Sakharov).

    For Laplace-type modes on a weakly curved background,
      δω_n / ω_n ∼ α_n R / ω_n² ,
      ΔE ∼ ∑_n (1/2) δω_n  ⇒  Γ_ind ⊃ −R/(16π G_ind),
    with 1/G_ind fixed by the spectral sum ∑ α_n/(2ω_n).

    Zeta / heat-kernel transcription (DERIVED identity):
      Tr e^{−tΔ} ∼ (4πt)^{−d/2} ∑_{k≥0} a_k t^k,
      a_1[bulk] = ∫_M √g (1/6 − ξ) R · 1 + …   (scalar)
    so the cutoff² piece of 1/G is proportional to a_1 / (16π²).
    Interacting 3-loop content lives in higher a_n and in matter vertices;
    the packet computes that via the Wick/Schwinger ledger (bulk).
    """
    xi, R = s.symbols("xi R")
    a1_bulk_scalar = (s.Rational(1, 6) - xi) * R
    at_conf = s.simplify(a1_bulk_scalar.subs(xi, s.Rational(1, 6)))
    assert at_conf == 0
    return {
        "form": "∑_n α_n/(2ω_n) ↔ zeta/HK a1 + interactions",
        "a1_bulk_scalar": a1_bulk_scalar,
        "a1_bulk_at_xi_1_6": 0,
        "label": "derived (standard HK); bulk 3-loop = ledger",
    }


def vassilevich_boundary_a1():
    """Manifold-with-boundary a1 (Vassilevich hep-th/0306138).

    Bulk:   a1 ⊃ ∫_M √g (1/6 − ξ) R
    Surface (Dirichlet scalar, standard normalization):
            a1^∂ ⊃ ∫_∂M √h [ (1/3 − ξ) ?  /  ξ-dependent + pure-K ] · K
    Robin BC replace the surface coefficient (scheme / BC dependent).

    At ξ=1/6 the ξ-dependent bulk vanishes; surface still has an extrinsic-
    curvature piece whose coefficient depends on BC (Dirichlet vs Robin).

    Newton matching in the packet: weakly curved R^4 / flat with compactly
    supported metric probe, NO geometric boundary. Stokes ⇒ surface integrals
    of a compactly supported curvature bump cancel. Therefore the GHY-like /
    extrinsic-K contribution to the *integrated* ∫√g R residue is 0.
    """
    xi = s.symbols("xi")
    # Standard Dirichlet surface linear-K density (Vassilevich Table 1 style):
    #   (1/3 − ξ) is the common textbook combination for mixed BC families;
    #   pure Dirichlet often quoted as K/3 independent of ξ after S→∞ limit.
    # We keep both as ASSUMED-cited forms and evaluate at ξ=1/6.
    dirichlet_K_coeff = s.Rational(1, 3)  # pure-K piece after Dirichlet limit
    robin_xi_piece = s.Rational(1, 3) - xi
    at = {
        "dirichlet_K_at_any_xi": dirichlet_K_coeff,
        "robin_xi_piece_at_1_6": s.simplify(robin_xi_piece.subs(xi, s.Rational(1, 6))),
        "geometric_boundary_in_Newton_probe": False,
        "integrated_GHY_contribution_to_EH": 0,
        "label_surface_formula": "assumed (Vassilevich BC table)",
        "label_Newton_vanishing": "derived (no ∂M / Stokes on compact support)",
    }
    assert at["robin_xi_piece_at_1_6"] == s.Rational(1, 6)
    assert at["integrated_GHY_contribution_to_EH"] == 0
    return at


def box_R_and_total_derivatives():
    """□R and total-derivative densities.

    In d=4, √g □R = ∂_μ(√g ∇^μ R) (up to the usual curvature identities for
    the Einstein-Hilbert variation sector). On a closed manifold, or for a
    compactly supported metric probe with falloff,
      ∫ √g □R = 0.
    A local effective Lagrangian may still contain c_□ □R with arbitrary
    finite c_□ (scheme / anomaly / improvement ambiguity), but that density
    does NOT contribute to the integrated Newton residue that the 1/24 test
    matches.

    DERIVED: δ_□R = 0 for the gap-relevant observable.
    """
    return {
        "integral_Box_R": 0,
        "local_coeff_free": True,
        "contributes_to_integrated_EH": False,
        "delta_simple_pole": 0.0,
        "label": "derived (Stokes / total derivative)",
    }


def vertex_contact_dictionary():
    """Coincident interaction-point / Kähler / metric-matching contacts.

    Feynman-graph coincidences (two or three vertices collapsing) produce
    local operators. In BPHZ/R* those UV/IR forests are subtracted into CTs;
    the *finite* remainder after subtraction is scheme-dependent
    (choice of finite counterterms / pure-metric R contact).

    The packet already:
      * includes Z_z and two-loop matter CT (derived from flat RG),
      * uses soft dictionary s_Y = −a_Y R/6 (ξ=1/6 HK),
      * identifies a finite pure-metric R contact as able to absorb the
        entire simple residual (FFFF_GAP_PHYSICAL.md).

    An EXTRA vertex-localized δ^{(4)}(x−y) R contact is therefore:
      — the same scheme class as the pure-metric contact,
      — NOT a new bulk Wick topology,
      — coefficient FREE unless fixed by an independent matching condition.

    Classical ξ R |X|² at ξ=1/6 is already in the action; its one-loop a1
    vanishes (no tadpole). It does not by itself predict a 3-loop finite
    contact of size GAP.
    """
    return {
        "classical_xi_R_phi2": "in action; a1→0 at ξ=1/6",
        "BPHZ_forests": "already in ledger CT / R*",
        "extra_finite_R_contact_coeff": None,  # free / scheme
        "same_class_as": "pure-metric contact (GAP_PHYSICAL)",
        "derived_extra_3loop_contact": 0.0,  # none derived
        "label": (
            "derived: forests accounted; "
            "assumed: any leftover finite contact is scheme"
        ),
    }


# ---------------------------------------------------------------------------
# (2) Compute candidate contributions
# ---------------------------------------------------------------------------

def compute_surface_contributions():
    """Assemble every candidate surface/contact piece vs the bulk ledger."""
    sak = sakharov_spectral_sum()
    bdy = vassilevich_boundary_a1()
    box = box_R_and_total_derivatives()
    vtx = vertex_contact_dictionary()

    pieces = {
        "bulk_a1_xi_1_6": {
            "value": 0.0,
            "label": "derived",
            "note": "scalar tadpole linear-R vanishes",
        },
        "GHY_extrinsic_K_integrated": {
            "value": float(bdy["integrated_GHY_contribution_to_EH"]),
            "label": "derived",
            "note": "no geometric ∂M in Newton probe",
        },
        "Box_R_integrated": {
            "value": box["delta_simple_pole"],
            "label": "derived",
            "note": "total derivative; Stokes",
        },
        "vertex_delta_contact_derived": {
            "value": vtx["derived_extra_3loop_contact"],
            "label": "derived (none fixed)",
            "note": "finite contact coeff is scheme / free",
        },
        "xi_1_6_extra_beyond_soft_dictionary": {
            "value": 0.0,
            "label": "derived",
            "note": "s_Y=−a_Y R/6 already is the ξ=1/6 dictionary",
        },
    }
    total_derived = sum(p["value"] for p in pieces.values())
    return {
        "sakharov": sak,
        "boundary": bdy,
        "box": box,
        "vertex": vtx,
        "pieces": pieces,
        "total_derived_surface": total_derived,
    }


# ---------------------------------------------------------------------------
# (3) Compare to gap
# ---------------------------------------------------------------------------

def compare_to_gap(total_derived: float):
    """Ratio predicted/target. Derived surface sum vs GAP."""
    if GAP == 0:
        ratio = float("nan")
    else:
        ratio = total_derived / GAP
    # Numerical coincidences (NOT derived as surface terms)
    ln2_over_72 = log(2) / 72.0
    coincidences = {
        "ln2/72": ln2_over_72,
        "1/104": 1.0 / 104.0,
        "7/720": 7.0 / 720.0,
    }
    coin_rows = []
    for name, val in coincidences.items():
        coin_rows.append(
            {
                "name": name,
                "value": val,
                "ratio_to_gap": val / GAP,
                "within_1sigma": abs(val - GAP) <= GAP_ERR,
                "label": "NOT derived — numerical proximity only",
            }
        )
    return {
        "GAP": GAP,
        "GAP_err": GAP_ERR,
        "predicted_derived": total_derived,
        "ratio": ratio,
        "within_uncertainty": abs(total_derived - GAP) <= GAP_ERR,
        "free_parameter": (
            "yes — any absorbing finite R contact is a free scheme coeff; "
            "derived surface/GHY/□R pieces are parameter-free and =0"
        ),
        "coincidences": coin_rows,
    }


# ---------------------------------------------------------------------------
# (4) Optional: Jacobson / entanglement-entropy cross-check
# ---------------------------------------------------------------------------

def entanglement_entropy_crosscheck():
    """Solodukhin / Jacobson area-law ↔ 1/G for the WZ field content.

    Free nonminimal complex scalar (Solodukhin):
      S_div ⊃ A/(48π) · (1/6 − ξ) / ε² + …
    At ξ=1/6 the quadratic area divergence vanishes — same as a1_bulk=0.
    Weyl fermion contributes a standard (ξ-independent) area piece already
    encoded in the one-loop HK baseline of the packet; it is not a new
    3-loop surface prediction.

    Contact-term ambiguity: for nonminimal scalars and gauge fields the
    map EE → 1/G has a finite contact ambiguity of the same scheme class
    as a pure-metric R contact (Jacobson; Solodukhin reviews). That
    ambiguity can absorb a finite matching residual but does NOT predict
    the number +0.00965.

    Full interacting 3-loop EE for WZ = same hardness as the bulk ledger;
    not attempted. Feasible subset: free-field consistency + ambiguity.
    """
    xi = s.Rational(1, 6)
    scalar_area_quad = s.simplify(s.Rational(1, 6) - xi)
    assert scalar_area_quad == 0
    # One complex scalar (WZ) + one Weyl fermion.
    # Fermion area-law UV piece ∝ N_f / (area ε²) with N_f=1 (Weyl);
    # its contribution to induced 1/G is the standard one-loop spinor HK
    # already in the packet baseline — not a 3-loop gap explanation.
    return {
        "scalar_EE_quad_at_xi_1_6": 0,
        "matches_a1_vanishing": True,
        "fermion_EE": "one-loop HK baseline (already in packet); not 3-loop gap",
        "contact_ambiguity": True,
        "ambiguity_predicts_GAP": False,
        "interacting_3loop_EE": "not attempted (same hardness as bulk ledger)",
        "label": (
            "derived: free EE quad ∝ (1/6−ξ) → 0; "
            "assumed: EE↔1/G finite contact ambiguity exists; "
            "derived: that ambiguity does not fix GAP"
        ),
    }


def higher_poles_unaffected():
    """Surface/contact finite pieces cannot reopen triple/double cancelation."""
    return {
        "triple_ok": True,
        "double_ok": True,
        "reason": (
            "GHY/□R integrate to 0; finite local R contact is simple/finite "
            "matching only (same as #16 ξv²)"
        ),
        "label": "derived",
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=== Surface / boundary / contact terms vs Newton gap ===\n")
    print("No tuning to 1/24. Dan: mode phase/freq/amp at interaction surfaces.\n")

    print("-- (1) Sakharov spectral sum; bulk vs boundary --")
    contrib = compute_surface_contributions()
    sak = contrib["sakharov"]
    print(f"  spectral form: {sak['form']}")
    print(f"  a1_bulk scalar at ξ=1/6: {sak['a1_bulk_at_xi_1_6']}  [{sak['label']}]")
    bdy = contrib["boundary"]
    print(
        f"  geometric ∂M in Newton probe? {bdy['geometric_boundary_in_Newton_probe']}"
    )
    print(
        f"  integrated GHY → EH residue: {bdy['integrated_GHY_contribution_to_EH']}"
        f"  [{bdy['label_Newton_vanishing']}]"
    )
    print(
        f"  Robin (1/3−ξ)|_ξ=1/6 density (if ∂M existed): "
        f"{bdy['robin_xi_piece_at_1_6']}  [{bdy['label_surface_formula']}]"
    )
    box = contrib["box"]
    print(
        f"  ∫□R: {box['integral_Box_R']}; contributes to EH? "
        f"{box['contributes_to_integrated_EH']}  [{box['label']}]"
    )
    vtx = contrib["vertex"]
    print(f"  vertex δ-contact: {vtx['label']}")
    print(f"  same class as: {vtx['same_class_as']}")

    print("\n-- (2) Candidate finite simple-pole pieces --")
    for name, p in contrib["pieces"].items():
        print(f"  {name:40} = {p['value']:+.6e}  [{p['label']}]  ({p['note']})")
    total = contrib["total_derived_surface"]
    print(f"  TOTAL derived surface/contact     = {total:+.6e}")

    print("\n-- (3) vs GAP +0.00965 ± 0.0024 --")
    cmp_ = compare_to_gap(total)
    print(f"  GAP = {cmp_['GAP']:+.8f} ± {cmp_['GAP_err']}")
    print(f"  predicted (derived surface sum) = {cmp_['predicted_derived']:+.8f}")
    print(f"  ratio predicted/GAP = {cmp_['ratio']:.6f}")
    print(f"  within ±1σ? {cmp_['within_uncertainty']}")
    print(f"  free parameter? {cmp_['free_parameter']}")
    print("  numerical coincidences (NOT derived as surface terms):")
    for row in cmp_["coincidences"]:
        print(
            f"    {row['name']:10} = {row['value']:.8f}  "
            f"ratio={row['ratio_to_gap']:.4f}  "
            f"within_1σ={row['within_1sigma']}  [{row['label']}]"
        )

    print("\n-- (4) Optional EE / Jacobson cross-check --")
    ee = entanglement_entropy_crosscheck()
    print(f"  scalar EE quad ∝ (1/6−ξ) at ξ=1/6: {ee['scalar_EE_quad_at_xi_1_6']}")
    print(f"  matches a1 vanishing? {ee['matches_a1_vanishing']}")
    print(f"  fermion EE: {ee['fermion_EE']}")
    print(f"  contact ambiguity exists? {ee['contact_ambiguity']}")
    print(f"  ambiguity predicts GAP? {ee['ambiguity_predicts_GAP']}")
    print(f"  3-loop interacting EE: {ee['interacting_3loop_EE']}")
    print(f"  [{ee['label']}]")

    print("\n-- Higher poles --")
    hp = higher_poles_unaffected()
    print(
        f"  triple_ok={hp['triple_ok']}, double_ok={hp['double_ok']}  [{hp['label']}]"
    )
    print(f"  {hp['reason']}")

    print("\nVerdict")
    print("-------")
    print(
        "(1) Sakharov/zeta: bulk a1(ξ=1/6)=0; 3-loop bulk = existing ledger.\n"
        "    Boundary/GHY and □R integrate to 0 for the Newton probe.\n"
        "    Vertex δ-contacts = finite scheme contacts already named in\n"
        "    GAP_PHYSICAL — not a dropped bulk topology.\n"
        "(2) Derived extra simple-pole from surface/□R/ξ/vertex = 0.\n"
        "(3) ratio predicted/GAP = 0. Absorbing contact is scheme-dependent\n"
        "    (free finite coeff), not derived. ln2/72≈GAP is coincidence only.\n"
        "(4) EE: free scalar area quad vanishes at ξ=1/6 (consistent); EE\n"
        "    contact ambiguity does not predict +0.00965.\n"
        "Plain: surface/boundary/contact analysis does NOT explain the gap\n"
        "without invoking a free scheme contact of size GAP."
    )
    print("\n=== SURFACE TERM COMPLETE ===")
    return {
        "total_derived": total,
        "ratio": cmp_["ratio"],
        "GAP": GAP,
    }


if __name__ == "__main__":
    main()
