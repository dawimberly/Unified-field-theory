#!/usr/bin/env python3
"""Gauge multiplet + Kabat / Donnelly–Wall contact vs the Newton gap.

Dan's surface/interaction-point idea via gauge fields. New-idea rung stacked
on ξ-running (#21). Extend the WZ model by an N=1 vector multiplet (U(1) or
SU(N)) and ask whether the known gauge contact / edge-mode contribution to
1/G can explain GAP=+0.00965±0.0024.

(1) Gauge-sector induced-R / a1 including Kabat contact (Solodukhin) and
    Donnelly–Wall edge-mode interpretation; coefficients from field content.
(2) SUSY: does the gaugino cancel the contact? Net N=1 vector multiplet?
    Express vs h; use N=2 relation ğ=h when a ratio is needed.
(3) Compare to GAP; check whether ln2 (Z2 edge entropy) arises; ln2/72?
(4) Higher poles.

No tuning to 1/24. Derived vs assumed labeled.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_gauge_contact.py

Requires: sympy
"""

from __future__ import annotations

from fractions import Fraction as F
from math import log, pi

import sympy as s

K = 16.0 * pi**2
C_LEDGER = 0.05131719
C_TARGET = 1.0 / 24.0
GAP = C_LEDGER - C_TARGET
GAP_ERR = 0.0024
LN2_OVER_72 = log(2) / 72.0


# ---------------------------------------------------------------------------
# (1) Kabat contact + Donnelly–Wall; heat-kernel a1
# ---------------------------------------------------------------------------

def kabat_maxwell_a1():
    """Solodukhin / Kabat: Maxwell + ghosts in d dimensions.

    a1^(1) = (d/6 − 1) ∫R for the vector operator (Tr X = R),
    ghosts: −2 × (minimal scalar a1 = ∫R/6),
    ⇒ a1/∫R = d/6 − 1 − 2/6 = D_s/6 − c_Kabat,
    with D_s(d=4)=d−2=2 and c_Kabat=1.

    At d=4: a1/∫R = 2/6 − 1 = −2/3.
    """
    d = s.symbols("d", positive=True)
    a1_vector_op = d / 6 - 1
    a1_ghosts = -2 * (s.Rational(1, 6))
    a1_total = s.simplify(a1_vector_op + a1_ghosts)
    D_s = d - 2
    c_kabat = 1
    assert s.simplify(a1_total - (D_s / 6 - c_kabat)) == 0
    at4 = s.simplify(a1_total.subs(d, 4))
    assert at4 == -s.Rational(2, 3)
    return {
        "a1_over_R_d": a1_total,
        "a1_over_R_d4": F(-2, 3),
        "D_s_d4": 2,
        "c_Kabat": c_kabat,
        "label_coeff": "derived (Solodukhin/Kabat heat kernel)",
    }


def donnelly_wall_prescription():
    """Donnelly–Wall: Kabat contact = entanglement of electromagnetic edge modes
    (normal electric flux on the entangling surface). Leading divergence matches
    one negative scalar degree of freedom confined to the surface.

    Prescription choice: use Donnelly–Wall edge-mode interpretation for the
    statistical meaning of Kabat's term. The coefficient c=1 itself is
    heat-kernel derived; the edge-mode story is the chosen interpretation.
    """
    return {
        "prescription": "Donnelly–Wall edge modes (≡ Kabat contact)",
        "surface_dof": "one negative scalar on entangling surface",
        "label_coeff": "derived (c_Kabat=1 from HK)",
        "label_interpretation": "assumed (DW edge-mode statistical reading)",
        "refs": "Kabat 1995; Donnelly–Wall 2012/2015; Solodukhin 2015",
    }


def gaugino_a1():
    """Majorana gaugino: no Kabat contact (Solodukhin: contact only for
    bosonic gauge mediators). On-shell dofs D=2 ⇒ physical contribution
    +D/6 = +1/3 to the 1/G density (same sign rule as bosons for physical
    dofs; fermions have no extra contact).
    """
    D = 2
    c_contact = 0
    a1 = F(D, 6) - c_contact
    assert a1 == F(1, 3)
    return {
        "a1_over_R": a1,
        "c_contact": 0,
        "label": "derived (Solodukhin: spin-1/2 has no contact)",
    }


def n1_vector_multiplet_a1(N=1, group="U(1)"):
    """Net a1/∫R for one N=1 vector multiplet × dim(G).

    U(1): dim=1. SU(N): dim=N²−1.
    Net = (D_V/6 − c_Kabat) + D_λ/6 = (2/6 − 1) + 2/6 = −1/3 per U(1).
    """
    kab = kabat_maxwell_a1()
    gau = gaugino_a1()
    per_u1 = kab["a1_over_R_d4"] + gau["a1_over_R"]
    assert per_u1 == F(-1, 3)
    if group == "U(1)":
        dim_G = 1
    elif group == "SU(N)":
        dim_G = N * N - 1
    else:
        raise ValueError(group)
    return {
        "group": group,
        "N": N,
        "dim_G": dim_G,
        "a1_per_u1": per_u1,
        "a1_total": per_u1 * dim_G,
        "contact_canceled_by_gaugino": False,
        "label": "derived (field content + Solodukhin a1)",
    }


# ---------------------------------------------------------------------------
# (2) SUSY cancelation; g² vs h; N=2 relation
# ---------------------------------------------------------------------------

def susy_cancelation_check():
    """Does SUSY cancel the Kabat contact against the gaugino?

    No: contact is gauge-boson / edge-mode specific; gaugino has c=0.
    Net multiplet a1 = −1/3 ≠ 0 (per U(1)).
    (cf. hep-ph/9504259: vector multiplets alone induce negative Newton
    constant — same qualitative conclusion, different a1 normalization.)
    """
    net = n1_vector_multiplet_a1()
    return {
        "contact_canceled": False,
        "net_a1_U1": net["a1_per_u1"],
        "remainder": net["a1_per_u1"],
        "label": "derived: gaugino has no contact; net ≠ 0",
    }


def channel_and_n2():
    """Map a1 → FFZ ledger channel (units h² z/K).

    Heat-kernel a1 renormalizes 1/G with a *quadratic* UV piece ∼ Λ² a1/(16π²),
    mass-independent at leading order. The FFZ three-loop simple pole is
    ∝ h² z/(K ε) — a different operator / RG structure.

    In dim-reg MS for massless gauge fields, that Λ² piece is absent from the
    h²z/K ledger. Hence the *derived* contribution to GAP's channel is 0.

    Assumed conversions (for honesty checks, not predictions):
      C1: replace Λ² → z  ⇒  F_R ∼ −a1 z/K  ⇒  coeff in h²z/K units = −a1/h²
      C2: N=2 relation ğ ≡ g²/K = h (y∼g for adjoint Yukawa) still leaves
          the free-field piece O(z/K)=O(h⁰·z/K), i.e. −a1/h² in ledger units.
    Neither conversion is forced by the packet; both are ASSUMED.
    """
    a1 = F(-1, 3)
    # Derived: MS ledger channel
    delta_c_MS = 0.0
    # Assumed C1/C2 at h=1: −a1/h² = 1/3
    h = s.symbols("h", positive=True)
    delta_c_assumed = -a1 / h**2
    at_h1 = float(delta_c_assumed.subs(h, 1))
    n2 = {
        "relation": "ğ = g²/K = h  (N=2 adjoint Yukawa y∼g)",
        "fixes_g_over_h": True,
        "still_needs_Lambda_to_z": True,
        "label": "assumed (N=2 identification); does not remove Λ²↔z choice",
    }
    return {
        "delta_c_MS_ledger": delta_c_MS,
        "delta_c_assumed_expr": delta_c_assumed,
        "delta_c_assumed_h1": at_h1,
        "n2": n2,
        "label_MS": "derived: free a1 is Λ², not h²z/K",
        "label_assumed": "assumed: Λ²→z conversion",
    }


# ---------------------------------------------------------------------------
# (3) vs GAP; ln2
# ---------------------------------------------------------------------------

def compare_to_gap(delta_c_derived: float, delta_c_assumed_h1: float):
    ratio_der = delta_c_derived / GAP
    ratio_ass = delta_c_assumed_h1 / GAP
    # ln2 from edge modes?
    # Continuum U(1) DW edge modes: continuous flux entropy — NOT ln2.
    # Discrete Z2 gauge theory: topological EE contains ln2, but that is NOT
    # the N=1 vector multiplet.
    ln2_from_gauge = False
    ln2_72_reproduced = False
    return {
        "GAP": GAP,
        "predicted_derived_MS": delta_c_derived,
        "ratio_derived": ratio_der,
        "assumed_Lambda_to_z_h1": delta_c_assumed_h1,
        "ratio_assumed_h1": ratio_ass,
        "ln2_from_U1_edge_modes": ln2_from_gauge,
        "ln2_72": LN2_OVER_72,
        "ln2_72_reproduced": ln2_72_reproduced,
        "ln2_note": (
            "Continuum Maxwell/DW edge modes do not produce ln2; "
            "Z2 would, but is not N=1 SYM. ln2/72 is NOT reproduced."
        ),
    }


def higher_poles():
    """Free one-loop vector multiplet: no triple/double poles.

    Adding dynamical gauge fields to the interacting WZ ledger would introduce
    new graphs; that re-ledger is OUT OF SCOPE. Within the existing pure-WZ
    ledger, gauge is absent so poles unchanged. A free multiplet a1 is
    simple/quadratic only.
    """
    return {
        "triple_ok": True,
        "double_ok": True,
        "reason": (
            "free N=1 vector a1 is Λ²/simple only; pure-WZ ledger unchanged; "
            "full gauge+WZ re-ledger not attempted"
        ),
        "label": "derived (free-field order); assumed (no new interacting gauge graphs)",
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=== Gauge / Kabat / Donnelly–Wall contact vs Newton gap ===\n")
    print("No tuning to 1/24. Model extension: + N=1 vector multiplet.\n")

    print("-- (1) Kabat contact + DW edge modes --")
    kab = kabat_maxwell_a1()
    dw = donnelly_wall_prescription()
    print(f"  Maxwell+ghosts a1/∫R (d=4) = {kab['a1_over_R_d4']}  [{kab['label_coeff']}]")
    print(f"  c_Kabat = {kab['c_Kabat']}")
    print(f"  prescription: {dw['prescription']}")
    print(f"  [{dw['label_coeff']}; {dw['label_interpretation']}]")
    gau = gaugino_a1()
    print(f"  Majorana gaugino a1/∫R = {gau['a1_over_R']} (contact={gau['c_contact']})  [{gau['label']}]")
    for grp, N in (("U(1)", 1), ("SU(N)", 2), ("SU(N)", 3)):
        mp = n1_vector_multiplet_a1(N=N, group=grp)
        print(
            f"  {grp} N={N}: dimG={mp['dim_G']}, a1_total={mp['a1_total']}  [{mp['label']}]"
        )

    print("\n-- (2) SUSY cancelation; g² vs h --")
    susy = susy_cancelation_check()
    print(f"  contact canceled by gaugino? {susy['contact_canceled']}")
    print(f"  net a1 per U(1) = {susy['net_a1_U1']}  [{susy['label']}]")
    ch = channel_and_n2()
    print(f"  MS ledger δc (h²z/K channel) = {ch['delta_c_MS_ledger']}  [{ch['label_MS']}]")
    print(f"  assumed Λ²→z at h=1: δc = {ch['delta_c_assumed_h1']:+.6f}  [{ch['label_assumed']}]")
    print(f"  N=2 ğ=h: {ch['n2']['relation']}  [{ch['n2']['label']}]")

    print("\n-- (3) vs GAP; ln2 --")
    cmp_ = compare_to_gap(ch["delta_c_MS_ledger"], ch["delta_c_assumed_h1"])
    print(f"  GAP = {cmp_['GAP']:+.8f} ± {GAP_ERR}")
    print(f"  predicted (derived MS) = {cmp_['predicted_derived_MS']:+.8f}")
    print(f"  ratio derived/GAP = {cmp_['ratio_derived']:.6f}")
    print(
        f"  assumed Λ²→z @h=1 = {cmp_['assumed_Lambda_to_z_h1']:+.6f}  "
        f"ratio={cmp_['ratio_assumed_h1']:.2f}  (NOT a prediction)"
    )
    print(f"  ln2 from U(1) edge modes? {cmp_['ln2_from_U1_edge_modes']}")
    print(f"  ln2/72 reproduced? {cmp_['ln2_72_reproduced']}  (ln2/72={cmp_['ln2_72']:.8f})")
    print(f"  {cmp_['ln2_note']}")

    print("\n-- (4) Higher poles --")
    hp = higher_poles()
    print(
        f"  triple_ok={hp['triple_ok']}, double_ok={hp['double_ok']}  [{hp['label']}]"
    )
    print(f"  {hp['reason']}")

    print("\nVerdict")
    print("-------")
    print(
        "(1) Kabat c=1 is HK-derived; DW edge modes are the chosen reading.\n"
        "    Maxwell+ghosts a1/R=−2/3; gaugino +1/3; net U(1) multiplet −1/3.\n"
        "(2) SUSY does NOT cancel the contact; net remainder −1/3 per U(1).\n"
        "    That remainder is a Λ²/1/G piece, not the h²z/K ledger channel.\n"
        "(3) Derived ratio to GAP = 0. Assumed Λ²→z @h=1 overshoots (ratio∼35).\n"
        "    ln2/72 is NOT reproduced from continuum gauge edge modes.\n"
        "(4) Free multiplet: no triple/double; pure-WZ ledger unchanged.\n"
        "Plain: gauge Kabat/DW contact does not explain the WZ Newton gap."
    )
    print("\n=== GAUGE CONTACT COMPLETE ===")
    return cmp_


if __name__ == "__main__":
    main()
