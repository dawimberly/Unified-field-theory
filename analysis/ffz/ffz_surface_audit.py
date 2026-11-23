#!/usr/bin/env python3
"""AG: Boundary / surface / total-derivative audit of the FFZ three-loop pole.

Attack GAP=+0.009651 (target 1/24 in h²z/K) using terms *internal* to the
original FFZ fields (WZ complex scalar + Weyl fermion; ξ=1/6; a_X=0, a_Y=-1).
External sectors gave δC=0 (#42/#54).

Question: were □R, GHY-type boundary, Euler density, or a2-type total
derivatives dropped in the three-loop computation? Compute each piece's
contribution to the simple pole C explicitly. Scheme-independent?

No tuning. DERIVED vs ASSUMED. Fields-only. Honest negatives.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_surface_audit.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import pi

C_LEDGER = 0.05131719
C_TARGET = 1.0 / 24.0
GAP = C_LEDGER - C_TARGET  # ≈ +0.00965052
GAP_ERR = 0.0024
K = 16.0 * pi**2


def ffz_field_content() -> dict:
    """Original FFZ Newton-test fields (no external gauge/GUT)."""
    return {
        "complex_scalar": {"xi": F(1, 6), "dofs": 2, "label": "derived (WZ)"},
        "weyl_fermion": {"dofs": 2, "label": "derived (WZ)"},
        "gauge": None,
        "prescription": "ξ=1/6, a_X=0, a_Y=-1, L=0, units h²z/K",
        "label": "derived (packet)",
    }


def three_loop_setup() -> dict:
    """What the bulk ledger actually computes."""
    return {
        "object": "3-loop simple pole of induced ∫√g R in h²z/K",
        "method": "Wick / Schwinger / soft dictionary on weakly curved R⁴",
        "geometric_boundary": False,
        "metric_probe": "compactly supported curvature bump",
        "C_ledger": C_LEDGER,
        "C_target": C_TARGET,
        "GAP": GAP,
        "label": "derived (ledger definition; #16–#23)",
    }


def term_box_R() -> dict:
    """□R total derivative.

    √g □R = ∂_μ(√g ∇^μ R) (up to standard identities). On closed M or for
    compactly supported probe: ∫√g □R = 0. Local finite c_□ is scheme /
    improvement free, but does not shift the *integrated* Newton residue C.
    """
    delta_C = 0.0
    return {
        "term": "□R",
        "dropped_from_ledger": "as local density may be omitted; integral = 0",
        "delta_C_simple_pole": delta_C,
        "scheme_independent": True,  # vanishing is SI
        "equals_GAP": False,
        "label": "derived (Stokes)",
    }


def term_ghy() -> dict:
    """Gibbons–Hawking–York / extrinsic-K surface term.

    Newton probe has no geometric ∂M. Stokes on compact support ⇒ surface
    integrals cancel. Vassilevich Dirichlet/Robin K densities exist only if
    a boundary is *introduced* (see AI); they are not present in the bulk
    three-loop setup.
    """
    return {
        "term": "GHY / extrinsic K",
        "geometric_boundary_in_setup": False,
        "delta_C_simple_pole": 0.0,
        "scheme_independent": True,  # vanishing without ∂M is SI
        "equals_GAP": False,
        "label": "derived (no ∂M)",
    }


def term_euler() -> dict:
    """Euler density E₄ (Gauss–Bonnet).

    ∫√g E₄ = 32π² χ(M) on closed 4-manifolds (topological). Metric variation
    of ∫E vanishes (Lovelock) — E does not source the Einstein tensor /
    linear-R Newton matching. Local a2 HK coefficient contains E; that
    structure is *not* the three-loop h²z/K simple pole of interacting WZ.

    Contribution to C: 0 (wrong operator / topological).
    Scheme: χ is SI, but irrelevant to GAP in h²z/K.
    """
    return {
        "term": "Euler density E₄",
        "integral": "32π² χ(M) (closed)",
        "varies_to_Einstein": False,
        "delta_C_simple_pole": 0.0,
        "scheme_independent_as_topology": True,
        "scheme_independent_into_h2zK": False,  # not in that channel
        "equals_GAP": False,
        "label": "derived (Gauss–Bonnet / Lovelock; channel split)",
    }


def term_a2_total_derivatives() -> dict:
    """Heat-kernel a2 total-derivative / improvement pieces.

    a2 ⊃ □R, E, Weyl², … . For free-field 1/G the UV-sensitive piece is a1
    (∝ Λ²), not a2. Interacting three-loop simple pole lives in matter
    Wick topologies already in the ledger — not in free a2 densities.

    Dropping free a2 TD pieces does not change C.
    """
    return {
        "term": "a2-type □R / E / improvements",
        "affects_free_1_over_G": "a1, not a2",
        "affects_3loop_h2zK": False,
        "delta_C_simple_pole": 0.0,
        "scheme_independent": True,  # zero contribution is SI
        "equals_GAP": False,
        "label": "derived (HK grading + channel)",
    }


def term_finite_vacuum_eh() -> dict:
    """Finite vacuum ∫√g R counterterm κ (#23).

    Allowed by curved-space renormalizability. Shifts C → C+κ.
    κ = −GAP removes the gap — but κ is *scheme*, not a derived surface
    term from the three-loop matter graphs.
    """
    kappa = -GAP
    return {
        "term": "finite vacuum EH κ",
        "delta_C_simple_pole": kappa,  # if chosen to absorb
        "forced_by_matter_graphs": False,
        "scheme_independent": False,
        "equals_minus_GAP": True,
        "kappa": kappa,
        "label": "derived allowed CT; finite κ ASSUMED/scheme (#23)",
        "plain": (
            f"κ={kappa:+.8f} absorbs GAP but is not a dropped bulk topology."
        ),
    }


def audit_table() -> dict:
    """Assemble all candidate dropped terms."""
    terms = [
        term_box_R(),
        term_ghy(),
        term_euler(),
        term_a2_total_derivatives(),
    ]
    derived_sum = sum(t["delta_C_simple_pole"] for t in terms)
    vac = term_finite_vacuum_eh()
    si_nonzero = [
        t["term"]
        for t in terms
        if t["delta_C_simple_pole"] != 0.0 and t.get("scheme_independent")
    ]
    return {
        "terms": terms,
        "derived_surface_TD_sum": derived_sum,
        "vacuum_kappa": vac,
        "SI_nonzero_into_C": si_nonzero,
        "closes_gap_without_scheme": False,
        "GAP": GAP,
        "GAP_after_derived": GAP - derived_sum,
        "honest_negative": derived_sum == 0.0,
        "plain": (
            "All derived □R / GHY / Euler / a2-TD contributions to the "
            f"h²z/K simple pole are 0. GAP stays {GAP:.6f}±{GAP_ERR}. "
            "Only scheme κ=−GAP absorbs it (#23) — not a dropped bulk term."
        ),
    }


def self_checks() -> None:
    a = audit_table()
    assert a["derived_surface_TD_sum"] == 0.0
    assert a["honest_negative"] is True
    assert a["closes_gap_without_scheme"] is False
    assert abs(a["GAP_after_derived"] - GAP) < 1e-12
    assert a["vacuum_kappa"]["kappa"] == -GAP
    assert term_euler()["delta_C_simple_pole"] == 0.0
    assert term_box_R()["delta_C_simple_pole"] == 0.0
    setup = three_loop_setup()
    assert setup["geometric_boundary"] is False
    assert abs(setup["GAP"] - GAP) < 1e-12


def main():
    print("=== AG: Surface / TD audit of FFZ three-loop simple pole ===\n")
    print("No tuning. Original FFZ fields only. Attack GAP.\n")

    fc = ffz_field_content()
    print("-- Field content --")
    print(f"  scalar ξ={fc['complex_scalar']['xi']}, gauge={fc['gauge']}")
    print(f"  [{fc['prescription']}]")

    setup = three_loop_setup()
    print("\n-- Three-loop setup --")
    print(f"  {setup['object']}")
    print(f"  geometric ∂M? {setup['geometric_boundary']}")
    print(f"  C_ledger={setup['C_ledger']:.8f}  target={setup['C_target']:.8f}")
    print(f"  GAP={setup['GAP']:+.8f} ± {GAP_ERR}")

    a = audit_table()
    print("\n-- Candidate dropped terms → δC --")
    for t in a["terms"]:
        print(
            f"  {t['term']:28} δC={t['delta_C_simple_pole']:+.6e}  "
            f"SI={t.get('scheme_independent', t.get('scheme_independent_into_h2zK'))}  "
            f"[{t['label']}]"
        )
    vac = a["vacuum_kappa"]
    print(
        f"  {vac['term']:28} δC={vac['delta_C_simple_pole']:+.6e}  "
        f"SI={vac['scheme_independent']}  [{vac['label']}]"
    )
    print(f"\n  derived TD/surface sum = {a['derived_surface_TD_sum']:+.6e}")
    print(f"  SI nonzero into C? {a['SI_nonzero_into_C']}")
    print(f"  GAP after derived = {a['GAP_after_derived']:+.8f}")
    print(f"\n  {a['plain']}")

    print("\nVerdict")
    print("-------")
    print(
        "□R, GHY, Euler, a2-TD: each δC=0 for the integrated h²z/K pole.\n"
        "None is a scheme-independent nonzero that equals GAP.\n"
        "κ=−0.00965052 absorbs (scheme) — already #23; not a bulk omission."
    )
    print("\n=== SURFACE AUDIT COMPLETE ===")
    return a


if __name__ == "__main__":
    main()
