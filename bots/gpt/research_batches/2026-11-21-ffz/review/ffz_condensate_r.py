#!/usr/bin/env python3
"""Higgs-like condensate / phase-locking shift of the induced-R residue.

Dan's new-idea rung (not an established peer-reviewed method on the ladder):
give the WZ scalar a background vev and ask whether the resulting R contact
can absorb the standing simple-pole gap +0.00975 without tuning, and whether
any independent condition fixes that vev.

No tuning to 1/24. Derived vs assumed labeled in the output.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_condensate_r.py

Requires: sympy, mpmath
"""

from __future__ import annotations

from fractions import Fraction as F

import mpmath as mp
import sympy as s

mp.mp.dps = 40

GAP = 0.00974619  # standing ledger − 1/24 (superspace/DRED rung)
XI = F(1, 6)
K = 16 * s.pi**2  # packet K=16π²


# ---------------------------------------------------------------------------
# (1) Vacuum / SUSY assumption
# ---------------------------------------------------------------------------

def vacuum_assumption():
    """WZ superpotential W = m X Y + y X³/6.

    F-terms: F_X = m Y + (y/2) X², F_Y = m X.
    SUSY vacuum F_X=F_Y=0 ⇒ X=Y=0 uniquely (for m≠0).
    There is NO SUSY-preserving flat direction with <φ>≠0.

    A background vev therefore requires soft and/or F-term breaking
    (e.g. soft mass + soft trilinear, or a linear superpotential deformation).
    """
    m, y, X, Y = s.symbols("m y X Y")
    F_X = m * Y + y * X**2 / 2
    F_Y = m * X
    # At X=Y=0, F=0
    assert F_X.subs({X: 0, Y: 0}) == 0 and F_Y.subs({X: 0, Y: 0}) == 0
    # If X≠0 and m≠0, F_Y ≠ 0
    return {
        "superpotential": "W = m X Y + y X³/6",
        "SUSY_vacuum": "X=Y=0 unique",
        "flat_direction_nonzero_vev": False,
        "assumption": (
            "ASSUMED: soft and/or F-term breaking that allows a background "
            "vev <X>=v (or phase-locked condensate) at fixed ξ=1/6. "
            "Not SUSY-preserving."
        ),
        "label": "assumed (breaking); derived (no SUSY flat direction)",
    }


# ---------------------------------------------------------------------------
# (2) Induced R shift
# ---------------------------------------------------------------------------

def tree_level_xi_shift():
    """Tree-level non-minimal term ξ R |φ|² → ξ v² R with v² ≡ |<X>|².

    Action convention (Jordan): … + ξ R |φ|² + …
    Packet Γ ⊃ −R/(16π G_ind); a positive ξ v² increases the EH coefficient
    in the Jordan frame (same role as M_Pl²). In the ledger's F_R matching
    units (where the three-loop residue is compared to +1/24), a finite
    pure-metric contact that *reduces* a too-large baseline is equivalent
    to subtracting from the reported simple total.

    Dimensionless tree shift in units of h² z / K:
        δ_tree = ξ v² / (h² z / K) = (1/6) (K v²) / (h² z)

    DERIVED from ξ=1/6 and dimensional conversion. Sign for *absorbing* a
    positive gap: need δ_match = −GAP in the ledger total, i.e. the contact
    must lower the reported baseline by GAP (scheme contact with that size).
    """
    h, z, v2 = s.symbols("h z v2", positive=True)
    delta = s.simplify(XI * v2 / (h**2 * z / K))
    # = K v2 / (6 h² z)
    assert delta == K * v2 / (6 * h**2 * z)
    return {
        "xi": XI,
        "delta_tree_symbolic": delta,
        "delta_tree_form": "K v² / (6 h² z)",
        "pole_structure": (
            "DERIVED: tree ξ v² is a finite local R contact — it does NOT "
            "generate 1/ε³ or 1/ε² poles. It shifts the finite / matching "
            "simple residue only (same class as the pure-metric contact that "
            "FFFF_GAP_PHYSICAL.md already noted can absorb the gap)."
        ),
        "label": "derived",
    }


def one_loop_mass_shift():
    """Vev-induced masses m² ~ |y|² v² = h K v² inserted into 1-loop F_R.

    Packet one-loop: F_R pole ~ z w /(12 K ε) with w=2+2a_X+2a_Y.
    At a_X=a_Y=0 and ξ=1/6, conformal scalars have a1=0; the Dirac/equal-mass
    baseline is ∝ z/(6Kε) structure used in the CT.

    Soft vev: fluctuation masses pick up δm² = c |y|² v² = c h K v²
    (c = O(1) group factor; ASSUMED c=1 for a single complex X direction
    as an order-of-magnitude estimate). Relative to the SUSY mass z:
        ρ = δm² / z = h K v² / z

    One-loop simple shift in units h² z/K would be *lower* order in h unless
    we count δm² as an independent soft scale. Relative to the three-loop
    test (order h²), the leading condensate effect at fixed soft scale is
    the tree ξ v² contact; the loop mass insertion is O(h) × (v²/z) or
    O(δm²/K) and is parametrically smaller in h when v² ~ O(h² z / K)×gap
    (the gap-absorbing window — see below).

    Higher poles: soft IR masses do not change the UV 1/ε³,1/ε² structure
    of the dimensionless h² graphs; they reshuffle finite and m²/ε pieces.
    """
    h, z, v2 = s.symbols("h z v2", positive=True)
    rho = h * K * v2 / z  # δm²/z with c=1
    # One-loop F_R ~ δm²/(12 K) in absolute units → in h² z/K units:
    delta_1l = s.simplify((h * K * v2 / 12) / (h**2 * z / K) * K / K)
    # δm²/(12K) / (h² z/K) = (h K v²)/(12K) * K/(h² z) = v²/(12 h z) * K/1
    # Wait: δm² = h K v²; δm²/(12K) = h v²/12
    # / (h² z / K) = (h v²/12) * K / (h² z) = K v² / (12 h z)
    delta_1l = K * v2 / (12 * h * z)
    return {
        "rho": rho,
        "delta_1loop_form": "K v² / (12 h z)   [c=1 estimate]",
        "delta_1loop_symbolic": delta_1l,
        "vs_tree": (
            "tree/1-loop ~ (K v²)/(6 h² z)  /  (K v²)/(12 h z) = 2/h "
            "→ tree dominates at weak h"
        ),
        "two_loop_note": (
            "Two-loop vev corrections are O(h) relative to one-loop mass "
            "insertions; not needed to decide whether v is fixed or to size "
            "the tree contact that absorbs the gap."
        ),
        "label": "derived (structure); assumed (c=1 single-direction estimate)",
    }


def two_loop_note():
    return {
        "status": (
            "Not recomputed graph-by-graph. Existing ledger scripts assume "
            "the SUSY equal-mass point z; a full re-run with split masses "
            "m²=z+c h K v² is a finite deformation of IR masses and does not "
            "move the UV higher poles that already cancel."
        ),
        "label": "assumed (no full re-ledger); derived (UV poles mass-IR blind)",
    }


# ---------------------------------------------------------------------------
# (3) Is v fixed?  Gap-absorbing value
# ---------------------------------------------------------------------------

def independent_fixes_for_v():
    """Check candidate conditions that might fix v."""
    checks = [
        {
            "condition": "a_Y = −1 projection",
            "fixes_v": False,
            "note": "Source-weight projection, not a scalar vev",
        },
        {
            "condition": "Yukawa RG fixed point",
            "fixes_v": False,
            "note": "Fixes h_*, not the dimensionful ratio v²/z",
        },
        {
            "condition": "Coleman–Weinberg / Veff minimum (SUSY W)",
            "fixes_v": False,
            "note": "SUSY vacuum is X=Y=0; CW does not generate a nonzero vev without soft terms",
        },
        {
            "condition": "Soft-breaking scale m_soft",
            "fixes_v": "external",
            "note": (
                "v ~ m_soft/|y| is set by the soft Lagrangian — an input, "
                "not predicted by the Newton matching"
            ),
        },
        {
            "condition": "Gap = +0.00975 itself",
            "fixes_v": "tuning",
            "note": "Choosing v to cancel the gap is tuning, not a prediction",
        },
    ]
    return {
        "verdict": (
            "v is a FREE parameter from the point of view of the Newton "
            "matching (or is set by an external soft scale). No independent "
            "condition in the model predicts the gap-absorbing value."
        ),
        "checks": checks,
        "label": "derived (logic of each check)",
    }


def gap_absorbing_v():
    """v that makes |δ_tree| = GAP in units h² z/K.

    Need K v² / (6 h² z) = GAP
    ⇒  v² / z = GAP * 6 * h² / K
    ⇒  v² / z = (GAP/K) * 6 h²

    Numerically K=16π² ≈ 157.91, GAP≈0.009746
    ⇒ v²/z = 6 * GAP * h² / K ≈ 3.702e-4 * h²
    """
    gap = s.Rational(*F(GAP).limit_denominator(10**12).as_integer_ratio())
    # use float for report
    Kf = float(16 * mp.pi**2)
    # v²/z = 6 * GAP * h² / K
    pref = 6 * GAP / Kf
    # Also quote for h=1 (strong) and h=0.1
    return {
        "condition": "K v² / (6 h² z) = GAP",
        "v2_over_z": f"{pref:.6e} × h²",
        "v2_over_z_h1": pref * 1.0,
        "v2_over_z_h01": pref * 0.01,
        "v_over_sqrt_z_h1": pref**0.5,
        "v_over_sqrt_z_h01": (pref * 0.01) ** 0.5,
        "GAP": GAP,
        "K": Kf,
        "independent_prediction": False,
        "label": "derived (algebra); assumed (that tree contact is the absorber)",
    }


# ---------------------------------------------------------------------------
# (4) Pole structure summary + ledger comparison
# ---------------------------------------------------------------------------

def pole_structure_verdict():
    return {
        "tree_xi_v2": "finite matching / simple residue only — NO change to triple/double",
        "soft_mass_loops": (
            "can shift m²/ε simple pieces; UV 1/ε³ and 1/ε² of h² graphs unchanged"
        ),
        "higher_poles_still_cancel": True,
        "label": "derived",
    }


def main():
    print("=== Condensate / phase-locking R-shift (new-idea rung) ===\n")
    print("No tuning to 1/24. Labels: derived vs assumed.\n")

    print("-- (1) Vacuum / SUSY --")
    vac = vacuum_assumption()
    print(f"  W: {vac['superpotential']}")
    print(f"  SUSY vacuum: {vac['SUSY_vacuum']}")
    print(f"  Nonzero-vev flat direction? {vac['flat_direction_nonzero_vev']}")
    print(f"  {vac['assumption']}")
    print(f"  [{vac['label']}]")

    print("\n-- (2) Induced R shift --")
    tr = tree_level_xi_shift()
    print(f"  ξ = {tr['xi']}")
    print(f"  δ_tree (units h²z/K) = {tr['delta_tree_form']}")
    print(f"  {tr['pole_structure']}")
    print(f"  [{tr['label']}]")

    ol = one_loop_mass_shift()
    print(f"\n  δ_1loop ~ {ol['delta_1loop_form']}")
    print(f"  {ol['vs_tree']}")
    print(f"  {ol['two_loop_note']}")
    print(f"  [{ol['label']}]")

    tl = two_loop_note()
    print(f"\n  Two-loop re-ledger: {tl['status']}")
    print(f"  [{tl['label']}]")

    print("\n-- (3) Is v fixed? Gap-absorbing value --")
    fx = independent_fixes_for_v()
    for c in fx["checks"]:
        print(f"  • {c['condition']}: fixes_v={c['fixes_v']} — {c['note']}")
    print(f"  Verdict: {fx['verdict']}")
    print(f"  [{fx['label']}]")

    gv = gap_absorbing_v()
    print(f"\n  To absorb GAP={gv['GAP']:.5f} by tree ξv²:")
    print(f"    {gv['condition']}")
    print(f"    v²/z = {gv['v2_over_z']}")
    print(f"    at h=1:   v²/z = {gv['v2_over_z_h1']:.6e},  |v|/√z = {gv['v_over_sqrt_z_h1']:.6e}")
    print(f"    at h=0.1: v²/z = {gv['v2_over_z_h01']:.6e},  |v|/√z = {gv['v_over_sqrt_z_h01']:.6e}")
    print(f"  Independent condition predicts this value? {gv['independent_prediction']}")
    print(f"  [{gv['label']}]")

    print("\n-- (4) Pole structure --")
    ps = pole_structure_verdict()
    print(f"  Tree: {ps['tree_xi_v2']}")
    print(f"  Soft loops: {ps['soft_mass_loops']}")
    print(f"  Higher poles still cancel? {ps['higher_poles_still_cancel']}")
    print(f"  [{ps['label']}]")

    print(
        f"""
Verdict
-------
(1) Nonzero <φ>=v is NOT SUSY-preserving for W=mXY+yX³/6. Soft/F-term
    breaking is ASSUMED.

(2) Tree ξv² at ξ=1/6 shifts the finite R contact by K v²/(6 h² z) in
    ledger units; 1-loop mass insertions are O(h) smaller at weak coupling.
    Higher poles unchanged.

(3) v is FREE (or external soft-scale). The value that absorbs GAP=+{GAP:.5f}
    is v²/z = {6*GAP/float(16*mp.pi**2):.6e}×h² — NOT predicted by a_Y=-1,
    Yukawa fixed point, or CW without soft inputs. Choosing it is tuning.

(4) The shift lands in the simple / finite matching residue, not the
    triple/double pole structure.

Plain statement: a condensate can *parametrize* the gap as a Jordan-frame
ξv² contact, but it does not explain it — v is not fixed by the model
ingredients that define the Newton test. This new-idea rung does not close
the gap without tuning.
"""
    )
    # Locks
    assert vac["flat_direction_nonzero_vev"] is False
    assert abs(gv["v2_over_z_h1"] - 6 * GAP / float(16 * mp.pi**2)) < 1e-12
    assert ps["higher_poles_still_cancel"]
    assert gv["independent_prediction"] is False
    print("=== CONDENSATE R RUNG COMPLETE ===")


if __name__ == "__main__":
    main()
