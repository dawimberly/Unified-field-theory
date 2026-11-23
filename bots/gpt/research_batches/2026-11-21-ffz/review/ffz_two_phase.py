#!/usr/bin/env python3
"""Two-phase mass/gravity: light = E/B only; matter = locked mass↔gravity pair.

Stacked on three-phase (#28). Avoids #28's photon–metric clash by keeping
universal metric coupling and restricting any *extra* gravity-phase piece to
trace sourcing (T^μ_μ ≠ 0 only for massive matter).

(1) Formalize: de Broglie/Compton clock, zitterbewegung, Dirac L–R mass term,
    Higgs phase, dilaton/trace; EP as m_i = m_g locking?
(2) Bounds on extra trace-coupled piece (Cassini-safe metric; MICROSCOPE/Eöt-Wash).
(3) Induced-gravity a1 sign: Proca vs massless; chiral/massive edge modes.
(4) Side-by-side vs three-phase (#28).

No tuning. DERIVED vs ASSUMED. Honest negatives.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_two_phase.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import pi

# From #27 / #28
A1_SM_XI0 = F(-43, 12)  # UV massless gauge count
A1_VECTOR_MASSLESS = F(-2, 3)
A1_VECTOR_PROCA = F(1, 2)
A1_WEYL = F(1, 12)
A1_SCALAR_MIN = F(1, 6)

CASSINI_GAMMA_MINUS_1 = 2.3e-5
# MICROSCOPE (Touboul et al. 2022): |η(Ti,Pt)| ≲ 10^{-15}
MICROSCOPE_ETA = 1.0e-15
# Eöt-Wash / torsion-balance composition EP (order; published)
EOTWASH_ETA = 1.0e-13


# ---------------------------------------------------------------------------
# (1) Formalization
# ---------------------------------------------------------------------------

def de_broglie_compton_clock() -> dict:
    """Mass as phase rotation: ψ ∼ e^{-i m c² t / ℏ}.

    Compton frequency ω_C = m c² / ℏ. This is a *single* U(1) phase clock for
    a massive particle — not by itself a locked pair with gravity.
    """
    return {
        "omega_C": "m c² / ℏ",
        "is_phase_rotation": True,
        "locks_to_gravity_derived": False,
        "label": "derived (de Broglie); gravity locking assumed/not shown",
    }


def zitterbewegung_and_dirac_LR() -> dict:
    """Dirac mass = L–R chirality mixing: m(ψ̄_L ψ_R + h.c.).

    Genuine 2-component locked structure (two Weyl spinors ↔ one Dirac).
    Zitterbewegung: interference of positive/negative energy at ∼ℏ/(2mc).
    This locks *chiral* phases, not a mass↔gravity phase pair.
    """
    return {
        "Dirac_mass_term": "m(ψ̄_L ψ_R + h.c.)",
        "locked_pair": "Weyl_L ↔ Weyl_R",
        "is_mass_gravity_pair": False,
        "zitter_scale": "ℏ/(2 m c)",
        "label_LR": "derived (Dirac)",
        "label_as_gravity": "assumed (not implied by L–R mixing)",
    }


def higgs_condensate_phase() -> dict:
    """Higgs: mass from vev; Goldstone eaten → longitudinal gauge modes.

    Condensate phase is the electroweak order parameter, not a gravitational
    phase. Couples to curvature only via ξ_H R |H|² (nonminimal; ASSUMED ξ_H).
    """
    return {
        "mass_from_vev": True,
        "gravity_phase_derived": False,
        "xi_H_channel": "ξ_H R |H|² (assumed ξ_H; a1 scalar weight)",
        "label": "derived (Higgs mechanism); gravity-phase reading assumed",
    }


def ep_as_locking() -> dict:
    """Equivalence principle: m_inertial = m_gravitational.

    In GR this is geometric (universal free fall / geodesic principle) —
    DERIVED as a structural feature of the metric theory, not as a dynamical
    phase-lock between a mass clock and a gravity wave phase.

    Interpreting m_i = m_g as 'two-phase locking' is ASSUMED language.
    Empirically η = 2(a1−a2)/(a1+a2) is consistent with 0 at ∼10^{-15}.
    """
    return {
        "m_i_equals_m_g_in_GR": True,
        "derived_as_geometry": True,
        "derived_as_two_phase_lock": False,
        "MICROSCOPE_eta": MICROSCOPE_ETA,
        "label_GR": "derived (metric universal coupling)",
        "label_two_phase_reading": "assumed interpretive overlay",
    }


def trace_dilaton_channel() -> dict:
    """Gravity coupled to trace / conformal anomaly (dilaton, scalar-tensor).

    T^μ_μ = 0 for classical Maxwell ⇒ no dilaton source from pure light.
    Massive matter / conformal anomaly ⇒ nonzero trace ⇒ dilaton on.
    Full metric still sourced by all T^{μν} (including EM energy).
    """
    return {
        "T_trace_EM": 0,
        "T_trace_massive_NR": "∼ −ρ c² ≠ 0",
        "dilaton_source_mass_only": True,
        "metric_source_all_Tmn": True,
        "label": "derived (Maxwell traceless; Einstein eq. uses full Tmn)",
    }


def two_phase_formalization() -> dict:
    return {
        "light": "2-phase E/B Maxwell wave; T^μ_μ=0; no dilaton source",
        "matter": "Dirac L–R locked pair (mass); Compton clock; Higgs vev",
        "gravity_phase_on": "nonzero T^μ_μ (dilaton/trace channel) — derived for that channel only",
        "mass_gravity_locked_pair_derived": False,
        "EP_locking_derived_as_phases": False,
        "honest": (
            "Dirac L–R and Compton clocks are real 2-phase/phase-rotation "
            "structures for mass. They do not derive a locked mass↔gravity "
            "wave pair. EP m_i=m_g is geometric in GR, not a two-phase lock. "
            "The only derived mass-only gravity knob is extra trace/dilaton "
            "sourcing."
        ),
        "de_broglie": de_broglie_compton_clock(),
        "dirac": zitterbewegung_and_dirac_LR(),
        "higgs": higgs_condensate_phase(),
        "EP": ep_as_locking(),
        "trace": trace_dilaton_channel(),
    }


# ---------------------------------------------------------------------------
# (2) Avoid #28 clash: universal metric + mass-only extra trace piece
# ---------------------------------------------------------------------------

def cassini_safe_version() -> dict:
    """Universal g_μν (light bends) + optional extra α_T φ T^μ_μ.

    Cassini constrains the PPN γ of the *metric* sector; with universal
    metric coupling, γ=1 is intact. The extra dilaton is a fifth force /
    composition-dependent acceleration constrained by EP tests.
    """
    return {
        "metric_universal": True,
        "light_bending_OK": True,
        "Shapiro_Cassini_OK": True,
        "Cassini_|gamma-1|": CASSINI_GAMMA_MINUS_1,
        "extra_trace_coupling": "α_T φ T^μ_μ  (ASSUMED strength α_T)",
        "label_split": "derived (channel split); α_T assumed free",
    }


def fifth_force_bounds() -> dict:
    """Allowed extra trace-coupled piece.

    For a long-range dilaton, EP violation η ∼ Δα × (composition contrast).
    MICROSCOPE: |η| ≲ 1×10^{-15} ⇒ |α_T|-scale fifth-force strength
    bounded at that order for composition-dependent couplings (schematic;
    precise map depends on dilaton charges — ASSUMED order-of-magnitude).

    Eöt-Wash: η ≲ 10^{-13} (weaker than MICROSCOPE for some compositions).
    Cassini also bounds |γ−1| for Brans–Dicke–like scalars: ω_BD ≳ 4×10^4
    ⇒ fractional scalar admixture ≲ 1/(2ω) ∼ 10^{-5} (same ballpark as γ).
    """
    # Schematic: |α_eff| ≲ MICROSCOPE_ETA for composition-dependent piece
    alpha_T_max_EP = MICROSCOPE_ETA
    # Brans–Dicke-like universal scalar: from Cassini γ−1 = −1/(ω+2)
    # |γ−1|<2.3e-5 ⇒ ω+2 > 1/2.3e-5 ⇒ ω ≳ 4.3e4
    omega_BD_min = 1.0 / CASSINI_GAMMA_MINUS_1 - 2.0
    scalar_admixture_max = 1.0 / (2.0 * omega_BD_min) if omega_BD_min > 0 else None
    return {
        "MICROSCOPE_eta": MICROSCOPE_ETA,
        "EotWash_eta": EOTWASH_ETA,
        "alpha_T_max_composition_schematic": alpha_T_max_EP,
        "omega_BD_min_Cassini": omega_BD_min,
        "scalar_admixture_max_universal": scalar_admixture_max,
        "closes_Newton_gap": False,
        "label_bounds": "assumed (published MICROSCOPE/Cassini/Eöt-Wash)",
        "label_gap": "derived: fifth-force α_T is free ≠ h²z/K ledger",
    }


# ---------------------------------------------------------------------------
# (3) Induced-gravity a1 sign
# ---------------------------------------------------------------------------

def a1_uv_sm() -> F:
    """UV SM (Λ ≫ v): 4 scalars + 45 Weyl + 12 massless vectors."""
    return 4 * A1_SCALAR_MIN + 45 * A1_WEYL + 12 * A1_VECTOR_MASSLESS


def a1_ir_proca_replacement() -> dict:
    """IR count: W±,Z as Proca (+1/2); γ+gluons massless (−2/3).

    ASSUMED as a counting exercise only — Sakharov induced G is UV.
    δ per converted vector: +1/2 − (−2/3) = +7/6.
    Three massive (W+,W−,Z): Δ = 7/2.
    a1_IR = a1_UV + 7/2 = −43/12 + 42/12 = −1/12  (still negative).
    """
    n_proca = 3
    n_massless = 9  # 8 gluons + γ
    a1 = (
        4 * A1_SCALAR_MIN
        + 45 * A1_WEYL
        + n_proca * A1_VECTOR_PROCA
        + n_massless * A1_VECTOR_MASSLESS
    )
    delta = n_proca * (A1_VECTOR_PROCA - A1_VECTOR_MASSLESS)
    assert delta == F(7, 2)
    assert a1 == A1_SM_XI0 + delta
    assert a1 == F(-1, 12)
    return {
        "a1_UV": A1_SM_XI0,
        "a1_IR_Proca": a1,
        "delta": delta,
        "sign_OK_IR": a1 > 0,
        "sign_OK_UV": A1_SM_XI0 > 0,
        "label": "derived arithmetic; IR count assumed irrelevant for UV Sakharov",
    }


def chiral_edge_mode_test() -> dict:
    """Can Dirac L–R / massive edge modes flip a1 without tuning?

    Donnelly–Wall: massless U(1) edge modes ↔ Kabat contact.
    Massive Proca: no residual gauge ⇒ no Kabat (already +1/2 vs −2/3).
    Dirac L–R locking is bulk mass; it does not remove photon/gluon contacts.

    Two-phase 'chiral-mixing boundary' that shares or cancels massless
    Kabat contacts is NOT derived from DW or Dirac mass.
    Shared-Kabat from #28 (Z3) gave +53/12 but was assumed — same status
    for any two-phase rewrite of that assumption.
    """
    # If all 12 vectors were Proca (no massless gauge — not SM): 
    a1_all_proca = 4 * A1_SCALAR_MIN + 45 * A1_WEYL + 12 * A1_VECTOR_PROCA
    # 2/3 + 15/4 + 6 = 0.667+3.75+6 = 10.417 = 125/12
    assert a1_all_proca == F(2, 3) + F(15, 4) + 6
    assert a1_all_proca == F(125, 12)
    return {
        "Dirac_LR_removes_photon_Kabat": False,
        "Proca_removes_own_Kabat": True,  # derived
        "a1_all_Proca_fantasy": a1_all_proca,
        "a1_all_Proca_is_SM": False,
        "two_phase_boundary_derived_flip": False,
        "shared_Kabat_Z3_from_28": F(53, 12),
        "shared_Kabat_label": "assumed (#28); not promoted by two-phase",
        "plain": (
            "IR Proca replacement leaves a1=−1/12 (still wrong sign). "
            "Only dropping all massless gauge contacts (not SM) or assumed "
            "shared-Kabat flips the sign. Two-phase chiral mixing does not "
            "derive a boundary cancelation for photons/gluons."
        ),
    }


def induced_sign_summary() -> dict:
    ir = a1_ir_proca_replacement()
    edge = chiral_edge_mode_test()
    return {
        "UV_a1": a1_uv_sm(),
        "UV_sign_OK": False,
        "IR_Proca": ir,
        "edge": edge,
        "derived_flip": False,
    }


# ---------------------------------------------------------------------------
# (4) Two-phase vs three-phase
# ---------------------------------------------------------------------------

def compare_two_vs_three() -> dict:
    return {
        "three_phase_28": {
            "claim": "E, B, gravity as 120° phases of one wave",
            "derived": False,
            "photon_issue": "decoupling gravity for photons clashes with Cassini",
            "a1_flip": "shared-Kabat +53/12 assumed",
        },
        "two_phase_this": {
            "claim": "light=E/B; matter=locked mass↔gravity via T^μ_μ",
            "derived_pieces": (
                "Dirac L–R lock; Compton clock; T^μ_μ=0 for EM; "
                "universal metric + optional dilaton"
            ),
            "not_derived": "mass↔gravity as phase-locked pair; EP as phase lock",
            "photon_issue": "avoided if metric universal (Cassini OK)",
            "a1_flip": "IR Proca → −1/12 still negative; no derived flip",
            "extra_trace_bound": f"η ≲ {MICROSCOPE_ETA:g} (MICROSCOPE)",
        },
        "improvement": (
            "Two-phase avoids the light-bending clash by construction when "
            "metric coupling stays universal. It still does not derive a "
            "mass–gravity phase lock or flip SM a1 without tuning."
        ),
        "gap_h2zK": "neither version touches the FFZ h²z/K ledger gap",
    }


# ---------------------------------------------------------------------------
# Self-checks
# ---------------------------------------------------------------------------

def self_checks() -> None:
    assert a1_uv_sm() == A1_SM_XI0 == F(-43, 12)
    ir = a1_ir_proca_replacement()
    assert ir["a1_IR_Proca"] == F(-1, 12)
    assert ir["sign_OK_IR"] is False
    edge = chiral_edge_mode_test()
    assert edge["two_phase_boundary_derived_flip"] is False
    assert edge["a1_all_Proca_fantasy"] == F(125, 12)
    tp = two_phase_formalization()
    assert tp["mass_gravity_locked_pair_derived"] is False
    assert tp["trace"]["T_trace_EM"] == 0
    b = fifth_force_bounds()
    assert b["closes_Newton_gap"] is False
    assert b["MICROSCOPE_eta"] == MICROSCOPE_ETA
    c = cassini_safe_version()
    assert c["light_bending_OK"] is True
    cmp_ = compare_two_vs_three()
    assert "avoided" in cmp_["two_phase_this"]["photon_issue"]


def main():
    print("=== Two-phase mass/gravity (Cassini-safe) ===\n")
    print("No tuning. Stacked on #28. Universal metric + mass-only trace extra.\n")

    print("-- (1) Formalization --")
    tp = two_phase_formalization()
    print(f"  Light: {tp['light']}")
    print(f"  Matter: {tp['matter']}")
    print(f"  Gravity phase on: {tp['gravity_phase_on']}")
    print(f"  Mass↔gravity locked pair derived? {tp['mass_gravity_locked_pair_derived']}")
    print(f"  EP as two-phase lock derived? {tp['EP_locking_derived_as_phases']}")
    d = tp["dirac"]
    print(f"  Dirac L–R: {d['Dirac_mass_term']} → locked {d['locked_pair']} [{d['label_LR']}]")
    print(f"    as mass↔gravity? {d['is_mass_gravity_pair']} [{d['label_as_gravity']}]")
    print(f"  Compton clock: {tp['de_broglie']['omega_C']} [{tp['de_broglie']['label']}]")
    print(f"  EP m_i=m_g geometric? {tp['EP']['derived_as_geometry']} [{tp['EP']['label_GR']}]")
    print(f"  Honest: {tp['honest']}")

    print("\n-- (2) Cassini-safe + fifth-force bounds --")
    c = cassini_safe_version()
    print(f"  Metric universal / light bending OK? {c['metric_universal']}/{c['light_bending_OK']}")
    print(f"  Extra: {c['extra_trace_coupling']} [{c['label_split']}]")
    b = fifth_force_bounds()
    print(f"  MICROSCOPE |η| ≲ {b['MICROSCOPE_eta']:.0e}")
    print(f"  Eöt-Wash |η| ≲ {b['EotWash_eta']:.0e}")
    print(f"  Cassini ⇒ ω_BD ≳ {b['omega_BD_min_Cassini']:.1e}; universal scalar admixture ≲ {b['scalar_admixture_max_universal']:.2e}")
    print(f"  Closes h²z/K gap? {b['closes_Newton_gap']} [{b['label_gap']}]")

    print("\n-- (3) Induced a1 sign --")
    sig = induced_sign_summary()
    print(f"  UV SM a1 = {sig['UV_a1']}  sign_OK={sig['UV_sign_OK']}")
    ir = sig["IR_Proca"]
    print(f"  IR Proca (3×W/Z): a1 = {ir['a1_IR_Proca']}  Δ={ir['delta']}  sign_OK={ir['sign_OK_IR']} [{ir['label']}]")
    e = sig["edge"]
    print(f"  Dirac L–R removes photon Kabat? {e['Dirac_LR_removes_photon_Kabat']}")
    print(f"  Proca removes own Kabat? {e['Proca_removes_own_Kabat']}")
    print(f"  All-Proca fantasy a1 = {e['a1_all_Proca_fantasy']} (not SM)")
    print(f"  Two-phase boundary derived flip? {e['two_phase_boundary_derived_flip']}")
    print(f"  #28 shared-Kabat {e['shared_Kabat_Z3_from_28']}: [{e['shared_Kabat_label']}]")
    print(f"  Plain: {e['plain']}")

    print("\n-- (4) Two-phase vs three-phase --")
    cmp_ = compare_two_vs_three()
    print(f"  Three-phase: derived={cmp_['three_phase_28']['derived']}; {cmp_['three_phase_28']['photon_issue']}")
    print(f"  Two-phase: photon issue {cmp_['two_phase_this']['photon_issue']}")
    print(f"  Two-phase a1: {cmp_['two_phase_this']['a1_flip']}")
    print(f"  Improvement: {cmp_['improvement']}")
    print(f"  Gap: {cmp_['gap_h2zK']}")

    print("\n-- Self-checks --")
    self_checks()
    print("  all asserts PASS")

    print(
        f"""
Verdict
-------
Dirac L–R and Compton clocks formalize mass as a 2-phase/phase structure;
they do not derive a mass↔gravity locked wave pair. EP m_i=m_g is geometric
in GR. Cassini-safe version: universal metric + optional trace dilaton with
|η| ≲ {MICROSCOPE_ETA:.0e}. UV a1=−43/12; IR Proca a1=−1/12 — still wrong
sign. No derived two-phase boundary flip (shared-Kabat remains assumed).
Better than three-phase on photons; still negative for a1 and the ledger gap.
"""
    )
    print("=== TWO-PHASE COMPLETE ===")


if __name__ == "__main__":
    main()
