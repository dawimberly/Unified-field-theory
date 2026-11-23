#!/usr/bin/env python3
"""Three-phase idea (E, B, gravity) + #22 vs #27 Majorana a1 reconciliation.

Stacked on SM-roster (#27).

(0) Reconcile N=1 vector multiplet a1:
      #22 used Solodukhin EE dof weight for the gaugino (+D/6=+1/3) → net −1/3
      #27 used Adler–Visser Sakharov weight (Majorana≡Weyl → +1/12) → net −7/12
    Derive which applies to induced 1/G (Sakharov) vs geometric entropy (EE).

(1) Formalize Dan's three-phase idea via GEM, Kaluza–Klein, Weyl duality, Z3.
    Check T^μ_μ=0 (EM traceless) vs light bending / EP.

(2) Z3 / boundary contact: can it flip SM a1=−43/12 to positive without tuning?

(3) Testable predictions vs Cassini |γ−1|<2.3e-5 and light deflection.

No tuning. DERIVED vs ASSUMED. Honest negatives.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_three_phase.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import pi, sqrt

# Standing SM a1 from #27
A1_SM_XI0 = F(-43, 12)
A1_SM_XI16 = F(-17, 4)
C_LEDGER = 0.05131719
C_TARGET = 1.0 / 24.0
GAP = C_LEDGER - C_TARGET

# Cassini Shapiro delay (Bertotti, Iess, Tortora 2003)
CASSINI_GAMMA_MINUS_1 = 2.3e-5  # |γ−1| < 2.3×10^{-5}


# ---------------------------------------------------------------------------
# (0) Majorana gaugino a1: #22 vs #27
# ---------------------------------------------------------------------------

def majorana_a1_reconciliation() -> dict:
    """Derive the correct Majorana weight for each observable.

    Observable A — geometric entropy / EE area law (Solodukhin):
      Every on-shell dof contributes +1/6 to the universal area coefficient;
      gauge fields get an extra Kabat contact −1. Fermions have no contact.
      Majorana: D=2 ⇒ a1_EE = +2/6 = +1/3.
      Maxwell+ghosts: −2/3.  Net N=1 multiplet_EE = −2/3 + 1/3 = −1/3.
      This is what #22 computed. Correct for EE / Kabat contact bookkeeping.

    Observable B — Sakharov induced 1/G (Adler–Visser / one-loop EA):
      Bracket C_AV = N0(1−6ξ) + N_Dirac − 4 N1, with N_Dirac = N_Weyl/2.
      In a1=C_AV/6 units: Weyl = +1/12. A 4d Majorana has the same on-shell
      content as one Weyl ⇒ a1_Sakharov(Majorana) = +1/12.
      Maxwell+ghosts: −2/3 (same HK). Net N=1 multiplet_Sak = −2/3+1/12 = −7/12.
      This is what #27 used. Correct for induced-G / Λ² channel.

    Why they differ: EE counts dofs with weight +1/6 each (all spins same sign
    in the area law). Sakharov weights fermions by the Dirac index in units of
    (1−6ξ), which is +1 per Dirac = +1/6 in a1 units = +1/12 per Weyl —
    four times smaller than the EE Majorana weight +1/3.

    For induced gravity (this program's Sakharov κ=0 reading): #27 is correct;
    #22's gaugino +1/3 and net −1/3 are the EE numbers misapplied to 1/G.
    Qualitative #22 conclusions still hold (Kabat not canceled; wrong-sign
    multiplet; Λ² ≠ h²z/K gap channel).
    """
    a1_vector = F(-2, 3)  # same both observables
    a1_maj_EE = F(1, 3)  # Solodukhin D/6
    a1_maj_Sak = F(1, 12)  # Adler–Visser Weyl weight
    net_EE = a1_vector + a1_maj_EE  # −1/3
    net_Sak = a1_vector + a1_maj_Sak  # −7/12
    assert net_EE == F(-1, 3)
    assert net_Sak == F(-7, 12)
    # Factor: EE/Sak = (1/3)/(1/12) = 4
    assert a1_maj_EE / a1_maj_Sak == 4
    return {
        "a1_vector_both": a1_vector,
        "a1_Majorana_EE_Solodukhin": a1_maj_EE,
        "a1_Majorana_Sakharov": a1_maj_Sak,
        "net_N1_EE_PR22": net_EE,
        "net_N1_Sakharov_PR27": net_Sak,
        "PR22_correct_for": "geometric entropy / EE / Kabat bookkeeping",
        "PR27_correct_for": "Sakharov induced 1/G (Λ² channel)",
        "for_induced_G_use": net_Sak,
        "PR22_error_if_read_as_1G": (
            "gaugino +1/3 is EE weight; for induced G use +1/12 ⇒ net −7/12 "
            "not −1/3. Qualitative conclusions of #22 unchanged."
        ),
        "label": "derived (observable split; same Maxwell −2/3)",
    }


# ---------------------------------------------------------------------------
# (1) Three-phase formalization
# ---------------------------------------------------------------------------

def gem_mashhoon() -> dict:
    """Gravitoelectromagnetism (Mashhoon): linearized GR → E_g, B_g.

    GEM writes ∇·E_g = −4πG ρ, ∇×B_g = −4πG/c² j + (1/c²)∂E_g/∂t, etc.
    E_g and B_g are metric potentials (g00, g0i), not a third phase of the
    Maxwell wave. Photons source the GEM fields via T^{μν} (energy), not via
    a 'gravity phase of the EM wave'.
    """
    return {
        "framework": "GEM (Mashhoon)",
        "E_B_G_as_three_phases": False,
        "what_it_is": "linearized GR rewritten as gravito-E/B; separate from Maxwell",
        "massless_decoupling": False,
        "label": "derived (GEM identity); three-phase reading assumed/rejected",
    }


def kaluza_klein() -> dict:
    """5D KK: g_MN → g_μν + A_μ + φ. Gravity + EM + dilaton from geometry.

    Not three phases of one 4D wave. The scalar φ couples to Tr T in some
    frames; photons (T^μ_μ=0) decouple from φ but still couple to g_μν.
    """
    return {
        "framework": "Kaluza–Klein 5D",
        "E_B_G_as_three_phases": False,
        "photon_dilaton_decoupling": True,  # T^μ_μ=0
        "photon_metric_decoupling": False,  # still curve spacetime
        "label": "derived (KK reduction); three-phase reading not supported",
    }


def weyl_duality() -> dict:
    """Duality / conformal: E' + iB' = e^{iα}(E + iB). Mixes E↔B only.

    A third 'gravity phase' at α=120° is not generated by Maxwell duality.
    Conformal/Weyl rescalings act on the metric, not as a third EM phase.
    """
    # 120° roots of unity
    phases = [0, 2 * pi / 3, 4 * pi / 3]
    return {
        "framework": "Weyl/conformal + EM duality rotations",
        "duality_mixes": "E and B only",
        "third_gravito_phase_from_duality": False,
        "z3_angles_rad": phases,
        "z3_from_Maxwell": False,
        "label_duality": "derived (EM duality is U(1) on F+i*F)",
        "label_z3_gravito": "assumed (not derived from Maxwell/Einstein)",
    }


def trace_coupling_unrealized_photons() -> dict:
    """Scalar-tensor / dilaton coupling ∝ T^μ_μ.

    EM stress tensor is traceless (Maxwell conformal): T^μ_μ(EM) = 0
    ⇒ dilaton/trace coupling vanishes for pure EM.
    Massive matter: T^μ_μ ∼ −ρ c² (nonrelativistic) ≠ 0 ⇒ coupling on.

    This realizes 'unrealized for photons' ONLY for the *trace* channel.
    Clash: photons DO gravitate via the full T^{μν} in GR (light bending,
    Shapiro delay). Equivalence principle: null geodesics feel the metric
    sourced by all T^{μν}, including EM energy.
    """
    # Maxwell: T_μν = F_μα F_ν^α − (1/4) g_μν F² ⇒ Tr = 0 in d=4
    T_trace_EM = 0
    # Light bending: GR predicts Δθ = 4GM/(c²b); observed (Eddington, VLBI)
    # Cassini: |γ−1| < 2.3e-5 ⇒ PPN γ≈1 (photons feel full curvature)
    return {
        "T_trace_EM": T_trace_EM,
        "trace_channel_decouples_photons": True,
        "full_Tmn_decouples_photons": False,
        "clash_light_bending": True,
        "clash_Shapiro": True,
        "clash_EP": True,
        "light_bending_GR": "4GM/(c²b)",
        "Cassini_|gamma-1|": CASSINI_GAMMA_MINUS_1,
        "label_trace": "derived (Maxwell traceless in d=4)",
        "label_clash": "derived (GR + experiment: photons gravitate)",
        "three_phase_as_trace_only": (
            "ASSUMED identification of 'gravity phase' with trace/dilaton "
            "channel; then photons unrealized there — but that is NOT GR "
            "gravity, and fails light-bending / Cassini if it replaces g_μν "
            "coupling."
        ),
    }


def z3_su3_analog() -> dict:
    """SU(3)/Z3 120° structure as a three-component analog.

    Color Z3 and EM duality U(1) are different groups; no derived map that
    places (E, B, gravito) at cube roots of unity as components of one wave.
    """
    omega = complex(-0.5, sqrt(3) / 2)  # e^{2πi/3}
    assert abs(omega**3 - 1) < 1e-12
    assert abs(1 + omega + omega**2) < 1e-12
    return {
        "z3_roots_sum_to_zero": True,
        "derived_from_Einstein_Maxwell": False,
        "label": "assumed analogy only; not derived",
    }


def three_phase_formalization() -> dict:
    gem = gem_mashhoon()
    kk = kaluza_klein()
    dual = weyl_duality()
    tr = trace_coupling_unrealized_photons()
    z3 = z3_su3_analog()
    # Can E,B,gravito be 120° apart with gravito off for massless?
    arranged = False  # not derived from existing frameworks
    return {
        "GEM": gem,
        "KK": kk,
        "duality": dual,
        "trace": tr,
        "Z3": z3,
        "E_B_G_120_derived": arranged,
        "honest": (
            "No existing framework (GEM, KK, EM duality, SU(3)/Z3) derives "
            "E,B,gravity as three phases of one wave at 120°. The only "
            "honest 'unrealized for photons' piece is trace/dilaton "
            "coupling (T^μ_μ=0), which clashes with light bending if taken "
            "as the full gravitational coupling."
        ),
    }


# ---------------------------------------------------------------------------
# (2) Induced-gravity sign test with Z3 / boundary contact
# ---------------------------------------------------------------------------

def z3_boundary_sign_test() -> dict:
    """Does a Z3 three-phase boundary arrangement flip SM a1 to positive?

    Kabat contact lives on the entangling surface (one negative scalar per
    U(1) gauge field). A Z3 arrangement of three 'phases' meeting at a
    boundary is not a derived SM structure.

    Try the only natural discrete options without a free continuous knob:

    (A) Replace 12 Kabat contacts by 12/3=4 'shared' contacts if three
        phases share one surface mode per triple: δa1 = +8*(+1)? 
        Each massless vector contributes −2/3 = −1/6 (phys) −1/6 (phys)
        −1 (Kabat). If Z3 shares one Kabat among 3 vectors:
        contact per vector −1/3 instead of −1 ⇒ a1_vec = 2/6 − 1/3 = 0.
        Then SM a1 = 4/6 + 45/12 + 12*0 = 2/3 + 15/4 = 41/12 > 0.
        BUT: sharing Kabat by 3 is ASSUMED (not derived from Maxwell or Z3
        gauge theory). Donnelly–Wall: each U(1) has its own edge Hilbert
        space; SU(3) has dim=8 independent contacts, not 8/3.

    (B) Multiply fermion weight by 3 (Z3 copies): still not enough vs
        vectors without touching Kabat, and is tuning/doubling.

    (C) Free boundary coefficient κ_b · (contact): continuous knob = tuning.
    """
    # Baseline SM ξ=0
    a1_base = A1_SM_XI0
    assert a1_base == F(-43, 12)

    # (A) shared-Kabat assumption: a1_vec = 2/6 − 1/3 = 0
    a1_vec_shared = F(2, 6) - F(1, 3)
    assert a1_vec_shared == 0
    a1_A = 4 * F(1, 6) + 45 * F(1, 12) + 12 * a1_vec_shared
    assert a1_A == F(2, 3) + F(15, 4)
    assert a1_A == F(53, 12)  # +4.416... positive

    # (B) triple fermion copies, standard Kabat
    a1_B = 4 * F(1, 6) + 3 * 45 * F(1, 12) + 12 * F(-2, 3)
    # 2/3 + 135/12 − 8 = 2/3 + 45/4 − 8
    assert a1_B == F(2, 3) + F(135, 12) - 8
    assert a1_B == F(47, 12)  # still positive, but field-doubling = tuning

    # (C) free κ_b: a1 = a1_base + 12 * κ_b * (something); need κ_b to cancel
    # To flip sign: need δ ≥ 43/12 from modifying 12 contacts.
    # Full removal of all Kabat (vectors at Ds/6 = +1/3 each):
    a1_no_kabat = 4 * F(1, 6) + 45 * F(1, 12) + 12 * F(2, 6)  # Ds/6 only
    assert a1_no_kabat == F(2, 3) + F(15, 4) + 4  # 12×(1/3)=4
    assert a1_no_kabat == F(101, 12)

    return {
        "a1_SM_baseline": a1_base,
        "option_A_shared_Kabat_a1": a1_A,
        "option_A_sign_OK": a1_A > 0,
        "option_A_label": "assumed (Kabat shared by Z3 triples; not DW/HK)",
        "option_B_triple_fermions_a1": a1_B,
        "option_B_label": "assumed (field doubling = tuning)",
        "option_C_free_boundary": "continuous κ_b — tuning by definition",
        "derived_without_assumption_flips_sign": False,
        "Donnelly_Wall": (
            "each U(1) has its own edge modes; SU(3) ⇒ 8 contacts, not 8/3"
        ),
        "plain": (
            "No derived Z3/boundary rearrangement flips SM a1 to positive. "
            "Shared-Kabat would, but is an unjustified assumption (tuning "
            "of contact multiplicity)."
        ),
    }


# ---------------------------------------------------------------------------
# (3) Testable predictions and bounds
# ---------------------------------------------------------------------------

def experimental_bounds() -> dict:
    """Predictions if 'gravity phase' = trace coupling replacing GR for photons.

    Would predict: reduced/absent light bending and Shapiro delay relative
    to GR, while massive bodies still attract. Ruled out by:

      * Cassini: |γ−1| < 2.3×10^{-5} (Bertotti et al. 2003)
      * VLBI light deflection: γ−1 consistent with 0 at ~10^{-4}
      * Lunar laser ranging / EP tests: composition-independent free fall

    Gravitational birefringence (different propagation for EM polarizations
    from a preferred 'phase' structure): CMB/astro bounds typically
    Δc/c ≲ 10^{-32}–10^{-37} for birefringence over cosmological paths
    (order-of-magnitude literature; not re-derived here — ASSUMED cited
    scale). No derived FFZ birefringence angle from Z3.
    """
    # If photons decouple from gravity entirely: effective γ→0 (PPN),
    # vs GR γ=1. Cassini rules out |γ−1| > 2.3e-5.
    gamma_GR = 1.0
    gamma_no_photon_coupling = 0.0  # schematic: no null deflection
    ruled_out = abs(gamma_no_photon_coupling - gamma_GR) > CASSINI_GAMMA_MINUS_1
    # Fractional room for a small 'phase' admixture ε in photon–metric
    # coupling: |ε| ≲ 2.3e-5 if it shifts γ by ε.
    eps_max = CASSINI_GAMMA_MINUS_1
    return {
        "Cassini_|gamma-1|_max": CASSINI_GAMMA_MINUS_1,
        "gamma_GR": gamma_GR,
        "photon_decoupling_ruled_out": ruled_out,
        "max_photon_metric_admixture": eps_max,
        "light_bending_status": "must match GR at ≲2×10^{-5} (Cassini γ)",
        "birefringence_derived_from_Z3": False,
        "birefringence_note": (
            "No derived Δn from three-phase; astro birefringence bounds "
            "exist in literature (ASSUMED citation scale only)"
        ),
        "massive_vs_light_EP": (
            "Any differential coupling massive vs photons is a composition/"
            "EP violation; LLR/EP tests constrain η ≲ 10^{-13} (ASSUMED "
            "published bound, not re-derived)"
        ),
        "label_Cassini": "assumed (published Bertotti+2003)",
        "label_prediction": "derived: full photon decoupling falsified",
    }


# ---------------------------------------------------------------------------
# Self-checks
# ---------------------------------------------------------------------------

def self_checks() -> None:
    r = majorana_a1_reconciliation()
    assert r["net_N1_EE_PR22"] == F(-1, 3)
    assert r["net_N1_Sakharov_PR27"] == F(-7, 12)
    assert r["for_induced_G_use"] == F(-7, 12)
    tp = three_phase_formalization()
    assert tp["E_B_G_120_derived"] is False
    assert tp["trace"]["T_trace_EM"] == 0
    assert tp["trace"]["clash_light_bending"] is True
    z = z3_boundary_sign_test()
    assert z["a1_SM_baseline"] == F(-43, 12)
    assert z["option_A_shared_Kabat_a1"] == F(53, 12)
    assert z["derived_without_assumption_flips_sign"] is False
    ex = experimental_bounds()
    assert ex["photon_decoupling_ruled_out"] is True
    assert ex["Cassini_|gamma-1|_max"] == 2.3e-5
    assert abs(GAP - 0.00965052) < 1e-5


def main():
    print("=== Three-phase idea + Majorana a1 reconciliation ===\n")
    print("No tuning. Stacked on #27.\n")

    print("-- (0) #22 vs #27 Majorana a1 --")
    r = majorana_a1_reconciliation()
    print(f"  Maxwell+ghosts (both):     {r['a1_vector_both']}")
    print(f"  Majorana EE (Solodukhin):  {r['a1_Majorana_EE_Solodukhin']}  → net N=1 = {r['net_N1_EE_PR22']}  [#22]")
    print(f"  Majorana Sakharov:         {r['a1_Majorana_Sakharov']} → net N=1 = {r['net_N1_Sakharov_PR27']}  [#27]")
    print(f"  For induced 1/G use:       {r['for_induced_G_use']}")
    print(f"  #22 correct for: {r['PR22_correct_for']}")
    print(f"  #27 correct for: {r['PR27_correct_for']}")
    print(f"  [{r['label']}]")
    print(f"  Correction note: {r['PR22_error_if_read_as_1G']}")

    print("\n-- (1) Three-phase formalization --")
    tp = three_phase_formalization()
    print(f"  GEM: three phases? {tp['GEM']['E_B_G_as_three_phases']} — {tp['GEM']['what_it_is']}")
    print(f"  KK: three phases? {tp['KK']['E_B_G_as_three_phases']}; photon dilaton off? {tp['KK']['photon_dilaton_decoupling']}; metric off? {tp['KK']['photon_metric_decoupling']}")
    print(f"  Duality: gravito from E↔B? {tp['duality']['third_gravito_phase_from_duality']} [{tp['duality']['label_duality']}]")
    print(f"  Z3 from Einstein–Maxwell? {tp['Z3']['derived_from_Einstein_Maxwell']} [{tp['Z3']['label']}]")
    tr = tp["trace"]
    print(f"  T^μ_μ(EM)={tr['T_trace_EM']} ⇒ trace channel unrealized for photons? {tr['trace_channel_decouples_photons']} [{tr['label_trace']}]")
    print(f"  Clash light bending/Shapiro/EP? {tr['clash_light_bending']}/{tr['clash_Shapiro']}/{tr['clash_EP']} [{tr['label_clash']}]")
    print(f"  Cassini |γ−1|<{tr['Cassini_|gamma-1|']}")
    print(f"  E,B,G at 120° derived? {tp['E_B_G_120_derived']}")
    print(f"  Honest: {tp['honest']}")

    print("\n-- (2) Z3 / boundary sign test --")
    z = z3_boundary_sign_test()
    print(f"  SM a1 baseline: {z['a1_SM_baseline']}")
    print(f"  (A) shared Kabat: a1={z['option_A_shared_Kabat_a1']} sign_OK={z['option_A_sign_OK']} [{z['option_A_label']}]")
    print(f"  (B) triple fermions: a1={z['option_B_triple_fermions_a1']} [{z['option_B_label']}]")
    print(f"  (C) {z['option_C_free_boundary']}")
    print(f"  Derived flip without assumption? {z['derived_without_assumption_flips_sign']}")
    print(f"  DW: {z['Donnelly_Wall']}")
    print(f"  Plain: {z['plain']}")

    print("\n-- (3) Predictions vs experiment --")
    ex = experimental_bounds()
    print(f"  Cassini |γ−1| < {ex['Cassini_|gamma-1|_max']}")
    print(f"  Full photon–gravity decoupling ruled out? {ex['photon_decoupling_ruled_out']}")
    print(f"  Max admixture ≲ {ex['max_photon_metric_admixture']}")
    print(f"  Z3 birefringence derived? {ex['birefringence_derived_from_Z3']}")
    print(f"  [{ex['label_Cassini']}; {ex['label_prediction']}]")

    print("\n-- Self-checks --")
    self_checks()
    print("  all asserts PASS")

    print(
        """
Verdict
-------
Majorana: for induced 1/G use +1/12 (Sakharov); #22's +1/3 is EE.
N=1 multiplet for Sakharov = −7/12 (correct #27); #22's −1/3 is EE-only.
Three-phase E/B/G at 120° is NOT derived from GEM, KK, or EM duality.
Trace coupling unrealizes photons (T^μ_μ=0) but clashes with light bending;
Cassini |γ−1|<2.3e-5 rules out replacing metric coupling. Z3/boundary
does not flip SM a1 without an assumed shared-Kabat or free κ_b (tuning).
"""
    )
    print("=== THREE-PHASE COMPLETE ===")


if __name__ == "__main__":
    main()
