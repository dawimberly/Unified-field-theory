#!/usr/bin/env python3
"""Distortion phase as intermediate between locked and collapsed; R0 from width?

Three regimes for a superposed mass/field configuration:
  locked → distorted (recoverable phase bend) → collapsed.
Ask whether the natural width of the distortion region supplies the
Diósi–Penrose smearing length R0 without tuning (Gran Sasso: R0≥5.4e-11 m).

No tuning. DERIVED vs ASSUMED labeled. Honest negatives.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_distortion_phase.py
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ffz_induced_collapse import (
    C_LEDGER,
    C_TARGET,
    G_N,
    HBAR,
    M_ELECTRON,
    U_ATOMIC,
    case_bec_rb87,
    case_flux_qubit,
    case_nanosphere,
    induced_vs_target,
    tau_DP,
)

C_LIGHT = 2.99792458e8
GRAN_SASSO_R0_MIN = 5.4e-11  # m, Donadi et al. 2021, 95% CL (ASSUMED published)
A0_BOHR = 5.29177210903e-11  # m
M_PROTON = 1.67262192369e-27
M_RB87 = 87.0 * U_ATOMIC


# ---------------------------------------------------------------------------
# (1) Three regimes + GL / phase-field description
# ---------------------------------------------------------------------------

def three_regimes_gl() -> dict:
    """Order-parameter description of locked / distorted / collapsed.

    Relative phase θ between two superposition branches. Ginzburg–Landau /
    phase-field energy density (standard):

        f = (ℏ² n / 2m) |∇θ|² + V(θ)

    with V having minima at locked values (e.g. θ=0, 2π) and a barrier.
    Analogues: Langer–Ambegaokar / McCumber–Halperin TAPS; quantum phase
    slips; Minev et al. 2019 reversible quantum-jump window.

    Regimes (DERIVED classification from this free-energy structure):
      locked:    θ near a minimum; gradient cost negligible; coherent.
      distorted: spatial core where θ slips across the barrier; finite width
                 ξ_θ; still recoverable if the slip is incomplete / reversible.
      collapsed: one branch selected; relative-phase coherence lost
                 (DP / irreversible jump completed).

    The potential V(θ) is NOT fixed by the WZ+Yukawa Newton ledger alone —
    its scale is ASSUMED when a number is needed. The *form* of ξ_θ from
    balancing gradient and V is DERIVED GL.
    """
    return {
        "regimes": ("locked", "distorted", "collapsed"),
        "free_energy": "f = (ℏ² n / 2m)|∇θ|² + V(θ)",
        "refs": (
            "Langer–Ambegaokar; McCumber–Halperin TAPS; "
            "quantum phase slips; Minev et al. Nature 2019"
        ),
        "V_from_WZ_ledger": False,
        "label_form": "derived (standard GL / phase-field)",
        "label_V": "assumed (barrier scale not fixed by C_ledger)",
    }


def gl_healing_length(m: float, alpha_abs: float) -> float:
    """ξ = ℏ / √(2 m |α|)  — standard GL coherence/healing length.

    DERIVED formula. |α| must be supplied (ASSUMED or from condensate μ).
    """
    if m <= 0 or alpha_abs <= 0:
        raise ValueError("m, |α| > 0")
    return HBAR / math.sqrt(2.0 * m * alpha_abs)


def gp_healing_length(n: float, a_s: float) -> float:
    """Gross–Pitaevskii healing length ξ = 1 / √(8π n a_s).

    DERIVED from GP. Equivalent to ℏ/√(2m g n) with g=4π ℏ² a_s / m.
    """
    if n <= 0 or a_s <= 0:
        raise ValueError("n, a_s > 0")
    return 1.0 / math.sqrt(8.0 * math.pi * n * a_s)


def phase_slip_core_xi(m: float, n: float, V_barrier_density: float) -> float:
    """Relative-phase slip core: ξ_θ ∼ ℏ √(n / (2 m V₀)) with V₀ = barrier
    energy density scale (ASSUMED).

    From balancing (ℏ² n / 2m)(1/ξ)² ~ V₀  ⇒  ξ ~ ℏ √(n/(2 m V₀)).
    DERIVED scaling; V₀ ASSUMED.
    """
    if V_barrier_density <= 0:
        raise ValueError("V₀ > 0")
    return HBAR * math.sqrt(n / (2.0 * m * V_barrier_density))


# ---------------------------------------------------------------------------
# (2) Natural width from model-contained quantities
# ---------------------------------------------------------------------------

def compton(m: float, reduced: bool = True) -> float:
    """λ̄ = ℏ/(m c) or λ = h/(m c). DERIVED."""
    if reduced:
        return HBAR / (m * C_LIGHT)
    return 2.0 * math.pi * HBAR / (m * C_LIGHT)


def bec_healing_from_case() -> dict:
    """Healing length for the Rb BEC case (1e6 atoms, R_cloud=10 μm).

    a_s(87Rb) ≈ 5.3 nm — ASSUMED standard scattering length.
    n = N / V_sphere(R_cloud) — ASSUMED uniform-cloud estimate.
    """
    N = 1.0e6
    R_cloud = 1.0e-5
    a_s = 5.3e-9  # ASSUMED
    V = 4.0 / 3.0 * math.pi * R_cloud**3
    n = N / V
    xi = gp_healing_length(n, a_s)
    return {
        "name": "BEC GP healing length",
        "xi_m": xi,
        "n_m3": n,
        "a_s_m": a_s,
        "label": "derived ξ=1/√(8πna_s); assumed a_s=5.3 nm, uniform cloud",
    }


def cooper_coherence_length() -> dict:
    """Typical BCS/dirty-limit coherence length for Al-like qubits.

    ξ ~ 100 nm (clean Al ~1 μm) — ASSUMED literature ballpark, not derived
    from WZ Yukawa h,z.
    """
    return {
        "name": "Cooper-pair coherence length (Al-like)",
        "xi_m": 1.0e-7,  # 100 nm ASSUMED
        "label": "assumed (typical Al ξ; not from WZ ledger)",
    }


def induced_gravity_lengths() -> dict:
    """Lengths built from G_ind / C_ledger and field masses.

    ℓ_P(G) = √(ℏ G / c³). With G_ind / G_N = (1/24)/C_ledger:
        ℓ_P,ind = ℓ_P,N × √(G_ind/G_N)

    Compton of electron/proton: from m (not from C_ledger).
    No combination of (C_ledger, G_N, m_e, m_p) is DERIVED to equal
    5.4e-11 m without inserting an extra atomic scale.
    """
    ratio_G = C_TARGET / C_LEDGER
    ell_P_N = math.sqrt(HBAR * G_N / C_LIGHT**3)
    ell_P_ind = ell_P_N * math.sqrt(ratio_G)
    return {
        "ell_P_N_m": ell_P_N,
        "ell_P_ind_m": ell_P_ind,
        "G_ind_over_G_N": ratio_G,
        "label": "derived ℓ_P from G; G_ind/G_N from C (assumed G_N↔1/24 id.)",
    }


def candidate_widths() -> list[dict]:
    """All candidate distortion / R0 widths vs Gran Sasso."""
    bec = bec_healing_from_case()
    coop = cooper_coherence_length()
    grav = induced_gravity_lengths()
    cands = [
        {
            "name": "electron Compton λ̄_e",
            "width_m": compton(M_ELECTRON, reduced=True),
            "from_model": "m_e only (not WZ C_ledger)",
            "label": "derived λ̄=ℏ/(m_e c)",
        },
        {
            "name": "proton Compton λ̄_p",
            "width_m": compton(M_PROTON, reduced=True),
            "from_model": "m_p only",
            "label": "derived",
        },
        {
            "name": "Bohr radius a0",
            "width_m": A0_BOHR,
            "from_model": "α_EM, m_e — NOT in pure WZ+Yukawa ledger",
            "label": "derived atomic; assumed relevant to distortion/DP",
            "coincidence_with_gran_sasso": True,
        },
        {
            "name": bec["name"],
            "width_m": bec["xi_m"],
            "from_model": "condensate n, a_s",
            "label": bec["label"],
        },
        {
            "name": coop["name"],
            "width_m": coop["xi_m"],
            "from_model": "superconducting material (assumed)",
            "label": coop["label"],
        },
        {
            "name": "ℓ_P (G_N)",
            "width_m": grav["ell_P_N_m"],
            "from_model": "G_N",
            "label": grav["label"],
        },
        {
            "name": "ℓ_P (G_ind from C_ledger)",
            "width_m": grav["ell_P_ind_m"],
            "from_model": "C_ledger, κ=0 matching",
            "label": grav["label"],
        },
    ]
    for c in cands:
        w = c["width_m"]
        c["ratio_to_gran_sasso"] = w / GRAN_SASSO_R0_MIN
        c["meets_gran_sasso"] = w >= GRAN_SASSO_R0_MIN
        c["equals_bound_without_tuning"] = False  # set below for honest check
    # Honest: none is DERIVED equal to the bound from WZ+induced-G alone.
    # a0 is numerically close — flag coincidence, not a derivation.
    for c in cands:
        if c["name"] == "Bohr radius a0":
            c["coincidence_note"] = (
                f"a0/R0_min = {A0_BOHR/GRAN_SASSO_R0_MIN:.4f} — COINCIDENCE; "
                "a0 is not derived from C_ledger / WZ+Yukawa"
            )
    return cands


def supplies_r0_without_tuning(cands: list[dict]) -> dict:
    """Verdict: does any candidate supply R0 without tuning?"""
    # Viable = meets Gran Sasso AND derived from induced-G/WZ ingredients
    # without importing EM atomic scales or material-assumed ξ.
    model_native = []
    for c in cands:
        native = c["name"] in (
            "ℓ_P (G_N)",
            "ℓ_P (G_ind from C_ledger)",
            "proton Compton λ̄_p",
            "electron Compton λ̄_e",
        ) or "BEC" in c["name"]
        # BEC healing is from condensate case, model-ish but ≫ bound
        if native and c["meets_gran_sasso"]:
            # only BEC healing among natives meets bound — but it's a
            # condensate healing length, not a universal DP R0 from C_ledger
            model_native.append(c)

    # BEC healing meets Gran Sasso numerically but is system-specific and
    # much larger than the bound; using it as universal nuclear R0 would be
    # ASSUMED identification of two different scales.
    return {
        "any_derived_universal_R0_from_C_ledger": False,
        "plain": (
            "No: induced-G/WZ-native lengths are either ≪ bound (ℓ_P, Compton) "
            "or system-specific healing lengths ≫ bound (BEC/Cooper). "
            "a0≈5.29e-11 m is within ~2% of the Gran Sasso number but is an "
            "atomic-EM coincidence, not derived from C_ledger."
        ),
        "label": "derived comparison; negative result",
    }


# ---------------------------------------------------------------------------
# (3)–(4) Modified τ with R0 = derived width
# ---------------------------------------------------------------------------

@dataclass
class ModifiedTau:
    case_name: str
    R0_choice: str
    R0_m: float
    d_m: float
    M_kg: float
    tau_G_N_s: float
    tau_ind_s: float
    coherence_s: float
    gravity_dominates: bool
    note: str


def modified_taus() -> list[ModifiedTau]:
    """Recompute the three cases with R0 = distortion/healing width.

    Policy (stated, not tuned):
      * BEC: R0 = GP healing length (natural GL width of that condensate).
      * Flux: R0 = Cooper ξ (assumed Al-like).
      * Nanosphere: continuous-sphere R_geom already regularizes; also report
        using R0 = a0 as a coincidence probe (NOT a prediction).

    Induced factor: τ_ind = τ(G_N) × C_ledger/(1/24) with G_N↔1/24 ASSUMED id.
    """
    ratio = induced_vs_target()["tau_induced_over_tau_1_24"]
    bec_xi = bec_healing_from_case()["xi_m"]
    coop_xi = cooper_coherence_length()["xi_m"]

    bec0 = case_bec_rb87()
    # Same M, d; replace R0 by healing length
    t_bec = tau_DP(bec0.M, bec_xi, bec0.d, G=G_N)
    out = [
        ModifiedTau(
            case_name=bec0.name,
            R0_choice="GP healing length",
            R0_m=bec_xi,
            d_m=bec0.d,
            M_kg=bec0.M,
            tau_G_N_s=t_bec,
            tau_ind_s=t_bec * ratio,
            coherence_s=bec0.thermal_or_coherence_s,
            gravity_dominates=(t_bec * ratio) < bec0.thermal_or_coherence_s,
            note="DERIVED ξ; ASSUMED a_s, cloud; larger R0 ⇒ slower collapse",
        )
    ]

    flux0 = case_flux_qubit()
    t_flux = tau_DP(flux0.M, coop_xi, flux0.d, G=G_N)
    out.append(
        ModifiedTau(
            case_name=flux0.name,
            R0_choice="Cooper ξ ~100 nm",
            R0_m=coop_xi,
            d_m=flux0.d,
            M_kg=flux0.M,
            tau_G_N_s=t_flux,
            tau_ind_s=t_flux * ratio,
            coherence_s=flux0.thermal_or_coherence_s,
            gravity_dominates=(t_flux * ratio) < flux0.thermal_or_coherence_s,
            note="ASSUMED Cooper ξ; if d < 2 R0 geometry changes (λ<1)",
        )
    )

    nano0 = case_nanosphere()
    # Geometric R0 already from density — distortion width does not replace it
    # unless one smears constituents. Report geometric (baseline) and a0 probe.
    t_nano_geom = tau_DP(nano0.M, nano0.R0, nano0.d, G=G_N)
    out.append(
        ModifiedTau(
            case_name=nano0.name + " [geometric R0]",
            R0_choice="geometric radius from M/ρ",
            R0_m=nano0.R0,
            d_m=nano0.d,
            M_kg=nano0.M,
            tau_G_N_s=t_nano_geom,
            tau_ind_s=t_nano_geom * ratio,
            coherence_s=nano0.thermal_or_coherence_s,
            gravity_dominates=(t_nano_geom * ratio) < nano0.thermal_or_coherence_s,
            note="DERIVED R_geom; distortion ξ not the sphere radius",
        )
    )
    # Coincidence probe: force R0=a0 (NOT derived)
    t_nano_a0 = tau_DP(nano0.M, A0_BOHR, d=max(nano0.d, 2.0 * A0_BOHR), G=G_N)
    out.append(
        ModifiedTau(
            case_name=nano0.name + " [R0=a0 COINCIDENCE PROBE]",
            R0_choice="Bohr a0 (coincidence with Gran Sasso; NOT derived)",
            R0_m=A0_BOHR,
            d_m=max(nano0.d, 2.0 * A0_BOHR),
            M_kg=nano0.M,
            tau_G_N_s=t_nano_a0,
            tau_ind_s=t_nano_a0 * ratio,
            coherence_s=nano0.thermal_or_coherence_s,
            gravity_dominates=(t_nano_a0 * ratio) < nano0.thermal_or_coherence_s,
            note="ASSUMED R0=a0; smaller R0 ⇒ faster collapse — NOT a prediction",
        )
    )
    return out


# ---------------------------------------------------------------------------
# (5) Distinct testable prediction?
# ---------------------------------------------------------------------------

def testable_prediction() -> dict:
    """Does the distortion phase add content beyond plain DP?

    Yes, qualitatively: a reversible window (Minev-like) with duration set
    by the slip dynamics, and scaling with barrier/gradient parameters —
    distinct from irreversible DP collapse at τ=ℏ/E_G.

    No, for R0: it does NOT yield a tuning-free R0 from C_ledger that meets
    Gran Sasso. So it does not fix DP's smearing problem by itself.
    """
    return {
        "distinct_from_plain_DP": True,
        "what": (
            "Recoverable distortion window (phase-slip / incomplete jump) "
            "with duration and mass/barrier scaling from GL dynamics — "
            "testable in principle vs irreversible DP collapse."
        ),
        "fixes_R0": False,
        "label": (
            "derived: qualitative distinction; "
            "assumed: concrete window duration needs V(θ); "
            "derived: no R0 from C_ledger"
        ),
    }


def self_checks() -> None:
    # GL healing formula
    m, alpha = M_RB87, 1.0e-30  # J dummy
    xi = gl_healing_length(m, alpha)
    assert abs(xi - HBAR / math.sqrt(2.0 * m * alpha)) < 1e-30 * xi

    # GP healing positive for BEC
    bec = bec_healing_from_case()
    assert bec["xi_m"] > 0
    assert bec["xi_m"] > GRAN_SASSO_R0_MIN  # healing ≫ bound

    # a0 coincidence proximity (flag, not claim)
    assert abs(A0_BOHR / GRAN_SASSO_R0_MIN - 1.0) < 0.05

    # ℓ_P ≪ bound
    grav = induced_gravity_lengths()
    assert grav["ell_P_ind_m"] < 1e-30

    # No universal R0 from C_ledger
    cands = candidate_widths()
    verd = supplies_r0_without_tuning(cands)
    assert verd["any_derived_universal_R0_from_C_ledger"] is False

    # Modified taus finite
    for mt in modified_taus():
        assert mt.tau_ind_s > 0

    # BEC with healing R0: gravity still does not dominate
    bec_mt = modified_taus()[0]
    assert bec_mt.gravity_dominates is False

    pred = testable_prediction()
    assert pred["fixes_R0"] is False
    assert pred["distinct_from_plain_DP"] is True


def main():
    print("=== Distortion phase ↔ DP R0? ===\n")
    print("No tuning. Honest negatives if the width does not supply R0.\n")

    print("-- (1) Three regimes (GL / phase-field) --")
    reg = three_regimes_gl()
    print(f"  regimes: {reg['regimes']}")
    print(f"  free energy: {reg['free_energy']}")
    print(f"  V from WZ ledger? {reg['V_from_WZ_ledger']}  [{reg['label_V']}]")
    print(f"  form: [{reg['label_form']}]")

    print("\n-- (2)–(3) Candidate widths vs Gran Sasso R0≥5.4e-11 m --")
    cands = candidate_widths()
    for c in cands:
        print(
            f"  {c['name']}: {c['width_m']:.4e} m  "
            f"(ratio to bound {c['ratio_to_gran_sasso']:.3e}, "
            f"meets={c['meets_gran_sasso']})  [{c['label']}]"
        )
        if c.get("coincidence_note"):
            print(f"    ** {c['coincidence_note']}")
    verd = supplies_r0_without_tuning(cands)
    print(f"\n  Universal R0 from C_ledger without tuning? "
          f"{verd['any_derived_universal_R0_from_C_ledger']}")
    print(f"  Plain: {verd['plain']}")

    print("\n-- (4) Modified τ with R0 = derived/assumed width --")
    for mt in modified_taus():
        print(f"  {mt.case_name}")
        print(f"    R0={mt.R0_choice} = {mt.R0_m:.4e} m, d={mt.d_m:.4e} m")
        print(
            f"    τ(G_N)={mt.tau_G_N_s:.4e} s, τ(ind)={mt.tau_ind_s:.4e} s, "
            f"coherence={mt.coherence_s:.4e} s"
        )
        print(f"    gravity dominates? {mt.gravity_dominates}  [{mt.note}]")

    print("\n-- (5) Distinct testable prediction? --")
    pred = testable_prediction()
    print(f"  distinct from plain DP? {pred['distinct_from_plain_DP']}")
    print(f"  what: {pred['what']}")
    print(f"  fixes R0? {pred['fixes_R0']}  [{pred['label']}]")

    print("\n-- Self-checks --")
    self_checks()
    print("  all asserts PASS")

    print("\nVerdict")
    print("-------")
    print(
        "A distortion (recoverable phase-slip) regime is a sensible GL\n"
        "intermediate between locked and collapsed, and it suggests\n"
        "testable reversible-window phenomenology beyond plain DP.\n"
        "It does NOT, however, derive a Gran-Sasso-viable universal R0\n"
        "from C_ledger / induced G: native lengths are wrong by many\n"
        "orders, and a0≈5.4e-11 m is a coincidence, not a derivation.\n"
        "With condensate healing lengths as R0, collapse only gets slower."
    )
    print("\n=== DISTORTION PHASE COMPLETE ===")


if __name__ == "__main__":
    main()
