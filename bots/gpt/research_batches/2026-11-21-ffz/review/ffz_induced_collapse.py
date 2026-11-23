#!/usr/bin/env python3
"""Penrose–Diósi gravity-induced collapse with induced G (Sakharov, κ=0).

Standard DP: τ = ℏ / E_G with E_G the gravitational self-energy of the
mass-density difference. Induced-gravity version: replace G → G_ind from
the Newton matching coefficient (ledger C≈0.0513 vs target 1/24).

No tuning. DERIVED vs ASSUMED labeled. Gap 0.0097 kept separate unless a
real link is derived (the C_ledger vs 1/24 ratio is a matching comparison,
not an explanation of the gap).

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_induced_collapse.py

Requires: numpy (stdlib math sufficient for core; numpy unused but ok)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# Physical constants (CODATA-ish)
G_N = 6.67430e-11  # m^3 kg^-1 s^-2
HBAR = 1.054571817e-34  # J s
U_ATOMIC = 1.66053906660e-27  # kg
M_ELECTRON = 9.1093837015e-31  # kg

# FFZ Newton matching coefficients (units h² z/K), κ=0 scheme
C_LEDGER = 0.05131719  # R* / standing ledger simple total
C_TARGET = 1.0 / 24.0  # 1/24 Y-cancelation / matching target
# NOTE: C_LEDGER - C_TARGET ≈ 0.00965 is the Newton gap. Using C_LEDGER
# here is a G_ind matching choice, NOT a derivation of that gap.


# ---------------------------------------------------------------------------
# (1) Standard Diósi–Penrose
# ---------------------------------------------------------------------------

def lambda_sep(d: float, R0: float) -> float:
    """λ = d / (2 R0)."""
    return d / (2.0 * R0)


def collapse_rate_Lambda(M: float, R0: float, d: float, G: float = G_N) -> float:
    """DP collapse rate Λ(d) for a uniform sphere (Wikipedia / standard).

    Λ = [6 G M² / (5 R0 ℏ)] × f(λ),  λ=d/(2R0)

    f(λ) = (5/3)λ² − (5/4)λ³ + (1/6)λ⁵    for 0 ≤ λ ≤ 1
    f(λ) = 1 − 5/(12λ)                     for λ ≥ 1

    DERIVED: textbook DP formula for uniform sphere.
    """
    if M <= 0 or R0 <= 0 or d < 0:
        raise ValueError("M, R0 > 0 and d ≥ 0 required")
    pre = 6.0 * G * M**2 / (5.0 * R0 * HBAR)
    lam = lambda_sep(d, R0)
    if lam <= 1.0:
        f = (5.0 / 3.0) * lam**2 - (5.0 / 4.0) * lam**3 + (1.0 / 6.0) * lam**5
    else:
        f = 1.0 - 5.0 / (12.0 * lam)
    return pre * f


def E_G(M: float, R0: float, d: float, G: float = G_N) -> float:
    """Gravitational self-energy of the density difference: E_G = ℏ Λ.

    DERIVED from τ = ℏ/E_G and Λ above.
    """
    return HBAR * collapse_rate_Lambda(M, R0, d, G=G)


def tau_DP(M: float, R0: float, d: float, G: float = G_N) -> float:
    """Collapse time τ = ℏ / E_G = 1/Λ."""
    return 1.0 / collapse_rate_Lambda(M, R0, d, G=G)


def tau_far(M: float, R0: float, G: float = G_N) -> float:
    """d ≫ R0 limit: τ = 5 R0 ℏ / (6 G M²). DERIVED."""
    return 5.0 * R0 * HBAR / (6.0 * G * M**2)


def gran_sasso_2020_note() -> dict:
    """Donadi et al., Nature Physics 2021 (Gran Sasso 2020 run).

    Ruled out the parameter-free DP version that takes R0 = nuclear
    positional spread (~few × 10^{-12} m for Ge). Lower bound
    R0 ≳ 5.4 × 10^{-11} m (95% CL) from non-observation of collapse-induced
    radiation.

    ASSUMED: we quote the published bound; not re-derived here.
    """
    return {
        "experiment": "Donadi et al., Nat. Phys. (2021); Gran Sasso",
        "parameter_free_R0": "~nuclear spread ~5e-12 m (Ge) — RULED OUT",
        "R0_lower_bound_m": 5.4e-11,
        "CL": "95%",
        "label": "assumed (published experimental result; not re-derived)",
    }


# ---------------------------------------------------------------------------
# (2) Induced-gravity version (κ=0, G → G_ind)
# ---------------------------------------------------------------------------

def g_ind_scaling(C: float, C_ref: float = C_TARGET) -> dict:
    """How G_ind scales with the Newton matching coefficient C.

    Packet: Γ_ind ⊃ −R/(16π G_ind). The simple-pole residue in units
    h²z/K is the coefficient C in F_R ~ C · h²z/(Kε). Matching identifies
    that residue with the EH coefficient, so

        1/G_ind ∝ C     (at this matching order / channel)
        G_ind   ∝ 1/C

    DERIVED: from Γ ⊃ −R/(16πG) and F_R ↔ C.
    ASSUMED: that the three-loop matching piece dominates the Newtonian
    G relevant for E_G (or that we compare only this piece's contribution).

    Therefore:
        E_G ∝ G  ⇒  E_G(C) / E_G(C_ref) = C_ref / C
        τ   ∝ 1/E_G  ⇒  τ(C) / τ(C_ref) = C / C_ref
    """
    if C <= 0 or C_ref <= 0:
        raise ValueError("C, C_ref must be positive")
    ratio_G = C_ref / C
    ratio_EG = ratio_G
    ratio_tau = C / C_ref
    return {
        "C": C,
        "C_ref": C_ref,
        "G_ind_over_G_ref": ratio_G,
        "E_G_over_E_G_ref": ratio_EG,
        "tau_over_tau_ref": ratio_tau,
        "label_scaling": "derived (1/G ∝ C)",
        "label_dominance": (
            "assumed: three-loop C dominates Newtonian G_ind for E_G, "
            "or comparison is only of this matching piece"
        ),
    }


def induced_vs_target() -> dict:
    """Ledger C vs 1/24: numerical τ ratio."""
    sc = g_ind_scaling(C_LEDGER, C_TARGET)
    return {
        **sc,
        "C_ledger": C_LEDGER,
        "C_1_over_24": C_TARGET,
        "tau_induced_over_tau_1_24": sc["tau_over_tau_ref"],
        "note_gap": (
            "C_ledger − 1/24 ≈ +0.00965 is the Newton gap. Using C_ledger "
            "here compares two matching choices for G_ind; it does NOT "
            "derive or explain the gap."
        ),
        "label": "derived ratio from C_ledger/C_target; gap link not derived",
    }


def induced_R0_smearing() -> dict:
    """Does induced gravity naturally supply the DP smearing length R0?

    DP needs a finite R0: point masses make E_G diverge. Candidates from
    induced gravity / Sakharov:
      * Planck length ℓ_P = √(ℏ G_ind/c³) ~ 10^{-35} m for G~G_N
      * Matter-loop Compton scale ħc/√z (model-dependent; particle-physics
        scales ≪ nuclear size if √z is GeV–TeV)

    Gran Sasso requires R0 ≳ 5.4×10^{-11} m. Nuclear size ~10^{-15}–10^{-14} m
    is already ruled out as the parameter-free choice.

    Plain answer: induced-gravity UV cutoffs do NOT naturally equal the
    R0 DP needs to survive experiment. Supplying R0 still requires an
    extra ASSUMED smearing (the known weak point of DP).
    """
    c_light = 2.99792458e8
    ell_P = math.sqrt(HBAR * G_N / c_light**3)
    return {
        "ell_P_G_N_m": ell_P,
        "nuclear_R0_m": 5e-15,
        "gran_sasso_R0_min_m": 5.4e-11,
        "naturally_supplies_viable_R0": False,
        "plain": (
            "No: ℓ_P and typical matter-loop scales are far below the "
            "Gran Sasso lower bound; induced gravity does not derive the "
            "smearing length DP needs."
        ),
        "label": "derived (scale comparison); R0 still assumed in DP",
    }


# ---------------------------------------------------------------------------
# (3) Condensate / phase-lock estimates
# ---------------------------------------------------------------------------

@dataclass
class CollapseCase:
    name: str
    M: float
    R0: float
    d: float
    tau_G_N: float
    tau_induced: float
    thermal_or_coherence_s: float
    coherence_note: str
    gravity_dominates: bool
    labels: str


def case_bec_rb87() -> CollapseCase:
    """(a) BEC of ~1e6 Rb-87 atoms split over ~1 mm.

    M = N m_Rb. R0 = cloud radius ~10 μm (ASSUMED typical BEC size).
    d = 1 mm. Observed BEC interferometry coherence: ~0.1–1 s typical
    (ASSUMED order-of-magnitude; can be shorter with noise).
    """
    N = 1.0e6
    m = 87.0 * U_ATOMIC
    M = N * m
    R0 = 1.0e-5  # 10 μm ASSUMED
    d = 1.0e-3
    tau_n = tau_DP(M, R0, d, G=G_N)
    ratio = induced_vs_target()["tau_induced_over_tau_1_24"]
    # Induced: scale G → G_ind relative to G_N using C_target as if G_N↔1/24
    # ASSUMED identification: G_N ↔ matching with C=1/24 for this comparison
    tau_ind = tau_n * ratio
    t_coh = 1.0  # seconds, ASSUMED order
    return CollapseCase(
        name="BEC ~1e6 Rb-87, d=1 mm, R0=10 μm",
        M=M,
        R0=R0,
        d=d,
        tau_G_N=tau_n,
        tau_induced=tau_ind,
        thermal_or_coherence_s=t_coh,
        coherence_note="ASSUMED ~1 s BEC interferometry coherence scale",
        gravity_dominates=tau_ind < t_coh,
        labels=(
            "DERIVED τ from DP formula; ASSUMED R0=10 μm cloud, "
            "coherence~1 s, G_N↔1/24 identification for induced ratio"
        ),
    )


def case_flux_qubit() -> CollapseCase:
    """(b) Superconducting flux-qubit / Cooper-pair condensate.

    Crude model: N_eff~1e9 Cooper pairs (ASSUMED), mass 2 m_e each,
    spatial separation ~1 μm (loop scale, ASSUMED), R0~1 μm.
    Real flux-qubit DP estimates are subtle (tiny mass displacement);
    this is an upper-bound-ish spatial-mass model. Coherence ~10–100 μs.
    """
    N_eff = 1.0e9
    M = N_eff * 2.0 * M_ELECTRON
    R0 = 1.0e-6
    d = 1.0e-6
    tau_n = tau_DP(M, R0, d, G=G_N)
    ratio = induced_vs_target()["tau_induced_over_tau_1_24"]
    tau_ind = tau_n * ratio
    t_coh = 5.0e-5  # 50 μs ASSUMED
    return CollapseCase(
        name="Flux qubit / Cooper pairs N_eff=1e9, d=R0=1 μm",
        M=M,
        R0=R0,
        d=d,
        tau_G_N=tau_n,
        tau_induced=tau_ind,
        thermal_or_coherence_s=t_coh,
        coherence_note="ASSUMED ~50 μs superconducting qubit coherence",
        gravity_dominates=tau_ind < t_coh,
        labels=(
            "DERIVED τ from DP; ASSUMED N_eff, R0, d, coherence; "
            "real flux mass displacement may be smaller ⇒ even larger τ"
        ),
    )


def case_nanosphere() -> CollapseCase:
    """(c) ~1e-14 kg nanosphere.

    R0 from silica-like density ρ~2e3 kg/m³ (ASSUMED).
    d = R0 (partial overlap) and also report d≫R0 limit.
    Optomechanical coherence: highly setup-dependent; ~1 ms ballpark
    at good isolation (ASSUMED comparison scale).
    """
    M = 1.0e-14
    rho = 2.0e3
    V = M / rho
    R0 = (3.0 * V / (4.0 * math.pi)) ** (1.0 / 3.0)
    d = R0  # λ = 1/2
    tau_n = tau_DP(M, R0, d, G=G_N)
    ratio = induced_vs_target()["tau_induced_over_tau_1_24"]
    tau_ind = tau_n * ratio
    t_coh = 1.0e-3  # 1 ms ASSUMED
    return CollapseCase(
        name=f"Nanosphere M=1e-14 kg, R0={R0:.3e} m, d=R0",
        M=M,
        R0=R0,
        d=d,
        tau_G_N=tau_n,
        tau_induced=tau_ind,
        thermal_or_coherence_s=t_coh,
        coherence_note="ASSUMED ~1 ms optomechanical coherence scale",
        gravity_dominates=tau_ind < t_coh,
        labels=(
            "DERIVED R0 from M/ρ and τ from DP; ASSUMED ρ=2e3, d=R0, "
            "coherence~1 ms"
        ),
    )


# ---------------------------------------------------------------------------
# (4) Phase framing
# ---------------------------------------------------------------------------

def phase_framing() -> dict:
    """Collapse criterion as accumulated relative phase E_G t / ℏ ~ 1.

    τ = ℏ/E_G is exactly the time when the gravitational relative phase
    between branches reaches O(1). This is a rewrite of the same criterion
    in Dan's phase/frequency/amplitude language: the two geometries accrue
    a relative phase at angular rate E_G/ℏ.

    Content: interpretive identification of collapse with gravitational
    phase mismatch; same quantitative prediction as τ=ℏ/E_G. No new
    dynamical law beyond DP unless one adds a concrete stochastic process
    (Diósi master equation) — that step is ASSUMED extra.
    """
    return {
        "criterion": "E_G t / ℏ ∼ 1  ⇔  t ∼ τ = ℏ/E_G",
        "rewrite_or_content": "rewrite of the same DP criterion",
        "extra_content": (
            "Identifies collapse with relative gravitational phase between "
            "branches (phase/frequency view). Does not by itself add new "
            "dynamics beyond DP; Diósi Lindblad dynamics is an extra assumption."
        ),
        "label": "derived equivalence; dynamics beyond τ ASSUMED",
    }


# ---------------------------------------------------------------------------
# Self-checks / test hooks
# ---------------------------------------------------------------------------

def self_checks() -> None:
    """Assertions for CI / unit test."""
    # Far-limit consistency: λ→∞ ⇒ τ → 5 R0 ℏ/(6 G M²)
    M, R0, d = 1.0e-12, 1.0e-6, 1.0  # λ ≫ 1
    t_far = tau_far(M, R0)
    t_dp = tau_DP(M, R0, d)
    assert abs(t_dp / t_far - 1.0) < 1e-6

    # Continuity at λ=1
    d1 = 2.0 * R0
    # evaluate both branches via formula internals
    pre = 6.0 * G_N * M**2 / (5.0 * R0 * HBAR)
    f_lo = 5.0 / 3.0 - 5.0 / 4.0 + 1.0 / 6.0  # λ=1
    f_hi = 1.0 - 5.0 / 12.0
    assert abs(f_lo - f_hi) < 1e-15
    assert abs(collapse_rate_Lambda(M, R0, d1) / (pre * f_hi) - 1.0) < 1e-12

    # Proton textbook ballpark: τ ~ 10^6–10^7 years for d ≫ R0~1 fm
    # (Wikipedia quotes ~10^6 yr; exact 5 R ℏ/(6 G m_p²) at R=1 fm is ~1.5×10^7 yr)
    m_p = 1.6726219e-27
    R_p = 1.0e-15
    t_p = tau_far(m_p, R_p)
    yr = 365.25 * 24 * 3600
    assert 1e5 * yr < t_p < 5e7 * yr

    # Scaling: τ ∝ 1/G, τ ∝ C
    sc = g_ind_scaling(C_LEDGER, C_TARGET)
    assert abs(sc["tau_over_tau_ref"] - C_LEDGER / C_TARGET) < 1e-12
    assert abs(sc["G_ind_over_G_ref"] - C_TARGET / C_LEDGER) < 1e-12

    # Induced R0 verdict
    assert induced_R0_smearing()["naturally_supplies_viable_R0"] is False

    # Phase rewrite
    assert "rewrite" in phase_framing()["rewrite_or_content"]

    # Cases finite and positive
    for case in (case_bec_rb87(), case_flux_qubit(), case_nanosphere()):
        assert case.tau_G_N > 0 and case.tau_induced > 0
        assert abs(case.tau_induced / case.tau_G_N - C_LEDGER / C_TARGET) < 1e-9


def main():
    print("=== Penrose–Diósi collapse with induced G (κ=0) ===\n")
    print("No tuning. Gap kept separate unless a real link is derived.\n")

    print("-- (1) Standard DP --")
    print("  τ = ℏ/E_G, E_G = ℏ Λ(d) for uniform sphere [DERIVED textbook]")
    gs = gran_sasso_2020_note()
    print(f"  Gran Sasso 2020: {gs['parameter_free_R0']}")
    print(f"  R0 ≳ {gs['R0_lower_bound_m']:.1e} m ({gs['CL']})  [{gs['label']}]")
    m_p, R_p = 1.6726219e-27, 1e-15
    print(
        f"  Proton d≫R0: τ = {tau_far(m_p, R_p)/year:.3e} yr "
        f"(textbook ~1e6–1e7 yr ballpark) [DERIVED]"
    )

    print("\n-- (2) Induced gravity G → G_ind --")
    inv = induced_vs_target()
    print(f"  C_ledger = {inv['C_ledger']:.8f}, 1/24 = {inv['C_1_over_24']:.8f}")
    print(f"  G_ind/G_(1/24) = C_(1/24)/C_ledger = {inv['G_ind_over_G_ref']:.6f}  [DERIVED]")
    print(f"  E_G ratio = same = {inv['E_G_over_E_G_ref']:.6f}  [DERIVED]")
    print(
        f"  τ(induced)/τ(1/24) = C_ledger/(1/24) = "
        f"{inv['tau_induced_over_tau_1_24']:.6f}  [DERIVED]"
    )
    print(f"  Dominance of 3-loop C: [{g_ind_scaling(C_LEDGER)['label_dominance']}]")
    print(f"  Gap note: {inv['note_gap']}")
    r0 = induced_R0_smearing()
    print(f"  ℓ_P(G_N) = {r0['ell_P_G_N_m']:.3e} m")
    print(f"  Naturally supplies viable R0? {r0['naturally_supplies_viable_R0']}")
    print(f"  Plain: {r0['plain']}  [{r0['label']}]")

    print("\n-- (3) Condensate / phase-lock cases --")
    for case in (case_bec_rb87(), case_flux_qubit(), case_nanosphere()):
        print(f"  {case.name}")
        print(f"    M={case.M:.4e} kg, R0={case.R0:.4e} m, d={case.d:.4e} m")
        print(f"    τ(G_N)={case.tau_G_N:.4e} s, τ(induced)={case.tau_induced:.4e} s")
        print(
            f"    coherence scale={case.thermal_or_coherence_s:.4e} s "
            f"({case.coherence_note})"
        )
        print(f"    gravity dominates? {case.gravity_dominates}")
        print(f"    [{case.labels}]")

    print("\n-- (4) Phase framing --")
    ph = phase_framing()
    print(f"  criterion: {ph['criterion']}")
    print(f"  rewrite or content: {ph['rewrite_or_content']}")
    print(f"  {ph['extra_content']}  [{ph['label']}]")

    print("\n-- Self-checks --")
    self_checks()
    print("  all asserts PASS")

    print("\nVerdict")
    print("-------")
    print(
        "Induced G from C_ledger≈0.0513 makes τ longer than the 1/24-matched\n"
        "value by C_ledger/(1/24)≈1.232 (weaker G ⇒ slower collapse).\n"
        "Induced UV scales do not supply a Gran-Sasso-viable R0.\n"
        "BEC and flux-qubit estimates: gravity collapse ≪ thermal/coherence\n"
        "scales. Nanosphere ~1e-14 kg: τ~ms scale — setup-dependent whether\n"
        "it can dominate. Phase E_G t/ℏ~1 is a rewrite of τ=ℏ/E_G."
    )
    print("\n=== INDUCED COLLAPSE COMPLETE ===")


# module-level year for proton print in main
year = 365.25 * 24 * 3600


if __name__ == "__main__":
    main()
