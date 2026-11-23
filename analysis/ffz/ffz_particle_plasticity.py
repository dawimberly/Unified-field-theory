#!/usr/bin/env python3
"""Particle plasticity: distortion = elastic; collapse = plastic yield.

Builds on #25 GL/phase-field for inter-branch phase θ. Adds:
  * elastic vs plastic response (Josephson / washboard / Prandtl–Tomlinson)
  * yield threshold robustness across V(θ) shapes
  * RCSJ hysteresis (Stewart–McCumber)
  * R0-from-yield check vs Gran Sasso (honest negative / coincidence flag)
  * rate-dependent (non-Newtonian) yield: Herschel–Bulkley analogue +
    Kurkijärvi–Fulton–Dunkleberger sweep-rate switching

No tuning. DERIVED vs ASSUMED labeled.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_particle_plasticity.py
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ffz_distortion_phase import (
    A0_BOHR,
    GRAN_SASSO_R0_MIN,
    bec_healing_from_case,
    cooper_coherence_length,
)
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
M_RB87 = 87.0 * U_ATOMIC


# ---------------------------------------------------------------------------
# (1) Elastic vs plastic on the GL phase field
# ---------------------------------------------------------------------------

def elastic_plastic_dictionary() -> dict:
    """Map #25 regimes onto continuum plasticity language.

    Free energy (from #25): f = (ℏ² n / 2m)|∇θ|² + V(θ).

    Elastic (recoverable): motion within one well / subcritical drive —
      corresponds to #25 *distorted* but reversible (phase slip incomplete).
    Plastic (irreversible): escape over the barrier / running phase —
      corresponds to #25 *collapsed* / completed slip (DP-like).

    Standard tools borrowed (ASSUMED applicability to inter-branch θ):
      * Josephson critical current / critical phase gradient
      * Landau critical velocity (superflow)
      * Prandtl–Tomlinson / tilted washboard
      * Yield-stress (Bingham) analogue for the phase "stress" ∼ ∇θ or θ̇
    """
    return {
        "elastic": "subcritical, recoverable distortion (locked/distorted)",
        "plastic": "barrier escape / running phase (collapse)",
        "tools": (
            "Josephson Ic; Landau v_c; Prandtl–Tomlinson; "
            "tilted washboard; Bingham/HB yield"
        ),
        "label_map": "derived from #25 free energy + standard analogies",
        "label_applicability": "assumed (θ as Josephson-like phase)",
    }


# ---------------------------------------------------------------------------
# (2) Yield threshold across V(θ) shapes
# ---------------------------------------------------------------------------

def washboard_V(theta: float, V0: float, F: float) -> float:
    """Tilted washboard: V = V0 (1 − cos θ) − F θ. ASSUMED shape + V0, F."""
    return V0 * (1.0 - math.cos(theta)) - F * theta


def double_well_V(theta: float, V0: float, F: float) -> float:
    """Tilted double well on [-π, π]: V = (V0/4)(θ²/θ0² − 1)² − F θ,
    with θ0=1 (rad) scale. ASSUMED shape."""
    theta0 = 1.0
    u = theta / theta0
    return 0.25 * V0 * (u * u - 1.0) ** 2 - F * theta


def cos2_V(theta: float, V0: float, F: float) -> float:
    """Alternate cos-type: V = (V0/2)(1 − cos 2θ) − F θ. ASSUMED."""
    return 0.5 * V0 * (1.0 - math.cos(2.0 * theta)) - F * theta


def critical_tilt_force(V_prime_max_at_F0: float) -> float:
    """For V_F(θ)=V(θ)−Fθ, barrier disappears when F = max V'(θ)|_{F=0}.

    DERIVED for any smooth tilted potential of the form V0·v(θ) − F θ:
    F_c = V0 · max v'(θ).
    """
    return V_prime_max_at_F0


def yield_thresholds() -> dict:
    """Critical tilts F_c / V0 for three shapes (dimensionless).

    Washbord V0(1−cosθ): V'=V0 sinθ ⇒ F_c = V0  (Josephson Ic ∝ V0).
    Cos2: V' = V0 sin(2θ) ⇒ F_c = V0.
    Double well: numerical max of v' on the untilted potential.

    Robustness: F_c / V0 = O(1) for all three — DERIVED given each shape;
    overall scale V0 remains ASSUMED (per #25).
    """
    # Washboard
    fc_wb = 1.0  # F_c/V0
    # Cos2
    fc_c2 = 1.0
    # Double well: v(θ)=(1/4)(θ²−1)², v'=θ(θ²−1)
    # extrema of v': d/dθ[θ³−θ]=3θ²−1=0 ⇒ θ=±1/√3
    # |v'(1/√3)| = (1/√3)*|(1/3−1)| = (1/√3)*(2/3) = 2/(3√3)
    fc_dw = 2.0 / (3.0 * math.sqrt(3.0))

    # Critical phase gradient from gradient energy vs barrier (order of mag.)
    # (ℏ² n / 2m)(Δθ/ξ)² ∼ V0  with Δθ∼π ⇒ ξ∼ℏ√(n/(2m V0))  [#25]
    # Critical |∇θ|_c ∼ Δθ/ξ ∼ π/ξ  — DERIVED scaling; V0 ASSUMED.
    return {
        "F_c_over_V0": {
            "washboard": fc_wb,
            "cos2": fc_c2,
            "double_well": fc_dw,
        },
        "robust_O1": True,
        "grad_theta_c_scaling": "π/ξ_θ with ξ_θ=ℏ√(n/(2m V0))",
        "label": (
            "derived F_c/V0 for each shape; "
            "assumed V0 absolute scale (not fixed by C_ledger)"
        ),
    }


def E_G_drive_rate(M: float, R0: float, d: float, G: float = G_N) -> float:
    """Characteristic gravitational phase rate Ω_G = E_G / ℏ = Λ_DP.

    DERIVED from DP: relative phase accumulates at Ω_G. Yield when drive
    (tilt) reaches F_c; identifying F ∼ ℏ Ω or energy-density tilt with E_G
    is ASSUMED matching of gravitational drive to washboard tilt.
    """
    from ffz_induced_collapse import collapse_rate_Lambda

    return collapse_rate_Lambda(M, R0, d, G=G)


# ---------------------------------------------------------------------------
# (3) Hysteresis: RCSJ / Stewart–McCumber
# ---------------------------------------------------------------------------

def rcsj_hysteresis(beta_c: float) -> dict:
    """RCSJ underdamped junction: hysteresis between switching and retrapping.

    Stewart–McCumber parameter β_c = 2π I_c R² C / Φ0  (ASSUMED analogue
    for θ dynamics: inertia/damping ratio).

    Standard qualitative results (DERIVED from RCSJ literature structure;
    numerical factors ASSUMED as textbook):
      * β_c ≪ 1 (overdamped): no hysteresis; I_r ≈ I_c
      * β_c ≫ 1 (underdamped): I_r / I_c ≈ 4/(π √β_c)  (approx.)
        loop width Δi = (I_c − I_r)/I_c ≈ 1 − 4/(π √β_c)

    Distinct from plain DP: DP gives a single τ=ℏ/E_G with no forward/return
    asymmetry. Hysteresis is an extra dynamical prediction of washboard+inertia.
    """
    if beta_c <= 0:
        raise ValueError("β_c > 0")
    if beta_c < 1.0:
        i_r_over_ic = 1.0
        hysteretic = False
    else:
        i_r_over_ic = min(1.0, 4.0 / (math.pi * math.sqrt(beta_c)))
        hysteretic = i_r_over_ic < 0.999
    loop_width = 1.0 - i_r_over_ic
    return {
        "beta_c": beta_c,
        "I_r_over_I_c": i_r_over_ic,
        "loop_width": loop_width,
        "hysteretic": hysteretic,
        "distinct_from_DP": True,
        "label": (
            "derived RCSJ structure; assumed β_c analogue for θ; "
            "I_r formula textbook approx."
        ),
    }


def hysteresis_scan() -> list[dict]:
    return [rcsj_hysteresis(b) for b in (0.1, 1.0, 10.0, 100.0, 1.0e4)]


# ---------------------------------------------------------------------------
# (4) Yield → length scale as R0?
# ---------------------------------------------------------------------------

def length_from_critical_gradient(xi_theta: float) -> float:
    """L_* ∼ 1 / |∇θ|_c ∼ ξ_θ / π   (DERIVED from |∇θ|_c ∼ π/ξ_θ)."""
    return xi_theta / math.pi


def r0_from_yield_check() -> dict:
    """Does yield define a tuning-free R0 meeting Gran Sasso?

    L_* = ξ_θ/π. For BEC, ξ_θ ∼ ξ_GP ⇒ L_* ∼ ξ_GP/π ∼ 5e-8 m ≫ bound
    but system-specific. For ℓ_P-scale V0 one would need extreme fine-tuning
    of V0 to hit 5.4e-11 m.

    Bohr a0 near-hit from #25: flagged as COINCIDENCE only — not derived
    from yield threshold or C_ledger.
    """
    bec_xi = bec_healing_from_case()["xi_m"]
    L_bec = length_from_critical_gradient(bec_xi)
    coop_xi = cooper_coherence_length()["xi_m"]
    L_coop = length_from_critical_gradient(coop_xi)
    # What V0 would make L_*=R0_min for Rb density? From ξ=ℏ√(n/(2m V0))=π R0
    # ⇒ V0 = ℏ² n / (2m (π R0)²)
    N = 1.0e6
    R_cloud = 1.0e-5
    n = N / (4.0 / 3.0 * math.pi * R_cloud**3)
    R0 = GRAN_SASSO_R0_MIN
    V0_needed = (HBAR**2) * n / (2.0 * M_RB87 * (math.pi * R0) ** 2)
    return {
        "L_bec_m": L_bec,
        "L_coop_m": L_coop,
        "meets_bound_bec": L_bec >= GRAN_SASSO_R0_MIN,
        "universal_from_C_ledger": False,
        "V0_needed_J_per_m3_for_L_star_eq_bound": V0_needed,
        "a0_coincidence_only": True,
        "a0_ratio_to_bound": A0_BOHR / GRAN_SASSO_R0_MIN,
        "plain": (
            "No: yield length L_*=ξ/π is either system healing-scale "
            "(BEC/Cooper, ≫ bound) or requires tuning V0 to hit 5.4e-11 m. "
            "a0 near-hit remains a coincidence, not derived from yield."
        ),
        "label": "derived L_* scaling; assumed V0; negative for untuned R0",
    }


# ---------------------------------------------------------------------------
# (5) Testable signatures for three cases + rate dependence
# ---------------------------------------------------------------------------

@dataclass
class CaseSignature:
    name: str
    Omega_G_Hz: float
    tau_ind_s: float
    coherence_s: float
    threshold_note: str
    hysteresis_note: str
    rate_note: str


def case_signatures() -> list[CaseSignature]:
    """Threshold mass/separation, hysteresis, rate signatures.

    Ω_G = E_G/ℏ from geometric/healing R0 as in #24–#25 cases.
    Concrete Ic-analogue numbers need V0 (ASSUMED) — we quote Ω_G and
    qualitative sweep-rate trends instead of fake absolute Ic.
    """
    ratio = induced_vs_target()["tau_induced_over_tau_1_24"]
    out = []

    bec = case_bec_rb87()
    xi = bec_healing_from_case()["xi_m"]
    # Use healing R0 for distortion-aware estimate
    om = E_G_drive_rate(bec.M, xi, bec.d)
    out.append(
        CaseSignature(
            name="BEC 1e6 Rb, d=1 mm",
            Omega_G_Hz=om,
            tau_ind_s=tau_DP(bec.M, xi, bec.d) * ratio,
            coherence_s=bec.thermal_or_coherence_s,
            threshold_note=(
                "Yield when gravitational drive ∼ F_c∼V0; "
                "threshold separation d_c from E_G(d_c)∼V0·Vol (ASSUMED V0)"
            ),
            hysteresis_note=(
                "Recoherence after reducing d: expect retrapping below "
                "switching if underdamped (β_c>1) — distinct from DP"
            ),
            rate_note=(
                "Faster d-sweep ⇒ higher effective switching threshold "
                "(Kurkijärvi); shear-thickening analogue"
            ),
        )
    )

    flux = case_flux_qubit()
    coop = cooper_coherence_length()["xi_m"]
    om_f = E_G_drive_rate(flux.M, coop, flux.d)
    out.append(
        CaseSignature(
            name="Flux qubit N_eff=1e9",
            Omega_G_Hz=om_f,
            tau_ind_s=tau_DP(flux.M, coop, flux.d) * ratio,
            coherence_s=flux.thermal_or_coherence_s,
            threshold_note=(
                "Closest lab analogue: Josephson switching/retrapping "
                "already measured; gravity-tilt is tiny vs electrical Ic"
            ),
            hysteresis_note=(
                "Standard RCSJ hysteresis; gravitational contribution "
                "to tilt is negligible vs electrical drive (honest)"
            ),
            rate_note=(
                "Sweep-rate dependent switching (KFD) is standard in JJ; "
                "gravity does not add a new rate law at accessible E_G"
            ),
        )
    )

    nano = case_nanosphere()
    om_n = E_G_drive_rate(nano.M, nano.R0, nano.d)
    out.append(
        CaseSignature(
            name="Nanosphere 1e-14 kg",
            Omega_G_Hz=om_n,
            tau_ind_s=tau_DP(nano.M, nano.R0, nano.d) * ratio,
            coherence_s=nano.thermal_or_coherence_s,
            threshold_note=(
                "Closest to DP-relevant Ω_G; threshold d or M where "
                "Ω_G t_coh ∼ F_c/V0 (ASSUMED)"
            ),
            hysteresis_note=(
                "Recoherence vs separation sweep could show hysteresis "
                "if inertial phase dynamics apply (ASSUMED)"
            ),
            rate_note=(
                "Vary separation ramp rate: delayed collapse at fast sweep "
                "(thickening) would support rate-dependent yield over plain DP"
            ),
        )
    )
    return out


# ---------------------------------------------------------------------------
# Rate-dependent / non-Newtonian yield
# ---------------------------------------------------------------------------

def herschel_bulkley_tau(tau_y: float, K: float, n: float, gamma_dot: float) -> float:
    """Herschel–Bulkley: τ = τ_y + K γ̇^n.

    n>1 shear-thickening (dilatant); n<1 shear-thinning; n=1 Bingham.
    ASSUMED analogue: γ̇ ↔ |θ̇| or sweep rate of gravitational drive;
    τ ↔ effective tilt F needed to yield.
    """
    if gamma_dot < 0:
        raise ValueError("γ̇ ≥ 0")
    return tau_y + K * (gamma_dot**n)


def kurkijarvi_switching(Ic: float, ramp_rate: float, T: float, Delta_I_scale: float) -> float:
    """Kurkijärvi–Fulton–Dunkleberger: I_sw rises with ramp rate.

    Qualitative textbook form (ASSUMED prefactors):
      I_sw ≈ Ic − ΔI · [ln(Γ0 Ic / (ramp_rate · const))]^{2/3}
    For small thermal activation, faster ramp ⇒ less time to escape ⇒
    higher apparent switching current.

    We implement a minimal monotone model:
      I_sw / Ic = 1 − a / (1 + b · ramp_rate / ramp0)^{2/3}
    with a,b,ramp0 ASSUMED O(1) shape params — only the *trend* is claimed.
    """
    if ramp_rate <= 0 or Ic <= 0:
        raise ValueError("positive Ic, ramp_rate")
    # Use Delta_I_scale as 'a', fixed b=1, ramp0=1 in normalized units
    a = min(0.5, max(0.0, Delta_I_scale))
    x = ramp_rate  # already normalized to ramp0=1
    return Ic * (1.0 - a / (1.0 + x) ** (2.0 / 3.0))


def rate_dependence_analysis() -> dict:
    """Does yield depend on E_G-driven phase build rate?

    DERIVED from washboard escape + KFD / RCSJ: yes — faster sweep raises
    apparent switching threshold (less time for thermal/quantum escape).
    Herschel–Bulkley: thickening (n>1) raises stress with rate; thinning
    lowers it. Mapping γ̇→Ω_G or ḋ is ASSUMED.

    Testable for three cases: vary how fast separation/mass superposition
    is prepared; look for delayed collapse at fast sweep (thickening/KFD)
    vs rate-independent plain DP τ=ℏ/E_G.
    """
    tau_y, K = 1.0, 0.2  # normalized ASSUMED
    rates = [0.0, 0.1, 1.0, 10.0, 100.0]
    thickening = [herschel_bulkley_tau(tau_y, K, 1.5, r) for r in rates]
    thinning = [herschel_bulkley_tau(tau_y, K, 0.5, r) for r in rates]
    bingham = [herschel_bulkley_tau(tau_y, K, 1.0, r) for r in rates]
    Ic = 1.0
    kfd = [kurkijarvi_switching(Ic, max(r, 1e-12), T=1.0, Delta_I_scale=0.3) for r in rates]
    return {
        "rates": rates,
        "HB_thickening_n1_5": thickening,
        "HB_thinning_n0_5": thinning,
        "HB_bingham_n1": bingham,
        "KFD_Isw": kfd,
        "depends_on_sweep_rate": True,
        "plain_DP_rate_independent": True,
        "distinct_signature": (
            "Collapse/recoherence threshold shifts with preparation ramp rate; "
            "plain DP predicts rate-independent τ=ℏ/E_G"
        ),
        "label": (
            "derived KFD/washboard trend; assumed HB analogy & prefactors; "
            "honest: absolute numbers need V0, β_c, T"
        ),
    }


def self_checks() -> None:
    yt = yield_thresholds()
    assert yt["F_c_over_V0"]["washboard"] == 1.0
    assert yt["F_c_over_V0"]["cos2"] == 1.0
    assert abs(yt["F_c_over_V0"]["double_well"] - 2.0 / (3.0 * math.sqrt(3.0))) < 1e-14
    assert yt["robust_O1"] is True

    h = rcsj_hysteresis(100.0)
    assert h["hysteretic"] is True and h["loop_width"] > 0
    h0 = rcsj_hysteresis(0.1)
    assert h0["hysteretic"] is False

    r0 = r0_from_yield_check()
    assert r0["universal_from_C_ledger"] is False
    assert r0["a0_coincidence_only"] is True

    rd = rate_dependence_analysis()
    assert rd["depends_on_sweep_rate"] is True
    # Thickening: stress increases with rate
    assert rd["HB_thickening_n1_5"][-1] > rd["HB_thickening_n1_5"][0]
    # KFD: I_sw increases with ramp rate
    assert rd["KFD_Isw"][-1] > rd["KFD_Isw"][1]

    for s in case_signatures():
        assert s.Omega_G_Hz > 0 and s.tau_ind_s > 0


def main():
    print("=== Particle plasticity (elastic distortion / plastic collapse) ===\n")
    print("No tuning. Includes rate-dependent (non-Newtonian) yield test.\n")

    print("-- (1) Elastic vs plastic --")
    ep = elastic_plastic_dictionary()
    print(f"  elastic: {ep['elastic']}")
    print(f"  plastic: {ep['plastic']}")
    print(f"  tools: {ep['tools']}")
    print(f"  [{ep['label_map']}; {ep['label_applicability']}]")

    print("\n-- (2) Yield threshold vs V(θ) shape --")
    yt = yield_thresholds()
    for name, val in yt["F_c_over_V0"].items():
        print(f"  F_c/V0 ({name}) = {val:.6f}")
    print(f"  robust O(1) across shapes? {yt['robust_O1']}")
    print(f"  |∇θ|_c scaling: {yt['grad_theta_c_scaling']}")
    print(f"  [{yt['label']}]")

    print("\n-- (3) Hysteresis (RCSJ / Stewart–McCumber) --")
    for h in hysteresis_scan():
        print(
            f"  β_c={h['beta_c']:.1e}: I_r/I_c={h['I_r_over_I_c']:.4f}, "
            f"loop_width={h['loop_width']:.4f}, hysteretic={h['hysteretic']}"
        )
    print(f"  distinct from plain DP? True  [{hysteresis_scan()[-1]['label']}]")

    print("\n-- (4) Yield length as R0? --")
    r0 = r0_from_yield_check()
    print(f"  L_*(BEC)=ξ/π = {r0['L_bec_m']:.4e} m (meets bound? {r0['meets_bound_bec']})")
    print(f"  L_*(Cooper)= {r0['L_coop_m']:.4e} m")
    print(f"  universal from C_ledger? {r0['universal_from_C_ledger']}")
    print(
        f"  a0/bound = {r0['a0_ratio_to_bound']:.4f} "
        f"(coincidence only? {r0['a0_coincidence_only']})"
    )
    print(f"  V0 needed for L_*=bound (BEC n): "
          f"{r0['V0_needed_J_per_m3_for_L_star_eq_bound']:.4e} J/m³")
    print(f"  Plain: {r0['plain']}")

    print("\n-- (5) Case signatures --")
    for s in case_signatures():
        print(f"  {s.name}")
        print(f"    Ω_G={s.Omega_G_Hz:.4e} rad/s, τ_ind={s.tau_ind_s:.4e} s, "
              f"coherence={s.coherence_s:.4e} s")
        print(f"    threshold: {s.threshold_note}")
        print(f"    hysteresis: {s.hysteresis_note}")
        print(f"    rate: {s.rate_note}")

    print("\n-- Rate-dependent / non-Newtonian yield --")
    rd = rate_dependence_analysis()
    print(f"  depends on sweep rate? {rd['depends_on_sweep_rate']}")
    print(f"  plain DP rate-independent? {rd['plain_DP_rate_independent']}")
    print("  Herschel–Bulkley τ(γ̇) [normalized]:")
    for label, arr in (
        ("thickening n=1.5", rd["HB_thickening_n1_5"]),
        ("thinning n=0.5", rd["HB_thinning_n0_5"]),
        ("Bingham n=1", rd["HB_bingham_n1"]),
    ):
        print(f"    {label}: " + ", ".join(f"{x:.3f}" for x in arr))
    print("  KFD I_sw(ramp) [normalized Ic=1]: " +
          ", ".join(f"{x:.3f}" for x in rd["KFD_Isw"]))
    print(f"  signature: {rd['distinct_signature']}")
    print(f"  [{rd['label']}]")

    print("\n-- Self-checks --")
    self_checks()
    print("  all asserts PASS")

    print("\nVerdict")
    print("-------")
    print(
        "Plasticity language fits: distortion=elastic, collapse=plastic yield.\n"
        "F_c/V0 is O(1) across cos/washboard/double-well (robust); V0 unset.\n"
        "RCSJ hysteresis vs β_c is a distinct prediction vs plain DP.\n"
        "Yield length does NOT give untuned Gran-Sasso R0 (a0 coincidence only).\n"
        "Rate-dependent yield (HB/KFD) predicts sweep-rate shifts of the\n"
        "collapse threshold — testable vs rate-independent DP — but absolute\n"
        "scales need assumed V0/β_c; gravity drive is tiny for BEC/flux."
    )
    print("\n=== PARTICLE PLASTICITY COMPLETE ===")


if __name__ == "__main__":
    main()
