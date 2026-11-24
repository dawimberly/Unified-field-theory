#!/usr/bin/env python3
"""Stage 14 — how to get and keep T_int of the Stage 12 cheapest
diamond at or below 4 K (ideally 1 K).

Same blackbody/Debye model as Stage 13, with a more careful
Kirchhoff radiative power. No parameters are fitted.
CITED / DERIVED / ASSUMED are tagged in comments and in the JSON.
"""
from __future__ import annotations
import json, os, sys
import numpy as np
from numpy import pi, sqrt, exp

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from stage12_lib import K_B, HBAR, C, EPS_DIAMOND, RHO_DIAMOND
from stage13_run import (
    CHEAP, EPS_P, THETA_D, M_C, IMCHI_PREFACTOR,
    im_chi_from_tand, heat_capacity,
)

# Physical extras
SIGMA_SB = 5.670374419e-8  # W/m^2/K^4, CITED CODATA
ZETA5 = 1.0369277551433699  # ζ(5)
H_PLANCK = 2 * pi * HBAR
M_HE = 4.002602 * 1.66053906660e-27  # kg
M_N2 = 28.014 * 1.66053906660e-27

M = CHEAP['m']
R = CHEAP['R']
V = (4.0 / 3.0) * pi * R ** 3
T_RUN = CHEAP['t_tot']  # 3 s
N_CARBON = M / M_C


def U_thermal(T, m=M):
    """Internal phonon energy from 0 to T in the Debye T^3 law (DERIVED).
    C = a T^3, U = a T^4 / 4 = C(T) T / 4.
    """
    if T <= 0:
        return 0.0
    return heat_capacity(T, m) * T / 4.0


def T_after_heat(T0, Q, m=M):
    """Final temperature after adding energy Q, Debye U∝T^4 (DERIVED)."""
    a = heat_capacity(1.0, m)  # = C/T^3
    u = U_thermal(T0, m) + Q
    if u <= 0:
        return 0.0
    return float((4.0 * u / a) ** 0.25)


def Q_between(T1, T2, m=M):
    """Energy needed to go from T1 to T2 (DERIVED)."""
    return U_thermal(T2, m) - U_thermal(T1, m)


# ---------------------------------------------------------------------------
# Radiative power, Rayleigh + Kirchhoff, frequency-independent Im χ
# P = (288 ζ(5)/π) ħ R^3 Imχ (kT/ħ)^5 / c^3     (DERIVED, Stage 14)
# ---------------------------------------------------------------------------
def P_rad(T, im_chi, R=R):
    """Emitted (or absorbed from a bath at T) radiative power, Watts."""
    if T <= 0 or im_chi <= 0:
        return 0.0
    pref = (288.0 * ZETA5 / pi) * HBAR * (R ** 3) * im_chi / (C ** 3)
    return float(pref * (K_B * T / HBAR) ** 5)


def radiative_tau_linear(T, im_chi, m=M, R=R):
    """Local linearised radiative time constant C / (dP/dT) (DERIVED).
    P∝T^5 so dP/dT = 5 P/T, τ = C T / (5 P).
    """
    P = P_rad(T, im_chi, R)
    if P <= 0:
        return np.inf
    return float(heat_capacity(T, m) * T / (5.0 * P))


def cooldown_time(T0, T_final, T_env, im_chi, m=M, R=R, max_t=1e12):
    """Time to fall from T0 to T_final with T_env bath (DERIVED ODE).
    C(T) dT/dt = -P(T) + P(T_env). Returns inf if T_final is unreachable.
    """
    if T0 <= T_final:
        return 0.0
    Penv = P_rad(T_env, im_chi, R)
    # If T_final is at or below the sink, require T_final > T_env
    if T_final <= T_env + 1e-9:
        return float('inf')

    # integrate dT / (Penv - P(T)) * C(T)  from T0 down to T_final
    # log-spaced temperatures
    Ts = np.geomspace(T0, T_final, 800)
    # trapezoid on τ = ∫ C dT / (Penv - P)
    Cs = np.array([heat_capacity(T, m) for T in Ts])
    dQdt = np.array([Penv - P_rad(T, im_chi, R) for T in Ts])  # negative when T>T_env
    # dt = C dT / (Penv - P); Penv-P < 0 for cooling, dT<0, dt>0
    if np.any(dQdt[:-1] >= 0) and T_final > T_env:
        # shouldn't happen for T > T_env
        pass
    integrand = Cs / dQdt  # both: C>0, dQdt<0, dT will be negative
    dT = np.diff(Ts)  # negative
    dt = integrand[:-1] * dT  # dT<0 and (Penv-P)<0 => dt>0
    if not np.all(np.isfinite(dt)):
        return float('inf')
    total = float(np.sum(dt))
    return total if total > 0 else float('inf')


# ---------------------------------------------------------------------------
# Mid-IR / 300 K leak using a measured 10.6 μm absorption (CITED)
# α ≈ 0.04 cm^{-1} = 4 m^{-1} at 10.6 μm for type IIa (Thomas et al. / typical)
# σ_abs ≈ V * α  in the optically thin geometric-path sense, or Rayleigh
# At 10 μm, kR ≈ 0.46, borderline Rayleigh. Use P ≈ 4 * σ_abs * σ_SB T^4
# with σ_abs = α * V  (mean-chord, optically thin) as a transparent-window estimate.
# ---------------------------------------------------------------------------
ALPHA_106UM_PER_M = 4.0  # 0.04 cm^{-1}, CITED mid-range of 0.033–0.07 cm^{-1}


def P_abs_300K_leak(view_factor, alpha=ALPHA_106UM_PER_M, R=R):
    """Power absorbed from a 300 K blackbody occupying `view_factor` of 4π (DERIVED)."""
    sigma_abs = alpha * V  # m^2, thin-limit
    flux = SIGMA_SB * 300.0 ** 4  # W/m^2, one-sided
    # isotropic bath: P = 4 π R^2 * (flux/4) * absorptivity, absorptivity ≈ σ_abs/(π R^2)
    # => P = σ_abs * flux. For a 4π bath it's σ_abs * 4 * (σT^4) if σ_abs is the
    # projected absorption cross section. In the thin-sphere limit the mean
    # absorbed power from a 4π blackbody is c u * σ_abs = 4 σT^4 * σ_abs
    # if σ_abs is the usual illumination cross section.
    # Our σ_abs = α V = α (4/3)π R^3, while projected area * (α * 4R/3) = π R^2 * α * 4R/3 = α V.
    # That's the illumination cross section. 4π bath => P = 4 σT^4 σ_abs.
    return float(view_factor * 4.0 * flux * sigma_abs)


# ---------------------------------------------------------------------------
# Optical absorption
# ---------------------------------------------------------------------------
# Lattice α at 1064 nm: 0.03 cm^{-1} = 3 m^{-1} on a CVD single crystal
# (CITED, simplified calorimetry, Element Six; type IIa often quoted <0.05 cm^{-1}).
# Best windows may be ~0.001 cm^{-1}; we also report that as optimistic.
ALPHA_1064_TYPICAL = 3.0    # m^{-1}   (= 0.03 cm^{-1})
ALPHA_1064_BEST = 0.1       # m^{-1}   (= 0.001 cm^{-1}), ASSUMED optimistic floor
# 532 nm: calorimetry on the same class of CVD gives higher α (the BNTU paper
# rose from 0.03 cm^{-1} at 1064 to 0.7 cm^{-1} at 364 nm). We take 0.1 cm^{-1}
# as a mid-visible order-of-magnitude (ASSUMED interpolation, flagged).
ALPHA_532_LATTICE = 10.0    # m^{-1}   (= 0.1 cm^{-1}), ASSUMED order of magnitude

# NV− absorption cross section at 532 nm: 3.1(0.8)×10^{-17} cm^2
# (CITED Wee / TPA paper, Phys. status or Opt. Express type measurements;
#  UAB compilation  σ_ZPL ≈ 4×10^{-17} cm^2). Use 3.1e-17 cm^2.
SIGMA_NV_532_M2 = 3.1e-17 * 1e-4  # cm^2 → m^2

# Mean fluorescence wavelength of NV ~ 700 nm; QE of bulk NV is high, ~0.7–1.
# Net heat fraction per absorbed 532 nm photon ≈ 1 - QE * (532/λ_fl).
# QE=0.8, λ_fl=700 nm => 1 - 0.8*532/700 = 0.39  (ASSUMED mid).
NV_QE = 0.8
NV_LAMBDA_FL_NM = 700.0
NV_HEAT_FRACTION = 1.0 - NV_QE * (532.0 / NV_LAMBDA_FL_NM)

E_PHOTON_532 = H_PLANCK * C / 532e-9
E_PHOTON_1064 = H_PLANCK * C / 1064e-9


def P_abs_lattice(I, alpha, R=R):
    """Optically thin lattice absorption P = α V I (DERIVED). I in W/m^2."""
    return float(alpha * V * I)


def intensity_from_power_waist(P, w):
    """Peak intensity of a Gaussian beam, I = 2P/(π w^2)."""
    return 2.0 * P / (pi * w ** 2)


# ---------------------------------------------------------------------------
# Free-molecular (Knudsen) gas cooling of a sphere (DERIVED kinetic theory)
# Collision rate Γ = n π R^2 <v>, <v> = sqrt(8 kT / (π m))
# Energy per accommodated monatomic collision: 2 k (T_p - T_g)
# P = Γ α 2 k (T_p - T_g) = 2 π R^2 α (p / T_g) <v> (T_p - T_g)
# Accommodation α_E is ASSUMED (0.5 default); He on cold surfaces can be smaller.
# ---------------------------------------------------------------------------
def v_mean(T, m_gas):
    return sqrt(8.0 * K_B * T / (pi * m_gas))


def gas_conductance(P, T_g, R=R, m_gas=M_HE, alpha=0.5):
    """G such that P_heat = G (T_p - T_g), in W/K (DERIVED)."""
    if P <= 0 or T_g <= 0:
        return 0.0
    return float(2.0 * pi * R ** 2 * alpha * (P / T_g) * v_mean(T_g, m_gas))


def gas_tau(P, T_g, T_p=None, m_gas=M_HE, alpha=0.5, m=M):
    """Linearised cooling time C(T_p)/G. T_p defaults to T_g."""
    if T_p is None:
        T_p = T_g
    G = gas_conductance(P, T_g, R, m_gas, alpha)
    if G <= 0:
        return np.inf
    return float(heat_capacity(T_p, m) / G)


# ---------------------------------------------------------------------------
# Cases for Im χ, same as Stage 13
# ---------------------------------------------------------------------------
CASES = dict(
    window=dict(label='window-grade CVD', tan_delta=1e-5,
                im_chi=float(im_chi_from_tand(1e-5))),
    millar22=dict(label='Millar 22 GHz cryo / cleaner ND', tan_delta=1e-3,
                  im_chi=float(im_chi_from_tand(1e-3))),
    dirty=dict(label='dirty nanodiamond / sp2 shell', tan_delta=1e-2,
               im_chi=float(im_chi_from_tand(1e-2))),
    stage12=dict(label='Stage 12 assumed Im χ=0.01', tan_delta=float(0.01 / IMCHI_PREFACTOR),
                 im_chi=0.01),
)


def main():
    os.chdir(ROOT)

    # --- heat capacity table
    T_list = [0.5, 1.0, 2.0, 4.0, 10.0, 20.0, 77.0, 300.0]
    cv_table = []
    for T in T_list:
        cv_table.append(dict(
            T=T, C_J_per_K=float(heat_capacity(T)),
            U_from_0_J=float(U_thermal(T)),
        ))

    Q_1_to_4 = float(Q_between(1.0, 4.0))
    Q_0_to_4 = float(U_thermal(4.0))
    Q_1_to_1p1 = float(Q_between(1.0, 1.1))
    # energy → T from 1 K and from 4 K
    energies = np.geomspace(1e-21, 1e-9, 49)  # 0.001 aJ to 1 nJ
    T_from_1K = [T_after_heat(1.0, Q) for Q in energies]
    T_from_4K = [T_after_heat(4.0, Q) for Q in energies]
    T_from_0K = [T_after_heat(0.0, Q) for Q in energies]

    # photons
    n532_1_to_4_full = Q_1_to_4 / E_PHOTON_532
    n532_1_to_4_net = Q_1_to_4 / (E_PHOTON_532 * NV_HEAT_FRACTION)
    n1064_1_to_4_full = Q_1_to_4 / E_PHOTON_1064

    # impulsive budget during the 3 s run (radiation negligible for clean)
    P_budget_impulse_1to4 = Q_1_to_4 / T_RUN  # W, if starting at 1 K, no cooling
    P_budget_impulse_4to5 = Q_between(4.0, 5.0) / T_RUN

    # --- radiative powers and times
    rad = {}
    cool = {}
    Pss = {}  # steady-state extra power that holds T_int at T, given T_env
    for key, c in CASES.items():
        im = c['im_chi']
        rad[key] = dict(
            im_chi=im, tan_delta=c['tan_delta'], label=c['label'],
            P_em_1K=P_rad(1.0, im), P_em_4K=P_rad(4.0, im),
            P_em_20K=P_rad(20.0, im), P_em_100K=P_rad(100.0, im),
            P_em_300K=P_rad(300.0, im),
            tau_lin_1K=radiative_tau_linear(1.0, im),
            tau_lin_4K=radiative_tau_linear(4.0, im),
            tau_lin_100K=radiative_tau_linear(100.0, im),
            tau_lin_300K=radiative_tau_linear(300.0, im),
        )
        cool[key] = {}
        for T0 in (300, 100, 77, 20, 10):
            for Tenv, Tfin in ((4.0, 5.0), (4.0, 4.2), (1.0, 4.0), (1.0, 5.0)):
                tag = f'from{int(T0)}K_env{Tenv:g}_to{Tfin:g}'
                cool[key][tag] = cooldown_time(T0, Tfin, Tenv, im)
        # extra power to sit at T_int with T_env
        Pss[key] = dict(
            P_extra_Tint4_Tenv1=P_rad(4.0, im) - P_rad(1.0, im),
            P_extra_Tint4_Tenv4=0.0,  # no margin
            P_extra_Tint5_Tenv4=P_rad(5.0, im) - P_rad(4.0, im),
        )

    # cooldown curves (T vs t) for the figure
    t_grid = np.geomspace(1e-3, 1e8, 400)
    curves = {}
    for key, c in CASES.items():
        im = c['im_chi']
        Tenv = 4.0
        # precompute a T(t) by integrating from 300 K
        Ts = np.geomspace(300.0, 4.05, 1000)
        ts = [0.0]
        for i in range(1, len(Ts)):
            dt = cooldown_time(Ts[i - 1], Ts[i], Tenv, im)
            ts.append(ts[-1] + (dt if np.isfinite(dt) else 1e99))
        curves[key] = dict(T=Ts.tolist(), t_s=ts)

    # --- gas cooling vs pressure
    P_gas = np.geomspace(1e-16, 1e5, 85)  # Pa
    gas = []
    for p in P_gas:
        gas.append(dict(
            P_Pa=float(p),
            tau_He_4K_s=gas_tau(p, 4.0, 4.0, M_HE, 0.5),
            tau_He_4K_alpha01=gas_tau(p, 4.0, 4.0, M_HE, 0.1),
            tau_He_300K_s=gas_tau(p, 300.0, 300.0, M_HE, 0.5),
            tau_N2_4K_s=gas_tau(p, 4.0, 4.0, M_N2, 0.5),
        ))
    # pressures where tau = 3 s and 1 s at 4 K He, α=0.5
    def p_for_tau(target, T=4.0, m_gas=M_HE, alpha=0.5):
        # G = C/tau = 2π R^2 α (p/T) <v> => p = (C/tau) * T / (2π R^2 α <v>)
        G = heat_capacity(T) / target
        return float(G * T / (2.0 * pi * R ** 2 * alpha * v_mean(T, m_gas)))

    p_tau_3s = p_for_tau(3.0)
    p_tau_1s = p_for_tau(1.0)
    p_tau_0p1s = p_for_tau(0.1)

    # --- optical absorption examples
    w_trap = 1.0e-6  # 1 μm waist, typical tweezer, ASSUMED for estimates
    examples = []
    for name, P, w, alpha, lam in [
        ('tweezer_100mW_1064_typical_alpha', 0.10, w_trap, ALPHA_1064_TYPICAL, 1064e-9),
        ('tweezer_100mW_1064_best_alpha', 0.10, w_trap, ALPHA_1064_BEST, 1064e-9),
        ('tweezer_1mW_1064_typical_alpha', 1e-3, w_trap, ALPHA_1064_TYPICAL, 1064e-9),
        ('alignment_1uW_1064_typical', 1e-6, w_trap, ALPHA_1064_TYPICAL, 1064e-9),
        ('NV_read_100uW_532_lattice_only', 1e-4, w_trap, ALPHA_532_LATTICE, 532e-9),
    ]:
        I = intensity_from_power_waist(P, w)
        Pabs = P_abs_lattice(I, alpha)
        examples.append(dict(
            name=name, P_beam_W=P, waist_m=w, I_W_m2=I, alpha_per_m=alpha,
            P_abs_W=Pabs, T_ss_if_Tenv4_impulse_in_3s=T_after_heat(4.0, Pabs * T_RUN),
            ratio_to_1to4_budget=Pabs / P_budget_impulse_1to4 if P_budget_impulse_1to4 else None,
        ))

    # single-NV 532 nm absorption at given intensity
    I_100uW = intensity_from_power_waist(1e-4, w_trap)
    P_nv_one = SIGMA_NV_532_M2 * I_100uW  # W, absorbed (before QE)
    P_nv_one_heat = P_nv_one * NV_HEAT_FRACTION
    # pulse lengths that deposit Q_1_to_4 as net heat, one NV, that intensity
    t_pulse_budget = Q_1_to_4 / P_nv_one_heat if P_nv_one_heat else np.inf

    # 300 K leak
    P300_full = P_abs_300K_leak(1.0)
    # view factor that deposits Q_1_to_4 in 3 s, and that equals P_em(4K) window-grade
    vf_impulse = (Q_1_to_4 / T_RUN) / P300_full
    vf_ss_window = rad['window']['P_em_4K'] / P300_full

    # MW / eddy qualitative numbers
    # SG: ΔB ~ ∇B Δx ≈ 3.74e4 * 3.85e-6 ≈ 0.14 T in 0.25 s (Stage 12 cheapest)
    dBdt_sg = 3.74e4 * 3.852e-6 / 0.25
    eddy = dict(
        dBdt_SG_T_per_s=float(dBdt_sg),
        note='For insulating type IIa diamond, eddy current dissipation from this slow Stern–Gerlach ramp is negligible (DERIVED: even a 2 nm graphite shell at σ~10^5 S/m gives ~10^{-29} W). A graphitic shell plus 2.87 GHz MW drive is the only plausible eddy risk and is not quantified here without a field geometry. Dielectric MW heating at tanδ=10^{-5} with a modest E also sits below the aW budget for short π-pulses; tanδ=10^{-2} and a high-E near-field antenna would need a dedicated estimate (ASSUMED not binding if pulses are short and the diamond is clean).',
    )

    literature = dict(
        hoang_li_2016='Hoang, Ahn, Bang, Li, Nat. Commun. 7, 12250 (2016): optically levitated NV nanodiamond, ESR thermometry. T from ~300 K at 1 atm to >450 K at 31 Torr under the trapping laser. Heating plus 1550/1064 nm ionization of NV−.',
        rahman_frangeskou_2016='Rahman, Frangeskou, et al., Sci. Rep. 6, 21633 (2016): optically levitated NV nanodiamonds burn in air; in N2 they graphitize below ≈10 mbar and vanish by ≈2 mbar. Internal T ≈ 750–800 K at 300–380 mW of 1064 nm while gas cooling still works. Not melting.',
        frangeskou_2018='Frangeskou, Rahman, Williams, Barker, Morley, New J. Phys. 20, 043016 (2018), arXiv:1608.04724: milled high-purity CVD nanodiamonds showed no measurable 1064 nm heating up to ~300 mW at 3–4 mbar; commercial HPHT NV powder has been reported >800 K. Nitrogen and surface graphitic carbon, not the lattice, dominate NIR absorption.',
        delord_hetet_2017='Delord, Nicolas, Bodini, Hétet, Appl. Phys. Lett. 111, 013101? / 111, (2017), arXiv:1701.06407: micron diamonds in a Paul trap (no IR trap laser). Green-laser heating alone: T up to 470 K at 700 μW, 1 mbar; 390 K at 40 μW, 0.2 mbar. Attributed to nitrogen in commercial HPHT powder. Particles lost from the trap below ~10^{-2} mbar under the green beam.',
        pettit_2017='Pettit et al., J. Opt. Soc. Am. B 34, C31 (2017): coherent control of a single NV in an optically levitated nanodiamond. After an initial low-vacuum exposure (cleaning), T stayed near room temperature at 2.5 kPa across the powers they tested. Still an optical trap, still gas-cooled, not XHV.',
        rahman_barker_2017='Rahman & Barker, Nat. Photonics 11, 634 (2017): anti-Stokes refrigeration of levitated Yb3+:YLF nanocrystals to 130 K. Demonstrated, but not diamond, and not below ~130 K.',
        kern_2017='Kern, Jeske, Lau, Greentree, Jelezko, Twamley, Phys. Rev. B 95, 235306 (2017), arXiv:1701.08505: theory of optical cryocooling of NV- or SiV-doped microdiamonds. Predict |ΔT|>10 K below room temperature at <100 mW in vacuum for NV (high QE); SiV less promising because of lower quantum efficiency. Not demonstrated, and not a 4 K protocol.',
        hsu_2016='Hsu, Ji, Lewandowski, D’Urso, Sci. Rep. 6, 30125 (2016): magneto-gravitational trap for diamagnetic diamond nanocrystal clusters, stable from atmosphere to high vacuum, no intrinsic trap heating. Feedback cooling of COM (not T_int) to <1 K in one degree of freedom. Loaded by spray at atmospheric pressure.',
        elahi_2024='Elahi, Schut, Dana, Grinin, Bose, Mazumdar, Geraci, arXiv:2411.02325 (2024/25): diamagnetic microchip traps designed for the QGEM protocol (mm-scale Z-trap for a weak Stern–Gerlach axis, 100 Hz / 500–800 Hz). Static/quasi-static currents, intended to avoid optical heating. A design paper, not a demonstration of a cold 10^{-15} Pa particle.',
        alpha_1064='Calorimetry on single-crystal CVD (Element Six): α ≈ 0.03 cm^{-1} at 1064 nm, rising toward the UV (0.7 cm^{-1} at 364 nm). Type IIa compilations often quote α < 0.05 cm^{-1} at 1064 nm. Mid-IR type IIa window: α ≈ 0.033–0.07 cm^{-1} at 10.6 μm (Thomas et al. / Crystran compilations).',
        sigma_NV='NV− absorption cross section at 532 nm: 3.1(0.8)×10^{-17} cm^2 (two-photon-referenced one-photon value). ZPL values around 4×10^{-17} cm^2 also appear. We use 3.1×10^{-17} cm^2.',
    )

    protocol = dict(
        feasible_outline=[
            'Material: high-purity type IIa / optical-grade CVD, oxygen-terminated (no graphitic shell), as few NV centres as the witness allows (ideally one). This is the only way Im χ and 532/1064 absorption stay in the window-grade / single-NV regime.',
            'Trap: all-magnetic (magneto-gravitational, Paul, or a QGEM microchip trap). No 1064 nm tweezer during pump-down or the 3 s run. Hsu 2016 demonstrated laser-free diamagnetic trapping of diamond from atmosphere to high vacuum; Elahi 2024 is the QGEM-shaped design. Not yet demonstrated at 10^{-15} Pa or at 4 K with a 5.6×10^{-15} kg particle.',
            'Environment: a 4 K (ideally 1 K) cryostat with tight radiation shields. A 300 K view factor of 10^{-7}–10^{-4} (depending on mid-IR α) already exceeds the 3 s energy budget. No room-temperature windows in line of sight.',
            'Load and thermalise with cold helium buffer gas at 4 K (or 1 K) while the magnetic trap is on. Gas cooling is fast (milliseconds to seconds) down to ~10^{-4}–10^{-6} Pa at 4 K, and is demonstrated physics, though not in this exact apparatus. Then pump toward 10^{-15} Pa with lasers off. The particle stays at T_int ≈ T_env because it is already there, not because it radiates.',
            'Initialise the NV spin while gas cooling still works (or after contact/loading, before the final pump). A 100 μW, 1 μm-waist 532 nm beam on a single NV deposits the entire 1 K→4 K budget in tens of microseconds of net heat at vacuum. Do not initialise in XHV and then wait — radiation will not dump the heat in time.',
            'During the 3 s interferometry: no 532, no 1064, only the slow Stern–Gerlach ramps and, if needed, short microwave π-pulses. Read the spin after recombination. Heating during a short final readout does not undo the already-written entanglement, but the pulse must still be short enough to preserve contrast.',
        ],
        demonstrated=[
            'Optical-trap heating and burning/graphitization of NV nanodiamonds in vacuum (Rahman 2016, Hoang/Li 2016).',
            'Green-laser heating of micron HPHT diamonds in a Paul trap with no IR tweezer (Delord/Hétet 2017).',
            'Reduced NIR heating in milled high-purity CVD nanodiamonds, still at a few mbar (Frangeskou 2018).',
            'Diamagnetic / magneto-gravitational trapping of diamond from air to high vacuum, no intrinsic heating (Hsu 2016).',
            'Anti-Stokes refrigeration of levitated Yb:YLF to 130 K (Rahman & Barker 2017) — not diamond, not 4 K.',
            'Buffer-gas cooling of levitated particles in general, and 4 K cryostats in general.',
        ],
        not_demonstrated=[
            'A 5.6×10^{-15} kg, 0.73 μm diamond with one (or few) NV, oxygen-terminated, in a magnetic trap at 10^{-15} Pa and T_int ≤ 4 K.',
            'Pump-down from buffer-gas pressures to 10^{-15} Pa around a trapped microdiamond in a 4 K apparatus.',
            'Stern–Gerlach splitting and Humpty-Dumpty recombination of such a trapped particle at the Stage 12 field precision.',
            'Optical cryocooling of diamond to (or anywhere near) 4 K. Kern et al. 2017 is theory for modest cooling below room temperature.',
            'An emissivity coating that speeds radiative cooling enough to matter at 4 K without destroying the PPT witness (Stage 13: Im χ ≳ 10^{-3} already squeezes the 4 K point).',
        ],
        rejected=[
            'Optical trapping or even a 1 μW alignment beam at 1064 nm during the XHV run: lattice α V I is 10^6–10^{11} times the 3 s energy budget, even for optimistic α.',
            'Relying on radiative cooling after a laser pulse: window-grade τ at 1–4 K is days; dirty-surface τ is minutes, still much longer than 3 s, and “dirty” already squeezes the witness at 4 K.',
            'Raising THz emissivity with a coating to speed cooling: the same Im χ enters the Stage 13 decoherence budget and kills C near 4–6 K.',
            'Laser refrigeration as the primary 4 K path: not demonstrated in diamond, requires a cooling laser (another absorption source), and existing theory aims at tens of kelvin below ambient, not 1–4 K.',
        ],
    )

    out = dict(
        particle=dict(m_kg=M, R_m=R, V_m3=float(V), N_carbon=float(N_CARBON),
                      t_run_s=T_RUN, theta_D_K=THETA_D,
                      C_1K=float(heat_capacity(1.0)), C_4K=float(heat_capacity(4.0))),
        heat_capacity_table=cv_table,
        energy_budget=dict(
            Q_1K_to_4K_J=Q_1_to_4,
            Q_1K_to_4K_aJ=Q_1_to_4 * 1e18,
            Q_0_to_4K_J=Q_0_to_4,
            Q_1K_to_1p1K_J=Q_1_to_1p1,
            Q_4K_to_5K_J=float(Q_between(4.0, 5.0)),
            P_budget_impulse_1to4_W=P_budget_impulse_1to4,
            P_budget_impulse_1to4_aW=P_budget_impulse_1to4 * 1e18,
            n_photons_532_full_abs_1to4=float(n532_1_to_4_full),
            n_photons_532_net_heat_1to4=float(n532_1_to_4_net),
            n_photons_1064_full_abs_1to4=float(n1064_1_to_4_full),
            E_photon_532_J=float(E_PHOTON_532),
            E_photon_1064_J=float(E_PHOTON_1064),
            NV_heat_fraction=NV_HEAT_FRACTION,
            note='Debye U=C T/4. From 1 K to 4 K takes only ~5.5 aJ. That is ~15 fully absorbed 532 nm photons, or ~40 absorbed 532 nm photons after a Stokes-and-QE correction of 0.39. Continuous allowed power over 3 s, with no cooling, is ~2 aW.',
        ),
        T_vs_energy=dict(
            Q_J=energies.tolist(),
            T_from_0K=T_from_0K,
            T_from_1K=T_from_1K,
            T_from_4K=T_from_4K,
        ),
        radiative=rad,
        cooldown=cool,
        cooldown_curves=curves,
        steady_state_power=Pss,
        gas=dict(
            points=gas,
            p_He_4K_tau_3s_Pa=p_tau_3s,
            p_He_4K_tau_1s_Pa=p_tau_1s,
            p_He_4K_tau_0p1s_Pa=p_tau_0p1s,
            alpha_assumed=0.5,
            formula='P = 2 π R^2 α (p/T_g) <v> (T_p-T_g), <v>=sqrt(8kT/πm), monatomic 2k per collision (DERIVED).',
        ),
        optical_examples=examples,
        NV_readout=dict(
            sigma_m2=SIGMA_NV_532_M2,
            I_100uW_1um_waist_W_m2=I_100uW,
            P_abs_one_NV_W=float(P_nv_one),
            P_heat_one_NV_W=float(P_nv_one_heat),
            t_pulse_for_1K_to_4K_s=float(t_pulse_budget),
            t_pulse_for_1K_to_4K_us=float(t_pulse_budget * 1e6),
            note='100 μW, 1 μm waist, one NV, heat fraction 0.39. The 1→4 K budget is spent in a few tens of microseconds of illumination in vacuum. A typical millisecond ODMR pulse is far over budget if it happens in XHV before the hold.',
        ),
        leak_300K=dict(
            alpha_10p6um_per_m=ALPHA_106UM_PER_M,
            P_abs_4pi_300K_W=float(P300_full),
            view_factor_for_3s_1to4_budget=float(vf_impulse),
            view_factor_equal_to_window_Pem_4K=float(vf_ss_window),
            note='A full 300 K bath, even with type IIa mid-IR absorption, deposits far more than 5 aJ in 3 s. Tight 4 K (or 1 K) radiation shields are mandatory; a room-temperature window in line of sight is not a small correction.',
        ),
        eddy_and_MW=eddy,
        literature=literature,
        protocol=protocol,
        verdict_one_liner='The 0.73 μm diamond has a 1 K heat capacity of 8.5×10^{-20} J/K; 5.5 aJ takes it from 1 K to 4 K. No optical beam can be on during the XHV run. The only feasible path is: start already cold in a 4 K (or 1 K) magnetically trapped, well-shielded cryostat, gas-cool, pump with lasers off, initialise before the final pump, read out after the hold. Radiative self-cooling and laser refrigeration will not get you there in time.',
    )
    path = os.path.join(ROOT, 'results', 'stage14.json')
    with open(path, 'w') as f:
        json.dump(out, f, indent=2)
    print('wrote', path)
    print('C(1K)=%.3e J/K  C(4K)=%.3e J/K' % (heat_capacity(1), heat_capacity(4)))
    print('Q 1→4 K = %.3e J = %.2f aJ' % (Q_1_to_4, Q_1_to_4*1e18))
    print('P budget 3s = %.2f aW' % (P_budget_impulse_1to4*1e18))
    print('532 nm photons full/net = %.1f / %.1f' % (n532_1_to_4_full, n532_1_to_4_net))
    print('P_em window 4K = %.3e W; dirty 4K = %.3e W' % (rad['window']['P_em_4K'], rad['dirty']['P_em_4K']))
    print('tau_lin window 4K = %.3e s; dirty 4K = %.3e s' % (rad['window']['tau_lin_4K'], rad['dirty']['tau_lin_4K']))
    print('cool window 300→5K @4K env: %.3e s' % cool['window']['from300K_env4_to5'])
    print('cool dirty  300→5K @4K env: %.3e s' % cool['dirty']['from300K_env4_to5'])
    print('He 4K P for tau=3s: %.3e Pa' % p_tau_3s)
    print('300K vf for budget: %.3e; P300 full=%.3e W' % (vf_impulse, P300_full))
    print('NV pulse budget: %.2f us' % (t_pulse_budget*1e6))
    for e in examples:
        print(' ', e['name'], 'Pabs=%.3e W  ratio=%.3e' % (e['P_abs_W'], e['ratio_to_1to4_budget']))


if __name__ == '__main__':
    main()
