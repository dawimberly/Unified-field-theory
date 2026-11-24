#!/usr/bin/env python3
"""Stage 11: gravity-mediated entanglement (Bose et al. / Marletto–Vedral).

Reproduces the 2017 numbers, maps the (mass, separation, time) requirement,
computes a decoherence budget, and records the interpretation plus the
Oppenheim classical-gravity tradeoff. No fitted parameters.
Every stored claim is tagged DERIVED, ASSUMED or CITED.
"""
from __future__ import annotations
import json, os
from math import pi, sqrt, exp, factorial, sin, cos, log10
import numpy as np
from scipy.special import zeta

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
os.chdir(ROOT)
os.makedirs('results', exist_ok=True)
os.makedirs('figures', exist_ok=True)

# CODATA 2018 / 2022 (CITED)
G = 6.67430e-11          # m^3 kg^-1 s^-2
HBAR = 1.054571817e-34   # J s
HPLANCK = 6.62607015e-34
K_B = 1.380649e-23
C = 299792458.0
U_DA = 1.66053906660e-27  # kg (dalton)
EPS0 = 8.8541878128e-12
L_P = 1.616255e-35       # Planck length
M_P = 2.176434e-8        # Planck mass
M_NUCLEON = 1.67262192369e-27

# Material / geometry defaults (CITED Bose et al. unless noted)
RHO_DIAMOND = 3510.0     # kg/m^3, typical diamond (ASSUMED; Bose used R~1 um with m~1e-14, slightly inconsistent)
EPS_DIAMOND = 5.7        # Bose
CHI_CP = ((EPS_DIAMOND - 1) / (EPS_DIAMOND + 2)) ** 2  # 0.3725

# Residual gas for collisions: N2 at a stated gas temperature (ASSUMED species)
M_N2 = 28.0134 * U_DA


def sphere_radius(m, rho=RHO_DIAMOND):
    return (3.0 * m / (4.0 * pi * rho)) ** (1.0 / 3.0)


def grav_phase(m1, m2, t, r):
    """Newtonian phase φ = V t / ħ with V = -G m1 m2 / r (DERIVED). Attractive => negative."""
    return - (G * m1 * m2 * t) / (HBAR * r)


def grav_phase_mag(m1, m2, t, r):
    return (G * m1 * m2 * t) / (HBAR * r)


def differential_phases(m1, m2, t, d, dx):
    """Bose geometry: centres d apart, each split by Δx along the joining axis.

    Distances: LL=RR=d, LR=d+Δx, RL=d-Δx.
    Returns subtracted phases Δφ_LR = φ_LR-φ_LL, Δφ_RL = φ_RL-φ_LL,
    and the entanglement combination Δφ_LR+Δφ_RL, using V=-Gmm/r.
    """
    phi_ll = grav_phase(m1, m2, t, d)
    phi_lr = grav_phase(m1, m2, t, d + dx)
    phi_rl = grav_phase(m1, m2, t, d - dx)
    dphi_lr = phi_lr - phi_ll
    dphi_rl = phi_rl - phi_ll
    return dict(phi_ll=phi_ll, phi_rr=phi_ll, phi_lr=phi_lr, phi_rl=phi_rl,
                dphi_lr=dphi_lr, dphi_rl=dphi_rl,
                dphi_sum=dphi_lr + dphi_rl,
                dphi_sum_analytic=(G * m1 * m2 * t / HBAR) * 2 * dx**2 / (d * (d**2 - dx**2)) * (-1.0),
                dphi_sum_bose_approx_mag=(G * m1 * m2 * t / HBAR / d) * (dx / d) ** 2)


def concurrence_bmv(dphi_lr, dphi_rl):
    """Pure-state concurrence of (|00>+e^{iα}|01>+e^{iβ}|10>+|11>)/2.
    C = |sin((α+β)/2)|  (DERIVED; entangled iff α+β != 2π n)."""
    return abs(sin(0.5 * (dphi_lr + dphi_rl)))


def pauli_correlators(alpha, beta):
    """⟨σi⊗σj⟩ on (|00>+e^{iα}|01>+e^{iβ}|10>+|11>)/2 (DERIVED)."""
    sxsz = 0.5 * (cos(beta) - cos(alpha))
    sysz = 0.5 * (sin(alpha) + sin(beta))
    sxsx = 0.5 * (1.0 + cos(alpha - beta))
    sysy = 0.5 * (-1.0 + cos(alpha - beta))
    szsz = 0.0
    W_bose = abs(sxsz - sysz)  # written definition W=|⟨σx⊗σz⟩-⟨σy⊗σz⟩|
    W_separable_max = sqrt(2.0)  # max of |⟨(σx-σy)⊗σz⟩| on product states
    return dict(sxsz=sxsz, sysz=sysz, sxsx=sxsx, sysy=sysy, szsz=szsz,
                W_written=W_bose, W_separable_max=W_separable_max,
                entangled_by_written_threshold_Wgt1=bool(W_bose > 1.0),
                written_threshold_valid=False,
                note=('W=|⟨(σx-σy)⊗σz⟩|. Product states reach √2, so W>1 does not '
                      'prove entanglement. Valid pure-state test: concurrence > 0, i.e. '
                      'Δφ_LR+Δφ_RL != 2πn.'))


def casimir_polder_energy(R, r, eps=EPS_DIAMOND):
    """Retarded CP energy of two dielectric spheres (CITED standard SI form).
    U = - [23 ħ c / (4π)] [R^3 (ε-1)/(ε+2)]^2 / r^7
    (equiv. form used by Bose up to typography of 4πε0 factors)."""
    chi = ((eps - 1) / (eps + 2)) ** 2
    return - (23.0 * HBAR * C / (4.0 * pi)) * (R**6) * chi / r**7


def gravity_energy(m1, m2, r):
    return - G * m1 * m2 / r


def cp_over_grav(m, R, r, eps=EPS_DIAMOND):
    ug = gravity_energy(m, m, r)
    uc = casimir_polder_energy(R, r, eps)
    return abs(uc / ug) if ug != 0 else np.inf


def r_min_casimir(m, R, ratio=0.1, eps=EPS_DIAMOND):
    """Closest approach at which |U_CP|/|U_G| = ratio (DERIVED). Scales as r_min ∝ R for fixed density."""
    chi = ((eps - 1) / (eps + 2)) ** 2
    # |U_CP|/|U_G| = 23 ħ c R^6 chi / (4π G m^2 r^6) = ratio
    # r^6 = 23 ħ c R^6 chi / (4π G m^2 ratio)
    num = 23.0 * HBAR * C * (R**6) * chi
    den = 4.0 * pi * G * (m**2) * ratio
    return (num / den) ** (1.0 / 6.0)


def vbar_gas(T, m_gas=M_N2):
    return sqrt(8.0 * K_B * T / (pi * m_gas))


def gamma_coll_resolution(P, T, R, m_gas=M_N2, sigma_factor=1.0):
    """Resolution-limited collisional decoherence (Δx ≫ thermal de Broglie of the gas).
    Each collision fully distinguishes the paths: Γ = n ⟨v⟩ σ, σ = sigma_factor * π R².
    sigma_factor=1 is the geometric cross section (ASSUMED; 2 is sometimes used for diffraction).
    CITED structure: Hornberger, Sipe & Arndt, PRA 70, 053608 (2004); Bose uses this regime.
    """
    n = P / (K_B * T)
    sigma = sigma_factor * pi * R**2
    return n * vbar_gas(T, m_gas) * sigma


def n_bb_photons(T):
    """Blackbody photon number density, 2 ζ(3) (kT/ħc)^3 / π² (CITED)."""
    return 2.0 * zeta(3) * (K_B * T / (HBAR * C))**3 / pi**2


def rayleigh_sigma(R, k, eps=EPS_DIAMOND):
    alpha = R**3 * (eps - 1) / (eps + 2)
    return (8.0 * pi / 3.0) * (k**4) * alpha**2


def gamma_bb_scattering_rate(T, R, eps=EPS_DIAMOND):
    """Order-of-magnitude thermal-photon scattering rate: n_γ σ_R(k_typ) c,
    with k_typ = 2.821 kT/ħc (Wien-ish for energy). Good to a factor of a few (ASSUMED k_typ).
    A proper frequency integral changes the O(1) prefactor, not the T^7 R^6 scaling in the Rayleigh regime."""
    if T <= 0:
        return 0.0
    k_typ = 2.821 * K_B * T / (HBAR * C)
    return n_bb_photons(T) * rayleigh_sigma(R, k_typ, eps) * C


def gamma_bb_scatt(T, R, dx, eps=EPS_DIAMOND):
    """Photon-scattering decoherence from the Rayleigh scattering rate (DERIVED interpolation).

    Γ_scatt ≈ n_γ(T) σ_R(k_typ) c is a trustworthy order-of-magnitude rate (factor of a few
    from the choice of k_typ and from frequency integration). Decoherence cannot exceed the
    scattering rate. In the long-wavelength regime Δx ≪ λ_th = ħc/kT, each photon only
    distinguishes the paths by a factor ~ (Δx/λ_th)^2:
        Γ = Γ_scatt * (Δx/λ_th)^2 / (1 + (Δx/λ_th)^2).
    This avoids published-prefactor ambiguity in the localization parameter Λ; Λ-formulas
    in the literature differ by large O(1) or 4π factors and are easy to mis-normalize.
    Emission from a hot internal temperature is the same rate with T -> T_int (ASSUMED
    the object is a black/grey body at T_int; diamond is actually a weak IR emitter, so
    this is conservative if T_int ≈ T_env).
    """
    if T <= 0:
        return 0.0
    lam = HBAR * C / (K_B * T)
    rate = gamma_bb_scattering_rate(T, R, eps)
    x = dx / lam
    return rate * (x * x) / (1.0 + x * x)


def mean_thermal_wavelength_gas(T, m_gas=M_N2):
    return HPLANCK / sqrt(2.0 * pi * m_gas * K_B * T)


# ---------------------------------------------------------------------------
# 1. Reproduce Bose et al. PRL 119, 240401 (2017), arXiv:1707.06050
# ---------------------------------------------------------------------------
bose = dict(
    m=1e-14, d=450e-6, dx=250e-6, tau=2.5, tau_acc=0.5,
    R=1e-6, T_internal=0.15, P=1e-15,
    closest=200e-6,  # d - Δx
    refs='Bose et al., PRL 119, 240401 (2017), arXiv:1707.06050; '
         'Marletto & Vedral, PRL 119, 240402 (2017), arXiv:1707.06036',
    quoted=dict(dphi_LR=-0.2, dphi_RL=0.7, W=1.16,
                CP_fraction_of_gravity=0.1, tau_coll_order_s=3.5),
)
ph = differential_phases(bose['m'], bose['m'], bose['tau'], bose['d'], bose['dx'])
# magnitudes (Bose eq. 4 uses +Gmmτ/ħr)
ph_mag = dict(phi_d=grav_phase_mag(bose['m'], bose['m'], bose['tau'], bose['d']),
              phi_closest=grav_phase_mag(bose['m'], bose['m'], bose['tau'], bose['closest']),
              phi_farthest=grav_phase_mag(bose['m'], bose['m'], bose['tau'], bose['d'] + bose['dx']))
# If one uses unsubtracted magnitudes with a minus on the far arm: (-φ_LR, +φ_RL) ~ (-0.23, +0.79)
W_derived = pauli_correlators(ph['dphi_lr'], ph['dphi_rl'])
W_quoted = pauli_correlators(bose['quoted']['dphi_LR'], bose['quoted']['dphi_RL'])
C_derived = concurrence_bmv(ph['dphi_lr'], ph['dphi_rl'])
C_quoted = concurrence_bmv(bose['quoted']['dphi_LR'], bose['quoted']['dphi_RL'])

# intro estimate: m=1e-14, d=100 um, τ=1 s (Bose opening paragraph)
phi_intro = grav_phase_mag(1e-14, 1e-14, 1.0, 100e-6)

# Casimir at closest approach
R_bose = bose['R']
m_from_R = 4.0 / 3.0 * pi * R_bose**3 * RHO_DIAMOND
cp_ratio_bose_m = cp_over_grav(bose['m'], R_bose, bose['closest'])
cp_ratio_density = cp_over_grav(m_from_R, R_bose, bose['closest'])
r_min_01 = r_min_casimir(bose['m'], R_bose, 0.1)
r_min_01_dens = r_min_casimir(m_from_R, R_bose, 0.1)

# Collisions: show T_gas = 300 K reproduces ~3.5 s; T_gas = 0.15 K is harsher at fixed P
g300 = gamma_coll_resolution(bose['P'], 300.0, R_bose, M_N2, 1.0)
g015 = gamma_coll_resolution(bose['P'], 0.15, R_bose, M_N2, 1.0)
# Helium at 0.15 K (cryogenic residual)
M_HE = 4.002602 * U_DA
g015_He = gamma_coll_resolution(bose['P'], 0.15, R_bose, M_HE, 1.0)

# Blackbody at T_env
dx = bose['dx']
bb = {T: dict(gamma_scatt_quadratic_or_sat=gamma_bb_scatt(T, R_bose, dx),
              gamma_scatt_rate=gamma_bb_scattering_rate(T, R_bose),
              lambda_th_m=HBAR * C / (K_B * T))
      for T in (0.15, 1.0, 4.0, 20.0, 77.0, 300.0)}

# Magnetic gradient needed for Δx (Bose eq. 7): Δx = (1/2) (g μB / m) ∂xB  τ_acc²
MU_B = 9.2740100783e-24
g_e = 2.0
# invert: dB/dx = 2 m Δx / (g μB τ_acc²)
dBdx_needed = 2.0 * bose['m'] * bose['dx'] / (g_e * MU_B * bose['tau_acc']**2)
# Bose quotes ~1e6 T/m. Check.
grav_accel_closest = G * bose['m'] / bose['closest']**2

# ---------------------------------------------------------------------------
# 2. Parameter-space maps
# ---------------------------------------------------------------------------
def geometry_for_mass(m, rho=RHO_DIAMOND, cp_ratio=0.1, dx_mode='bose_ratio'):
    """Casimir-safe geometry. dx_mode='bose_ratio' keeps Δx/(d-Δx)=250/200, d/(d-Δx)=450/200.
    dx_mode='fixed_dx' keeps Δx=250 μm and sets d = r_min+Δx (not always possible if r_min large)."""
    R = sphere_radius(m, rho)
    rmin = r_min_casimir(m, R, cp_ratio)
    if dx_mode == 'bose_ratio':
        closest = rmin
        d = 450.0 / 200.0 * closest
        dx = 250.0 / 200.0 * closest
    else:
        dx = 250e-6
        closest = rmin
        d = closest + dx
    return dict(m=m, R=R, rmin=rmin, d=d, dx=dx, closest=closest)


def t_for_phase(m, phase_target, kind='entangle_sum_mag', rho=RHO_DIAMOND, cp_ratio=0.1):
    """Interaction time needed so the chosen phase reaches phase_target (DERIVED)."""
    geo = geometry_for_mass(m, rho, cp_ratio)
    d, dx = geo['d'], geo['dx']
    if kind == 'simple':
        # Δφ = G m² t / (ħ d)
        return phase_target * HBAR * d / (G * m * m), geo
    if kind == 'closest':
        return phase_target * HBAR * geo['closest'] / (G * m * m), geo
    # |Δφ_LR + Δφ_RL| = (G m² t / ħ) * 2 Δx² / (d (d²-Δx²))
    pre = (G * m * m / HBAR) * 2.0 * dx**2 / (d * (d**2 - dx**2))
    return phase_target / pre, geo


# grid of masses
masses = np.logspace(-24, -10, 281)
records = []
for m in masses:
    rec = {'m': float(m)}
    for target in (0.01, 1.0):
        for kind in ('simple', 'entangle_sum_mag'):
            t, geo = t_for_phase(m, target, kind)
            rec[f't_{kind}_{target}'] = float(t)
            if kind == 'entangle_sum_mag':
                rec.update({f'geo_{k}': float(geo[k]) for k in ('R', 'd', 'dx', 'closest', 'rmin')})
    # required pressure so Γ_coll * t_0.01 = 0.1, T_gas=300 K, geometric σ
    t01, geo = t_for_phase(m, 0.01, 'entangle_sum_mag')
    # Γ = n <v> π R² = P/(kT) * v * π R²  => P = 0.1 / t * kT / (v π R²)
    v = vbar_gas(300.0)
    P_req_300 = (0.1 / t01) * K_B * 300.0 / (v * pi * geo['R']**2) if t01 > 0 else np.nan
    rec['P_max_for_t01_Tgas300_Pa'] = float(P_req_300)
    # T_env such that Γ_bb(T) * t_0.01 = 0.1 (binary search)
    lo, hi = 1e-3, 400.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if gamma_bb_scatt(mid, geo['R'], geo['dx']) * t01 > 0.1:
            hi = mid
        else:
            lo = mid
    rec['T_env_max_for_t01_K'] = float(0.5 * (lo + hi))
    records.append(rec)

# Current experimental markers (CITED)
# Fein et al., Nature Physics 15, 1242 (2019): oligoporphyrins >25 kDa, up to ~2000 atoms, λdB down to 53 fm.
# Pedalino et al., arXiv:2507.21211 (Nature 2026 per journal record s41586-025-09917-9): Na clusters, mean 172 kDa.
# Delić et al., Science 367, 892 (2020): ~150 nm-diameter silica, ground-state cooling, not a spatial cat.
# Tebbenjohanns et al., Nature / arXiv:2103.03853 (2021): silica, measured diameter 106±5 nm, m ~ 1 fg = 1e-18 kg.
experiments = {
    'Fein_2019_25kDa_interference': dict(
        m=25000 * U_DA, dx=266e-9,  # optical grating period scale (ASSUMED typical KDTLI period; delocalization is at least this)
        t=2.0 / 300.0,  # 2 m interferometer, v ~ 300 m/s (CITED order, Fein)
        kind='matter-wave interference',
        cite='Fein et al., Nat. Phys. 15, 1242 (2019)'),
    'Pedalino_2025_172kDa_Na_interference': dict(
        m=172000 * U_DA, dx=133e-9,  # 133 nm gratings (CITED arXiv:2507.21211)
        t=None,  # v ~ 160 m/s; length not taken as a single number here
        v=160.0,
        kind='matter-wave interference',
        cite='Pedalino, Ramírez-Galindo & Ferstl, arXiv:2507.21211'),
    'Delic_2020_silica_ground_state': dict(
        m=None, diameter=150e-9, rho=2200.0,  # ~150 nm diameter silica
        omega=2 * pi * 305e3,  # CITED Ωx/2π ≈ 305 kHz
        kind='levitated ground-state cooling (NOT a spatial cat)',
        cite='Delić et al., Science 367, 892 (2020), arXiv:1911.04406'),
    'Tebbenjohanns_2021_silica_ground_state': dict(
        m=1e-18, diameter=106e-9,  # m∼1 fg = 1e-18 kg, measured 106±5 nm
        kind='levitated ground-state cooling (NOT a spatial cat)',
        cite='Tebbenjohanns et al., arXiv:2103.03853; Nature 595, 373 (2021)'),
    'Westphal_2021_classical_gravity': dict(
        m_source=92.1e-6, m_test=90.7e-6,  # 92.1 mg and 90.7 mg gold
        kind='classical gravity measurement (no superposition)',
        cite='Westphal, Hepach, Pfaff & Aspelmeyer, Nature 591, 225 (2021), arXiv:2009.09546'),
    'Bose_2017_target': dict(
        m=1e-14, dx=250e-6, t=2.5, d=450e-6,
        kind='QGEM / BMV proposal target',
        cite='Bose et al., PRL 119, 240401 (2017)'),
}
# fill derived masses / zpf
RHO_SIO2 = 2200.0
experiments['Delic_2020_silica_ground_state']['m'] = (4/3)*pi*(0.5*150e-9)**3 * RHO_SIO2
experiments['Delic_2020_silica_ground_state']['x_zpf'] = sqrt(
    HBAR / (2 * experiments['Delic_2020_silica_ground_state']['m'] * experiments['Delic_2020_silica_ground_state']['omega']))
experiments['Tebbenjohanns_2021_silica_ground_state']['x_zpf_if_300kHz'] = sqrt(
    HBAR / (2 * 1e-18 * 2 * pi * 300e3))

# Phase two such objects would accumulate at a Casimir-safe Bose-ratio geometry, using a generous t
# of 10 ms for molecular interferometers (order of flight) and 1 s for levitated (trap lifetime, still not a cat).
phase_at_records = {}
for name, ex in experiments.items():
    if 'm' not in ex or ex['m'] is None or ex['m'] >= 1e-8:
        continue
    m = ex['m']
    if 'interference' in ex['kind']:
        t_use = ex.get('t') or 0.01
        # even if we give them the Bose time, report both
    elif 'cooling' in ex['kind']:
        t_use = 1.0
    else:
        t_use = ex.get('t', 1.0)
    t_need_01, geo = t_for_phase(m, 0.01, 'entangle_sum_mag')
    t_need_1, _ = t_for_phase(m, 1.0, 'entangle_sum_mag')
    ph_rec = differential_phases(m, m, t_use, geo['d'], geo['dx'])
    phase_at_records[name] = dict(
        m=m, t_used=t_use, t_for_0p01=t_need_01, t_for_1=t_need_1,
        dphi_sum_mag=abs(ph_rec['dphi_sum']),
        gap_orders_vs_0p01=log10(0.01 / abs(ph_rec['dphi_sum'])) if ph_rec['dphi_sum'] != 0 else None,
        geometry=geo, kind=ex['kind'], cite=ex['cite'])

# Orders of magnitude: Bose target vs records
m_bose = 1e-14
m_fein = 25000 * U_DA
m_na = 172000 * U_DA
m_teb = 1e-18
m_delic = experiments['Delic_2020_silica_ground_state']['m']
gap = dict(
    mass_orders_vs_Fein=log10(m_bose / m_fein),
    mass_orders_vs_Pedalino172kDa=log10(m_bose / m_na),
    mass_orders_vs_Tebbenjohanns_cooled=log10(m_bose / m_teb),
    mass_orders_vs_Delic_cooled=log10(m_bose / m_delic),
    note=('Mass gap vs matter-wave interference is ~8 orders (Fein 25 kDa) or ~7.5 (Na 172 kDa). '
          'The cooled-nanoparticle mass gap is only ~4 orders, but those particles are in a ground-state '
          'oscillator state with Δx ~ picometres, not a 250 μm cat. Phase scales as m² t / d, so the '
          'phase gap is much larger than the mass gap (see phase_at_records).'),
)

# Required P, T at the Bose point
P_bose_for_10pct = (0.1 / 3.5) * K_B * 300.0 / (vbar_gas(300.0) * pi * R_bose**2)
# invert Γ_bb(T)*3.5s = 0.1
lo, hi = 1e-3, 400.0
for _ in range(60):
    mid = 0.5 * (lo + hi)
    if gamma_bb_scatt(mid, R_bose, dx) * 3.5 > 0.1:
        hi = mid
    else:
        lo = mid
T_env_bose_10pct = 0.5 * (lo + hi)

# ---------------------------------------------------------------------------
# 4. Oppenheim postquantum / CQ gravity: decoherence–diffusion tradeoff
#    CITED: Oppenheim, Sparaciari, Šoda & Weller-Davies, Nat. Commun. 14, 7910
#    (2023), arXiv:2203.01982; updated upper bounds Janse et al., PRR 6, 033076
#    (2024), arXiv:2403.08912.
# ---------------------------------------------------------------------------
def oppenheim_bounds(sigma_a=1e-7, T_avg=100.0, N=1e26, rN=1e-15,
                     Vb=1e15, mN=None,
                     M_lambda=1e-24, V_lambda=1e-25, R_lambda=1e-9, lam=10.0):
    """Reproduce Oppenheim et al. conservative numbers (DERIVED from their inputs).
    mN: if None, use gold-nucleus mass so that N~1e26 would be inconsistent;
    they used N~1e26 (nucleon-like) in one place and got 10^{-1} kg s in another,
    which matches using nuclear mass in the discrete upper bound. We compute both.
    """
    out = {}
    # Ultra-local continuous: D2(x,y)=D2 δ(x-y)
    # upper: D2 <= σa² N rN⁴ T / (Vb G²)
    upper_cont = (sigma_a**2 * N * rN**4 * T_avg) / (Vb * G**2)
    lower_cont = (M_lambda**2) / (V_lambda * lam)
    out['ultralocal_continuous'] = dict(
        units='kg^2 s m^{-3}',
        upper_Cavendish=upper_cont, lower_interferometry=lower_cont,
        ratio_lower_over_upper=lower_cont / upper_cont,
        already_excluded=bool(lower_cont > upper_cont),
        inputs=dict(sigma_a=sigma_a, T_avg=T_avg, N=N, rN=rN, Vb=Vb,
                    M_lambda=M_lambda, V_lambda=V_lambda, lam=lam))
    # Ultra-local discrete: (lP^3/mP) D2 squeezed between
    # upper = σa² N T rN⁴ / (mN G²), lower = M_λ / λ
    def discrete(mN):
        upper = (sigma_a**2 * N * T_avg * rN**4) / (mN * G**2)
        lower = M_lambda / lam
        return dict(units='kg s  [this is (l_P^3/m_P) D2]',
                    upper=upper, lower=lower, gap_orders=log10(upper / lower) if upper > lower else None,
                    already_excluded=bool(lower > upper), mN=mN)
    out['ultralocal_discrete_nucleon_mN'] = discrete(M_NUCLEON)
    # gold nucleus ~ 197 u
    out['ultralocal_discrete_gold_nucleus_mN'] = discrete(197 * U_DA)
    # Non-local continuous (Diósi-Penrose-like kernel): lP² D2
    # upper = σa² N T rN³ / G², lower = M_λ² / (R_λ λ)
    upper_nl = (sigma_a**2 * N * T_avg * rN**3) / (G**2)
    lower_nl = (M_lambda**2) / (R_lambda * lam)
    out['nonlocal_continuous'] = dict(
        units='kg^2 s m^{-1}  [this is l_P^2 D2]',
        upper_Cavendish=upper_nl, lower_interferometry=lower_nl,
        gap_orders=log10(upper_nl / lower_nl) if upper_nl > lower_nl else None,
        already_excluded=bool(lower_nl > upper_nl),
        note='Oppenheim quotes 10^{-9} >= lP^2 D2 >= 10^{-40} kg s m^{-1} with these inputs.')
    return out

opp_original = oppenheim_bounds()
# Update lower bound with Pedalino 172 kDa, conservative: if they held coherence for ~1 ms
# (flight), M=2.86e-22 kg. Volume: Na cluster of 7000 atoms, r ~ a few nm.
# ASSUMED V_λ: 5 nm x 5 nm x (delocalization ~ 1 um) as a generous wavepacket volume, analogous
# to Oppenheim's fullerene estimate (1 nm x 1 nm x 100 nm). This is a rough estimate, flagged.
m_na = 172000 * U_DA
# 172 kDa Na: N_atoms = 172000/23 ≈ 7478, r_Wigner-Seitz Na ~ 0.21 nm, R ~ 4 nm
R_na = 4e-9
V_na_oppenheim_style = (4e-9) * (4e-9) * (1e-6)  # 1 um delocalization ASSUMED
# Flight time: not given as a single T; use 1 ms conservative (shorter => weaker lower bound? 
# λ < 1/T_coh, so longer coherence gives a tighter (higher) lower bound on D2? 
# lower D2 >= M^2 / (V λ) and λ < 1/T so lower >= M^2 T / V, longer coherence raises the floor.
# Conservative (weaker) lower bound: shorter T, larger V.
opp_na_conservative = oppenheim_bounds(M_lambda=m_na, V_lambda=V_na_oppenheim_style,
                                       R_lambda=R_na, lam=1.0 / 0.001)
# Fein 25 kDa, 2m / 300 m/s ~ 7 ms, R~1 nm, delocalization ASSUMED 1 um
opp_fein = oppenheim_bounds(M_lambda=25000 * U_DA, V_lambda=(1e-9)*(1e-9)*(1e-6),
                            R_lambda=1e-9, lam=1.0 / 7e-3)

# Janse et al. 2024: FOM improved by ~15 orders vs Cavendish for absolute Earth-bound
# acceleration noise (Gisler nanowire). We scale only the upper bounds by 1e-15, with a caveat.
opp_janse_note = (
    'Janse et al., Phys. Rev. Research 6, 033076 (2024), arXiv:2403.08912: the Cavendish '
    'figure of merit FOM=σa² N ΔT used by Oppenheim is 1e14 m²/s³; Gisler et al. nanowire '
    '(m_eff=9.3e-15 kg) improves FOM by ~15 orders (absolute accel. noise on Earth). '
    'LISA Pathfinder and atom interferometry look even better but are differential / relative '
    'measurements, and atom interferometers are also the experiments that set the lower bound, '
    'so using them on both sides is questionable (their caveat). We therefore only scale the '
    'upper bounds by 1e-15 as a CITED estimate for the best uncontroversial tightening, and '
    'do not claim the remaining models are ruled out.')

# ---------------------------------------------------------------------------
# Assemble JSON
# ---------------------------------------------------------------------------
out = dict(
    title='Stage 11: gravity-mediated entanglement (BMV / QGEM)',
    constants=dict(G=G, hbar=HBAR, k_B=K_B, c=C, u=U_DA, l_P=L_P, m_P=M_P,
                   source='CODATA 2018/2022'),
    section1_bose_reproduction=dict(
        inputs_CITED=bose,
        phases_DERIVED=dict(**{k: float(v) for k, v in ph.items()}, **{k: float(v) for k, v in ph_mag.items()}),
        intro_phi_mag_100um_1s=float(phi_intro),
        witness=dict(
            derived_phases=dict(concurrence=float(C_derived), **{k: (float(v) if not isinstance(v, (bool, str)) else v)
                                                                for k, v in W_derived.items()}),
            quoted_phases=dict(concurrence=float(C_quoted), **{k: (float(v) if not isinstance(v, (bool, str)) else v)
                                                              for k, v in W_quoted.items()}),
            discrepancy=('Bose quotes Δφ_LR∼−0.2, Δφ_RL∼0.7, W∼1.16. With V=−Gmm/r and the stated '
                         'geometry we get Δφ_LR=+0.126, Δφ_RL=−0.440 (signs from attraction). The magnitudes '
                         'of the unsubtracted phases are φ(d)=0.352, φ_closest=0.791, φ_farthest=0.226, so '
                         'their −0.2 and 0.7 look like rounded unsubtracted φ_LR and φ_RL (0.226 and 0.791) '
                         'with an unexplained minus on the far arm, not the subtracted Δφ of their eq. (2). '
                         'The written witness W=|⟨σx⊗σz⟩−⟨σy⊗σz⟩| equals 0.107 (derived phases) or 0.330 '
                         '(their quoted phases), not 1.16. We could not reconstruct 1.16 from the written '
                         'definition. Independently, W>1 is not a valid entanglement test: product states '
                         'reach W=√2. A valid test for this family is concurrence C=|sin((Δφ_LR+Δφ_RL)/2)|, '
                         'which is 0.156 (derived) or 0.247 (quoted phases). Both are entangled, just not '
                         'by the W>1 argument as written.')),
        casimir=dict(
            formula='U_CP=-(23 ħ c /4π) [R^3 (ε-1)/(ε+2)]^2 / r^7  (SI, retarded two-sphere)',
            ratio_at_200um_R1um_m1e14=float(cp_ratio_bose_m),
            ratio_if_m_from_diamond_density=float(cp_ratio_density),
            m_of_1um_diamond_kg=float(m_from_R),
            r_min_for_ratio_0p1_m=float(r_min_01),
            r_min_density_consistent_m=float(r_min_01_dens),
            note=('Bose says CP is a tenth of gravity at 200 μm. We get 0.050 with their (m,R) pair, '
                  'or 0.023 if m is taken from diamond density at R=1 μm. Same order; the 0.1 is a '
                  'rounded estimate and their (m,R) pair is not density-consistent (1 μm of diamond '
                  'is 1.47e-14 kg).')),
        gradient_and_drop=dict(
            dBdx_T_per_m=float(dBdx_needed),
            bose_quoted_T_per_m=1e6,
            match_orders=True,
            grav_accel_at_closest_m_s2=float(grav_accel_closest),
            bose_quoted=1e-16,
            drop_time_total_s=3.5,
            drop_height_half_g_t2_m=0.5 * 9.81 * 3.5**2,
            note='Bose quotes 1e6 T/m and g_grav~1e-16 m/s^2; both match. 3.5 s free fall is 60 m, so drop tower or space.'),
        decoherence_at_bose_point=dict(
            P_Pa=1e-15, T_internal_K=0.15,
            tau_coll_N2_Tgas300K_s=float(1.0 / g300),
            tau_coll_N2_Tgas0p15K_s=float(1.0 / g015),
            tau_coll_He_Tgas0p15K_s=float(1.0 / g015_He),
            note_coll=('At P=1e-15 Pa with geometric πR² collisions on a 1 μm sphere, the decoherence '
                       'time is 2.8e3 s for N2 at 300 K and 62 s for N2 at 0.15 K (DERIVED). Both are '
                       'comfortably longer than the 3.5 s drop: Γ t = 0.0013 (300 K) or 0.056 (0.15 K). '
                       'Bose’s “same order as 3.5 s” is therefore pessimistic by 1–3 orders at this P. '
                       'Γ t = 1 at 300 K would require P ≈ 8e-13 Pa, still very good vacuum but not XHV.'),
            blackbody_by_Tenv=bb,
            T_env_for_Gamma_t_eq_0p1_over_3p5s_K=float(T_env_bose_10pct),
            P_max_Gamma_t_eq_0p1_over_3p5s_Tgas300_Pa=float((0.1 / 3.5) / (g300 / bose['P'])),
            bb_note=('At 0.15 K, λ_th≈1.5 cm ≫ 250 μm, so photon scattering is in the soft quadratic '
                     'regime and Γ t is ~10^{-17}, negligible, as Bose says. Room-temperature radiation '
                     'would destroy the superposition in nanoseconds (resolution-limited scattering rate '
                     '~10^{11}/s on a 1 μm dielectric). Internal temperature must also be cryogenic, or '
                     'the object’s own thermal emission does the same. The environment can actually be '
                     'several kelvin before photon scattering competes with 3.5 s (see T_env_max).'),
        ),
    ),
    section2_parameter_space=dict(
        geometry='For each mass, R from diamond density 3510 kg/m^3; closest approach from |U_CP|/|U_G|=0.1; '
                 'Δx/closest=250/200, d/closest=450/200 (Bose ratios, ASSUMED shape).',
        masses=dict(min=float(masses[0]), max=float(masses[-1]), n=len(masses)),
        grid_file_keys='see parameter_grid (downsampled in this json: every 20th point; full curves are in the figure)',
        parameter_grid=[records[i] for i in range(0, len(records), 20)],
        experiments=experiments,
        phase_at_records=phase_at_records,
        gap_orders=gap,
        bose_point_requirements=dict(
            P_max_10pct_loss_3p5s_Tgas300_Pa=float((0.1 / 3.5) / (g300 / bose['P'])),
            T_env_max_10pct_loss_3p5s_K=float(T_env_bose_10pct),
        ),
    ),
    section3_status=dict(
        experiment_done=False,
        as_of='2026-11',
        groups_CITED=[
            'Bose / UCL and the QGEM white-paper consortium (arXiv:2509.01586, 2025): building traps for micron diamonds / nanodiamonds; first spatial superpositions expected at 10–1000 fm, then push to 10–100 μm; two-diamond setups and EM screening in parallel. “Significant steps … within the next decade” with sufficient resources. Underground lab may be needed for vibration / GGN.',
            'Folman (Ben-Gurion): Stern–Gerlach interferometry with atoms already demonstrated (e.g. Science Advances 7, eabg2879 (2021)); contributing to the spin-based diamond pathway (white paper coauthor).',
            'Aspelmeyer (Vienna / IQOQI): levitated optomechanics, ground-state cooling (Delić 2020); smallest source-mass gravity measurement, 90.7 mg test / 92.1 mg source gold spheres (Westphal et al., Nature 591, 225 (2021)). Public timeline ~15 years for a viable entanglement experiment (New Scientist, 12 May 2025), described as a self-imposed deadline.',
            'Geraci (now Northwestern / formerly UNR): optically levitated nanoparticles and short-range gravity; Bose 2017 and white-paper coauthor.',
            'Arndt (Vienna): heaviest matter-wave interference — 25 kDa molecules (Fein 2019) and 172 kDa sodium clusters (Pedalino et al., arXiv:2507.21211). Complementary mass path, still many orders below QGEM.',
            'Others on the white paper: Ulbricht (Southampton), Barker (UCL), Kim (Imperial), Milburn, Kent, Purdy, Dutt, and several theory groups.',
        ],
        alternatives_CITED=[
            'Krisnanda, Tham, Paternostro & Paterek, arXiv:1906.08808: optomechanical (always-trapped) masses, or masses released from traps. Trapped version wants LIGO-like mirrors plus squeezing or very long coherence; released masses can be lighter. No macroscopic cat is required in their Gaussian-state analysis.',
            'Pedernales & Plenio, arXiv:2311.04745 (and earlier Pedernales/Plenio nanodiamond work): force sensitivity of delocalised mechanics is set by spatial delocalisation; cat states and momentum-squeezed Gaussians entangle at the same rate if equally delocalised. This is the “mass / delocalisation determines the rate” alternative to spin-SG cats, not a different gravitational mechanism.',
            'Casimir screening (van de Kamp, Marshman, Bose & Mazumdar, arXiv:2006.06931; white paper): a conducting plate can suppress EM backgrounds and may allow a smaller d, lowering the mass or time demand.',
        ],
        timeline_CITED=dict(
            white_paper_2025='Significant steps toward QGEM possible in the next decade (i.e. mid-2030s), not a claim that the witness will be measured by then. First milestone is any spatial superposition of a nanoparticle, even 10–1000 fm.',
            aspelmeyer_2025='~15 years from May 2025, i.e. around 2040, for a viable experiment; self-imposed. He also said a strictly spacelike-separated, “foolproof” version is not realistic on any human timescale (New Scientist, 12 May 2025).',
            honest_read=('No group has a 10^{-14} kg object in a ~250 μm superposition for seconds, which is the actual requirement. The mass gap in interference is 7–8 orders; the phase gap is larger because Δφ ∝ m². A completed BMV test before the late 2030s would be a surprise. Collapse-model tests with a single large superposition can come earlier and already matter.'),
        ),
        ggn_CITED='Toroš, Marshman, Mazumdar & Bose, Phys. Rev. Research 3, 023178 (2021), arXiv:2007.15029: gravity-gradient noise is not a fundamental blocker. Seismic/atmospheric are small; human/vehicle/aircraft can be kept out by access restrictions of about 2 m / 10 m / 60 m. Underground helps. This is a site and isolation problem, not a no-go.',
    ),
    section4_interpretation=dict(
        positive_result=dict(
            locc=('Marletto–Vedral / Bose: if two quantum systems become entangled and the only interaction is gravitational, and if that interaction is mediated (not action at a distance), then the mediator cannot be a classical channel, because LOCC cannot create entanglement. That is the clean core.'),
            what_it_shows=('It shows that gravity can entangle quantum matter, so gravity is not a purely classical communication channel between quantum systems. In linearized GR language, Bose et al. argue this needs off-diagonal amplitudes in the gravitational-field coherent-state basis.'),
            what_it_does_not_show=[
                'It does not detect gravitons as particles, and does not reach Planck-scale quantum gravity.',
                'Fragkos, Kopp & Pikovski, AVS Quantum Sci. 4, 045601 (2023), arXiv:2206.00558: the entanglement is generated by the Newtonian interaction (a constraint / Poisson-gauge term), which can be written as a non-local stress-energy coupling with no local mediator. Relativity does not force a unique local-channel reading. A positive result would still be the first test of how a spatial superposition sources a Newtonian field, which has never been done.',
                'Anastopoulos & Hu: the 2018 arXiv comment (1804.11315) was withdrawn (arXiv lists it as withdrawn; copyright, not a retraction of the scientific point). The later paper “Gravity, Quantum Fields and Quantum Information”, arXiv:2202.02789, argues that in the weak-field Newtonian limit the potential is slaved to the mass density, so seeing entanglement from that term does not by itself prove that the true radiative degrees of freedom (gravitons) are quantum.',
                'A laboratory BMV protocol with d~200 μm and t~2.5 s is not spacelike separated (d/c ~ 0.7 ps). A classical or quantum non-gravitational process that can act over that time is not forbidden by causality; it has to be excluded by shielding and control experiments. Aspelmeyer’s “foolproof” spacelike version is far harder.',
            ],
        ),
        theories_killed=dict(
            positive_entanglement=[
                'Fully classical gravitational channels (gravity as LOCC) and classical-only Newton potential acting as a c-number field with no quantum coherence.',
                'Objective collapse models that would have already decohered a 10^{-14} kg / 250 μm superposition on the 1 s timescale (Diósi–Penrose and similar), if such a superposition is first created at all. Aspelmeyer, C. R. Physique 27, 1 (2026) argues a positive GIE result rules out current gravitational-collapse models once lab constraints are included.',
                'Oppenheim-type postquantum classical gravity is not automatically killed by a GIE detection in every formulation; the debate is live (see critiques and later CQ papers). A detected entanglement with coherence surviving the CQ tradeoff would put severe pressure on those models, but we do not claim a theorem here.',
            ],
            negative_or_null=[
                'A null result at insufficient mass/time/vacuum does not kill quantum gravity; it is just a failed attempt.',
                'A null result in a regime where a quantum Newtonian interaction should have entangled the masses, after EM backgrounds are demonstrated to be under control, would be a genuine crisis for standard quantized linearized gravity and would favour collapse models or a classical gravitational channel.',
            ],
        ),
        oppenheim=dict(
            statement=('Oppenheim et al., Nat. Commun. 14, 7910 (2023), arXiv:2203.01982: any consistent classical-gravity / quantum-matter coupling has a tradeoff. Long quantum coherence requires spacetime diffusion (noise in the metric / Newtonian potential). Precision mass/acceleration measurements cap the diffusion from above; matter-wave coherence caps the decoherence rate and so floors the diffusion. The two together squeeze the theory.'),
            original_bounds_recomputed=opp_original,
            lower_bound_with_Fein2019=opp_fein,
            lower_bound_with_Pedalino172kDa_ASSUMED_volume=opp_na_conservative,
            janse2024_cite=opp_janse_note,
            verdict=('Ultra-local continuous noise (D2 δ(x−y) everywhere) is already excluded by Oppenheim’s own conservative Cavendish-vs-fullerene squeeze: we recompute lower/upper ≈ 10^{-24} / 10^{-41} kg² s m^{-3}, a 17-order overlap. The ultra-local discrete and non-local continuous (Diósi–Penrose-like) models are not excluded on those conservative numbers (gaps of tens of orders). Janse et al. 2024 tighten the upper side by ~15 orders with modern absolute force-noise data, which shrinks but does not close those gaps if one refuses differential measurements. This is the nearer-term test: it uses existing torsion, nanowire, and interferometry data, and improves every time a quieter force sensor or a heavier coherent superposition is published. It does not require two 10^{-14} kg cats.'),
        ),
    ),
    section5_amateur_contributions=dict(
        useful=[
            'Independent decoherence budgets and design optimisation: map (m, R, d, Δx, t, P, T_gas, T_env, T_int, gas species, dielectric function) to Δφ, concurrence, and Γ t, including Casimir-safe d and optional conductive shields (van de Kamp et al.). This is finite-dimensional numerical work and is exactly what proposals need checked by a second calculation.',
            'Recompute and stress-test published feasibility numbers. We already found that the headline W∼1.16 does not follow from the written witness and quoted phases, and that Bose’s CP “tenth” is ~0.05 with the standard SI formula. An open, cited errata/clarification note would be useful.',
            'Oppenheim squeeze with current data: update the lower bound with each new interference mass (25 kDa, 172 kDa, and whatever comes next) and the upper bound with each new low-force-noise measurement, being explicit about absolute vs differential and about the wavepacket volume. Janse et al. already started this; it is ongoing and checkable.',
            'Open problems listed in the 2025 QGEM white paper (arXiv:2509.01586) that are calculation- or simulation-heavy: EM screening / patch potentials, Casimir with a real finite-conductivity plate between two diamonds, surface-spin noise on nanodiamond (they quote 0.01 spins/nm² ⇒ >10^4 unpaired spins on a nanodiamond), gravity-gradient and vibration transfer functions, motional dynamical decoupling sequences, and whether dark-matter or collapse-model decoherence would dominate before gravity entangles.',
        ],
        not_useful=[
            'Inventing a new unified-field picture without a number that one of these experiments can kill.',
            'Trying to simulate Planck-scale quantum gravity or a first-principles graviton S-matrix for this apparatus: the entanglement, if it happens, is a low-energy Newtonian / linearized-GR effect.',
        ],
    ),
)

# store a denser grid separately would bloat; the figure uses the full `records` in memory
json.dump(out, open('results/stage11.json', 'w'), indent=1, default=float)

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0))

# --- left: required time vs mass, Casimir-safe Bose-ratio geometry
ax = axes[0]
m_grid = np.array([r['m'] for r in records])
t01 = np.array([r['t_entangle_sum_mag_0.01'] for r in records])
t1 = np.array([r['t_entangle_sum_mag_1.0'] for r in records])
tsim01 = np.array([r['t_simple_0.01'] for r in records])
ax.loglog(m_grid, t1, color='C0', lw=2, label=r'$|\Delta\phi_{LR}+\Delta\phi_{RL}|=1$ rad')
ax.loglog(m_grid, t01, color='C0', lw=2, ls='--', label=r'$|\Delta\phi_{LR}+\Delta\phi_{RL}|=0.01$ rad')
ax.loglog(m_grid, tsim01, color='C0', lw=1, ls=':', label=r'simple $G m^2 t/(\hbar d)=0.01$')
# Bose point
ax.scatter([1e-14], [2.5], marker='*', s=160, color='C3', zorder=5, label='Bose 2017 target (2.5 s, $10^{-14}$ kg)')
ax.axhline(2.5, color='C3', lw=0.6, ls=':')
# records: flight / hold times
ax.axvspan(m_fein, m_na, color='0.7', alpha=0.35, label='matter-wave interference (25–172 kDa)')
ax.axvline(m_fein, color='0.4', lw=0.8, ls='--')
ax.axvline(m_na, color='0.4', lw=0.8, ls='--')
ax.axvspan(m_teb, m_delic, color='C2', alpha=0.25, label='levitated ground-state cooling (~$10^{-18}$ kg)')
ax.set_xlabel('mass of each object  [kg]')
ax.set_ylabel('interaction time needed  [s]')
ax.set_xlim(1e-24, 1e-11)
ax.set_ylim(1e-6, 1e16)
ax.set_title('Casimir-safe geometry (diamond, $U_{CP}/U_G=0.1$)\nBose path ratios  (DERIVED)')
ax.legend(fontsize=7, loc='upper right', framealpha=0.92)
ax.text(m_fein, 3e-6, 'Fein\n25 kDa', fontsize=7, ha='right', va='bottom', color='0.2')
ax.text(m_na, 3e-6, 'Na\n172 kDa', fontsize=7, ha='left', va='bottom', color='0.2')
ax.text(3e-18, 3e-6, 'cooled\nNP', fontsize=7, ha='center', va='bottom', color='C2')
# annotate the 0.01 rad crossing at 2.5 s (the 1 rad crossing is essentially the Bose star)
i = int(np.argmin(np.abs(t01 - 2.5)))
ax.annotate(f"0.01 rad at 2.5 s\nm={m_grid[i]:.1e} kg",
            xy=(m_grid[i], 2.5), xytext=(2e-17, 4e3),
            fontsize=7, arrowprops=dict(arrowstyle='->', color='0.3'), color='0.2')

# --- right: mass vs superposition size, records vs target
ax = axes[1]
# shade the Bose-like target region
ax.scatter([1e-14], [250e-6], marker='*', s=200, color='C3', zorder=6, label='Bose 2017 QGEM target')
ax.annotate('m = $10^{-14}$ kg\n$\\Delta x$ = 250 μm\n$t$ = 2.5 s', xy=(1e-14, 250e-6),
            xytext=(3e-17, 3e-3), fontsize=8, color='C3',
            arrowprops=dict(arrowstyle='->', color='C3'))
# interference records
ax.scatter([m_fein], [266e-9], marker='o', s=70, color='k', zorder=5, label='Fein 2019 interference (25 kDa)')
ax.scatter([m_na], [133e-9], marker='s', s=70, color='k', zorder=5, label='Pedalino 2025 interference (172 kDa)')
# cooled: zpf, not a cat — plot as open triangles
ax.scatter([m_delic], [experiments['Delic_2020_silica_ground_state']['x_zpf']],
           marker='^', s=70, facecolors='none', edgecolors='C2', linewidths=1.5, zorder=5,
           label='Delić 2020 ground-state (zpf, not a cat)')
ax.scatter([m_teb], [experiments['Tebbenjohanns_2021_silica_ground_state']['x_zpf_if_300kHz']],
           marker='v', s=70, facecolors='none', edgecolors='C2', linewidths=1.5, zorder=5,
           label='Tebbenjohanns 2021 ground-state (zpf, not a cat)')
# gap arrows
ax.annotate('', xy=(1e-14, 250e-6), xytext=(m_na, 133e-9),
            arrowprops=dict(arrowstyle='->', color='C3', lw=1.2, ls='--'))
ax.text(3e-19, 3e-6,
        f"mass gap vs 172 kDa: {log10(m_bose/m_na):.1f} orders\n"
        f"mass gap vs cooled NP: {log10(m_bose/m_teb):.1f} orders\n"
        f"but those NPs are not 250 μm cats\n"
        f"phase $\\propto m^2 t/d$: gap is larger still",
        fontsize=7, color='0.15', bbox=dict(facecolor='white', edgecolor='0.8', alpha=0.9))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-24, 1e-11)
ax.set_ylim(1e-15, 1e-2)
ax.set_xlabel('mass  [kg]')
ax.set_ylabel('spatial delocalisation $\\Delta x$  [m]')
ax.set_title('Where we are vs what BMV needs')
ax.legend(fontsize=6.5, loc='lower right')
ax.axhline(250e-6, color='C3', lw=0.5, ls=':')

fig.suptitle('Stage 11 — Gravity-mediated entanglement: the experimental gap', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig('figures/stage11_bmv.png', dpi=120)
print('wrote results/stage11.json and figures/stage11_bmv.png')

# brief stdout summary
print('--- Bose phases ---')
print('phi(d), phi_closest, phi_far:', ph_mag)
print('dphi_LR, dphi_RL, sum:', ph['dphi_lr'], ph['dphi_rl'], ph['dphi_sum'])
print('C derived/quoted, W derived/quoted:', C_derived, C_quoted, W_derived['W_written'], W_quoted['W_written'])
print('CP ratio (m=1e-14,R=1um,r=200um):', cp_ratio_bose_m)
print('dB/dx T/m:', dBdx_needed, 'g_grav:', grav_accel_closest)
print('tau_coll 300K / 0.15K:', 1/g300, 1/g015)
print('T_env max 10% / 3.5s:', T_env_bose_10pct)
print('mass gaps (orders) Fein, Na, Teb:', gap['mass_orders_vs_Fein'], gap['mass_orders_vs_Pedalino172kDa'], gap['mass_orders_vs_Tebbenjohanns_cooled'])
print('Oppenheim continuous excluded?', opp_original['ultralocal_continuous']['already_excluded'],
      opp_original['ultralocal_continuous']['lower_interferometry'],
      opp_original['ultralocal_continuous']['upper_Cavendish'])
print('discrete gold-nucleus gap orders:', opp_original['ultralocal_discrete_gold_nucleus_mN']['gap_orders'])
print('nonlocal gap orders:', opp_original['nonlocal_continuous']['gap_orders'])
