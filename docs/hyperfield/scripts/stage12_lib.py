"""Stage 12 physics library: BMV/QGEM geometries, witness, decoherence budget.

Every formula is tagged in comments as DERIVED, ASSUMED or CITED.
Units are SI throughout.
"""
from __future__ import annotations
from math import pi, sqrt, sin, cos, exp, log10, factorial
from dataclasses import dataclass, asdict, field
import numpy as np
from scipy.special import zeta

# --- constants (CODATA 2018/2022, CITED) ---
G = 6.67430e-11
HBAR = 1.054571817e-34
HPLANCK = 6.62607015e-34
K_B = 1.380649e-23
C = 299792458.0
EPS0 = 8.8541878128e-12
U_DA = 1.66053906660e-27
MU_B = 9.2740100783e-24
E_CHARGE = 1.602176634e-19
G_EARTH = 9.81
G_SPIN = 2.0  # electron g-factor, NV ms=±1 uses ±g μB; we use g=2 as Bose

# diamond (CITED Bose / van de Kamp unless noted)
RHO_DIAMOND = 3510.0
EPS_DIAMOND = 5.7
CHI_DIAMOND = (EPS_DIAMOND - 1) / (EPS_DIAMOND + 2)  # 0.610
# Im((ε-1)/(ε+2)) for IR absorption/emission. Diamond is a weak far-IR absorber.
# VDK say phonon emission is negligible below T_int=4 K. We keep this as an
# ASSUMED conservative value for sensitivity scans; 0 is the optimistic limit.
IM_CHI_DEFAULT = 0.01

M_N2 = 28.0134 * U_DA
M_HE = 4.002602 * U_DA

# 1 e·cm in C m (Schut dipole unit)
E_CM = E_CHARGE * 0.01  # 1.602e-21 C m


def sphere_radius(m, rho=RHO_DIAMOND):
    return (3.0 * m / (4.0 * pi * rho)) ** (1.0 / 3.0)


# ---------------------------------------------------------------------------
# 1. Entangling phases
# ---------------------------------------------------------------------------
def phi_newton(m, t, r):
    """Magnitude G m^2 t / (ħ r). Sign of V=-Gmm/r is handled in combinations."""
    return G * m * m * t / (HBAR * r)


def phases_linear(m, t, d, dx):
    """Linear Bose/VDK geometry.

    Convention (CITED Schut 2024 / VDK): d is the *closest* (inner-inner) separation
    during the hold, Δx is the superposition width of each mass.
    Distances: inner-inner = d, same-side (LL=RR) = d+Δx, outer-outer = d+2Δx.

    Returns φ1, φ2, Φ=φ1+φ2 in the Schut convention:
      φ1 = Gmm t/ħ (1/d - 1/(d+Δx))
      φ2 = Gmm t/ħ (1/(d+2Δx) - 1/(d+Δx))
    so Φ = φ1+φ2 = Gmm t/ħ [1/d + 1/(d+2Δx) - 2/(d+Δx)].
    This equals Bose's Δφ_LR+Δφ_RL up to the overall sign from V=-Gmm/r.
    """
    if d <= 0 or dx <= 0 or d + 2 * dx <= 0:
        return dict(phi1=np.nan, phi2=np.nan, Phi=np.nan, ok=False)
    a = phi_newton(m, t, d)
    b = phi_newton(m, t, d + dx)
    c = phi_newton(m, t, d + 2 * dx)
    phi1 = a - b
    phi2 = c - b
    return dict(phi1=phi1, phi2=phi2, Phi=phi1 + phi2, ok=True,
                distances=dict(inner=d, mid=d + dx, outer=d + 2 * dx))


def phases_linear_centers(m, t, D, dx):
    """Linear geometry with D = centre-to-centre (Bose 2017). Closest = D-Δx."""
    return phases_linear(m, t, D - dx, dx)


def phases_parallel(m, t, d, dx):
    """Parallel geometry (CITED Schut et al.).

    Masses separated by d horizontally; each split by Δx vertically (same direction).
    Same-spin distances = d; opposite-spin = sqrt(d^2+Δx^2).
      φ1 = φ2 = Gmm t/ħ [1/sqrt(d^2+Δx^2) - 1/d]   (negative, cross paths farther)
      Φ = φ1+φ2
    Schut's trap-paper '|ω_ent|=0.01 Hz' is |φ1|/t, not |Φ|/t.
    """
    if d <= 0 or dx < 0:
        return dict(phi1=np.nan, phi2=np.nan, Phi=np.nan, ok=False)
    a = phi_newton(m, t, d)
    b = phi_newton(m, t, sqrt(d * d + dx * dx))
    phi1 = b - a  # 1/sqrt - 1/d, negative
    return dict(phi1=phi1, phi2=phi1, Phi=2 * phi1, ok=True,
                distances=dict(same=d, cross=sqrt(d * d + dx * dx)))


# ---------------------------------------------------------------------------
# 2. Witness, concurrence, shots
# ---------------------------------------------------------------------------
def dephased_ppt_lmin(Phi, gamma_t, gamma_d_t=0.0):
    """Smallest PPT eigenvalue of the QGEM two-qubit state (DERIVED from
    Schut et al., PRR 6, 013199 (2024), appendix eqs. 32-33).

    Φ = φ1+φ2. Local decoherence (which-path scattering) is modelled as the
    standard QGEM damping e^{-γ t} on single-qubit off-diagonals and e^{-2γ t}
    on the two-body coherences. An extra common dephasing e^{-γ_d t} on the
    entangling coherence (path-jitter, GGN, dipole orientation noise) is
    included as a visibility factor on the sin term (ASSUMED, equivalent to
    a Gaussian average over a phase jitter with vis = e^{-γ_d t}).

    λ_min = (1/4)[ 1 - e^{-2γt} - 2 e^{-γt} vis |sin(Φ/2)| ]
    Entangled iff λ_min < 0 (two-qubit PPT = NPT).
    """
    vis = exp(-gamma_d_t)
    e1 = exp(-gamma_t)
    e2 = exp(-2.0 * gamma_t)
    return 0.25 * (1.0 - e2 - 2.0 * e1 * vis * abs(sin(0.5 * Phi)))


def concurrence_dephased(Phi, gamma_t, gamma_d_t=0.0):
    """Wootters concurrence of the same two-qubit family (DERIVED).

    For this X-state, C = max(0, e^{-γt} vis |sin(Φ/2)| - (1-e^{-2γt})/2 )
    which vanishes exactly when λ_min ≥ 0. Check: γ=0, vis=1 => C=|sin(Φ/2)|
    matching Stage 11's pure-state result.
    """
    vis = exp(-gamma_d_t)
    e1 = exp(-gamma_t)
    e2 = exp(-2.0 * gamma_t)
    return max(0.0, e1 * vis * abs(sin(0.5 * Phi)) - 0.5 * (1.0 - e2))


def ppt_negative(Phi, gamma_t, gamma_d_t=0.0):
    return dephased_ppt_lmin(Phi, gamma_t, gamma_d_t) < 0.0


def n_shots_witness_5sigma(Phi, gamma_t, gamma_d_t=0.0):
    """Number of shots for a 5σ detection of λ_min < 0 (DERIVED, conservative).

    We treat λ_min as estimated with a Bernoulli/correlator error of order
    1/sqrt(N) (each Pauli correlator has variance ≤ 1). Need sqrt(N) |λ_min| ≥ 5,
    i.e. N ≥ 25 / λ_min^2, when λ_min < 0. Returns inf if not entangled.
    """
    lam = dephased_ppt_lmin(Phi, gamma_t, gamma_d_t)
    if lam >= 0:
        return np.inf
    return 25.0 / (lam * lam)


def n_shots_phase_5sigma(Phi, n_target=None):
    """CITED Schut, Grinin et al., PRR 5, 043170 (2023): Φ_eff ≥ 5/sqrt(N).
    This is a phase-detection criterion, not a PPT-witness 5σ. Looser.
    Returns N needed so that |Φ| ≥ 5/sqrt(N), or the 5σ phase for a given N.
    """
    if n_target is not None:
        return 5.0 / sqrt(n_target)
    if abs(Phi) <= 0:
        return np.inf
    return 25.0 / (Phi * Phi)


# ---------------------------------------------------------------------------
# 3. Casimir, dipole, patch
# ---------------------------------------------------------------------------
def U_cp_sphere_sphere(R, r, eps=EPS_DIAMOND):
    """Retarded CP energy of two dielectric spheres (SI, CITED Bose/VDK eq. 4)."""
    chi = ((eps - 1) / (eps + 2)) ** 2
    return - (23.0 * HBAR * C / (4.0 * pi)) * (R ** 6) * chi / r ** 7


def r_min_cp_unshielded(rho=RHO_DIAMOND, eps=EPS_DIAMOND, ratio=0.1):
    """Closest approach with |U_CP|/|U_G| = ratio. Mass-independent at fixed density
    (DERIVED; VDK eq. 5 give 157 μm for ratio=0.1, diamond)."""
    chi = (eps - 1) / (eps + 2)
    pref = (1.0 / ratio) * 23.0 * HBAR * C / (4.0 * pi * G)  # VDK's '10' is 1/0.1
    geo = (3.0 / (4.0 * pi * rho) * chi) ** 2
    return (pref * geo) ** (1.0 / 6.0)


def F_cp_sphere_plate(R, z, eps=EPS_DIAMOND, eta_cond=1.0):
    """Retarded Casimir-Polder *force* (attractive) between a dielectric sphere
    and a conducting plate (CITED VDK eq. 13 / Ford 1998):
        F = - (3 ħ c / 2π) ((ε-1)/(ε+2)) R^3 / z^5
    eta_cond ≤ 1 is an ASSUMED finite-conductivity reduction (1 = perfect).
    """
    chi = (eps - 1) / (eps + 2)
    return - eta_cond * (3.0 * HBAR * C / (2.0 * pi)) * chi * (R ** 3) / z ** 5


def F_dipole_plate(p, z, theta=0.0):
    """Image-dipole force, sphere-plate (CITED Schut 2023 eq. 16):
        F = + 1/(4πϵ0) * 3 p^2 / (16 z^4) * [1 + cos^2 θ]
    Attractive for the usual alignment; we return the signed magnitude toward the plate
    (positive = attractive). p in C m, z in m.
    """
    return (1.0 / (4.0 * pi * EPS0)) * (3.0 * p * p / (16.0 * z ** 4)) * (1.0 + cos(theta) ** 2)


def p_schut_default():
    """Permanent + induced dipole used by Schut et al.: p = 10^{-2} e cm (CITED,
    possibly an overestimate for sub-micron diamonds; they say so)."""
    return 1e-2 * E_CM


# ---------------------------------------------------------------------------
# 4. Decoherence rates
# ---------------------------------------------------------------------------
def v_thermal(T, m_gas=M_N2):
    return sqrt(8.0 * K_B * T / (pi * m_gas))


def gamma_air_geometric(P, T, R, m_gas=M_N2):
    """Resolution-limited collisions with σ = π R^2, ⟨v⟩ = sqrt(8kT/πm) (DERIVED)."""
    n = P / (K_B * T)
    return n * pi * R * R * v_thermal(T, m_gas)


def gamma_air_vdk(P, T, R, m_gas=M_N2):
    """Short-wavelength air-molecule rate used in the QGEM papers
    (CITED VDK appendix B6 / Romero-Isart-type):
        Γ = (16π/3) n R^2 sqrt(2π k T / m)
    This is 8π/3 ≈ 8.38 times the geometric π R^2 ⟨v⟩ rate (DERIVED factor).
    """
    n = P / (K_B * T)
    return (16.0 * pi / 3.0) * n * R * R * sqrt(2.0 * pi * K_B * T / m_gas)


def lambda_bb_scatter(T, R, eps=EPS_DIAMOND):
    """Long-wavelength blackbody *scattering* localization (CITED VDK B7).
    Λ_sc has units 1/(s m^2); Γ = Λ (Δx)^2 in the Δx ≪ λ_th regime.
    Dimensionally there is no extra 1/ħ in the prefactor (ħ sits in kT/ħc).
    """
    chi = ((eps - 1) / (eps + 2)) ** 2
    pref = factorial(8) * zeta(9) * 8.0 / (9.0 * pi)
    return pref * C * (R ** 6) * (K_B * T / (HBAR * C)) ** 9 * chi


def lambda_bb_abs_em(T, R, im_chi=IM_CHI_DEFAULT):
    """Absorption (T=T_env) or emission (T=T_int) localization (CITED VDK B7):
        Λ = 16 π^5 c R^3 / 189 * (kT/ħc)^6 * Im((ε-1)/(ε+2))
    """
    pref = 16.0 * pi ** 5 * C / 189.0
    return pref * (R ** 3) * (K_B * T / (HBAR * C)) ** 6 * im_chi


def gamma_bb(T_env, T_int, R, dx, im_chi=IM_CHI_DEFAULT, eps=EPS_DIAMOND):
    """Photon decoherence rate, interpolating long-wavelength and resolution-limited.

    Γ = (Λ_sc+Λ_abs+Λ_em) Δx^2 / (1+(Δx/λ_th)^2)   with λ_th=ħc/kT_env
    so the rate never exceeds the order of the long-wavelength estimate at Δx~λ_th.
    (DERIVED interpolation; the cap is conservative at large Δx.)
    """
    if T_env <= 0 and T_int <= 0:
        return 0.0
    Lam = 0.0
    if T_env > 0:
        Lam += lambda_bb_scatter(T_env, R, eps)
        Lam += lambda_bb_abs_em(T_env, R, im_chi)
    if T_int > 0:
        Lam += lambda_bb_abs_em(T_int, R, im_chi)
    Tref = max(T_env, 1e-6)
    lam = HBAR * C / (K_B * Tref)
    return Lam * dx * dx / (1.0 + (dx / lam) ** 2)


# ---------------------------------------------------------------------------
# 5. Stern-Gerlach, Humpty-Dumpty, drop height
# ---------------------------------------------------------------------------
def accel_sg(m, dBdx, g=G_SPIN):
    """Spin-dependent acceleration a = (g μB / m) ∂xB (CITED Bose eq. 7)."""
    return g * MU_B * dBdx / m


def dBdx_for_split(m, dx, tau_acc, g=G_SPIN):
    """Invert Δx = (1/2) a τ_acc^2  =>  ∂xB (CITED Bose)."""
    if tau_acc <= 0:
        return np.inf
    a = 2.0 * dx / (tau_acc ** 2)
    return a * m / (g * MU_B)


def drop_height(t_tot, g=G_EARTH):
    return 0.5 * g * t_tot ** 2


def wavepacket_sigma(m, f_trap=1e6, T_com=None):
    """Initial COM width. Ground state of a trap at f_trap (Hz), or thermal
    if T_com is given: σ = sqrt(kT / m ω^2) (ASSUMED initial state)."""
    omega = 2.0 * pi * f_trap
    sig_gs = sqrt(HBAR / (2.0 * m * omega))
    if T_com is None:
        return sig_gs, 'ground-state'
    return max(sig_gs, sqrt(K_B * T_com / (m * omega ** 2))), 'thermal'


def humpty_rel_gradient(dx, sigma0):
    """Required relative gradient mismatch δ(∂B)/∂B ≲ 2 σ_0 / Δx (DERIVED).

    Recombination overlap is fixed by the *initial* wavepacket (Humpty-Dumpty;
    CITED Folman et al., Sci. Adv. 7, eabg2879 (2021); Englert-Schwinger-Scully).
    Position mismatch δx ≈ (δa/a) Δx must stay below σ_0.
    """
    if dx <= 0:
        return np.inf
    return 2.0 * sigma0 / dx


# ---------------------------------------------------------------------------
# 6. GGN / acceleration noise (order-of-magnitude from Toroš et al.)
# ---------------------------------------------------------------------------
def gamma_ggn_toros_rule(Phi_target=0.01):
    """CITED Toroš, Marshman, Mazumdar & Bose, PRR 3, 023178 (2021):
    they require Γ_GGN ≪ Φ_eff ~ 0.01 and argue this is met with access
    restrictions ~2 m (people), 10 m (vehicles), 60 m (aircraft), optionally
    underground. We treat residual GGN as a dephasing rate that can be
    designed below 0.01 / 10 = 1e-3 Hz (ASSUMED, following their ≪).
    Returns that residual rate, not a no-go.
    """
    return 1e-3  # Hz, residual after mitigation, ASSUMED


# ---------------------------------------------------------------------------
# 7. NV spin T2
# ---------------------------------------------------------------------------
# CITED Bose 2017 / white paper arXiv:2509.01586:
# electronic NV T2 ~ 1 s below 77 K with dynamical decoupling (target, not a
# demonstrated 1 s in a falling microdiamond). Nuclear T2 is much longer
# (they swap to nuclear spin for the hold).
T2_ELECTRONIC_S = 1.0     # ASSUMED target
T2_NUCLEAR_S = 100.0      # ASSUMED, conservative vs 'very long'


# ---------------------------------------------------------------------------
# 8. A single experimental point
# ---------------------------------------------------------------------------
@dataclass
class Point:
    name: str
    geometry: str          # 'linear' or 'parallel'
    shielded: bool
    trapped: bool
    m: float
    d: float               # closest (linear) or horizontal (parallel) separation
    dx: float
    t_int: float
    tau_acc: float = 0.5
    T_env: float = 1.0
    T_int_temp: float = 1.0
    T_gas: float = 4.0
    P: float = 1e-15
    p_dipole: float = 1e-2 * E_CM
    z_plate: float = None  # sphere-to-plate distance (defaults to d/2 - W/2)
    W_plate: float = 1e-6
    eta_cond: float = 1.0
    im_chi: float = IM_CHI_DEFAULT
    f_trap: float = 1e6
    T_com: float = None
    m_gas: float = M_N2
    notes: str = ''


def plate_gap(pt: Point):
    """Distance from the inner surface of the sphere to the plate face."""
    if pt.z_plate is not None:
        return pt.z_plate
    # Linear shielded: plate in the middle, closest approach d, plate thickness W
    # gap from sphere surface to plate: d/2 - W/2 - 0 (centre-to-plate is d/2)
    # For parallel shielded, d is the horizontal centre-to-centre through the plate
    # so centre-to-plate is d/2, gap = d/2 - W/2 - R? Force formulas use centre-to-plate
    # as z in VDK (x in their eq. 13 is sphere-centre to plate). Use centre-to-plate.
    return 0.5 * pt.d - 0.5 * pt.W_plate


def evaluate(pt: Point, air_model='vdk'):
    """Full budget + witness for one point. Returns a plain dict."""
    R = sphere_radius(pt.m)
    t_tot = pt.t_int + 2.0 * pt.tau_acc
    if pt.geometry == 'linear':
        ph = phases_linear(pt.m, pt.t_int, pt.d, pt.dx)
    elif pt.geometry == 'parallel':
        ph = phases_parallel(pt.m, pt.t_int, pt.d, pt.dx)
    else:
        raise ValueError(pt.geometry)
    Phi = ph['Phi']
    # magnetic gradient needed to reach Δx in tau_acc
    dBdx = dBdx_for_split(pt.m, pt.dx, pt.tau_acc)
    a_mag = accel_sg(pt.m, dBdx) if np.isfinite(dBdx) else np.nan
    # Casimir / dipole / gravity at the relevant closest EM distance
    if pt.shielded:
        z = plate_gap(pt)
        Fcp = abs(F_cp_sphere_plate(R, z, eta_cond=pt.eta_cond)) if z > 0 else np.inf
        Fd = F_dipole_plate(pt.p_dipole, z) if z > 0 else np.inf
        a_cp = Fcp / pt.m if pt.m > 0 else np.inf
        a_dip = Fd / pt.m if pt.m > 0 else np.inf
        # gravity between the two masses (not screened)
        if pt.geometry == 'linear':
            a_grav = G * pt.m / pt.d ** 2
        else:
            a_grav = G * pt.m / pt.d ** 2
        # VDK constraint: a_ca / a_mag ≤ 0.1 (free-fall recombination)
        ratio_cp_mag = a_cp / a_mag if a_mag else np.inf
        ratio_dip_mag = a_dip / a_mag if a_mag else np.inf
        cp_ss_ratio = 0.0  # screened
        r_cp = z
    else:
        # unshielded: sphere-sphere CP vs gravity at closest
        if pt.geometry == 'linear':
            r_cp = pt.d
        else:
            r_cp = pt.d
        Ucp = abs(U_cp_sphere_sphere(R, r_cp))
        Ug = G * pt.m * pt.m / r_cp
        cp_ss_ratio = Ucp / Ug if Ug else np.inf
        z = None
        a_cp = a_dip = ratio_cp_mag = ratio_dip_mag = 0.0
        Fcp = Fd = 0.0
        a_grav = G * pt.m / r_cp ** 2

    # Decoherence
    if air_model == 'vdk':
        g_air = gamma_air_vdk(pt.P, pt.T_gas, R, pt.m_gas)
    else:
        g_air = gamma_air_geometric(pt.P, pt.T_gas, R, pt.m_gas)
    g_bb = gamma_bb(pt.T_env, pt.T_int_temp, R, pt.dx, pt.im_chi)
    g_ggn = gamma_ggn_toros_rule()
    # spin: if hold is on nuclear spin, only the 2*tau_acc window needs electronic T2
    g_spin_elec = 1.0 / T2_ELECTRONIC_S
    g_spin_nuc = 1.0 / T2_NUCLEAR_S
    # path decoherence integrated: air and GGN act for t_tot; blackbody uses a
    # representative Δx (hold value) for t_tot — slightly conservative during split.
    gamma_path = g_air + g_bb  # 1/s, resolution-limited air + BB
    gamma_t = gamma_path * t_tot
    # residual GGN + dipole-orientation / placement jitter as extra dephasing
    # Dipole dephasing: Schut 2023 found worst-case (uncontrolled p=1e-2 e cm)
    # can demand sub-fm placement. We convert a fractional path jitter
    # δz/z into a phase jitter. ASSUMED: with 'good' dipole control, take
    # γ_d from GGN residual only; also report an uncontrolled-dipole flag.
    gamma_d_t = g_ggn * t_tot

    lam = dephased_ppt_lmin(Phi, gamma_t, gamma_d_t)
    C = concurrence_dephased(Phi, gamma_t, gamma_d_t)
    N5 = n_shots_witness_5sigma(Phi, gamma_t, gamma_d_t)
    N5_phase = n_shots_phase_5sigma(Phi)

    sig0_gs, _ = wavepacket_sigma(pt.m, pt.f_trap, None)
    sig0_th, _ = wavepacket_sigma(pt.m, pt.f_trap, 0.001)  # 1 mK thermal, if used
    hd_gs = humpty_rel_gradient(pt.dx, sig0_gs)
    hd_1mK = humpty_rel_gradient(pt.dx, sig0_th)

    hdrop = drop_height(t_tot) if not pt.trapped else 0.0

    # geometric clash
    clash = False
    if pt.geometry == 'linear' and pt.d <= 2 * R:
        clash = True
    if pt.shielded and z is not None and z <= R:
        clash = True
    if pt.geometry == 'parallel' and pt.d <= 2 * R + (pt.W_plate if pt.shielded else 0):
        clash = True

    return dict(
        name=pt.name, geometry=pt.geometry, shielded=pt.shielded, trapped=pt.trapped,
        m=pt.m, R=R, d=pt.d, dx=pt.dx, t_int=pt.t_int, tau_acc=pt.tau_acc, t_tot=t_tot,
        T_env=pt.T_env, T_int=pt.T_int_temp, T_gas=pt.T_gas, P=pt.P,
        phi1=ph['phi1'], phi2=ph['phi2'], Phi=Phi,
        dBdx=dBdx, a_mag=a_mag, drop_height_m=hdrop,
        z_plate=z, a_cp=a_cp, a_dip=a_dip, a_grav=a_grav,
        ratio_cp_mag=ratio_cp_mag, ratio_dip_mag=ratio_dip_mag,
        cp_ss_ratio=cp_ss_ratio,
        gamma_air=g_air, gamma_bb=g_bb, gamma_ggn=g_ggn,
        gamma_path=gamma_path, gamma_t=gamma_t, gamma_d_t=gamma_d_t,
        lambda_min=lam, concurrence=C, entangled=bool(lam < 0 and not clash),
        N_shots_5sig_witness=N5, N_shots_5sig_phase=N5_phase,
        sigma0_gs=sig0_gs, hd_rel_gs=hd_gs, sigma0_1mK=sig0_th, hd_rel_1mK=hd_1mK,
        clash=clash,
        # spin windows
        t_electronic=2.0 * pt.tau_acc, t_nuclear=pt.t_int,
        electronic_ok=bool(2.0 * pt.tau_acc <= T2_ELECTRONIC_S),
        nuclear_ok=bool(pt.t_int <= T2_NUCLEAR_S),
    )


def required_pressure(pt: Point, gamma_air_budget, air_model='vdk'):
    """P such that Γ_air * t_tot = gamma_air_budget (DERIVED)."""
    R = sphere_radius(pt.m)
    t_tot = pt.t_int + 2.0 * pt.tau_acc
    if t_tot <= 0 or gamma_air_budget <= 0:
        return np.inf
    # Γ = k * P, k from the unit-pressure rate
    pt1 = Point(**{**asdict(pt), 'P': 1.0})
    g1 = (gamma_air_vdk if air_model == 'vdk' else gamma_air_geometric)(
        1.0, pt.T_gas, R, pt.m_gas)
    return gamma_air_budget / (g1 * t_tot)
