#!/usr/bin/env python3
"""Stage 16 — heavier trapped test masses for the gravity-entanglement (BMV/QGEM) test.

Question: the Stage 12 cheapest workable point (diamond, 5.6e-15 kg, d = 11 um,
dx = 3.85 um, 3 s) is 7.3 orders of magnitude heavier than the heaviest object
ever put in a spatial superposition. Do heavier test masses (bigger diamonds,
superconducting Meissner-levitated spheres, silica, dense metal spheres released
in space) close that gap once the full decoherence + thermal budget is applied?

Rules: no fitted parameters. Every input is tagged CITED (with source), ASSUMED
(a design choice, kept identical across candidates unless stated) or DERIVED.
Physics re-uses scripts/stage12_lib.py (phases, PPT witness, gas/BB rates,
Casimir) unchanged; Stage 16 adds:
  * long-wavelength gas decoherence (Romero-Isart 2011) for sub-nm dx,
  * an uncorrelated acceleration-noise channel (generalises Grossardt 2020),
  * conductor (superconductor / metal) polarisabilities,
  * per-material Debye heat capacity and a witness-derived T_int ceiling,
  * a joint (mass, delocalisation, time) state-of-the-art ledger.
Units SI.
"""
from __future__ import annotations
import json, os, sys
from math import pi, sqrt, log10, factorial
import numpy as np
from scipy.special import zeta

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from stage12_lib import (G, HBAR, K_B, C, U_DA, MU_B, G_SPIN,
                         phases_parallel, phases_linear, dephased_ppt_lmin,
                         concurrence_dephased, gamma_air_vdk, T2_NUCLEAR_S,
                         T2_ELECTRONIC_S)

M_HE = 4.002602 * U_DA
N_SHOTS_MAX = 1e4          # ASSUMED: "healthy" 5-sigma PPT witness (Stage 12 cheap point used 5.2e3)
ACC_ALLOC = 0.10           # ASSUMED: fraction of the |Phi|/2 entanglement budget given to acceleration/force noise
TAU_ACC = 0.25             # ASSUMED, Stage 12 cheap point: 0.25 s split + 0.25 s recombine windows
F_TRAP = 100.0             # ASSUMED representative trap frequency (Hz). CITED ranges: diamagnetic chip 62-785 Hz
                           # (Elahi et al. arXiv:2411.02325); SC traps 30-160 Hz (Gutierrez Latorre 2023), <=240 Hz (Hofer 2023)
S_GAP = 5.0e-6 - 0.7259e-6 # ASSUMED, Stage 12: sphere-surface to plate-face gap held at its cheap-point value 4.27 um
W_PLATE = 1.0e-6           # ASSUMED, Stage 12 plate thickness
DBDX_MAX = 1e7             # ASSUMED, Stage 12 cut on the Stern-Gerlach gradient (T/m)
RATIO_CUT = 0.1            # ASSUMED, Stage 12 style 10% cut for Casimir vs trap / vs gravity
GAS_ALLOC = 0.40           # ASSUMED, Stage 12 rule: residual gas may use 40% of what is left of |Phi|/2
P_SOA = 1e-15              # Pa, the XHV already assumed in QGEM designs (Stage 12 ledger; inferred, not gauge-measured)
N_SOA = P_SOA / (K_B * 4.0)  # DERIVED: scored as a number density (1.8e7 m^-3, i.e. 1e-15 Pa at 4 K), because a cryopumped
                             # vacuum fixes n, and P = n k T is not comparable between a 4 K and a 15 mK can

# ---------------------------------------------------------------------------
# Materials (CITED unless marked)
# ---------------------------------------------------------------------------
MATERIALS = dict(
    diamond=dict(rho=3510.0, A=12.011, thetaD=2200.0, conductor=False, eps=5.7,
                 # window-grade CVD tan(delta)=1e-5 -> Im chi = 0.288 tan(delta) (Stage 13, DERIVED/CITED)
                 im_chi=0.288 * 1e-5, Tcap=None,
                 src='rho, eps: Stage 12 (CITED Bose/van de Kamp); theta_D 2200 K: Stage 13; tan(delta)=1e-5 window-grade CVD: Stage 13 (Garin et al.; Diamond Materials spec)'),
    lead=dict(rho=11340.0, A=207.2, thetaD=105.0, conductor=True, eps=None,
              im_chi=1e-6,  # ASSUMED: sub-gap surface-resistance absorption of a superconductor; negligible at 15 mK anyway
              Tc=7.2, Tcap=7.2 / 2,  # ASSUMED: keep T_int <= Tc/2 to stay deep in the Meissner state
              src='rho 11340 kg/m3, Tc 7.2 K: standard handbook; theta_D 105 K: Van Sciver, Helium Cryogenics, Table 2.1 (https://www.jlab.org/sites/default/files/magnet-group/materials/low_temperature_materials_properties.pdf); Hofer 2023 used 90/10 PbSn, rho 10900'),
    niobium=dict(rho=8570.0, A=92.906, thetaD=275.0, conductor=True, eps=None,
                 im_chi=1e-6, Tc=9.2, Tcap=9.2 / 2,
                 src='rho 8570, Tc 9.2 K: Pino et al. QST 3, 025001 (2018) App. A; theta_D 275 K: Van Sciver Table 2.1'),
    silica=dict(rho=2200.0, A=60.08 / 3.0, thetaD=470.0, conductor=False, eps=3.8,
                im_chi=3 * 3.8 * 1e-3 / (3.8 + 2) ** 2,  # ASSUMED tan(delta)=1e-3 at mm-wave for fused silica
                Tcap=None,
                src='rho 2200 and eps 3.8: standard fused-silica values (ASSUMED as such); theta_D ~470 K ASSUMED; tan(delta)=1e-3 ASSUMED'),
    osmium=dict(rho=22590.0, A=190.23, thetaD=500.0, conductor=True, eps=None,
                im_chi=1e-3,  # ASSUMED normal-metal THz absorption, order-of-magnitude
                Tcap=None,
                src='rho 22590: Krisnanda et al. npj QI 6, 12 (2020) (arXiv:1906.08808); theta_D ~500 K ASSUMED; Im chi 1e-3 ASSUMED'),
)


def radius(m, rho):
    return (3.0 * m / (4.0 * pi * rho)) ** (1.0 / 3.0)


def z_zpf(m, f=F_TRAP):
    return sqrt(HBAR / (2.0 * m * 2.0 * pi * f))


# ---------------------------------------------------------------------------
# Decoherence
# ---------------------------------------------------------------------------
def v_mean(T, ma):
    return sqrt(8.0 * K_B * T / (pi * ma))


def gas_lambda_LW(P, T, R, ma=M_HE):
    """CITED Romero-Isart, PRA 84, 052121 (2011), eq. 26 (arXiv:1110.4495):
    Lambda_air = 8 sqrt(2 pi) m_a vbar P R^2 / (3 sqrt 3 hbar^2)."""
    return 8.0 * sqrt(2.0 * pi) * ma * v_mean(T, ma) * P * R * R / (3.0 * sqrt(3.0) * HBAR ** 2)


def gas_rate(P, T, R, dx, ma=M_HE):
    """DERIVED interpolation: saturated (short-wavelength) rate = Stage 12 van de Kamp
    formula; long-wavelength Lambda dx^2 for dx << lambda_eff, with
    lambda_eff^2 = Gamma_sat / Lambda_LW. Gamma = Gamma_sat u/(1+u), u = dx^2/lambda_eff^2."""
    g_sat = gamma_air_vdk(P, T, R, ma)
    lam2 = g_sat / gas_lambda_LW(P, T, R, ma)
    u = dx * dx / lam2
    return g_sat * u / (1.0 + u), g_sat, sqrt(lam2)


def epstein_gamma(P, T, R, m, ma=M_HE, delta=1.0):
    """Epstein drag damping rate (CITED kinetic theory; delta=1 specular .. 1.39 diffuse)."""
    n = P / (K_B * T)
    return (4.0 * pi / 3.0) * delta * R * R * n * ma * v_mean(T, ma) / m


def bb_rate(T_env, T_int, R, dx, mat):
    """Blackbody localisation (CITED Romero-Isart 2011 / van de Kamp 2020 / Pino 2018 eqs. 4-5),
    same interpolation as Stage 12.  Conductors: Re chi_E = 1; scattering gets the
    perfect-conductor factor 1.25 = (10/3)/(8/3) from the electric+magnetic Rayleigh
    cross-section of a perfectly conducting sphere (DERIVED, Jackson)."""
    if mat['conductor']:
        chi2 = 1.25
    else:
        e = mat['eps']
        chi2 = ((e - 1) / (e + 2)) ** 2
    im = mat['im_chi']
    kc = lambda T: K_B * T / (HBAR * C)
    lam_sc = factorial(8) * zeta(9) * 8.0 / (9.0 * pi) * C * R ** 6 * kc(T_env) ** 9 * chi2
    ae = lambda T: 16.0 * pi ** 5 * C * R ** 3 / 189.0 * kc(T) ** 6 * im
    Lam = lam_sc + ae(T_env) + ae(T_int)
    lth = HBAR * C / (K_B * T_env)
    return Lam * dx * dx / (1.0 + (dx / lth) ** 2)


# ---------------------------------------------------------------------------
# Casimir (CITED Stage 12 formulae + PFA)
# ---------------------------------------------------------------------------
def casimir_plate_stiffness(R, s, mat):
    """Spring constant of sphere-plate Casimir attraction, surface gap s.
    R < s: retarded Casimir-Polder F = (3 hbar c/2pi) chi R^3/z^5, z = centre-plate (Stage 12).
    R >= s: PFA, perfect conductors F = pi^3 hbar c R/(360 s^3) (CITED; dielectric scaled by chi).
    Returns k = |dF/dz| (DERIVED)."""
    chi = 1.0 if mat['conductor'] else (mat['eps'] - 1) / (mat['eps'] + 2)
    if R < s:
        z = R + s
        F = (3 * HBAR * C / (2 * pi)) * chi * R ** 3 / z ** 5
        return 5 * F / z
    F = chi * pi ** 3 * HBAR * C * R / (360 * s ** 3)
    return 3 * F / s


def casimir_ss_curvature_ratio(R, D, m):
    """Unshielded sphere-sphere (released, linear geometry). CITED Krisnanda 2020:
    PFA energy E = -f0 (pi^3/1440) hbar c R / g^2, g = D - 2R; f0 = 1 (perfect-conductor
    upper bound, ASSUMED). Ratio of entangling curvatures |E''|/|U_G''| with
    U_G = -G m^2/D (DERIVED) = 3 A D^3 / (G m^2 g^4), A = pi^3 hbar c R / 1440."""
    g = D - 2 * R
    if g <= 0:
        return np.inf
    A = pi ** 3 * HBAR * C * R / 1440.0
    return 3 * A * D ** 3 / (G * m * m * g ** 4)


# ---------------------------------------------------------------------------
# Heat (Debye phonons only; for metals this is a LOWER bound on C, so budgets are conservative)
# ---------------------------------------------------------------------------
def U_debye(T, m, mat):
    N = m / (mat['A'] * U_DA)
    return (3.0 * pi ** 4 / 5.0) * N * K_B * T ** 4 / mat['thetaD'] ** 3


# ---------------------------------------------------------------------------
# Candidates (platform environment = ASSUMED from the cited experiment/proposal)
# ---------------------------------------------------------------------------
CANDIDATES = [
    dict(key='diamond_spin', label='(a) diamond + NV spin, diamagnetic trap, plate',
         mat='diamond', geom='trap_plate', split='spin', T_env=1.0, T_int=1.0, T_gas=4.0, P=1e-15,
         t_list=[0.5, 1.0, 2.5, 5.0], m_range=(1e-16, 1e-9),
         env_src='T_env=T_int=1 K, T_gas=4 K, P=1e-15 Pa: Stage 12 defaults (QGEM-design XHV)'),
    dict(key='lead_meissner', label='(b) Pb sphere, Meissner trap, SC plate',
         mat='lead', geom='trap_plate', split='deloc', T_env=0.015, T_int=0.015, T_gas=0.015, P=1e-15,
         t_list=[0.5, 1.0, 2.5, 5.0], m_range=(1e-16, 1e-6),
         env_src='15 mK: Hofer et al. PRL 131, 043603 (2023) dilution fridge; P=1e-15 Pa (=1e-17 mbar): Pino et al. 2018 (inferred cryo-XHV, Gabrielse)'),
    dict(key='niobium_meissner', label='(b\') Nb sphere, Meissner chip (Pino 2018 env.)',
         mat='niobium', geom='trap_plate', split='deloc', T_env=0.05, T_int=0.05, T_gas=0.05, P=1e-15,
         t_list=[0.5, 1.0, 2.5, 5.0], m_range=(1e-16, 1e-6),
         env_src='T=50 mK, P=1e-17 mbar: Pino et al. QST 3, 025001 (2018) App. A'),
    dict(key='silica_trap', label='(c) silica, magneto-grav./diamagnetic trap, plate',
         mat='silica', geom='trap_plate', split='deloc', T_env=1.0, T_int=1.0, T_gas=4.0, P=1e-15,
         t_list=[0.5, 1.0, 2.5, 5.0], m_range=(1e-16, 1e-9),
         env_src='same environment as (a)'),
    dict(key='osmium_space', label='(c\') Os spheres released in space (Krisnanda/MAQRO-style)',
         mat='osmium', geom='free_unshielded', split='deloc', T_env=20.0, T_int=20.0, T_gas=20.0, P=1e-13,
         t_list=[1.0, 2.5, 10.0, 30.0, 100.0], m_range=(1e-14, 1e-5),
         env_src='T<=20 K, P<=1e-15 mbar (=1e-13 Pa), ~100 s: MAQRO BPS2023 white paper arXiv:2202.01535 requirements'),
    dict(key='osmium_space_cold', label='(c\'\') Os released in space at 2.7 K (Krisnanda optimistic)',
         mat='osmium', geom='free_unshielded', split='deloc', T_env=2.7, T_int=2.7, T_gas=2.7, P=None,
         t_list=[1.0, 2.5, 10.0, 30.0, 100.0], m_range=(1e-14, 1e-5),
         env_src='T = 2.7 K (CMB-limited passive shield) as assumed by Krisnanda et al. 2020 for space; better than the MAQRO 20 K requirement'),
]

# ---------------------------------------------------------------------------
# State-of-the-art ledger (m, dx, t), CITED
# ---------------------------------------------------------------------------
SOA = [
    dict(key='pedalino2025', label='Pedalino 2025 cluster interference', m=2.86e-22, dx=133e-9, t=10e-3, com=True,
         src='Stage 12 ledger: 172 kDa, 133 nm grating, ~10 ms flight'),
    dict(key='rossi2025', label='Rossi 2025 delocalised nanosphere', m=1.2e-18, dx=73e-12, t=2.6e-5, com=True,
         src='Rossi et al., PRL (2025), arXiv:2408.01264: 1.2 fg silica, coherence length 73+-34 pm (x3.5 of 21 pm); '
             't = 2 x pi/(2 Omega_e) + pi/(2 Omega_m) with 18 us period, r = 2.45 -> 26 us (DERIVED)'),
    dict(key='bild2023', label='Bild 2023 16 ug HBAR cat (internal mode)', m=16.2e-9, dx=2.1e-18, t=80e-6, com=False,
         src='Bild et al., Science 380, 274 (2023), arXiv:2211.00449: M_eff 16.2 ug, delocalised over 2.1e-18 m; '
             '~80 us (Physics Today 2023). NOT a centre-of-mass displacement'),
]

# demonstrated acceleration/force noise (CITED / DERIVED)
ACC_SOA = dict(
    lpf=dict(val=1.74e-15, src='LISA Pathfinder differential acceleration noise 1.74 fm s^-2/sqrt(Hz) above 2 mHz, '
             'Armano et al. PRL 120, 061101 (2018), https://sci.esa.int/web/lisa-pathfinder/-/59964-armano-et-al (2 kg masses, space, mHz band)'),
)
# Fuchs et al. Sci. Adv. 10, eadk2949 (2024): 0.43 mg, gamma/2pi = 2.9 uHz, T < 100 mK, force noise 0.5 fN/rtHz (measured)
FUCHS = dict(m=0.43e-6, gamma=2 * pi * 2.9e-6, T=0.1, F_meas=0.5e-15)
FUCHS['a_thermal'] = sqrt(4 * FUCHS['gamma'] * K_B * FUCHS['T'] / FUCHS['m'])  # DERIVED, upper T
FUCHS['a_meas'] = FUCHS['F_meas'] / FUCHS['m']
HOFER = dict(m=5.6e-9, Q=2.6e7, f=240.0, T=0.015)  # CITED Hofer 2023 (Q and f quoted as maxima, not necessarily simultaneous)
HOFER['gamma'] = 2 * pi * HOFER['f'] / HOFER['Q']


# ---------------------------------------------------------------------------
# Geometry and a single evaluation
# ---------------------------------------------------------------------------
def phase_parallel_stable(m, t, d, dx):
    """|Phi| for the Schut parallel geometry, cancellation-free (DERIVED; equals stage12_lib.phases_parallel):
    |Phi| = 2 (G m^2 t / hbar d) [1 - 1/sqrt(1+u)],  u = (dx/d)^2,  1 - 1/sqrt(1+u) = u / (sqrt(1+u)(1+sqrt(1+u)))."""
    u = (dx / d) ** 2
    q = sqrt(1.0 + u)
    return 2.0 * G * m * m * t / (HBAR * d) * u / (q * (1.0 + q))


def phase_linear_stable(m, t, d, dx):
    """|Phi| for the linear geometry, closest inner separation d (DERIVED, exact identity):
    1/d + 1/(d+2dx) - 2/(d+dx) = 2 dx^2 / (d (d+dx)(d+2dx))."""
    return G * m * m * t / HBAR * 2.0 * dx * dx / (d * (d + dx) * (d + 2 * dx))


def geometry(m, cand, mat):
    """Returns centre separation D, closest phase distance, and Casimir diagnostic."""
    R = radius(m, mat['rho'])
    if cand['geom'] == 'trap_plate':
        D = 2.0 * (R + S_GAP) + W_PLATE        # ASSUMED: Stage 12 gap rule, scaled with R
        k_cas = casimir_plate_stiffness(R, S_GAP, mat)
        ratio = k_cas / (m * (2 * pi * F_TRAP) ** 2)
        return dict(R=R, D=D, cas_ratio=ratio, kD=D / R)
    # free, unshielded: smallest D/R on a fixed ladder that passes the Casimir curvature cut (rule, not a fit)
    for kD in np.geomspace(2.05, 200.0, 400):
        D = kD * R
        r = casimir_ss_curvature_ratio(R, D, m)
        if r <= RATIO_CUT:
            return dict(R=R, D=D, cas_ratio=r, kD=kD)
    return dict(R=R, D=np.nan, cas_ratio=np.inf, kD=np.nan)


def evaluate(m, dx, t_int, cand, T_int=None, P=None):
    mat = MATERIALS[cand['mat']]
    g = geometry(m, cand, mat)
    R, D = g['R'], g['D']
    t_tot = t_int + 2 * TAU_ACC
    if not np.isfinite(D) or dx >= 0.5 * D:
        return None
    if cand['geom'] == 'trap_plate':
        Phi = phase_parallel_stable(m, t_int, D, dx)
    else:
        Phi = phase_linear_stable(m, t_int, D - dx, dx)
    T_i = cand['T_int'] if T_int is None else T_int
    g_bb = bb_rate(cand['T_env'], T_i, R, dx, mat)
    gamma_d_t = ACC_ALLOC * Phi / 2.0
    # pressure: Stage 12 rule (DERIVED): gas gets GAS_ALLOC of what is left of |Phi|/2 after BB and noise.
    # Gamma_gas is linear in P and lambda_eff is P-independent, so P_req is explicit.
    g1, _, lam_eff = gas_rate(1.0, cand['T_gas'], R, dx)
    left = Phi / 2.0 - g_bb * t_tot - gamma_d_t
    P_req = GAS_ALLOC * left / (g1 * t_tot) if left > 0 else 0.0
    P_use = P_req if P is None else P
    g_gas = g1 * P_use
    g_sat = gamma_air_vdk(P_use, cand['T_gas'], R, M_HE)
    gamma_t = (g_gas + g_bb) * t_tot
    lam = dephased_ppt_lmin(Phi, gamma_t, gamma_d_t)
    N5 = 25.0 / lam ** 2 if lam < 0 else np.inf
    # acceleration-noise requirement that the ACC_ALLOC allocation implies (DERIVED):
    #   Gamma_acc = m^2 S_a dx^2 / (4 hbar^2)   (one-sided white S_a, uncorrelated between masses)
    #   Gamma_acc t_tot = ACC_ALLOC |Phi|/2  =>  S_a = 2 ACC_ALLOC hbar^2 |Phi| / (m^2 dx^2 t_tot)
    S_a = 2 * ACC_ALLOC * HBAR ** 2 * Phi / (m * m * dx * dx * t_tot)
    gamma_m_max = m * S_a / (4 * K_B * cand['T_env'])     # thermal damping bath at T_env (DERIVED, FDT)
    n_req = P_use / (K_B * cand['T_gas'])
    out = dict(m=m, R=R, D=D, kD=g['kD'], dx=dx, t_int=t_int, t_tot=t_tot, Phi=Phi,
               gamma_gas=g_gas, gamma_gas_sat=g_sat, lambda_gas_eff=lam_eff, gamma_bb=g_bb,
               gamma_t=gamma_t, gamma_d_t=gamma_d_t, lambda_min=lam, N5=N5,
               C=concurrence_dephased(Phi, gamma_t, gamma_d_t),
               cas_ratio=g['cas_ratio'], sqrt_Sa_req=sqrt(S_a), gamma_m_max=gamma_m_max,
               Q_min_100Hz=(2 * pi * F_TRAP / gamma_m_max) if gamma_m_max > 0 else np.inf,
               P_req=P_use, n_req=n_req,
               vac_gap=(log10(N_SOA / n_req) if P_use > 0 else np.inf),
               zpf=z_zpf(m), deloc_factor=dx / (2 * z_zpf(m)))
    ok = (N5 <= N_SHOTS_MAX) and P_use > 0 and (g['cas_ratio'] <= RATIO_CUT) and dx < 0.5 * D
    if cand['split'] == 'spin':
        dBdx = 2 * dx * m / (G_SPIN * MU_B * TAU_ACC ** 2)   # Stage 12 rule
        out['dBdx'] = dBdx
        out['HD_rel'] = 2 * z_zpf(m) / dx                    # Humpty-Dumpty, ground state of F_TRAP
        ok = ok and dBdx <= DBDX_MAX and t_int <= T2_NUCLEAR_S and 2 * TAU_ACC <= T2_ELECTRONIC_S
    out['ok'] = bool(ok)
    return out


DX_GRID = np.geomspace(1e-18, 1e-3, 1501)


def min_dx_point(m, t_int, cand):
    """Smallest superposition that gives a healthy witness under all cuts (rule: smallest dx)."""
    for dx in DX_GRID:
        r = evaluate(m, dx, t_int, cand)
        if r is not None and r['ok']:
            return r
    return None


def diagnose(m, cand):
    """Why no dx works at this mass (first failing cut at the dx that maximises -lambda_min)."""
    best = None
    for t in cand['t_list']:
        for dx in DX_GRID[::10]:
            r = evaluate(m, dx, t, cand)
            if r is None:
                continue
            if best is None or r['lambda_min'] < best['lambda_min']:
                best = r
    if best is None:
        return 'geometry (dx >= D/2 or no Casimir-safe D)'
    why = []
    if best['N5'] > N_SHOTS_MAX: why.append('shots(N5=%.2g)' % best['N5'])
    if best['cas_ratio'] > RATIO_CUT: why.append('casimir')
    if best.get('dBdx', 0) > DBDX_MAX: why.append('gradient')
    if best['dx'] >= 0.5 * best['D']: why.append('dx>D/2')
    if best['P_req'] <= 0: why.append('no gas budget')
    return ','.join(why) or 'other'


def T_int_ceiling(pt, cand):
    """Highest internal temperature that still keeps N5 <= N_SHOTS_MAX, everything else fixed (DERIVED);
    capped by the material limit (Tc/2 for superconductors, ASSUMED) and 300 K."""
    mat = MATERIALS[cand['mat']]
    cap = mat['Tcap'] if mat['Tcap'] else 300.0
    f = lambda T: evaluate(pt['m'], pt['dx'], pt['t_int'], cand, T_int=T, P=pt['P_req'])['N5'] <= N_SHOTS_MAX
    lo = cand['T_int']
    if f(cap):
        return cap, 'material/300K cap'
    hi = cap
    for _ in range(80):
        mid = sqrt(lo * hi)
        if f(mid):
            lo = mid
        else:
            hi = mid
    return lo, 'witness (BB emission)'


def heat_budget(pt, cand, T_hi):
    mat = MATERIALS[cand['mat']]
    Q = U_debye(T_hi, pt['m'], mat) - U_debye(cand['T_int'], pt['m'], mat)
    return Q, Q / pt['t_tot']


def joint_gap(m, dx, t, com_only=True):
    best = None
    for s in SOA:
        if com_only and not s['com']:
            continue
        parts = dict(m=max(0.0, log10(m / s['m'])), dx=max(0.0, log10(dx / s['dx'])), t=max(0.0, log10(t / s['t'])))
        tot = sum(parts.values())
        if best is None or tot < best[0]:
            best = (tot, s['key'], parts)
    return best


# ---------------------------------------------------------------------------
def main():
    res = dict(stage=16, rules=dict(N_SHOTS_MAX=N_SHOTS_MAX, ACC_ALLOC=ACC_ALLOC, TAU_ACC=TAU_ACC, F_TRAP=F_TRAP,
                                    S_GAP=S_GAP, W_PLATE=W_PLATE, DBDX_MAX=DBDX_MAX, RATIO_CUT=RATIO_CUT,
                                    tags='all ASSUMED design rules carried over from Stage 12 or stated here'),
               materials=MATERIALS, soa=SOA, acc_soa=ACC_SOA,
               fuchs=FUCHS, hofer=HOFER, candidates={})

    # ---- cross-checks (DERIVED) ----
    chk = {}
    chk['stable_vs_stage12_parallel'] = phase_parallel_stable(5.6234e-15, 2.5, 11e-6, 3.852e-6) / abs(phases_parallel(5.6234e-15, 2.5, 11e-6, 3.852e-6)['Phi'])
    chk['stable_vs_stage12_linear'] = phase_linear_stable(1e-14, 2.5, 200e-6, 250e-6) / abs(phases_linear(1e-14, 2.5, 200e-6, 250e-6)['Phi'])
    # 1. Grossardt 2020: S_a << 2 G hbar / d^3 (our Gamma_acc with alloc -> 1) at d = 157 um (unshielded diamond)
    d157 = 156.9e-6
    chk['grossardt_limit_157um_fm'] = sqrt(2 * G * HBAR / d157 ** 3) * 1e15
    chk['grossardt_paper_fm'] = 0.07
    chk['plate_limit_11um_fm'] = sqrt(2 * G * HBAR / 11e-6 ** 3) * 1e15
    # 2. Pino 2018: 1/gamma_a = 1.6 s for R = 1 um, He, 50 mK, 1e-15 Pa
    chk['pino_gas_coh_time_s'] = 1.0 / gamma_air_vdk(1e-15, 0.05, 1e-6, M_HE)
    # 3. LW gas Lambda vs Epstein damping x FDT: hbar^2 Lambda / (m gamma kT) should be O(1)
    R_, m_ = 1e-6, 1e-14
    chk['LW_gas_vs_epstein_FDT_ratio'] = (HBAR ** 2 * gas_lambda_LW(1e-15, 4.0, R_)) / (
        m_ * epstein_gamma(1e-15, 4.0, R_, m_) * K_B * 4.0)
    # 4. Krisnanda 2020 picked 100 ug Os at D = 3R; our Casimir curvature ratio there
    R_os = radius(1e-7, 22590.0)
    chk['krisnanda_casimir_ratio_100ug_3R'] = casimir_ss_curvature_ratio(R_os, 3 * R_os, 1e-7)
    chk['krisnanda_paper_ratio'] = 1e-2
    chk['min_mass_for_3R_ratio_0p1_kg'] = sqrt(81 * pi ** 3 * HBAR * C / (1440 * G * 0.1))
    # 5. Fuchs 2024: thermal acceleration noise implied by their damping at T = 100 mK
    chk['fuchs_thermal_acc_fm'] = FUCHS['a_thermal'] * 1e15
    chk['fuchs_measured_acc_fm'] = FUCHS['a_meas'] * 1e15
    chk['hofer_gamma_s'] = HOFER['gamma']
    chk['hofer_thermal_acc_fm'] = sqrt(4 * HOFER['gamma'] * K_B * HOFER['T'] / HOFER['m']) * 1e15
    chk['lambda_eff_He_4K_nm'] = gas_rate(1e-15, 4.0, 1e-6, 1e-9)[2] * 1e9
    chk['lambda_eff_He_15mK_nm'] = gas_rate(1e-15, 0.015, 1e-6, 1e-9)[2] * 1e9
    chk['lambda_th_He_4K_nm'] = 2 * pi * HBAR / sqrt(2 * pi * M_HE * K_B * 4.0) * 1e9
    res['cross_checks'] = chk

    # ---- baseline: Stage 12 cheapest point re-run under Stage 16 rules ----
    cand_a = CANDIDATES[0]
    # (i) exact Stage 12 reproduction with stage12_lib (N2 'air', P = 2.4e-15 Pa, GGN 1e-3 Hz)
    import stage12_lib as s12
    p12 = s12.Point(name='cheap', geometry='parallel', shielded=True, trapped=True, m=5.6234e-15, d=11e-6,
                    dx=3.852e-6, t_int=2.5, tau_acc=0.25, P=2.40e-15, z_plate=5e-6)
    e12 = s12.evaluate(p12)
    res['stage12_reproduction'] = dict(Phi=e12['Phi'], C=e12['concurrence'], N5=e12['N_shots_5sig_witness'],
                                       gamma_air_N2=e12['gamma_air'],
                                       gamma_air_He_same_P=gamma_air_vdk(2.40e-15, 4.0, e12['R'], M_HE))
    # (ii) same point, Stage 16 rules (He residual gas, pressure from the 40% rule, 10% noise allocation)
    base = evaluate(5.6234e-15, 3.852e-6, 2.5, cand_a)
    base_mat = MATERIALS['diamond']
    Tc_b, why_b = T_int_ceiling(base, cand_a)
    Qb, Pb = heat_budget(base, cand_a, 4.0)
    Qbw, Pbw = heat_budget(base, cand_a, Tc_b)
    base.update(P_req=base['P_req'], T_int_ceiling=Tc_b, T_int_ceiling_reason=why_b,
                heat_budget_1to4K_aW=Pb * 1e18, heat_budget_to_ceiling_aW=Pbw * 1e18,
                joint_gap=joint_gap(base['m'], base['dx'], base['t_tot']),
                mass_gap_vs_interference=log10(base['m'] / SOA[0]['m']))
    res['baseline'] = base
    base_joint = base['joint_gap'][0]

    # ---- scans ----
    best_overall = None
    for cand in CANDIDATES:
        mat = MATERIALS[cand['mat']]
        m_grid = np.geomspace(cand['m_range'][0], cand['m_range'][1], 31)
        rows = []
        for m in m_grid:
            per_t = []
            for t in cand['t_list']:
                p = min_dx_point(m, t, cand)
                if p is None:
                    continue
                jg = joint_gap(p['m'], p['dx'], p['t_tot'])
                jg_all = joint_gap(p['m'], p['dx'], p['t_tot'], com_only=False)
                p['joint_gap'] = jg[0]; p['joint_gap_ref'] = jg[1]; p['joint_gap_parts'] = jg[2]
                p['joint_gap_incl_bild'] = jg_all[0]
                p['acc_gap_vs_LPF'] = max(0.0, log10(ACC_SOA['lpf']['val'] / p['sqrt_Sa_req']))
                p['mass_gap_vs_interference'] = log10(p['m'] / SOA[0]['m'])
                p['vac_gap_pos'] = max(0.0, p['vac_gap'])
                p['score'] = p['joint_gap'] + p['acc_gap_vs_LPF'] + p['vac_gap_pos']
                per_t.append(p)
            if not per_t:
                rows.append(dict(m=m, feasible=False, why=diagnose(m, cand)))
                continue
            # best t for this mass = smallest joint gap (rule), ties -> smaller acc gap
            p = min(per_t, key=lambda q: (round(q['score'], 6), q['joint_gap']))
            Tcl, why = T_int_ceiling(p, cand)
            Q, P = heat_budget(p, cand, Tcl)
            p.update(feasible=True, T_int_ceiling=Tcl, T_int_ceiling_reason=why,
                     heat_budget_J=Q, heat_budget_W=P, heat_budget_aW=P * 1e18,
                     heat_budget_vs_stage14=P / 1.8161706e-18,
                     per_t=[dict(t_int=q['t_int'], dx=q['dx'], joint_gap=q['joint_gap'], N5=q['N5'],
                                 sqrt_Sa_req=q['sqrt_Sa_req']) for q in per_t],
                     t25=[dict(dx=q['dx'], joint_gap=q['joint_gap'], sqrt_Sa_req=q['sqrt_Sa_req'])
                          for q in per_t if abs(q['t_int'] - 2.5) < 1e-9])
            rows.append(p)
        feas = [r for r in rows if r.get('feasible')]
        best = min(feas, key=lambda q: (q['score'], q['joint_gap'])) if feas else None
        res['candidates'][cand['key']] = dict(label=cand['label'], material=cand['mat'], env_src=cand['env_src'],
                                              geom=cand['geom'], split=cand['split'], rows=rows, best=best,
                                              T_env=cand['T_env'], P=cand['P'])
        if best is not None:
            score = best['score']
            if best_overall is None or score < best_overall[0]:
                best_overall = (score, cand['key'])
    res['best_overall'] = dict(key=best_overall[1], score_joint_plus_acc=best_overall[0],
                               baseline_joint_gap=base_joint,
                               baseline_score=base_joint + max(0.0, log10(ACC_SOA['lpf']['val'] / base['sqrt_Sa_req'])) + max(0.0, base['vac_gap']))

    # ---- scaling law (DERIVED): for R >> gap, D -> 2R and the requirement is mass-independent in zpf units ----
    def deloc_large_R(rho, Phi=0.5, t=2.5, f=F_TRAP):
        # parallel, small dx: |Phi| = G m^2 t dx^2/(hbar D^3), D^3 = 6 m/(pi rho)
        # dx/(2 zpf) = 0.5 sqrt(12 Phi omega/(pi rho G t))
        return 0.5 * sqrt(12 * Phi * 2 * pi * f / (pi * rho * G * t))
    res['scaling'] = {k: dict(deloc_factor_large_R=deloc_large_R(v['rho']),
                              m_dx2_large_R=6 * 0.5 * HBAR / (pi * v['rho'] * G * 2.5))
                      for k, v in MATERIALS.items()}
    res['scaling']['note'] = ('Large-R limit, parallel geometry, Phi=0.5, t=2.5 s, f=100 Hz: required dx/(2 z_zpf) and m dx^2 are '
                              'mass-independent; only density helps (DERIVED).')

    with open(os.path.join(ROOT, 'results', 'stage16.json'), 'w') as fh:
        json.dump(res, fh, indent=1, default=lambda o: float(o) if isinstance(o, (np.floating,)) else str(o))

    # ---- print summary ----
    print('cross-checks:')
    for k, v in chk.items():
        print('  %-38s %.4g' % (k, v))
    b = base
    print('stage12 reproduction:', res['stage12_reproduction'])
    print('\nbaseline (Stage 12 cheap point, Stage 16 rules): Phi=%.3f N5=%.3g C=%.3f gas=%.3g bb=%.3g sqrtSa=%.3g fm  '
          'P_req=%.3g gamma_m_max=%.3g  Tceil=%.2f K (%s)  heat(1-4K)=%.2f aW, to ceiling=%.3g aW  joint=%.2f via %s %s'
          % (b['Phi'], b['N5'], b['C'], b['gamma_gas'], b['gamma_bb'], b['sqrt_Sa_req'] * 1e15, b['P_req'], b['gamma_m_max'],
             b['T_int_ceiling'], b['T_int_ceiling_reason'], b['heat_budget_1to4K_aW'], b['heat_budget_to_ceiling_aW'],
             b['joint_gap'][0], b['joint_gap'][1], b['joint_gap'][2]))
    for key, c in res['candidates'].items():
        print('\n==', c['label'])
        for r in c['rows']:
            if not r.get('feasible'):
                print('  m=%.2e  infeasible: %s' % (r['m'], r['why']))
                continue
            print('  m=%.2e R=%.2e D=%.2e t=%.1f dx=%.3e Phi=%.3f N5=%.2g gas=%.2g bb=%.2g joint=%.2f[%s] mgap=%.2f '
                  'sqrtSa=%.3g fm accgap=%.2f P=%.2g vac=%.2f gmax=%.2g deloc=%.3g Tceil=%.3g heat=%.3g aW (x%.2g) score=%.2f%s'
                  % (r['m'], r['R'], r['D'], r['t_int'], r['dx'], r['Phi'], r['N5'], r['gamma_gas'], r['gamma_bb'],
                     r['joint_gap'], r['joint_gap_ref'], r['mass_gap_vs_interference'], r['sqrt_Sa_req'] * 1e15,
                     r['acc_gap_vs_LPF'], r['P_req'], r['vac_gap'], r['gamma_m_max'], r['deloc_factor'], r['T_int_ceiling'], r['heat_budget_aW'],
                     r['heat_budget_vs_stage14'], r['score'], '  dBdx=%.2g' % r['dBdx'] if 'dBdx' in r else ''))
        bb = c['best']
        if bb:
            print('  BEST: m=%.3e dx=%.3e t=%.1f joint=%.2f acc_gap=%.2f' % (bb['m'], bb['dx'], bb['t_int'], bb['joint_gap'], bb['acc_gap_vs_LPF']))
    print('\nbest overall:', res['best_overall'])
    print('scaling:', {k: (round(v['deloc_factor_large_R']), '%.2e' % v['m_dx2_large_R']) for k, v in res['scaling'].items() if k != 'note'})


if __name__ == '__main__':
    main()
