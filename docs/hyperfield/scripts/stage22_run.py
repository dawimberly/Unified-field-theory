#!/usr/bin/env python3
"""Stage 22 - designing the experiment that closes the non-local Oppenheim window.

Window (Stage 15/17): l_P^2 D2 lower bound from Pedalino 2025 (172 kDa Na, 1.47e-29) vs upper bound from Gisler 2022
(FOM 0.298 m^2 s^-3 -> 6.69e-26); open by 3.66 orders. Units: kg^2 s m^-1 (Stage 15 labels said kg s m^-1; erratum).

Lower bound (Oppenheim et al. Nat Commun 14, 7910 (2023), Eq. 47, published version):
    l_P^2 D2 >= N_lam M^2 / (R_lam lam),  R_lam = (r_N^2 dx)^(1/3), lam = 1/tau, r_N = 1e-15 m    (as Stage 15)
Talbot-Lau (spatial or time domain): L/L_T = xi, L_T = m v d^2/h, flight tau = 2L/v = 2 xi d^2 M / h   (v drops out)
=> closing condition  2 xi M^4 d^(5/3) / (h m_a r_N^(2/3)) >= U_ceiling                                  (DERIVED)

Everything else is a constraint evaluated along that curve. The interference model is Pedalino et al. Methods Eqs. 1-9
(scripts/stage22_tl.py), validated here against their measured V = 0.10 (172 kDa, P = 62/15.2/68 mW).
"""
import json, os, sys
import numpy as np
from math import pi, sqrt, log10, exp
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import stage22_tl as TL

H = TL.H; HBAR = TL.HBAR; C = TL.C; EPS0 = TL.EPS0; U = TL.U
KB = 1.380649e-23; G = 6.6743e-11; GRAV = 9.81; EV = 1.602176634e-19
R_N = 1e-15
S15 = json.load(open(os.path.join(ROOT, 'results/stage15.json')))
CL = S15['closing']['nonlocal_DP_like']
U_CEIL = CL['from_below_coherence']['need_lower_at_least']        # 6.69e-26  (Gisler 2022)
L_FLOOR = CL['from_below_coherence']['current_lower']             # 1.47e-29  (Pedalino 2025)
FOM_GISLER = 0.298; FOM_NEED = CL['from_above_force_noise']['FOM_needed']
K_NL = R_N ** 3 / G ** 2
M_A = dict(Na=22.98976928 * U, Au=196.966570 * U, Si=28.0855 * U, SiO2=20.03 * U)
OMEGA_E = 7.2921e-5; LAT_VIENNA = 48.2

def nl_bound(M, m_a, tau, dx):
    return (M / m_a) * M ** 2 * tau / (R_N ** 2 * dx) ** (1 / 3)

def M_close(xi, d, m_a, order_factor=2.0):
    """tau = order_factor * xi * d^2 M / h (2 for symmetric TL: tau = 2L/v)."""
    return (U_CEIL * h_ma_r(m_a, d) / (order_factor * xi * d ** 2)) ** 0.25

def h_ma_r(m_a, d):
    return H * m_a * (R_N ** 2 * d) ** (1 / 3)

def orders_left(val):
    return max(0.0, log10(U_CEIL / val))

# ---------------------------------------------------------------- 0. validation of the interferometer model
M172 = 172e3 * U
Vq, S0q, _ = TL.visibility(M172, 160, 10, 0.983, 133e-9, P1=0.062, P2=0.0152, P3=0.068, dm_frac=0.135, quantum=True)
Vc, _, _ = TL.visibility(M172, 160, 10, 0.983, 133e-9, P1=0.062, P2=0.0152, P3=0.068, dm_frac=0.135, quantum=False)
validation = dict(V_quantum_x0p78=0.78 * Vq, V_classical_x0p78=0.78 * Vc, measured='0.10 +- 0.01 (and 0.08 +- 0.01)',
                  transmission_S0=S0q, inputs='m 143-197 kDa uniform, v 160 +- 10 m/s, P1/P2/P3 = 62/15.2/68 mW, '
                  'vertical waist w_y = 600 um (CITED H = 575-620 um), global factor 0.78 (CITED)')
floor_check = nl_bound(M172, M_A['Na'], 2 * 0.983 / 160, 133e-9)

# ---------------------------------------------------------------- 1. closing curve
closing_curve = {}
for lab, d, sp, of in [('MUSCLE_266nm_Na', 133e-9, 'Na', 2.0), ('355nm_Na', 177.5e-9, 'Na', 2.0),
                       ('MUSCLE_266nm_Au', 133e-9, 'Au', 2.0), ('LUMI_d266nm_Na', 266e-9, 'Na', 2.0),
                       ('OTIMA_157nm_Au_2ndorder', 78.8e-9, 'Au', 2.0)]:
    rows = []
    for xi in (0.5, 0.8, 1.0, 1.5, 2.0):
        M = M_close(xi, d, M_A[sp], of)
        tau = of * xi * d ** 2 * M / H
        rows.append(dict(xi=xi, M_kDa=M / U / 1e3, tau_s=tau, T_arm_s=tau / 2,
                         v_for_L_0p983_m=0.983 / (tau / 2), L_arm_for_v140_m=140 * tau / 2,
                         sag_mm=1e3 * 0.5 * GRAV * (tau / 2) ** 2, n_atoms=M / M_A[sp]))
    closing_curve[lab] = dict(d_nm=d * 1e9, species=sp, rows=rows)
# fixed-L MUSCLE geometry: closing mass vs velocity and the quantum (xi >= xi_min) line
fixedL = []
for v in (10, 15, 20, 25, 30, 40, 60, 100, 160):
    L = 0.983; Rl = (R_N ** 2 * 133e-9) ** (1 / 3)
    Mc = (U_CEIL * M_A['Na'] * Rl * v / (2 * L)) ** (1 / 3)
    xi = L * H / (Mc * v * (133e-9) ** 2)
    fixedL.append(dict(v=v, M_close_kDa=Mc / U / 1e3, xi_at_M_close=xi))

# ---------------------------------------------------------------- 2. physical constraints along the curve
def drude_eps(omega, a, gamma_size_A=1.0):
    wp = 5.71 * EV / HBAR; gam = 0.0276 * EV / HBAR + gamma_size_A * 1.07e6 / a   # CITED Na Drude; A=1 ASSUMED
    return 1 - wp ** 2 / (omega ** 2 + 1j * gam * omega)

def thermal_exponent(a, T, d, tau, mode='emit'):
    """Long-wavelength (k d << 1) decoherence exponent on the first fringe harmonic of a symmetric TL:
    E = (d^2 tau / 18) * Int R(omega) k^2 domega, R = sigma(omega) omega^2 nbar / (pi^2 c^2)  (Kirchhoff).
    sigma = electric (4 pi k a^3 Im[(eps-1)/(eps+2)]) + magnetic ((2 pi/15) k^3 a^5 Im eps) dipole absorption
    for 'emit'/'absorb'; Rayleigh (8 pi/3) k^4 a^6 |(eps-1)/(eps+2)|^2 x 2 for 'scatter'."""
    om = np.logspace(11, 15.5, 4000)
    k = om / C; eps = drude_eps(om, a)
    chi = (eps - 1) / (eps + 2)
    x = HBAR * om / (KB * T); nbar = 1 / np.expm1(np.clip(x, 1e-9, 700))
    if mode == 'scatter':
        sig = 2 * (8 * pi / 3) * k ** 4 * a ** 6 * np.abs(chi) ** 2
    else:
        sig = 4 * pi * k * a ** 3 * np.imag(chi) + (2 * pi / 15) * k ** 3 * a ** 5 * np.imag(eps)
    R = sig * om ** 2 * nbar / (pi ** 2 * C ** 2)
    return (d ** 2 * tau / 18) * np.trapezoid(R * k ** 2, om), np.trapezoid(R, om)

def c6_cluster_gas(a, alpha_g, E_g_eV):
    """C6 = (3 hbar/2) a^3 alpha_g wM wg/(wM+wg) (single-oscillator Casimir-Polder; polarizability volumes) [J m^6]."""
    wM = 5.71 / sqrt(3) * EV / HBAR; wg = E_g_eV * EV / HBAR
    return 1.5 * HBAR * a ** 3 * alpha_g * wM * wg / (wM + wg) * 4 * pi * EPS0 / (4 * pi * EPS0)

GASES = dict(H2=dict(m=2.016 * U, alpha=0.805e-30, E=15.4, r=0.14e-9), N2=dict(m=28.0 * U, alpha=1.74e-30, E=15.6, r=0.18e-9),
             He=dict(m=4.0026 * U, alpha=0.205e-30, E=24.6, r=0.13e-9))

def collision(a, tau, gas='H2', T=300.0, budget=0.1):
    gz = GASES[gas]; vbar = sqrt(8 * KB * T / (pi * gz['m']))
    c6 = c6_cluster_gas(a, gz['alpha'], gz['E'])
    s_vdw = 8.083 * (c6 / (HBAR * vbar)) ** 0.4          # total cross section, r^-6 potential (Massey-Mohr)
    s_geo = pi * (a + gz['r']) ** 2
    s = max(s_vdw, s_geo)
    n_req = budget / (s * vbar * tau)
    return dict(C6_J_m6=c6, sigma_vdW_m2=s_vdw, sigma_geo_m2=s_geo, sigma_used_m2=s, v_gas=vbar,
                n_req_m3=n_req, p_req_mbar=n_req * KB * T / 100, J_req_m2s=n_req * vbar / 4,
                exponent_at_MUSCLE_9em9_mbar=(9e-7 / (KB * T)) * s * vbar * tau)

def material_grating_rpi(a, v, b, f=0.6):
    """Point-dipole vdW (Lifshitz, single Mie oscillator): C3 = hbar a^3 wM f / 8; phase pi at r_pi=(C3 b/(pi hbar v))^(1/3)."""
    wM = 5.71 / sqrt(3) * EV / HBAR
    C3 = HBAR * a ** 3 * wM * f / 8
    return C3, (C3 * b / (pi * HBAR * v)) ** (1 / 3)

def cluster_radius(M, sp='Na'):
    rs = dict(Na=2.08e-10, Au=1.59e-10)[sp]          # Wigner-Seitz radii (CITED textbook values)
    return rs * (M / M_A[sp]) ** (1 / 3)

def gravity_rotation(T_arm, v, d, a_y=0.6e-3, lat=LAT_VIENNA):
    s = 0.5 * GRAV * T_arm ** 2
    Om_v = OMEGA_E * np.sin(np.radians(lat))
    dvv_accept = a_y / (4 * s)                         # sagitta spread 2 s dv/v <= a_y/2
    theta_comp = Om_v * v / GRAV                       # tilt that cancels the 1st-order dv dispersion of Coriolis
    L = v * T_arm
    def resid(dvv): return Om_v * L ** 2 / v * dvv ** 2       # 2nd-order Coriolis residual [m]
    def tilt_tol(dvv): return d * v ** 2 / (20 * GRAV * L ** 2 * dvv)
    vert_kvec_dvv = d / (8 * GRAV * T_arm ** 2)        # if grating vector were vertical
    return dict(sag_m=s, dv_over_v_max_vertical_accept=dvv_accept, theta_comp_rad=theta_comp,
                coriolis_shift_uncomp_m=2 * Om_v * v * T_arm ** 2, coriolis_resid_at_1pct_m=resid(0.01),
                tilt_tolerance_rad_at_1pct=tilt_tol(0.01), dv_over_v_max_if_kvec_vertical=vert_kvec_dvv,
                a_y_m=a_y)

# ---------------------------------------------------------------- 3. design scan (ion route + neutral long baseline)
D = 133e-9; SP = 'Na'
B_DARK = 15.0                      # CITED lower end of 15-100 cps dark counts (Pedalino Methods)
ALPHA_REL = 0.5 / 4.5              # CITED alpha_266 uncertainty
SIG_REL = 0.2                      # ASSUMED sigma_ion uncertainty after x8 mass extrapolation of the CITED fit
SRC_ION_RATE = 10e-12 / EV         # 10 pA mass-selected (CITED 'tens of pA', Palmer/Pratontep 2005 via Kemper 2006)
ETA_NEUT = 0.3; ETA_DET = 0.1; T_TRANS = 10.0; A_COLL = 0.5e-3     # ASSUMED
DM = 0.04                          # CITED M/dM ~ 25 (lateral-TOF mass filter)

def design_point(xi, L, route='ion', dvv=0.01, n_grid=(0.5, 1, 2, 3, 4, 6), p_grid=(0.5, 1, 1.5, 2, 3, 4, 5, 6, 8),
                 margin=1.0, T_int=10.0, gas='H2', p_mbar=1e-10):
    M = M_close(xi, D, M_A[SP]) * margin ** 0.25
    T_arm = xi * D ** 2 * M / H; v = L / T_arm; tau = 2 * T_arm; a = cluster_radius(M)
    dvv_req = dvv; dvv = min(dvv, gravity_rotation(T_arm, v, D)['dv_over_v_max_vertical_accept'])
    th = {m: thermal_exponent(a, T, D, tau, m)[0] for m, T in (('emit', T_int), ('absorb', 300.0), ('scatter', 300.0))}
    col = collision(a, tau, gas)
    E_col = (p_mbar * 100 / (KB * 300)) * col['sigma_used_m2'] * col['v_gas'] * tau
    deco = exp(-sum(th.values()) - E_col)
    best = None
    for n0 in n_grid:
        for ph in p_grid:
            kw = dict(m_kg=M, v_mean=v, dv_sigma=dvv * v, L=L, d=D, n0_13=n0, phi2=-ph, dm_frac=DM)
            Vq, S0, P = TL.visibility(quantum=True, **kw)
            Vc, _, _ = TL.visibility(quantum=False, **kw)
            VQ = [Vq]; VC = [Vc]
            for asc in (1 - ALPHA_REL, 1 + ALPHA_REL):
                for ssc in (1 - SIG_REL, 1 + SIG_REL):
                    VQ.append(TL.visibility(quantum=True, alpha_scale=asc, sigma_scale=ssc, **kw)[0])
                    VC.append(TL.visibility(quantum=False, alpha_scale=asc, sigma_scale=ssc, **kw)[0])
            f = 0.78 * deco
            dV = f * abs(Vq - Vc)
            sys_sep = f * min(abs(x - y) for x in VQ for y in VC)   # worst case under alpha +-11% (CITED) x sigma_ion +-20% (ASSUMED)
            if route == 'ion':
                svt = sqrt(KB * T_TRANS / M); fax = min(1.0, (A_COLL / (2 * L)) * v / (sqrt(2 * pi) * svt))
                s_in = SRC_ION_RATE * ETA_NEUT * fax ** 2 * ETA_DET
            else:
                s_in = route                       # numeric rate of clusters entering G1 (counted) [1/s]
            s = s_in * S0 * exp(-E_col)
            T5 = 50 * (s + B_DARK) / (s ** 2 * sys_sep ** 2) if sys_sep > 0 else float('inf')
            rec = dict(n0_13=n0, phi0_2=ph, P1_P3_W=P[0], P2_W=P[1], Vq=f * Vq, Vc=f * Vc, dV=dV, dV_worst_alpha=sys_sep,
                       S0=S0, rate_s=s, T5sigma_s=T5)
            if best is None or T5 < best['T5sigma_s']:
                best = rec
    gr = gravity_rotation(T_arm, v, D)
    mg = material_grating_rpi(a, v, 10e-9)
    if route == 'ion':
        svt = sqrt(KB * T_TRANS / M); fax = min(1.0, (A_COLL / (2 * L)) * v / (sqrt(2 * pi) * svt))
        best['f_geo'] = fax ** 2; best['s_in'] = SRC_ION_RATE * ETA_NEUT * fax ** 2 * ETA_DET
    return dict(xi=xi, L_arm_m=L, M_kDa=M / U / 1e3, n_atoms=M / M_A[SP], radius_nm=a * 1e9, v_m_s=v, T_arm_s=T_arm,
                tau_s=tau, bound=nl_bound(M, M_A[SP], tau, D), dv_over_v=dvv, dv_over_v_requested=dvv_req,
                coriolis_resid_over_d=gr['coriolis_resid_at_1pct_m'] * (dvv / 0.01) ** 2 / D, thermal_exponents=th,
                collision=col, collision_exponent=E_col, p_mbar=p_mbar, gravity=gr,
                material_grating_10nm=dict(C3_J_m3=mg[0], r_pi_nm=mg[1] * 1e9, half_period_nm=D / 2 * 1e9),
                kinetic_energy_eV=0.5 * M * v ** 2 / EV, stray_potential_tol_V=M * v ** 2 * dvv / EV, best=best)

designs = {}
for xi in (0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5):
    designs[f'ion_L0.983_xi{xi}'] = design_point(xi, 0.983, 'ion')
designs['ion_L0.983_xi0.6_margin10'] = design_point(0.6, 0.983, 'ion', margin=10.0)
designs['ion_L0.983_xi0.6_T300K_int'] = design_point(0.6, 0.983, 'ion', T_int=300.0)
designs['ion_L0.983_xi0.6_p9e-9'] = design_point(0.6, 0.983, 'ion', p_mbar=9e-9)

# neutral long-baseline (no slowing): v = 140 m/s CITED; clusters at G1 per second ASSUMED from Pedalino inference
# Inferred MUSCLE count rate: V = 0.10 +- 0.01 => N_s >= 2/0.01^2 = 2e4 counts per scan; scan ~ 18 steps x 4 s (CITED step
# 15 nm, dwell <= 4 s; 2 periods ASSUMED) => ~2.8e2 cps detected for the 143-197 kDa window at S0 ~ 0.12.
INFER_172 = 2e4 / 72.0
neutral = {}
for xi in (0.8, 1.0):
    T_arm = xi * D ** 2 * M_close(xi, D, M_A[SP]) / H
    L = 140 * T_arm
    for frac_label, frac in (('flux_ratio_1e-2', 1e-2), ('flux_ratio_1e-3', 1e-3)):
        # entering rate scaled from 172 kDa inference: / S0(172) x mass-window ratio (ASSUMED) x velocity selection
        sel = 0.0125 / 0.0714 * 0.8                # Delta v/v accepted 1.25% of a 7.1% (10/140) Gaussian sigma (DERIVED, approx)
        s_in = INFER_172 / validation['transmission_S0'] * frac * sel
        neutral[f'xi{xi}_{frac_label}'] = design_point(xi, L, s_in, dvv=0.0125)
        neutral[f'xi{xi}_{frac_label}']['s_in_assumed'] = s_in

# ---------------------------------------------------------------- 4. route comparison (same formula, same r_N)
def route(name, M_amu, sp, tau, dx, status, flag=None, note=''):
    M = M_amu * U; b = nl_bound(M, M_A[sp], tau, dx)
    return dict(name=name, M_amu=M_amu, species=sp, tau_s=tau, dx_m=dx, bound=b, log10_bound=log10(b),
                orders_left=log10(U_CEIL / b), status=status, flag=flag, note=note)

routes = [
    route('Pedalino 2025 (floor)', 172e3, 'Na', 2 * 0.983 / 160, 133e-9, 'demonstrated'),
    route('MUSCLE 1 MDa @25 m/s, tau=2L/v (Stage 15 convention)', 1e6, 'Na', 2 * 0.983 / 25, 133e-9, 'proposed (Pedalino Fig. 3c)'),
    route('MUSCLE 1 MDa @25 m/s, tau=L/v (as Stage 17 used)', 1e6, 'Na', 0.983 / 25, 133e-9, 'proposed', flag='Stage 17 convention slip'),
    route('LUMI-type d=266 nm, 1 MDa Na at xi=1 (needs material G1/G3)', 1e6, 'Na', 2 * 1.0 * (266e-9) ** 2 * 1e6 * U / H, 266e-9, 'proposed',
          flag='material gratings excluded for metal clusters (vdW)'),
    route('OTIMA Au 1e6 amu, 2nd Talbot order (tau = 4 T_T = 60 ms)', 1e6, 'Au', 4 * 1e6 * U * (78.8e-9) ** 2 / H, 78.8e-9, 'proposed (Nimmrichter 2011)'),
    route('Bateman/Ulbricht 1e6 amu Si, release->screen 286 ms, dx 150 nm', 1e6, 'Si', 0.286, 150e-9, 'proposed (2014), not realized'),
    route('Bateman/Ulbricht, grating->screen only 126 ms', 1e6, 'Si', 0.126, 150e-9, 'proposed', note='conservative tau'),
    route('Neumeier/Romero-Isart 1e8 amu, 1 ms, 5 nm', 1e8, 'SiO2', 1e-3, 5e-9, 'proposed', flag='dx << particle size (convention-sensitive)'),
    route('Rossi 2025 delocalization 73 pm, 26.5 us', 7.2266e8, 'SiO2', 26.5e-6, 73e-12, 'demonstrated', flag='dx << particle size (convention-sensitive)'),
]
# closing masses for the alternatives
alt_close = dict(
    OTIMA_Au_2nd_order_kDa=M_close(2.0, 78.8e-9, M_A['Au'], 2.0) / U / 1e3,
    Bateman_Si_tau_286ms_amu=((U_CEIL * M_A['Si'] * (R_N ** 2 * 150e-9) ** (1 / 3) / 0.286) ** (1 / 3)) / U,
)

# ---------------------------------------------------------------- 5. ceiling side
def fom_thermal(T, f0, Q, m_a):
    return 4 * KB * T * 2 * pi * f0 / (Q * m_a)
ceiling = dict(
    FOM_needed=FOM_NEED, FOM_Gisler=FOM_GISLER, improvement_orders=log10(FOM_GISLER / FOM_NEED),
    required_T_omega_over_Q_K_per_s=dict(SiN_20u=FOM_NEED * M_A['SiO2'] / (4 * KB), NdFeB_63p6u=FOM_NEED * 63.6 * U / (4 * KB),
                                        Au=FOM_NEED * M_A['Au'] / (4 * KB), W=FOM_NEED * 183.84 * U / (4 * KB)),
    candidates=[
        dict(name='LISA Pathfinder 2018 (differential, L1)', FOM=1.78e-5, status='demonstrated', counts='excluded by Stage 15 rule',
             orders_beyond_needed=log10(FOM_NEED / 1.78e-5)),
        dict(name='Asenbaum 2017 atom gradiometer (differential)', FOM=2.41e-11, status='demonstrated', counts='excluded by Stage 15 rule',
             orders_beyond_needed=log10(FOM_NEED / 2.41e-11)),
        dict(name='Amaral 2025 0.43 mg maglev, measured', FOM=0.881, status='demonstrated', counts='absolute',
             orders_short=log10(0.881 / FOM_NEED)),
        dict(name='Amaral 2025 maglev THERMAL floor at T=20 mK (ASSUMED), Q=9.3e6, 26.7 Hz (CITED)',
             FOM=fom_thermal(0.020, 26.7, 9.3e6, 63.6 * U), status='not reached (vibration-limited)', counts='absolute'),
        dict(name='same maglev, T=20 mK, Q=1e8 (ASSUMED)', FOM=fom_thermal(0.020, 26.7, 1e8, 63.6 * U), status='projection', counts='absolute'),
        dict(name='room-T torsion pendulum, f0=1 mHz, Q=5e5, Al (ASSUMED class values)', FOM=fom_thermal(300, 1e-3, 5e5, 26.98 * U),
             status='thermal floor, class estimate', counts='absolute'),
        dict(name='cryogenic torsion pendulum 4 K, f0=1 mHz, Q=1e7 (ASSUMED)', FOM=fom_thermal(4, 1e-3, 1e7, 26.98 * U),
             status='projection', counts='absolute'),
        dict(name='Stage 15 projection proj_thermal_maglev: FOM = Gisler/1e4 (ASSUMED)', FOM=2.98e-5,
             status='projection', counts='absolute'),
    ])
for cnd in ceiling['candidates']:
    cnd['closes'] = cnd['FOM'] <= FOM_NEED
    cnd['margin_orders'] = log10(FOM_NEED / cnd['FOM'])

# ---------------------------------------------------------------- 6. model re-check (Oppenheim Table 1 / Eq. 47)
recheck = dict(
    published_eq47='sigma_a^2 N r_N^3 dT / G^2 >= l_P^2 D2 >= N_lam M_lam^2 / (R_lam lam), R_lam ~ V_lam^(1/3) ~ 1e-12 m; '
                   'lower bound quoted 1e-35 kg^2 s m^-1 (Nat Commun 14, 7910, Methods) - CITED, fetched this stage',
    arxiv_v1_eq47='arXiv:2203.01982v1 reads l_P^2 D2 >= M_lam^2/(R_lam lam) with R_lam ~ 1e-9 m and lower bound 1e-40 '
                  '(no N_lam factor): superseded by the journal version',
    stage15_used='N M^2 tau / R with R = (r_N^2 dx)^(1/3), r_N = 1e-15 m, upper K_NL * FOM, K_NL = r_N^3/G^2: matches the '
                 'published Eq. 47 exactly (DERIVED check: Gerlich-2011 inputs give 2.2e-34 vs Table-1 rounding 1e-35).',
    units_erratum='Stage 15/17 json label l_P^2 D2 in kg s m^-1; correct is kg^2 s m^-1 (published text). Numbers unchanged.',
    reproduce_floor=floor_check,
    newer_work=[
        dict(ref='Oppenheim & Sajjad, arXiv:2605.05375 (v2 Sep 2026)', effect='covariant, dimensionless scalar/tensor couplings; '
             'LVK tensor ceiling ~3.3e-67, LISA Pathfinder scalar; matter-wave floor ~1e-71. Different parametrization: does not '
             'tighten the Table-1 non-local l_P^2 D2 window.'),
        dict(ref='Oppenheim, arXiv:2511.10738 (Nov 2025)', effect='stochastic Klein-Gordon; theory, no new bound on l_P^2 D2.'),
        dict(ref='Hirotani & Matsumura, arXiv:2603.29230 (v2 Jun 2026)', effect='geodesic-deviation strain spectra; bound on '
             'D0^org of the original model from GW detectors; different kernel/parameter.'),
        dict(ref='Angeli et al., arXiv:2501.13030 (v3 Jul 2026)', effect='LOCC-gravity diffusion protocol with mK torsion pendulum; '
             'proposal; a different parameter.'),
        dict(ref='Zenodo 10.5281/zenodo.24089422 (2026, unrefereed)', effect='argues V is set by a free smearing length l, not the '
             'wavepacket; supports treating dx-dependence as convention (our flag), no new number on this window.'),
        dict(ref='hexagonmath.org 2611.00594v1 (Nov 2026; AI-generated per its provenance, not on arXiv under that ID)',
             effect='claims XENONnT excludes Gaussian covariant CQ gravity by >1e8; unrefereed; not used.'),
    ],
    verdict='No 2025-26 paper found that tightens the Table-1 non-local l_P^2 D2 window itself; the 3.66-order gap stands.')

# internal temperature at which emission alone reaches exponent 0.1 (xi = 0.8 design)
d06 = designs['ion_L0.983_xi0.8']; a06 = d06['radius_nm'] * 1e-9   # the cheapest design (checked below); a06 = d06['radius_nm'] * 1e-9
Ts = np.linspace(50, 3000, 300); ex = [thermal_exponent(a06, T, D, d06['tau_s'], 'emit')[0] for T in Ts]
T_int_0p1 = float(np.interp(0.1, ex, Ts)) if max(ex) > 0.1 else None
thermal_scan = dict(T_int_K=[100, 300, 600, 1000, 1500], emit_exponent=[thermal_exponent(a06, T, D, d06['tau_s'], 'emit')[0] for T in (100, 300, 600, 1000, 1500)],
                    T_int_for_exponent_0p1_K=T_int_0p1)
# cheapest point: all constraints pass, minimum statistical 5-sigma time (objective stated before looking)
def passes(dsg):
    b = dsg['best']
    return (dsg['coriolis_resid_over_d'] <= 0.1 and b['dV_worst_alpha'] >= 0.05
            and dsg['collision_exponent'] <= 0.1001 + 0.05 and sum(dsg['thermal_exponents'].values()) <= 0.1)
cands = {k: v for k, v in designs.items() if k.count('_') == 2}
for k, v in cands.items(): v['passes_all'] = passes(v)
ok = {k: v for k, v in cands.items() if v['passes_all']}
cheapest = min(ok, key=lambda k: ok[k]['best']['T5sigma_s']) if ok else None
out = dict(stage=22, cheapest=cheapest, thermal_scan=thermal_scan, U_ceiling=U_CEIL, floor=L_FLOOR, gap_orders=log10(U_CEIL / L_FLOOR), validation=validation,
           closing_curve=closing_curve, fixed_L_MUSCLE=fixedL, designs=designs, neutral_long_baseline=neutral,
           routes=routes, alt_closing_masses=alt_close, ceiling=ceiling, model_recheck=recheck,
           assumptions=dict(B_dark_cps=B_DARK, src_ion_rate_s=SRC_ION_RATE, eta_neut=ETA_NEUT, eta_det=ETA_DET,
                            T_transverse_K=T_TRANS, collimation_aperture_m=A_COLL, dm_frac=DM, a_y_m=0.6e-3,
                            discrimination='5 sigma: dV_worst(alpha +-11%, sigma_ion +-20%) >= 5 sqrt(2 (N_s+N_b))/N_s',
                            decoherence_applied_to_both_models=True, global_factor=0.78))
json.dump(out, open(os.path.join(ROOT, 'results/stage22.json'), 'w'), indent=1, default=float)

print('validation Vq=%.3f Vc=%.3f floor=%.3g' % (validation['V_quantum_x0p78'], validation['V_classical_x0p78'], floor_check))
for k, c in closing_curve.items():
    print(k, [(r['xi'], round(r['M_kDa']), round(r['tau_s'] * 1e3, 1), round(r['v_for_L_0p983_m'], 1), round(r['L_arm_for_v140_m'], 2)) for r in c['rows']])
print('fixedL', [(r['v'], round(r['M_close_kDa']), round(r['xi_at_M_close'], 2)) for r in fixedL])
for k, dsg in {**designs, **neutral}.items():
    b = dsg['best']
    print(k, 'M=%.0f kDa v=%.2f tau=%.1f ms Ek=%.2f eV th=%s Ecol=%.2g p_req=%.2g mbar(H2) sag=%.1f mm dvv_acc=%.3g tilt=%.2g rpi=%.0f nm | n0=%g phi=%g Vq=%.3f Vc=%.3f dVw=%.3f S0=%.3f s=%.3g T5=%.3g s' % (
        dsg['M_kDa'], dsg['v_m_s'], dsg['tau_s'] * 1e3, dsg['kinetic_energy_eV'], {a: '%.1e' % x for a, x in dsg['thermal_exponents'].items()},
        dsg['collision_exponent'], dsg['collision']['p_req_mbar'], dsg['gravity']['sag_m'] * 1e3, dsg['gravity']['dv_over_v_max_vertical_accept'],
        dsg['gravity']['tilt_tolerance_rad_at_1pct'], dsg['material_grating_10nm']['r_pi_nm'], b['n0_13'], b['phi0_2'], b['Vq'], b['Vc'], b['dV_worst_alpha'], b['S0'], b['rate_s'], b['T5sigma_s']))
for r in routes:
    print('%-75s bound=%.3g left=%+.2f %s' % (r['name'], r['bound'], r['orders_left'], r['flag'] or ''))
print(alt_close); print('cheapest', cheapest, thermal_scan); print({k:(v['passes_all'], round(v['dv_over_v'],4), round(v['coriolis_resid_over_d'],3), v['best'].get('f_geo'), v['best'].get('s_in')) for k,v in cands.items()})
for cnd in ceiling['candidates']:
    print('%-90s FOM=%.3g closes=%s margin=%+.2f' % (cnd['name'], cnd['FOM'], cnd['closes'], cnd['margin_orders']))
print(ceiling['required_T_omega_over_Q_K_per_s'])
