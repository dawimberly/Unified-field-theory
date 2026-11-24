"""Stage 18 - common-mode vs differential acceleration noise for the trapped BMV/QGEM point.

Question (Stage 16 follow-up): Stage 16 needs sqrt(S_a) <= 0.90 fm s^-2/rtHz per mass at the cheap
point (diamond, m = 5.6e-15 kg, dx = 3.85 um, D = 11 um, t_int = 2.5 s, 1 K).  Both masses share one
chip / cryostat / platform.  How much of the ground noise is common-mode, how much common-mode
rejection (CMRR) is realistic, and does that close the gap on the ground?

Labels: DERIVED (computed here), CITED (number from a verified source, URL given), ASSUMED (our choice).
No parameter is tuned to make anything pass.
"""
from __future__ import annotations
import json, os, sys
from math import pi, sqrt, log10
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
from stage12_lib import G, HBAR, K_B, dephased_ppt_lmin, gamma_air_vdk          # noqa: E402
from stage16_run import (CANDIDATES, evaluate, epstein_gamma, M_HE, TAU_ACC,     # noqa: E402
                         ACC_ALLOC, S_GAP, W_PLATE, MATERIALS)

G0 = 9.81
MU0 = 4e-7 * pi
E_CHG = 1.602176634e-19
E_CM = E_CHG * 1e-2

# ---------------------------------------------------------------------------
# 0. Baseline point (Stage 16 rules) and its which-path lever arm
# ---------------------------------------------------------------------------
CAND = CANDIDATES[0]                       # (a) diamond + NV spin, diamagnetic trap, plate, 1 K
BASE = evaluate(5.6234e-15, 3.852e-6, 2.5, CAND)
M, DX, D, R, PHI = BASE['m'], BASE['dx'], BASE['D'], BASE['R'], BASE['Phi']
T_INT = BASE['t_int']
KAPPA = M * DX / HBAR                      # rad per (m/s) of velocity change; DERIVED
Z_C = S_GAP + R                            # plate face to sphere centre (5.0 um), Stage 12 geometry
GT = BASE['gamma_t']                       # gas + blackbody exponent (local, uncorrelated), Stage 16

# Pino et al. 2018 point quoted by Stage 19 (needs 0.68 fm s^-2/rtHz at 3.3e-14 kg): used only as a
# second requirement line (CITED from Stage 19 of this repo; not recomputed here).
PINO_REQ = 0.68e-15


# ---------------------------------------------------------------------------
# 1. Exact two-qubit model with correlated Gaussian which-path phase noise (DERIVED)
# ---------------------------------------------------------------------------
def rho_noisy(Phi, gamma_t, sig1, sig2, c, s1=1, s2=1):
    """QGEM parallel-geometry state, basis |ab>, a,b in {L=0,R=1}; clean phases [0,-Phi/2,-Phi/2,0].
    Local which-path decoherence (gas, BB): e^{-gamma_t} per flipped qubit index (Stage 12/16 convention).
    Acceleration noise: branch R of mass k picks up theta_k = s_k * kappa * int a_k w dt, jointly Gaussian,
    sd sig_k, correlation c.  Exact average: E exp(i(da*th1+db*th2)) = exp(-(da^2 s1^2 + db^2 s2^2 + 2 c da db
    s1 s2 sig1 sig2... )/2) written with variances below."""
    ph = np.array([0.0, -Phi / 2, -Phi / 2, 0.0])
    psi = 0.5 * np.exp(1j * ph)
    r = np.outer(psi, psi.conj())
    idx = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for i, (a, b) in enumerate(idx):
        for j, (a2, b2) in enumerate(idx):
            da, db = a - a2, b - b2
            loc = np.exp(-gamma_t * (abs(da) + abs(db)))
            var = (da * s1 * sig1) ** 2 + (db * s2 * sig2) ** 2 + 2 * c * da * db * s1 * s2 * sig1 * sig2
            r[i, j] *= loc * np.exp(-0.5 * var)
    return r


def lmin(r):
    pt = r.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    return float(np.linalg.eigvalsh(pt).min())


def lmin_cs(Phi, gamma_t, var_c, var_u):
    """Common part (variance var_c, fully shared) + independent part (var_u each)."""
    v = var_c + var_u
    c = var_c / v if v > 0 else 0.0
    s = sqrt(v)
    return lmin(rho_noisy(Phi, gamma_t, s, s, c))


# requirement: noise may raise lambda_min (and hence N5 = 25/lambda^2) exactly as much as the Stage 16
# allocation did (ASSUMED rule carried over: same shot count, N5 = 5.0e3).
L_CLEAN = lmin(rho_noisy(PHI, GT, 0, 0, 0))
L_TARGET = BASE['lambda_min']
assert L_TARGET > L_CLEAN


def var_allowed(c=0.0, s2=1):
    f = lambda v: lmin(rho_noisy(PHI, GT, sqrt(v), sqrt(v), c, 1, s2)) - L_TARGET
    return brentq(f, 1e-12, 50.0)


# mapping variance <-> white acceleration PSD through the split profile w(t) (DERIVED):
# ramp up TAU_ACC, hold T_INT, ramp down TAU_ACC (linear ramps ASSUMED for the profile shape)
T_TOT = T_INT + 2 * TAU_ACC
def w_t(t):
    if t < TAU_ACC:
        return t / TAU_ACC
    if t < TAU_ACC + T_INT:
        return 1.0
    if t < T_TOT:
        return (T_TOT - t) / TAU_ACC
    return 0.0
INT_W2 = T_INT + 2 * TAU_ACC / 3.0          # int w^2 dt, DERIVED


def W2(f):
    """|W(f)|^2, W(f) = int w(t) e^{-2 pi i f t} dt, analytic for trapezoid (DERIVED)."""
    if f < 1e-9:
        return (T_INT + TAU_ACC) ** 2
    om = 2 * pi * f
    # trapezoid = (box of length T_INT+TAU) convolved with box of length TAU, /TAU
    L1 = T_INT + TAU_ACC
    s = lambda L: L * np.sinc(f * L)          # |FT| of unit box of length L
    return float((s(L1) * s(TAU_ACC) / TAU_ACC) ** 2)


def var_from_white(Sa):
    """Var theta for one-sided white S_a (m^2 s^-4 / Hz): kappa^2 S_a (1/2) int w^2 dt."""
    return KAPPA ** 2 * Sa * 0.5 * INT_W2


def white_req(var):
    return var / (KAPPA ** 2 * 0.5 * INT_W2)


F_CAL = 1e-3     # ASSUMED: phase drifts slower than ~100 shots x 10 s cycle are tracked and removed by calibration
F_MAX = 300.0


def white_equivalent(asd_fn, fcal=None):
    """Return the white ASD that gives the same Var theta as the coloured ASD asd_fn(f) (DERIVED)."""
    lf = np.linspace(np.log10(F_CAL if fcal is None else fcal), np.log10(F_MAX), 6000)
    f = 10 ** lf
    integrand = np.array([asd_fn(x) ** 2 * W2(x) for x in f]) * f * np.log(10)
    I = np.trapezoid(integrand, lf)
    return sqrt(I / (0.5 * INT_W2))


def line_white_equivalent(a0, f0):
    """Sinusoidal acceleration line, amplitude a0, random phase: Var theta = kappa^2 a0^2 |W(f0)|^2 / 2."""
    var = KAPPA ** 2 * a0 ** 2 * W2(f0) / 2.0
    return sqrt(white_req(var)), var


# ---------------------------------------------------------------------------
# 2. Source models: per-mass acceleration along the split axis x (parallel to the plate)
# ---------------------------------------------------------------------------
# Peterson (1993) NLNM, CITED: USGS OFR 93-322; table values checked against the GFZ NMSOP table
# (T=0.8 s: -169.2 dB, T=1.24 s: -163.7 dB) https://gfzpublic.gfz.de/rest/items/item_4017_4/component/file_4018/content
NLNM = [(0.10, -162.36, 5.64), (0.17, -166.7, 0.0), (0.40, -170.0, -8.30), (0.80, -166.4, 28.90),
        (1.24, -168.6, 52.48), (2.40, -159.98, 29.81), (4.30, -141.1, 0.0), (5.0, -71.36, -99.77),
        (6.0, -97.26, -66.49), (10.0, -132.18, -31.57), (12.0, -205.27, 36.16), (15.6, -37.65, -104.33),
        (21.9, -114.37, -47.10), (31.6, -160.58, -16.28), (45.0, -187.5, 0.0), (70.0, -216.47, 15.70),
        (101.0, -185.0, 0.0), (154.0, -168.34, -7.61), (328.0, -217.43, 11.90), (600.0, -258.28, 26.60),
        (10000.0, -346.88, 48.75)]


def nlnm_asd(f):
    P = 1.0 / f
    P = min(max(P, 0.10), 99999.0)
    row = NLNM[0]
    for r in NLNM:
        if P >= r[0]:
            row = r
    dB = row[1] + row[2] * log10(P)
    return sqrt(10 ** (dB / 10.0))


ISO_DISP_1HZ = 1e-11          # CITED aLIGO HAM/BSC-ISI residual 0.01 nm/rtHz at 1 Hz (sensor-noise limited)
                              # Matichard et al. CQG 32, 185003 (2015) https://doi.org/10.1088/0264-9381/32/18/185003
                              # quoted in RSI 96, 064504 (2025) https://pubs.aip.org/aip/rsi/article/96/6/064504/3350986
OMEGA_RL = 9.1e-11            # CITED ring laser G Wettzell resolution, rad/s/rtHz (upper bound on ambient rotation)
                              # https://syrte.obspm.fr/jsr/journees2005/pdf/s5_01_Kluegel.pdf (Table 1)
LPF_READ = 34.8e-15           # CITED LPF differential interferometer readout 34.8 fm/rtHz, PRL 116, 231101 (2016)
                              # https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.231101
EPS_AXIS = 1e-3               # ASSUMED: misalignment between split axis and the vertical (rad)
DI_I_WIRE = 1e-7              # ASSUMED: fractional current noise of a good chip-wire supply, /rtHz
DI_I_PERSIST = 1e-10          # ASSUMED: persistent-mode superconducting coil, /rtHz
SIGMA_AU_1K = 1e9             # ASSUMED: thin Au film conductivity at 1 K (S/m), residual-resistivity limited
B0 = 1.0                      # ASSUMED: field at the diamond in the diamagnetic trap (T); B dB/dz ~ 2e3 T^2/m needed
CHI_D = 2.2e-5                # CITED |chi| diamond (volume SI), standard handbook value
P_DIP = 1e-2 * E_CM           # CITED Schut 2023 default permanent dipole (Stage 12)
# electric-field noise: Labaziewicz et al. PRL 100, 013001 (2008), arXiv:0706.3763 (6 K surface traps, 75-150 um):
# n_dot = 2.1 quanta/s at ~1 MHz (88Sr+); S_E = 4 m hbar omega n_dot / q^2 (their eq.) -> DERIVED below.
# distance scaling beta = 4.0(2) at 5 K: Sedlacek et al. PRA 97, 020302 (2018) https://doi.org/10.1103/physreva.97.020302
M_SR = 87.9 * 1.66054e-27
S_E_REF = 4 * M_SR * HBAR * 2 * pi * 1e6 * 2.1 / E_CHG ** 2
Z_REF = 100e-6                # ASSUMED representative ion height for that data point (75-150 um traps)


def src_table():
    S = []
    req = sqrt(white_req(var_allowed(0.0)))  # placeholder, overwritten in main
    # 1 seismic, unisolated quiet site (NLNM)
    S.append(dict(key='seismic_nlnm', label='seismic, quiet site (NLNM), no isolation', c=1.0,
                  kind='common', asd=nlnm_asd, label_src='CITED Peterson 1993 NLNM'))
    # 2 seismic behind aLIGO-class isolation, flat in acceleration at the 1-Hz value (shape ASSUMED)
    a_iso = (2 * pi) ** 2 * ISO_DISP_1HZ
    S.append(dict(key='seismic_aligo_iso', label='seismic behind aLIGO-class active isolation', c=1.0,
                  kind='common', asd=lambda f, a=a_iso: a, label_src='CITED Matichard 2015 (1e-11 m/rtHz @1 Hz)'))
    # 3 tilt-horizontal coupling g*theta, theta = Omega/(2 pi f)
    S.append(dict(key='tilt', label='tilt x g (ring-laser bound)', c=1.0, kind='common',
                  asd=lambda f: G0 * OMEGA_RL / (2 * pi * f), label_src='CITED G Wettzell resolution'))
    # 4 rotation: alpha * D/2 per mass, anti-correlated
    S.append(dict(key='rotation_diff', label='rotation alpha x D/2 (anti-correlated)', c=-1.0, kind='differential',
                  asd=lambda f: 2 * pi * f * OMEGA_RL * D / 2, label_src='CITED G Wettzell; DERIVED lever arm'))
    # 5 trap current noise via axis misalignment
    S.append(dict(key='trap_current_wire', label='trap-drive current noise (chip wires), eps g 2dI/I', c=1.0,
                  kind='common', asd=lambda f: EPS_AXIS * G0 * 2 * DI_I_WIRE, label_src='ASSUMED dI/I=1e-7, eps=1e-3'))
    S.append(dict(key='trap_current_persist', label='trap current noise (persistent SC coil)', c=1.0,
                  kind='common', asd=lambda f: EPS_AXIS * G0 * 2 * DI_I_PERSIST, label_src='ASSUMED dI/I=1e-10'))
    # 6 magnetic Johnson noise of a normal-metal plate (thin-film limit, Varpula & Poutanen J. Appl. Phys. 55,
    #   4015 (1984) form S_B = mu0^2 k T sigma t / (8 pi z^2), ASSUMED prefactor); F_x ~ (chi V/mu0) B0 dB/z
    V = 4 / 3 * pi * R ** 3
    SB = MU0 ** 2 * K_B * CAND['T_env'] * SIGMA_AU_1K * W_PLATE / (8 * pi * Z_C ** 2)
    aJ = CHI_D * V / MU0 * B0 * sqrt(SB) / Z_C / M
    S.append(dict(key='johnson_plate', label='magnetic Johnson noise, Au plate 1 K (B0=1 T)', c=0.0,
                  kind='differential', asd=lambda f, a=aJ: a, label_src='ASSUMED sigma=1e9 S/m, t=1 um; form CITED',
                  extra=dict(sqrt_SB=sqrt(SB))))
    # 7 electric-field (patch) noise on the permanent dipole: F ~ p dE/z, S_E scaled (z_ref/z)^4
    SE_white = S_E_REF * (Z_REF / Z_C) ** 4
    S.append(dict(key='efield_white', label='surface E-field noise on dipole, white from 1 MHz', c=0.0,
                  kind='differential', asd=lambda f, s=SE_white: P_DIP * sqrt(s) / Z_C / M,
                  label_src='CITED Labaziewicz 2008, Sedlacek 2018; extrapolation ASSUMED'))
    S.append(dict(key='efield_1overf', label='surface E-field noise on dipole, 1/f from 1 MHz', c=0.0,
                  kind='differential', asd=lambda f, s=SE_white: P_DIP * sqrt(s * 1e6 / f) / Z_C / M,
                  label_src='CITED Labaziewicz 2008 normalisation; 1/f to 1 Hz ASSUMED'))
    # 8 seismic Newtonian noise, |dg| ~ 2 pi G rho xi (Rayleigh-wave order of magnitude, ASSUMED prefactor)
    S.append(dict(key='newtonian_seismic', label='Newtonian noise from NLNM ground motion', c=1.0, kind='common',
                  asd=lambda f: 2 * pi * G * 2600.0 * nlnm_asd(f) / (2 * pi * f) ** 2,
                  label_src='ASSUMED 2 pi G rho xi, rho=2600'))
    return S, dict(aJ=aJ, SB=SB, SE_white=SE_white)


def ggn_person(r, m_p=70.0, dr=0.3, f0=1.0):
    """Person of 70 kg moving 0.3 m at ~1 Hz at distance r (ASSUMED).  Trap is fixed to the ground, so each
    test mass feels the full delta(G m/r^2) relative to its trap: a0 = 2 G m dr / r^3 (common); the
    differential (tidal) part is 6 G m dr D / r^4."""
    a_c = 2 * G * m_p * dr / r ** 3
    a_d = 6 * G * m_p * dr * D / r ** 4
    return a_c, a_d


def main():
    res = dict(stage=18, baseline=dict(m=M, dx=DX, D=D, R=R, Phi=PHI, t_int=T_INT, t_tot=T_TOT,
                                       kappa_rad_per_mps=KAPPA, z_centre=Z_C,
                                       stage16_sqrt_Sa_req=BASE['sqrt_Sa_req']))
    # --- 1. witness response to correlated noise --------------------------------------------------
    v0 = var_allowed(0.0)
    vp = var_allowed(1.0)
    vm = var_allowed(-1.0)
    req0, reqp, reqm = (sqrt(white_req(v)) for v in (v0, vp, vm))
    eps_list = [1e-4, 1e-3, 1e-2, 1e-1]
    mism = []
    for e in eps_list:      # trap-frequency mismatch: theta2 = (1+2 dw/w) theta1 for SG-in-trap (dx ~ 1/w^2)
        f = lambda v: lmin(rho_noisy(PHI, GT, sqrt(v), (1 + 2 * e) * sqrt(v), 1.0)) - L_TARGET
        vv = brentq(f, 1e-12, 50.0)
        mism.append(dict(dw_over_w=e, sqrt_Sa_allowed=sqrt(white_req(vv)), cmrr_vs_uncorr=sqrt(vv / v0)))
    # sweep sigma for the figure
    sig = np.geomspace(0.03, 5, 60)
    sweep = {k: [lmin(rho_noisy(PHI, GT, s, s, c)) for s in sig] for k, c in
             (('uncorrelated', 0.0), ('common', 1.0), ('anti', -1.0))}
    # allowed region: common ASD vs independent ASD
    sc_grid = np.geomspace(1e-17, 1e-14, 61)
    su_allowed = []
    for sc in sc_grid:
        vc = var_from_white(sc ** 2)
        g = lambda vu: lmin_cs(PHI, GT, vc, vu) - L_TARGET
        if g(0.0) >= 0:
            su_allowed.append(0.0)
        else:
            su_allowed.append(sqrt(white_req(brentq(g, 0.0, 50.0))))
    res['witness'] = dict(lambda_clean=L_CLEAN, lambda_target=L_TARGET, N5_target=25 / L_TARGET ** 2,
                          var_allowed=dict(uncorrelated=v0, common=vp, anti=vm),
                          sqrt_Sa_allowed=dict(uncorrelated=req0, common=reqp, anti=reqm),
                          intrinsic_cmrr=dict(common=sqrt(vp / v0), anti=sqrt(vm / v0)),
                          mismatch=mism,
                          sigma_sweep=dict(sigma=sig.tolist(), **sweep),
                          region=dict(sqrt_Sc=sc_grid.tolist(), sqrt_Su_allowed=su_allowed),
                          stage16_model_req=BASE['sqrt_Sa_req'])
    REQ = min(req0, reqp, reqm)
    # --- 2. sources --------------------------------------------------------------------------------
    S, aux = src_table()
    rows = []
    for s in S:
        we = white_equivalent(s['asd'])
        we01 = white_equivalent(s['asd'], 0.1)   # ASSUMED aggressive variant: drifts below 0.1 Hz tracked
        rows.append(dict(key=s['key'], label=s['label'], kind=s['kind'], c=s['c'], src=s['label_src'],
                         asd_1Hz=s['asd'](1.0), white_equiv=we, white_equiv_fcal01=we01,
                         orders_over_req_fcal01=log10(we01 / REQ), orders_over_req=log10(we / REQ),
                         orders_over_pino=log10(we / PINO_REQ)))
    # GGN people (line at 1 Hz)
    ggn = []
    for r in (2.0, 10.0, 30.0, 100.0, 300.0):
        a_c, a_d = ggn_person(r)
        we_c, _ = line_white_equivalent(a_c, 1.0)
        we_d, _ = line_white_equivalent(a_d, 1.0)
        ggn.append(dict(r=r, a_common=a_c, a_diff=a_d, white_equiv_common=we_c, white_equiv_diff=we_d,
                        orders_common=log10(we_c / REQ)))
    # exclusion radius for one person (DERIVED): white_equiv_common(r) = REQ
    r_excl = brentq(lambda r: line_white_equivalent(ggn_person(r)[0], 1.0)[0] - REQ, 1.0, 1e5)
    # pulse-tube line at 1.4 Hz (CITED typical PT frequency): allowed displacement amplitude (DERIVED)
    f_pt = 1.4
    a_allow = sqrt(2 * white_req(1.0) * 0 + 2 * v0 / (KAPPA ** 2 * W2(f_pt)))  # Var = kappa^2 a0^2 W2/2 = v_allowed
    x_allow = a_allow / (2 * pi * f_pt) ** 2
    # plate-vibration coupling through the image-dipole stiffness (Stage 12 a_dip = 0.25 m/s^2 at 5 um)
    a_dip = 0.25
    kd_over_m = 4 * a_dip / 5e-6
    dz_allow = REQ / (EPS_AXIS * kd_over_m)
    # reference-sensor (feed-forward) floors at 1 Hz -> white-equivalent
    lpf_asd = lambda f: (2 * pi * f) ** 2 * LPF_READ
    floors = dict(
        aligo_isolation_sensor=dict(val=white_equivalent(lambda f: (2 * pi) ** 2 * ISO_DISP_1HZ),
                                    src='aLIGO ISI GS13-limited, 1e-11 m/rtHz at 1 Hz (flat-acceleration shape ASSUMED)'),
        lpf_readout=dict(val=white_equivalent(lpf_asd),
                         src='LPF interferometer 34.8 fm/rtHz x (2 pi f)^2, space (PRL 116, 231101)'),
        lpf_mHz=dict(val=1.74e-15, src='LPF 1.74 fm s^-2/rtHz, mHz band only (Armano 2018)'),
        asenbaum_diff=dict(val=5.4e-11 * G0, src='Asenbaum PRL 125, 191101 (2020): eta 5.4e-11/rtHz -> 5.3e-10 m s^-2/rtHz differential'),
        mcguirk_diff=dict(val=4e-9 * G0, src='McGuirk et al. arXiv:physics/0105088: 4e-9 g/rtHz differential'),
    )
    for k in floors:
        floors[k]['orders_over_req'] = log10(floors[k]['val'] / REQ)
    # thermal channel
    gam_max = M * REQ ** 2 / (4 * K_B * CAND['T_env'])
    gam_E = epstein_gamma(BASE['P_req'], CAND['T_gas'], R, M, delta=1.39)
    gam_LW_equiv = M * gam_E * K_B * CAND['T_gas'] * DX ** 2 / HBAR ** 2   # would-be LW rate (not physical in SW)
    gamma_J = M * aux['aJ'] ** 2 / (4 * K_B * CAND['T_env'])
    # Johnson: sigma*t needed for aJ <= REQ (a ~ sqrt(T sigma t))
    sigt_max = SIGMA_AU_1K * W_PLATE * (REQ / aux['aJ']) ** 2
    # electric: z needed for white case (a ~ z^-3) and dipole needed at 5 um
    aE_w = [r for r in rows if r['key'] == 'efield_white'][0]['white_equiv']
    z_need = Z_C * (aE_w / REQ) ** (1 / 3)
    p_need = P_DIP * REQ / aE_w
    # charge: one stray electron at 5 um (white case)
    a_e1 = E_CHG * sqrt(aux['SE_white']) / M
    thermal = dict(gamma_max_1K=gam_max, Q_min_100Hz=2 * pi * 100 / gam_max,
                   gamma_epstein_at_Preq=gam_E, gas_LW_formula_rate=gam_LW_equiv,
                   gas_SW_rate_used=BASE['gamma_gas'],
                   gamma_eddy_johnson=gamma_J, orders_eddy_over_max=log10(gamma_J / gam_max),
                   sigma_t_max_Sm=sigt_max, sigma_t_assumed=SIGMA_AU_1K * W_PLATE)
    res['sources'] = rows
    res['ggn_person'] = dict(rows=ggn, exclusion_radius_m=r_excl)
    res['pulse_tube'] = dict(f=f_pt, a_allowed=a_allow, x_allowed=x_allow)
    res['plate_vibration'] = dict(kdip_over_m=kd_over_m, eps_axis=EPS_AXIS, dz_allowed=dz_allow)
    res['reference_floors'] = floors
    res['thermal'] = thermal
    res['efield'] = dict(S_E_ref=S_E_REF, z_ref=Z_REF, S_E_5um_white=aux['SE_white'], z_needed=z_need,
                         p_needed_ecm=p_need / E_CM, a_one_electron=a_e1)
    res['johnson'] = dict(sqrt_SB=sqrt(aux['SB']), a=aux['aJ'])
    res['REQ'] = REQ
    res['PINO_REQ'] = PINO_REQ
    # --- 3. bottom line ----------------------------------------------------------------------------
    common = [r for r in rows if r['kind'] == 'common']
    diff = [r for r in rows if r['kind'] == 'differential']
    best_common_ground = min(r['white_equiv'] for r in rows if r['key'] in ('seismic_aligo_iso',))
    tilt = [r for r in rows if r['key'] == 'tilt'][0]['white_equiv']
    # Scenario A: witness-intrinsic CMRR only (what the PPT test actually gives)
    scen = {}
    def tot(keys):
        return sqrt(sum([r['white_equiv'] ** 2 for r in rows if r['key'] in keys]))
    keysA = ['seismic_aligo_iso', 'tilt', 'rotation_diff', 'trap_current_persist', 'johnson_plate',
             'efield_white', 'newtonian_seismic']
    scen['A_ground_isolated_intrinsic_cmrr'] = dict(keys=keysA, total=tot(keysA))
    # Scenario B: perfect feed-forward of all common-mode down to the best ground inertial sensor floor
    keysB_diff = ['rotation_diff', 'johnson_plate', 'efield_white']
    scen['B_feedforward_ground_sensor'] = dict(total=sqrt(tot(keysB_diff) ** 2 + floors['aligo_isolation_sensor']['val'] ** 2))
    # Scenario C: space-grade reference (LPF readout) and SC/no-metal plate, E-field excluded (best case inertial)
    scen['C_lpf_reference_inertial_only'] = dict(total=sqrt(floors['lpf_readout']['val'] ** 2 +
                                                           [r for r in rows if r['key'] == 'rotation_diff'][0]['white_equiv'] ** 2))
    # Scenario D: ground, inertial only, best isolation, drifts below 0.1 Hz tracked, plate surface noise ignored
    keysD = ['seismic_aligo_iso', 'tilt', 'rotation_diff', 'newtonian_seismic']
    scen['D_ground_inertial_only_fcal01'] = dict(keys=keysD, total=sqrt(sum(r['white_equiv_fcal01'] ** 2 for r in rows if r['key'] in keysD)))
    for k in scen:
        scen[k]['orders_over_req'] = log10(scen[k]['total'] / REQ)
        scen[k]['orders_over_pino'] = log10(scen[k]['total'] / PINO_REQ)
    binding = max(rows, key=lambda r: r['white_equiv'] if r['key'] in keysA else 0)
    res['scenarios'] = scen
    res['binding_in_A'] = binding['key']
    with open(os.path.join(ROOT, 'results', 'stage18.json'), 'w') as fh:
        json.dump(res, fh, indent=1, default=float)
    # --- print --------------------------------------------------------------------------------------
    print('baseline m=%.3g dx=%.3g D=%.3g Phi=%.3f kappa=%.3g rad/(m/s) z=%.3g' % (M, DX, D, PHI, KAPPA, Z_C))
    print('lambda clean %.4f target %.4f (N5 %.0f)' % (L_CLEAN, L_TARGET, 25 / L_TARGET ** 2))
    print('allowed sigma^2: uncorr %.4f common %.4f anti %.4f ; sqrtSa req fm: %.3f %.3f %.3f (stage16 %.3f)' %
          (v0, vp, vm, req0 * 1e15, reqp * 1e15, reqm * 1e15, BASE['sqrt_Sa_req'] * 1e15))
    print('intrinsic witness CMRR (ASD): common %.3f anti %.3f' % (sqrt(vp / v0), sqrt(vm / v0)))
    for m_ in mism:
        print('  mismatch dw/w=%.0e: allowed %.3f fm, cmrr %.3f' % (m_['dw_over_w'], m_['sqrt_Sa_allowed'] * 1e15, m_['cmrr_vs_uncorr']))
    print('REQ used = %.3f fm' % (REQ * 1e15))
    for r in rows:
        print('%-22s %-12s asd(1Hz)=%.2e white-eq=%.2e  orders=%+.2f  (fcal0.1: %.2e %+.2f) [%s]' %
              (r['key'], r['kind'], r['asd_1Hz'], r['white_equiv'], r['orders_over_req'], r['white_equiv_fcal01'],
               r['orders_over_req_fcal01'], r['src']))
    for g in ggn:
        print('GGN person r=%5.0f m: a_c=%.2e (we %.2e, %+.2f) a_d=%.2e (we %.2e)' %
              (g['r'], g['a_common'], g['white_equiv_common'], g['orders_common'], g['a_diff'], g['white_equiv_diff']))
    print('person exclusion radius %.0f m' % r_excl)
    print('pulse tube 1.4 Hz allowed a0=%.2e m/s^2, x0=%.2e m' % (a_allow, x_allow))
    print('plate vibration allowed dz=%.2e m/rtHz' % dz_allow)
    for k, v in floors.items():
        print('floor %-24s %.2e  (%+.2f)  %s' % (k, v['val'], v['orders_over_req'], v['src']))
    print('thermal:', {k: ('%.3g' % v) for k, v in thermal.items()})
    print('efield:', {k: ('%.3g' % v) for k, v in res['efield'].items()})
    print('johnson:', res['johnson'])
    for k, v in scen.items():
        print('scenario %-34s total %.2e  orders %+.2f (Pino %+.2f)' % (k, v['total'], v['orders_over_req'], v['orders_over_pino']))
    print('binding in A:', binding['key'])


if __name__ == '__main__':
    main()
