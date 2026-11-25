#!/usr/bin/env python3
"""Stage 27 - signature of CQ / Diosi-type gravitational noise, and how far 'tuning' can push the Stage 26 shortfalls.
Builds on Stage 26: every probe sees N_E = A_E(l)(kappa + 1/(4 kappa)) (DP kernel; minimum at kappa = 1/2).
Labels: DERIVED (computed here), CITED (number from a source), ASSUMED (our choice, stated).
"""
import json, os, time
import numpy as np
from math import pi, sqrt, log10
from scipy.special import erf
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
H = 6.62607015e-34; HBAR = H / (2 * pi); G = 6.6743e-11; U = 1.66053906660e-27; KB = 1.380649e-23
GH = G * HBAR
t0 = time.time()
s26 = json.load(open(os.path.join(ROOT, 'results/stage26.json')))
E26 = np.array(s26['scan']['ell'])
def lu26(name, ell): return float(np.interp(np.log(ell), np.log(E26), np.log10(np.array(s26['scan']['u'][name], float))))

# ================================================================== 1. SIGNATURE (DERIVED)
# Normalised body integral J = I/(M_A M_B), Gamma_AB(kappa) = 4 pi G hbar J_AB (kappa + 1/(4 kappa)).
# self, plateau (spacing << l << R): J = 1/(3V)                      -> Gamma = (4 pi/3) G hbar / V
# self, body << l:                    J = 1/(24 pi^1.5 l^3)            -> Gamma = G hbar /(6 sqrt(pi) l^3)
# self, incoherent (l << spacing):    J = <m^2>/<m> /(24 pi^1.5 l^3 M) -> composition dependent
# cross, two bodies, gap >> l:  J_ij = -d_i d_j [1/(4 pi r)] = (delta_ij - 3 n_i n_j)/(4 pi d^3)   (shell theorem: exact for spheres)
#      longitudinal J = -1/(2 pi d^3) -> Gamma_AB = -2 G hbar/d^3 ; transverse +1/(4 pi d^3) -> +G hbar/d^3
# cross, two point-like bodies (R << l): J_ij = -d_i d_j phi, phi = erf(r/2l)/(4 pi r)
def phi_derivs(r, ell):
    a = 2 * ell
    e = erf(r / a); g = (2 / sqrt(pi)) * np.exp(-(r / a) ** 2) / a           # d/dr erf(r/a)
    phi = e / (4 * pi * r)
    d1 = (g / r - e / r ** 2) / (4 * pi)
    gp = -2 * r / a ** 2 * g
    d2 = (gp / r - 2 * g / r ** 2 + 2 * e / r ** 3) / (4 * pi)
    return phi, d1, d2
def rho_points(d, ell):
    """correlation coefficient of two identical point-like bodies smeared at l (independent of kappa)."""
    _, d1, d2 = phi_derivs(d, ell); J0 = 1 / (24 * pi ** 1.5 * ell ** 3)
    return -d2 / J0, -d1 / d / J0
def rho_spheres_plateau(d, V):  # equal spheres, spacing << l << R, gap >> l
    return -3 * V / (2 * pi * d ** 3), 3 * V / (4 * pi * d ** 3)

sig = {}
sig['volume_constant_Gamma_times_V_min'] = (4 * pi / 3) * GH          # m^5 s^-3, at kappa = 1/2 (= Diosi kappa=1)
sig['cross_longitudinal_min'] = 'Gamma_AB = -2 G hbar / d^3 (kappa=1/2); = Kafri-Taylor-Milburn mutual noise'
sig['spectrum'] = 'white (Markovian): flat from DC to at least the X-ray band (1e21 rad/s) if the X-ray bound is to apply; no model cutoff'
sig['T_P_dependence'] = 'none (no temperature, pressure, charge, magnetic or laser dependence)'
sig['composition'] = dict(plateau='independent of composition and mass at fixed volume (Gamma ~ 1/V)',
                          below_l='universal, body-independent G hbar/(6 sqrt(pi) l^3)',
                          incoherent_l_below_0p6A='~ <m^2>/<m> / M (heavy nuclei noisier per unit mass) - excluded region (l < 4.5 A)',
                          Au_over_Si_incoherent_equal_mass=(196.97) / (28.0855))
# scaling exponents of acceleration PSD for a sphere of radius R, density rho (ASSUMED standard background models)
sig['scaling_table'] = [
    dict(source='CQ / Diosi-type (plateau)', R=-3, rho=0, T=0, P=0, f=0, d_cross='-2 G hbar/d^3 (anti-correlated along d)'),
    dict(source='CQ / Diosi-type (body << l)', R=0, rho=0, T=0, P=0, f=0, d_cross='-> +1 for d << l'),
    dict(source='residual-gas Brownian (molecular regime)', R=-4, rho=-2, T=0.5, P=1, f=0, d_cross='0 (independent molecules)'),
    dict(source='blackbody radiation-pressure / recoil', R=-4, rho=-2, T=5, P=0, f=0, d_cross='~0'),
    dict(source='internal (structural) damping', R='geometry', rho='material', T=1, P=0, f=-1, d_cross='0'),
    dict(source='electrostatic patch / charge', R='surface/gap', rho=-2, T='weak', P=0, f=-1, d_cross='0 (housing-specific)'),
    dict(source='seismic / platform', R=0, rho=0, T=0, P=0, f='coloured', d_cross='~+1 (common mode)'),
    dict(source='laser shot / readout (displacement-white)', R=0, rho=0, T=0, P=0, f=4, d_cross='0'),
]

# correlation coefficients for real geometries
geo = []
V_lpf = 0.046 ** 3
for lab, d, V in (('LISA Pathfinder TM1-TM2 (46 mm cubes, 376 mm)', 0.376, V_lpf),
                  ('torsion dumbbell, two Au spheres 2.18e-4 kg, 40 mm apart (ASSUMED geometry)', 0.040, 2.18e-4 / 19300),
                  ('two levitated silica spheres R = 4.8 um (Monteiro-class), d = 3R', 3 * 4.83e-6, 4 / 3 * pi * 4.83e-6 ** 3),
                  ('two levitated spheres R = 1 um, d = 3 um', 3e-6, 4 / 3 * pi * 1e-6 ** 3),
                  ('two touching spheres (d = 2R)', 2.0, 4 / 3 * pi * 1.0 ** 3)):
    rl, rt = rho_spheres_plateau(d, V)
    geo.append(dict(geometry=lab, d=d, V=V, rho_long=rl, rho_trans=rt, Gamma_cross_min=2 * GH / d ** 3,
                    Gamma_self_min=(4 * pi / 3) * GH / V))
d_rho_1pct_over_R = (300 * 4 / 3 * pi / (2 * pi)) ** (1 / 3)   # |rho_long| = 0.01 for spheres: d^3 = 300 V/(2 pi)
# point-like bodies (R << l): rho vs d/l
dl = np.logspace(-2, 2, 400)
rpL, rpT = rho_points(dl, 1.0)
zero_cross = float(dl[np.where(np.diff(np.sign(rpL)))[0][0]])
pts = dict(d_over_l=dl.tolist(), rho_long=rpL.tolist(), rho_trans=rpT.tolist(), longitudinal_zero_crossing_d_over_l=zero_cross,
           far_field='rho_long -> -12 sqrt(pi) (l/d)^3', check_far=float(rpL[-1] / (-12 * sqrt(pi) / dl[-1] ** 3)))
dsph = np.logspace(np.log10(2.0), 3, 300)
sph = dict(d_over_R=dsph.tolist(), rho_long=(-3 * (4 / 3 * pi) / (2 * pi * dsph ** 3)).tolist(), rho_trans=(3 * (4 / 3 * pi) / (4 * pi * dsph ** 3)).tolist())

# ================================================================== 2. TUNING (DERIVED gains)
def gain_psd(T, B, eps):
    """excess-noise power estimation: statistical floor 1/sqrt(T B) on the PSD, systematic floor eps (background-model accuracy)."""
    stat = 0.5 * log10(T * B); sysg = log10(1 / eps)
    return dict(stat_orders=stat, sys_orders=sysg, realistic_orders=min(stat, sysg))
EPS = dict(conservative=0.10, optimistic=0.05)   # ASSUMED background-model accuracy on the PSD
tune = {}
base = dict(LPF=lu26('LISA Pathfinder 2018 (Armano)', 1e-5), MAQRO_baseline=lu26('MAQRO baseline 1e10 u silica, dx=R, 100 s', 1e-9),
            torsion=lu26('Westphal 2021 mm torsion balance', 1e-6), levitated=lu26('Monteiro 2020 levitated ng sphere', 1e-6),
            Gisler=lu26('Gisler 2022 SiN nanowire', 1e-9))
# (a) PSD integration - published durations / bands (ASSUMED where noted)
A_params = dict(LPF=dict(T=1.1e6, B=0.019, note='~13 d quiet run (Armano 2018), white band 1-20 mHz (ASSUMED band)'),
                torsion=dict(T=1e6, B=0.01, note='ASSUMED 12 d, 10 mHz band'),
                levitated=dict(T=1e5, B=10.0, note='ASSUMED 1 d, 10 Hz band'),
                Gisler=dict(T=1e4, B=1e3, note='ASSUMED'))
for k, p in A_params.items():
    for lab, e in EPS.items():
        g = gain_psd(p['T'], p['B'], e)
        tune.setdefault(k, {})['a_' + lab] = dict(**g, new_log10u=base[k] - g['realistic_orders'], note=p['note'])
# LPF public-data reanalysis: pressure-intercept fit of the white (Brownian) branch over ~500 d (Armano 2024, arXiv:2405.05207)
# S_white(t) = kappa_gas P(t) + S_0 ; the floor S_0 contains CQ. Brownian PSD fell by ~(5.2/1.74)^2 = 8.9 over the mission (CITED 2016 vs 2017 ASDs).
lpf_reanalysis = dict(psd_drop_factor=(5.2 / 1.74) ** 2, assumed_runs=20, assumed_per_run_psd_error=0.05,
                      intercept_error_fraction_of_final_psd=None)
# simple linear-fit propagation: points spread over x = P/P_final in [1, 8.9], equal relative errors 5% of each point
x = np.linspace(1, lpf_reanalysis['psd_drop_factor'], lpf_reanalysis['assumed_runs']); sigy = 0.05 * x
W = 1 / sigy ** 2; Sw, Sx, Sxx = W.sum(), (W * x).sum(), (W * x * x).sum()
sig_intercept = sqrt(Sxx / (Sw * Sxx - Sx ** 2))           # in units of final PSD (x=1 point = kappa*P_final ~ total)
lpf_reanalysis['intercept_error_fraction_of_final_psd'] = sig_intercept
lpf_reanalysis['gain_orders_2sigma'] = log10(1 / (2 * sig_intercept))
lpf_reanalysis['new_log10u'] = base['LPF'] - lpf_reanalysis['gain_orders_2sigma']
lpf_reanalysis['note'] = ('feasible with public Level-1..3 Delta-g + housekeeping (pressure proxy) from the LPF Legacy Archive; '
                          'cold run (274.8 K) already below the warm runs (Armano 2024)')
tune['LPF_public_reanalysis'] = lpf_reanalysis
# MAQRO: replace the conservative lambda <= 1/tau rule by the statistical contrast resolution
lamDP_maqro = (1 / 100.0) / 10 ** lu26('MAQRO baseline 1e10 u silica, dx=R, 100 s', 1e-9)   # DP(kappa=1) rate at the plateau
mq = {}
for lab, Ndet in (('one run, 1e4 detections (CITED MAQRO up to ~1e4 data points per run)', 1e4),
                  ('30 runs over the 2-yr mission (ASSUMED 50% duty)', 3e5)):
    Delta = 2 * sqrt(2 / Ndet)           # 2-sigma visibility resolution (ASSUMED Poisson, V ~ 1)
    lam_min = -np.log(1 - Delta) / 100.0
    newu = log10(lam_min / lamDP_maqro)
    # l-range excluded: u_MAQRO(l) scales with the Stage 26 curve; shift by the same gain
    curve = np.log10(np.array(s26['scan']['u']['MAQRO baseline 1e10 u silica, dx=R, 100 s'], float)) - (base['MAQRO_baseline'] - newu)
    ex = (curve < 0) & (E26 >= 4.5e-10)
    mq[lab] = dict(N_det=Ndet, visibility_resolution=Delta, lambda_min=lam_min, lambda_DP_kappa1=lamDP_maqro, new_log10u=newu,
                   gain_orders=base['MAQRO_baseline'] - newu,
                   excludes_l=(float(E26[ex][0]), float(E26[ex][-1])) if ex.any() else None)
tune['MAQRO_baseline_statistics'] = mq
# (c) cross-correlation of two sensors: threshold c_thr = 2 sqrt(Gamma_A Gamma_B / (T B)) on the correlated part c = 2 G hbar/d^3
cc = []
for lab, d, Gam, T, B in (('LISA Pathfinder TMs (hypothetical independent channels)', 0.376, 7.57e-31, 1.1e6, 0.019),
                          ('two Monteiro-class levitated spheres, d = 3R = 14.5 um', 3 * 4.83e-6, 4.34e-13, 1e5, 10.0),
                          ('two Gisler-class nanowires 10 um apart', 10e-6, 5.30e-13, 1e4, 1e3),
                          ('two Westphal-class torsion masses 40 mm apart', 0.04, 8.65e-18, 1e6, 0.01)):
    c = 2 * GH / d ** 3; thr = 2 * Gam / sqrt(T * B)
    cc.append(dict(setup=lab, d=d, Gamma_bg=Gam, T=T, B=B, correlated_signal_min=c, threshold_2sigma=thr, orders_short=log10(thr / c)))
tune['c_cross_correlation'] = cc
# (d) composition / volume differential
tune['d_composition'] = ('plateau (all l in the surviving window for bodies larger than l): CQ is composition-independent at fixed V '
                         '-> differential composition signal = 0 (use only as a veto). Volume differential: CQ ~ V^-1 vs gas ~ R^-4/rho^2; '
                         'separation limited by the same background-model accuracy as (a).')

# (b) modulation: DP / translation-invariant kernels have NO external handle (noise exists in vacuum).
# State-dependent (Oppenheim ultra-local discrete) kernel: field noise is sourced by matter -> a moving source mass modulates it.
# record noise <dm dm> = Dhat rho/4 ; Gamma_x on a body from a point source M_s on the x axis at r: (G^2 Dhat/4) M_s / r^4
mXe = 131.29 * U; GX = GH / (6 * sqrt(pi) * 4.5e-10 ** 3)               # XENONnT, per axis two-sided
lpf = dict(R=(3 * V_lpf / (4 * pi)) ** (1 / 3), rho=1.928 / V_lpf, Gamma=7.57e-31)
Dhat_max_LPF_self = 15 * lpf['R'] * lpf['Gamma'] / (pi * lpf['rho'] * G ** 2)      # Gamma_x = pi rho G^2 Dhat/(15 R)
a_LXe = (mXe / 2940.0) ** (1 / 3)
def Dhat_max_xray_env(ell):     # liquid-Xe neighbours source field noise on each nucleus, cutoff max(a, l)
    rc = max(a_LXe, ell); return GX / (G ** 2 / 4 * 2940.0 * (4 * pi / 3) / rc)
def Dhat_min_xray_heat(ell): return HBAR ** 2 / (ell ** 2 * mXe * GX)
ELL = np.logspace(-10, -4, 241)
dmin = np.array([Dhat_min_xray_heat(e) for e in ELL]); dmax = np.minimum(Dhat_max_LPF_self, [Dhat_max_xray_env(e) for e in ELL])
Wd = np.where(dmax > dmin, np.log10(dmax / dmin), np.nan); Wd = np.where(ELL <= 52e-6, Wd, np.nan)
l_min_disc = float(ELL[np.isfinite(Wd)][0])
# source-mass lock-in (ASSUMED source: 10 kg W at 5 cm; LPF-class sensor; T = 1e6 s, B = 0.019 Hz; 2-sigma on Delta S between states)
src = dict(M_s=10.0, r=0.05)
def Dhat_detectable(Gam, T, B): thr = 2 * 2 * Gam / sqrt(T * B); return 4 * thr * src['r'] ** 4 / (G ** 2 * src['M_s'])
D_det_space = Dhat_detectable(7.57e-31, 1.1e6, 0.019)
D_det_torsion = Dhat_detectable(8.65e-18, 1e6, 0.01) * (0.01 / 0.05) ** 4 * (10 / 1.0)    # torsion: 1 kg at 1 cm
cover = np.where(np.isfinite(Wd), np.log10(np.minimum(dmax, np.maximum(dmin, D_det_space)) / dmin), np.nan)
disc = dict(Dhat_max_LPF_self_kg=Dhat_max_LPF_self, a_LXe=a_LXe, l_min=l_min_disc,
            window_at={('%.0e' % e): float(np.interp(np.log(e), np.log(ELL), Wd)) for e in (3e-6, 1e-5, 5e-5)},
            source=src, Dhat_detectable_LPFclass_kg=D_det_space, Dhat_detectable_torsion_1kg_1cm_kg=D_det_torsion,
            fraction_of_window_covered_LPFclass=float(np.nanmin(np.where(np.isfinite(Wd), (np.log10(dmax) - np.log10(np.maximum(dmin, D_det_space))) / Wd, np.nan))),
            torsion_orders_short_at_window_top=float(log10(D_det_torsion / Dhat_max_LPF_self)),
            scan=dict(ell=ELL.tolist(), Dhat_min=dmin.tolist(), Dhat_max=dmax.tolist(), window=Wd.tolist()))
tune['b_modulation'] = dict(DP='no external handle: translation-invariant D0 -> noise identical in vacuum; only body swaps (V) or background modulation (P,T) -> same limit as (a)',
                            discrete=disc)

summary = dict(
    LPF=dict(stage26=base['LPF'], best=min(tune['LPF']['a_optimistic']['new_log10u'], lpf_reanalysis['new_log10u']), route='(a)/public reanalysis'),
    MAQRO_baseline=dict(stage26=base['MAQRO_baseline'], best=mq['30 runs over the 2-yr mission (ASSUMED 50% duty)']['new_log10u'], route='contrast statistics'),
    torsion=dict(stage26=base['torsion'], best=tune['torsion']['a_optimistic']['new_log10u'], route='(a)'),
    levitated=dict(stage26=base['levitated'], best=min(tune['levitated']['a_optimistic']['new_log10u'], cc[1]['orders_short']), route='(a) or (c) pair'),
)
out = dict(stage=27, signature=sig, geometries=geo, d_over_R_for_1pct=d_rho_1pct_over_R, points=pts, spheres=sph,
           tuning=tune, summary=summary,
           sources=dict(kryhin_sudhir='https://arxiv.org/abs/2309.09105 (PRL 134, 061501)', ktm='https://arxiv.org/abs/1401.0946',
                        altamirano='https://arxiv.org/abs/1612.07735', tilloy_diosi='https://arxiv.org/abs/1706.01856',
                        helou='https://arxiv.org/abs/1606.03637 (PRD 95, 084054)', lpf_kg_update='https://arxiv.org/abs/2411.17588',
                        armano2024='https://arxiv.org/abs/2405.05207', boileau='https://arxiv.org/abs/2204.03867',
                        angeli='https://arxiv.org/abs/2501.13030', lpf_archive='https://lpf.esac.esa.int/lpfsa/',
                        heasarc_lpf='https://heasarc.gsfc.nasa.gov/lpf/'),
           runtime_s=time.time() - t0)
json.dump(out, open(os.path.join(ROOT, 'results/stage27.json'), 'w'), indent=1, default=float)
print(json.dumps(dict(geo=geo, d1pct=d_rho_1pct_over_R, zero=zero_cross, far=pts['check_far'], summary=summary,
                      lpf=lpf_reanalysis, mq=mq, cc=cc, disc={k: v for k, v in disc.items() if k != 'scan'},
                      a={k: tune[k] for k in A_params}), indent=1, default=str))
