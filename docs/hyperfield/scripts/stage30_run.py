#!/usr/bin/env python3
"""Stage 30 - Big G vs test-mass volume: can CQ/DP gravitational noise (volume law, Stages 27/29) explain the ~500 ppm G spread?
Data: the 16 CODATA 2018 = CODATA 2022 input values (CODATA 2022, arXiv:2409.03787, Table XXX: "there is no new relevant input datum
for the 2022 adjustment").  Geometry: CITED where verified this stage, otherwise ASSUMED (recalled from the original papers, not re-checked).
Labels: DERIVED / CITED / ASSUMED.
"""
import json, os, csv, time
import numpy as np
from math import pi, sqrt, log10
from scipy import stats
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
t0 = time.time(); rng = np.random.default_rng(30)
HBAR = 1.054571817e-34; G0 = 6.67430e-11; GH = G0 * HBAR
G_CODATA, U_CODATA = 6.67430, 0.00015          # CITED CODATA 2018/2022, x1e-11, Birge expansion 3.9

def N_spheroid_axis(m):     # demagnetizing factor along the symmetry axis of a spheroid with c/a = m (Osborn 1945)
    if abs(m - 1) < 1e-9: return 1 / 3
    if m > 1: e = sqrt(1 - 1 / m ** 2); return (1 - e ** 2) / e ** 3 * (np.arctanh(e) - e)
    e = sqrt(1 / m ** 2 - 1); return (1 + e ** 2) / e ** 3 * (e - np.arctan(e))
# ---------------------------------------------------------------- data (G, u: CITED CODATA 2022 Table XXX)
# method classes: TB-dyn (time-of-swing), TB-comp (static/servo/compensation, incl. AAF), SUSP (suspended body / pendulum displacement),
# STAT (stationary body, beam balance), AI (atom interferometer)
# test-mass: per suspended BODY (a dumbbell = two bodies); shape 'sphere','cyl(h=d)','plate(t/L)','block'; geo_src CITED/ASSUMED
D = [
 dict(id='NIST-82', ref='Luther & Towler, PRL 48, 121 (1982)', G=6.67248, u=0.00043, method='TB-dyn',
      test='W dumbbell, 2 small spheres, 7 g total', mat='W', rho=19300, m_body=3.5e-3, shape='sphere', aspect=1.0,
      src='2 W spheres 10.5 kg', Ms=10.5, r=0.05, T_int=1e5, geo='ASSUMED'),
 dict(id='TR&D-96', ref='Karagioz & Izmailov, Meas. Tech. 39, 979 (1996)', G=6.6729, u=0.0005, method='TB-dyn',
      test='unknown (not found)', mat='?', rho=np.nan, m_body=np.nan, shape='?', aspect=np.nan,
      src='?', Ms=np.nan, r=np.nan, T_int=1e5, geo='UNKNOWN'),
 dict(id='LANL-97', ref='Bagley & Luther, PRL 78, 3047 (1997)', G=6.67398, u=0.00070, method='TB-dyn',
      test='W dumbbell (NIST-82 apparatus)', mat='W', rho=19300, m_body=3.5e-3, shape='sphere', aspect=1.0,
      src='2 W spheres 10.5 kg', Ms=10.5, r=0.05, T_int=1e5, geo='ASSUMED'),
 dict(id='UWash-00', ref='Gundlach & Merkowitz, PRL 85, 2869 (2000)', G=6.674255, u=0.000092, method='TB-comp',
      test='Pyrex plate 76 x 42 x 1.5 mm', mat='Pyrex', rho=2230, m_body=2230 * 0.076 * 0.042 * 0.0015, shape='plate', aspect=0.0015 / sqrt(0.076 * 0.042),
      src='4 stainless spheres ~8 kg', Ms=8.1, r=0.17, T_int=1e5, geo='ASSUMED'),
 dict(id='BIPM-01', ref='Quinn et al., PRL 87, 111101 (2001)', G=6.67559, u=0.00027, method='TB-comp',
      test='4 Cu-0.7%Te cylinders 1.2 kg, h=d', mat='CuTe', rho=8900, m_body=1.2, shape='cyl(h=d)', aspect=1.0,
      src='4 Cu-Te cylinders 12 kg', Ms=12.0, r=0.2, T_int=1e5, geo='CITED (Quinn/Speake BIPM summary)'),
 dict(id='UWup-02', ref='Kleinevoss et al. (2002), Meas. Sci. Technol. 10, 492 (1999) apparatus', G=6.67422, u=0.00098, method='SUSP',
      test='2 pendulum bobs ~3.2 kg (brass/Cu)', mat='brass', rho=8500, m_body=3.2, shape='cyl(h=d)', aspect=1.0,
      src='2 brass masses 576 kg', Ms=576.0, r=0.4, T_int=1e5, geo='ASSUMED'),
 dict(id='MSL-03', ref='Armstrong & Fitzgerald, PRL 91, 201101 (2003)', G=6.67387, u=0.00027, method='TB-comp',
      test='Cu test mass ~0.5 kg (shape not re-checked)', mat='Cu', rho=8960, m_body=0.5, shape='cyl(h=d)', aspect=1.0,
      src='Cu? source masses ~30 kg', Ms=30.0, r=0.2, T_int=1e5, geo='ASSUMED'),
 dict(id='HUST-05', ref='Hu, Guo & Luo, PRD 71, 127505 (2005)', G=6.67222, u=0.00087, method='TB-dyn',
      test='Cu sphere 32.26 g (dumbbell end)', mat='Cu', rho=8960, m_body=32.26e-3, shape='sphere', aspect=1.0,
      src='2 SS cylinders 6.25 kg, 100x100 mm', Ms=6.25, r=0.08, T_int=1e5, geo='CITED (Luo et al., Phil Trans A 2014 review, PMC4173272)'),
 dict(id='UZur-06', ref='Schlamminger et al., PRD 74, 082001 (2006)', G=6.67425, u=0.00012, method='STAT',
      test='2 Cu cylinders ~1.1 kg', mat='Cu', rho=8960, m_body=1.1, shape='cyl(h=d)', aspect=1.0,
      src='2 Hg tanks, 13.5 t total', Ms=13500.0, r=0.3, T_int=1e5, geo='ASSUMED'),
 dict(id='HUST-09', ref='Luo et al., PRL 102, 240801 (2009); Tu et al., PRD 82, 022001 (2010)', G=6.67349, u=0.00018, method='TB-dyn',
      test='gilded glass/quartz block 91.47 x 12.01 x 26.22 mm, 63.38 g', mat='quartz', rho=63.38e-3 / (0.09147 * 0.01201 * 0.02622),
      m_body=63.38e-3, shape='block', aspect=0.01201 / sqrt(0.09147 * 0.02622),
      src='2 SS316 spheres 778.18 g, 57.15 mm', Ms=0.77818, r=0.05, T_int=1e5, geo='CITED (PMC4173272)'),
 dict(id='BIPM-14', ref='Quinn et al., PRL 111, 101102 (2013); Phil Trans A 372, 20140032 (2014)', G=6.67554, u=0.00016, method='TB-comp',
      test='4 Cu-0.7%Te cylinders 1.2 kg, h=d', mat='CuTe', rho=8900, m_body=1.2, shape='cyl(h=d)', aspect=1.0,
      src='4 Cu-Te cylinders 12 kg', Ms=12.0, r=0.2, T_int=1e5, geo='CITED (as BIPM-01)'),
 dict(id='LENS-14', ref='Rosi et al., Nature 510, 518 (2014); Prevedelli et al., Phil Trans A 372, 20140030 (2014)', G=6.67191, u=0.00099, method='AI',
      test='87Rb atoms (point-like, R << l)', mat='Rb', rho=np.nan, m_body=87 * 1.66054e-27, shape='atom', aspect=np.nan,
      src='W cylinders ~516 kg', Ms=516.0, r=0.15, T_int=1e5, geo='ASSUMED (source); atoms CITED'),
 dict(id='UCI-14', ref='Newman et al., Phil Trans A 372, 20140025 (2014)', G=6.67435, u=0.00013, method='TB-dyn',
      test='fused silica plate 40.59 x 40.59 x 2.943 mm, 10.60 g', mat='SiO2', rho=10.5966e-3 / (0.0405902 * 0.04059 * 0.002943),
      m_body=10.5966e-3, shape='plate', aspect=0.002943 / 0.0405901,
      src='2 OFHC Cu rings 59 kg (OD 520.7, ID 312.42, 48.26 mm)', Ms=59.0, r=0.325, T_int=1e5, geo='CITED (Newman et al. 2014, Table 3)'),
 dict(id='HUSTT-18', ref='Li et al., Nature 560, 582 (2018), time-of-swing', G=6.674184, u=0.000078, method='TB-dyn',
      test='quartz block (taken as HUST-09-like)', mat='quartz', rho=2200, m_body=63e-3, shape='block', aspect=0.01201 / sqrt(0.09147 * 0.02622),
      src='SS spheres ~0.78 kg', Ms=0.78, r=0.05, T_int=1e5, geo='ASSUMED'),
 dict(id='HUSTA-18', ref='Li et al., Nature 560, 582 (2018), angular-acceleration feedback', G=6.674484, u=0.000077, method='TB-comp',
      test='quartz block/plate (taken as HUST-09-like)', mat='quartz', rho=2200, m_body=63e-3, shape='block', aspect=0.01201 / sqrt(0.09147 * 0.02622),
      src='SS spheres (AAF)', Ms=7.8, r=0.15, T_int=1e5, geo='ASSUMED'),
 dict(id='JILA-18', ref='Parks & Faller, PRL 105, 110801 (2010); revised 2019 (CODATA JILA-18)', G=6.67260, u=0.00025, method='SUSP',
      test='2 pendulum bobs ~0.78 kg (Cu)', mat='Cu', rho=8960, m_body=0.78, shape='cyl(h=d)', aspect=1.0,
      src='4 W cylinders ~120 kg', Ms=120.0, r=0.15, T_int=1e5, geo='ASSUMED'),
]
for d in D:
    d['ppm'] = (d['G'] - G_CODATA) / G_CODATA * 1e6; d['u_ppm'] = d['u'] / d['G'] * 1e6
    d['V'] = d['m_body'] / d['rho'] if np.isfinite(d['rho']) else np.nan
    d['V_eff'] = d['V']           # homogeneous bodies: V_eff = M^2/int rho^2 = V (coatings neglected)
    if d['shape'] in ('plate', 'block'):
        Nn = N_spheroid_axis(d['aspect']); d['N_normal'] = Nn; d['N_inplane'] = (1 - Nn) / 2
    elif d['shape'] in ('sphere', 'cyl(h=d)'):
        d['N_normal'] = d['N_inplane'] = 1 / 3      # cylinder h=d: N_axis = 0.311 (close to 1/3); ASSUMED 1/3
    else:
        d['N_normal'] = d['N_inplane'] = np.nan
    # relevant axis for the signal: plates/blocks are driven along their normal (torque from forces normal to the face) - ASSUMED
    d['N_rel'] = d['N_normal']
# weighted mean of the raw data (check against CODATA 6.67430 with the 3.9 Birge factor)
g = np.array([d['G'] for d in D]); u = np.array([d['u'] for d in D])
w = 1 / u ** 2; Gw = np.sum(w * g) / np.sum(w); chi2 = np.sum(((g - Gw) / u) ** 2)
out = dict(stage=30, codata=dict(G=G_CODATA, u=U_CODATA, note='CODATA 2022 Table XXX = CODATA 2018 input data (16 values, no new datum)',
                                 Gw_recomputed=Gw, chi2=chi2, dof=len(g) - 1, birge=sqrt(chi2 / (len(g) - 1)),
                                 spread_max_min_ppm=(g.max() - g.min()) / G_CODATA * 1e6))
# DerSimonian-Laird between-experiment scatter tau (ppm)
y = np.array([d['ppm'] for d in D]); s = np.array([d['u_ppm'] for d in D])
wi = 1 / s ** 2; ybar = np.sum(wi * y) / np.sum(wi); Q = np.sum(wi * (y - ybar) ** 2)
tau2 = max(0.0, (Q - (len(y) - 1)) / (np.sum(wi) - np.sum(wi ** 2) / np.sum(wi)))
out['scatter'] = dict(Q=Q, tau_ppm=sqrt(tau2), sd_unweighted_ppm=float(np.std(y, ddof=1)), median_u_ppm=float(np.median(s)))

# ---------------------------------------------------------------- 2. statistics
def wls(x, y, s, tau=0.0):
    W = 1 / (s ** 2 + tau ** 2); X = np.vstack([np.ones_like(x), x]).T
    C = np.linalg.inv(X.T @ (W[:, None] * X)); b = C @ (X.T @ (W * y)); return b, np.sqrt(np.diag(C))
def perm_p(x, y, s, tau, n=20000):
    b0 = wls(x, y, s, tau)[0][1]; cnt = 0
    for _ in range(n):
        if abs(wls(rng.permutation(x), y, s, tau)[0][1]) >= abs(b0): cnt += 1
    return (cnt + 1) / (n + 1)
def analyse(name, xfun, subset):
    idx = [i for i, d in enumerate(D) if subset(d) and np.isfinite(xfun(d))]
    x = np.array([xfun(D[i]) for i in idx]); yy = y[idx]; ss = s[idx]; z = (g[idx] - G_CODATA) / u[idx]
    # tau re-estimated on the subset
    wi_ = 1 / ss ** 2; yb = np.sum(wi_ * yy) / np.sum(wi_); Q_ = np.sum(wi_ * (yy - yb) ** 2)
    t2 = max(0.0, (Q_ - (len(yy) - 1)) / (np.sum(wi_) - np.sum(wi_ ** 2) / np.sum(wi_))); t = sqrt(t2)
    rs_ppm = stats.spearmanr(x, yy); rs_z = stats.spearmanr(x, z); kt = stats.kendalltau(x, yy)
    b, eb = wls(x, yy, ss, t); bn, ebn = wls(x, yy, ss, 0.0)
    jk = np.array([wls(np.delete(x, k), np.delete(yy, k), np.delete(ss, k), t)[0][1] for k in range(len(x))])
    return dict(variable=name, n=len(x), ids=[D[i]['id'] for i in idx],
                spearman_rho_ppm=rs_ppm.statistic, spearman_p_ppm=rs_ppm.pvalue, spearman_rho_z=rs_z.statistic, spearman_p_z=rs_z.pvalue,
                kendall_tau=kt.statistic, kendall_p=kt.pvalue,
                wls_slope_ppm_per_unit=b[1], wls_slope_err=eb[1], wls_tau_ppm=t, wls_t=b[1] / eb[1],
                wls_p_normal=2 * stats.norm.sf(abs(b[1] / eb[1])), wls_perm_p=perm_p(x, yy, ss, t, 5000),
                wls_no_tau_slope=bn[1], wls_no_tau_t=bn[1] / ebn[1],
                jackknife_slope_range=[float(jk.min()), float(jk.max())], jackknife_sign_changes=int(np.sum(np.sign(jk) != np.sign(b[1]))))
macro = lambda d: d['method'] != 'AI' and d['geo'] != 'UNKNOWN'
documented = lambda d: macro(d) and d['geo'].startswith('CITED')
VARS = dict(log10_Veff=lambda d: log10(d['V_eff']) if np.isfinite(d['V_eff']) else np.nan,
            density=lambda d: d['rho'],
            log10_mass=lambda d: log10(d['m_body']) if np.isfinite(d['m_body']) else np.nan,
            log10_CQ_Gamma_proxy_3N_over_V=lambda d: log10(3 * d['N_rel'] / d['V_eff']) if np.isfinite(d['V_eff']) else np.nan)
res = {k: analyse(k, f, macro) for k, f in VARS.items()}
res_doc = {k: analyse(k, f, documented) for k, f in VARS.items()}
# method: Kruskal-Wallis (ppm) on all 16; plus random-effects group chi2
groups = {}
for d in D: groups.setdefault(d['method'], []).append(d['ppm'])
kw = stats.kruskal(*groups.values())
grp = {}
for k, v in groups.items():
    ii = [i for i, d in enumerate(D) if d['method'] == k]; W_ = 1 / (s[ii] ** 2 + tau2)
    grp[k] = dict(n=len(ii), mean_ppm=float(np.sum(W_ * y[ii]) / np.sum(W_)), err_ppm=float(1 / sqrt(np.sum(W_))), members=[D[i]['id'] for i in ii])
Wall = 1 / (s ** 2 + tau2); mu = np.sum(Wall * y) / np.sum(Wall)
chi2_grp = sum(((v['mean_ppm'] - mu) / v['err_ppm']) ** 2 for v in grp.values()); dof_grp = len(grp) - 1
method = dict(kruskal_H=kw.statistic, kruskal_p=kw.pvalue, groups=grp, re_chi2=chi2_grp, re_dof=dof_grp, re_p=stats.chi2.sf(chi2_grp, dof_grp))
# look-elsewhere: 5 primary tests (4 continuous on macro subset via permutation WLS + method KW)
prim = {k: v['wls_perm_p'] for k, v in res.items()}; prim['method_KW'] = method['kruskal_p']
pmin = min(prim.values()); ntest = len(prim)
lee = dict(primary_p=prim, n_tests=ntest, min_p=pmin, bonferroni=min(1, pmin * ntest), sidak=1 - (1 - pmin) ** ntest,
           note='also tried: Spearman on ppm and z, documented-geometry subset - effective number of looks > 5, so Bonferroni(5) is lenient')
# sensitivity: perturb ASSUMED volumes by lognormal factor 3 (1 sigma), 2000 draws
sens = []
idx = [i for i, d in enumerate(D) if macro(d)]
for _ in range(2000):
    x = np.array([log10(D[i]['V_eff']) + (rng.normal(0, log10(3)) if D[i]['geo'].startswith('ASSUMED') else 0) for i in idx])
    sens.append(stats.spearmanr(x, y[idx]).pvalue)
sens = np.array(sens)
# same-volume pairs (direct test, no statistics)
pairs = [('BIPM-14', 'UZur-06'), ('BIPM-14', 'JILA-18'), ('HUSTT-18', 'HUSTA-18'), ('BIPM-01', 'BIPM-14'), ('UCI-14', 'UWash-00')]
byid = {d['id']: d for d in D}
pair_out = []
for a, b_ in pairs:
    A, B = byid[a], byid[b_]
    pair_out.append(dict(pair=f'{a} vs {b_}', V_ratio=A['V_eff'] / B['V_eff'], dG_ppm=A['ppm'] - B['ppm'],
                         sigma_comb_ppm=sqrt(A['u_ppm'] ** 2 + B['u_ppm'] ** 2), n_sigma=(A['ppm'] - B['ppm']) / sqrt(A['u_ppm'] ** 2 + B['u_ppm'] ** 2)))
# (a) the CQ prediction is excess SCATTER for small V, not a signed shift: test |residual| vs log V
idm = [i for i, d in enumerate(D) if macro(d)]
xv = np.array([log10(D[i]['V_eff']) for i in idm])
abs_test = dict(spearman_absppm=stats.spearmanr(xv, np.abs(y[idm])).__dict__ if False else None)
r1 = stats.spearmanr(xv, np.abs(y[idm])); r2 = stats.spearmanr(xv, np.abs((g[idm] - G_CODATA) / u[idm]))
abs_test = dict(spearman_rho_absppm=r1.statistic, p_absppm=r1.pvalue, spearman_rho_absz=r2.statistic, p_absz=r2.pvalue,
                note='CQ (if anything) predicts rho < 0 here (smaller V -> more scatter)')
# (b) method fixed effects + log V (does the volume trend survive within method classes?)
meths = sorted(set(D[i]['method'] for i in idm))
X = np.column_stack([xv] + [[1.0 if D[i]['method'] == m else 0.0 for i in idm] for m in meths])
Wf = 1 / (s[idm] ** 2 + tau2); C = np.linalg.pinv(X.T @ (Wf[:, None] * X)); bf = C @ (X.T @ (Wf * y[idm]))
within = dict(slope_ppm_per_decade=bf[0], err=sqrt(C[0, 0]), t=bf[0] / sqrt(C[0, 0]), p=2 * stats.norm.sf(abs(bf[0] / sqrt(C[0, 0]))), methods=meths)
# within the two big classes separately
for m in ('TB-dyn', 'TB-comp'):
    ii = [i for i in idm if D[i]['method'] == m]
    within[m] = dict(n=len(ii), spearman=stats.spearmanr([log10(D[i]['V_eff']) for i in ii], y[ii]).statistic,
                     p=stats.spearmanr([log10(D[i]['V_eff']) for i in ii], y[ii]).pvalue)
# association between method and volume (confound check)
mv = {m: float(np.median([log10(D[i]['V_eff']) for i in idm if D[i]['method'] == m])) for m in meths}
out['statistics_extra'] = dict(abs_residual_vs_logV=abs_test, logV_with_method_fixed_effects=within, median_log10V_by_method=mv)
out['statistics'] = dict(macro_subset=res, documented_subset=res_doc, method=method, look_elsewhere=lee,
                         sensitivity_assumed_volumes=dict(draws=len(sens), median_spearman_p=float(np.median(sens)),
                                                          frac_p_below_0p05=float(np.mean(sens < 0.05)), frac_p_below_0p01=float(np.mean(sens < 0.01))),
                         same_volume_pairs=pair_out)

# ---------------------------------------------------------------- 3. physics: what CQ noise can do to a G measurement
s28 = json.load(open(os.path.join(ROOT, 'results/stage28.json')))
GV_MAX = max(v['GammaV_bound'] for v in s28['bounds'].values())       # most permissive Stage 28 variant (95%)
GV_DP = 4 * pi / 3 * GH                                                 # Diosi kappa=1 / CQ kappa=1/2 floor
GAM_ATOM_MAX = GH / (6 * sqrt(pi) * 4.5e-10 ** 3)                       # Stage 26 X-ray ceiling per nucleus/atom (point-like regime)
D2_MAX = 1.169e-15                                                      # Stage 29 matter-sourced ceiling, kg s
REQ = 500e-6
phys = []
for d in D:
    if not np.isfinite(d['Ms']): continue
    a_sig = G0 * d['Ms'] / d['r'] ** 2           # ASSUMED order-of-magnitude signal acceleration on the test body
    T = d['T_int']                                # ASSUMED 1e5 s (1 day) - shorter than any real campaign -> conservative
    if d['method'] == 'AI':
        Tint = 0.16; Nat = 1e5; Nshot = 1e5       # ASSUMED Rosi 2014-like: T = 160 ms, 1e5 atoms/shot, 1e5 shots
        var_a_shot = 2 * GAM_ATOM_MAX / (3 * Tint) / Nat
        sig_a = sqrt(var_a_shot / Nshot); gam = GAM_ATOM_MAX; gv_note = 'point-like: Gamma <= G hbar/(6 sqrt(pi) R0^3) (X-ray, white assumed)'
    else:
        gam = 3 * d['N_rel'] * GV_MAX / d['V_eff']; sig_a = sqrt(gam / T); gv_note = 'Gamma = 3 N_rel (Gamma V)_max / V_eff'
    gam_dp = 3 * d['N_rel'] * GV_DP / d['V_eff'] if d['method'] != 'AI' else GH / (6 * sqrt(pi) * 1e-9 ** 3)
    gam_src = G0 ** 2 * D2_MAX * d['Ms'] / (4 * d['r'] ** 4)
    # rectification bound for time-of-swing: heating dE/dt = M Gamma/2, steady extra energy M Gamma Q/(2 w0) (ASSUMED Q=1e5, period 500 s,
    # amplitude A = 1 mrad, r_g = 3 cm, amplitude coefficient of the period c_A = 1 per rad^2, signal frequency shift 1e-3): all deliberately generous
    rect = np.nan
    if d['method'] == 'TB-dyn':
        Qf, w0, A, rg, cA, fsig = 1e5, 2 * pi / 500, 1e-3, 0.03, 1.0, 1e-3
        dA2 = gam * Qf / (w0 ** 3 * rg ** 2)           # = 2 dE/(I w0^2)
        rect = 2 * cA * dA2 / fsig
    phys.append(dict(id=d['id'], a_signal=a_sig, Gamma_max=gam, Gamma_DP_floor=gam_dp, Gamma_matter_sourced_max=gam_src, note=gv_note,
                     scatter_ppm_at_max=sig_a / a_sig * 1e6, scatter_ppm_DP_floor=(sqrt(gam_dp / T) / a_sig * 1e6) if d['method'] != 'AI' else np.nan,
                     matter_sourced_ppm=sqrt(gam_src / T) / a_sig * 1e6, rectification_bias_ppm=(rect * 1e6) if np.isfinite(rect) else np.nan,
                     orders_below_500ppm=log10(REQ * 1e6 / (sig_a / a_sig * 1e6))))
need = []
for d in D:
    if d['method'] == 'AI' or not np.isfinite(d['Ms']): continue
    a_sig = G0 * d['Ms'] / d['r'] ** 2; gam_need = (REQ * a_sig) ** 2 * d['T_int']
    need.append(dict(id=d['id'], Gamma_needed=gam_need, GammaV_needed=gam_need * d['V_eff'] / (3 * d['N_rel']),
                     orders_above_stage28=log10(gam_need * d['V_eff'] / (3 * d['N_rel']) / GV_MAX)))
out['physics'] = dict(GammaV_max_stage28=GV_MAX, GammaV_DP_floor=GV_DP, Gamma_atom_max=GAM_ATOM_MAX, Dhat_max=D2_MAX,
                      per_experiment=phys, needed_for_500ppm=need,
                      max_scatter_ppm=max(p['scatter_ppm_at_max'] for p in phys),
                      max_rectification_ppm=float(np.nanmax([p['rectification_bias_ppm'] for p in phys])),
                      mean_bias_mechanisms=[
                          'Linear response: white zero-mean force noise adds variance only; <F> unchanged. CQ models reproduce the Newtonian mean '
                          '(Tilloy-Diosi 2017 Eq. 10: unitary pair potential + double commutator; Kryhin-Sudhir Eq. 5: averaged Poisson equation sourced by <f>).',
                          'DP/CQ self-noise is independent of source-mass position (Stage 29), so any heating-induced shift is identical in both source '
                          'configurations and cancels in the near-far difference that defines G.',
                          'Second-order (rectification) via amplitude-dependent torsion period: bounded above (column rectification_bias_ppm) with generous inputs.',
                          'Matter-sourced (discrete) kernel: noise near the source mass differs between configurations, but enters as variance; bounded above '
                          '(matter_sourced_ppm).'])
out['oppenheim_prior'] = ('arXiv:2203.01982 Sec. V already floated this: "Since we have found that CQ theories predict an uncertainty in mass measurements '
                          'it is perhaps intriguing that different experiments to measure Newton\'s constant G yield results whose relative uncertainty differ '
                          'by as much as 5e-4 ... masses with smaller volume should yield larger uncertainty and this would be the effect to look for".')
# CSV
with open(os.path.join(ROOT, 'results/stage30_G_table.csv'), 'w', newline='') as fh:
    cols = ['id', 'ref', 'method', 'G', 'u', 'ppm', 'u_ppm', 'test', 'mat', 'rho', 'm_body', 'V_eff', 'shape', 'aspect', 'N_normal', 'N_inplane',
            'src', 'Ms', 'r', 'geo']
    wr = csv.writer(fh); wr.writerow(cols)
    for d in D: wr.writerow([d[c] for c in cols])
out['table'] = D
out['runtime_s'] = time.time() - t0
json.dump(out, open(os.path.join(ROOT, 'results/stage30.json'), 'w'), indent=1, default=float)
pr = lambda o: print(json.dumps(o, indent=1, default=float))
pr(out['codata']); pr(out['scatter'])
for k, v in res.items(): print('MACRO', k, {kk: v[kk] for kk in ('n', 'spearman_rho_ppm', 'spearman_p_ppm', 'spearman_p_z', 'wls_slope_ppm_per_unit', 'wls_slope_err', 'wls_tau_ppm', 'wls_perm_p', 'jackknife_slope_range', 'wls_no_tau_t')})
for k, v in res_doc.items(): print('DOC', k, {kk: v[kk] for kk in ('n', 'spearman_rho_ppm', 'spearman_p_ppm', 'wls_perm_p')})
pr(out['statistics_extra']); pr(method); pr(lee); pr(out['statistics']['sensitivity_assumed_volumes']); pr(pair_out)
for p in phys: print(p['id'], '%.2e ppm scatter(max)  DP %.2e  src %.2e  rect %.2e  orders %.1f' % (p['scatter_ppm_at_max'], p['scatter_ppm_DP_floor'], p['matter_sourced_ppm'], p['rectification_bias_ppm'], p['orders_below_500ppm']))
pr(need)
