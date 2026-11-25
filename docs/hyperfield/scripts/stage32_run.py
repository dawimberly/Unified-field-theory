#!/usr/bin/env python3
"""Stage 32 -- forecast: can a joint (emission angle x bond length) split of the Grundmann et al. 2020 H2 data tell
H0 (birth delay scales with R: travel time R cos(beta)/c + nondipole) from H1 (H0 + do-si-do release term d1*cos(beta), d1 R-independent)?

All inputs were fetched this session (see results/stage32_sources.csv). Labels: CITED, DIGITIZED, DERIVED, ASSUMED.
Run:  python scripts/stage32_run.py      (writes results/stage32.json, results/stage32_forecast.csv, results/stage32_curves.csv,
                                          results/stage32_sources.csv, figures/stage32_split_forecast.png)
"""
import os, sys, json, csv, hashlib
import numpy as np
from math import sqrt
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
import stage32_analyze_joint as J

D31 = os.path.join(ROOT, 'data', 'stage31'); D32 = os.path.join(ROOT, 'data', 'stage32'); RES = os.path.join(ROOT, 'results'); FIG = os.path.join(ROOT, 'figures')
S = {}
def src(key, ref, url, loc, quote, label='CITED', **vals):
    S[key] = dict(ref=ref, url=url, location=loc, quote=quote, label=label, values=vals)

def md5(p): return hashlib.md5(open(p, 'rb').read()).hexdigest()

# ------------------------------------------------------------------------------------------------ sources (all fetched 2026-11-24 CT)
src('grundmann_text', 'Grundmann et al., Science 370, 339 (2020), author version arXiv:2010.08298v1 (text saved as data/stage32/grundmann2020_arXiv2010.08298v1_text.txt)',
    'https://arxiv.org/pdf/2010.08298', 'main text p.5 Eqs.(1)-(2); p.5 KER relation; p.6 subset S',
    '"tau_b = cos(alpha_0) R/v_ph" (Eq. 2) ; "propagates with the phase velocity v_ph = E_e/p_e" ; "using left-handed circularly polarized photons with an energy of 800 eV" ; '
    '"In the reflection approximation, the internuclear distance R at the moment of photoabsorption is related to the kinetic energy release via KER = e^2/(4 pi eps0 R)" ; '
    '"we inspected the interference pattern of the fast electron for different angles beta ... in Fig. 3 for the subset S (see Fig. 2 B)"',
    Eph_eV=800.0)
src('grundmann_fig3_caption', 'Grundmann et al. 2020 (arXiv:2010.08298v1)', 'https://arxiv.org/pdf/2010.08298', 'p.11 Fig. 3 caption; p.10 Fig. 2 caption',
    '"fast electrons [Ee = (735+-15) eV] ... for the average internuclear distance of R = (0.74 +- 0.02) A [selected subset S as shown in Fig. 2 B]" ; '
    '"The maxima are obtained using Gaussian fits" ; Fig.2: "The subset S of the data is used for panel A and for the subsequent analysis of the birth time delay."',
    Ee_eV=735.0, dEe_eV=15.0, R_S_A=0.74, dR_S_A=0.02)
src('grundmann_eq3', 'Grundmann et al. 2020 (arXiv:2010.08298v1)', 'https://arxiv.org/pdf/2010.08298', 'p.7 Eq.(3) and following',
    '"tau_b\' = cos(beta) R/c" ; "For this case the expected birth time delay was +-247 zs for R = 0.74 A" ; red line: "the atomic nondipole effect tilts the electron angular '
    'distributions from each center slightly in forward direction ... This fact led to an additional angular shift of the interference pattern that slightly increased the slope" ; '
    '"more theoretical work including electron-electron correlation is needed to clarify the deviation"', tau_Rc_S_zs=247.0)
src('grundmann_sm', 'Grundmann et al. 2020, Supplementary Materials (in arXiv:2010.08298v1)', 'https://arxiv.org/pdf/2010.08298', 'SM Materials and Methods p.13-14',
    '"The roughly 300.000 ion-ion-electron coincidences displayed in Fig. 2 have been obtained during a total acquisition time of approximately 132 h." ; '
    '"The accuracy to which we could specify the zero point of the momentum distribution of the fast electron was approximately +-0.04 atomic units in all three directions" ; '
    '"we determined the maximum and minimum slope of the linear fit function in Fig. 3 C (including +-1 sigma estimation error) and indicated this range by the purple-shaded area" ; '
    '"The statistical errors ... are the standard deviations of the mean values from the Gaussian fits estimated as the square root of the respective diagonal element of the covariance matrix"',
    N_coinc=300000, t_h=132, dp0_au=0.04)
src('grundmann_energy', 'Grundmann et al. 2020 (arXiv:2010.08298v1)', 'https://arxiv.org/pdf/2010.08298', 'main text p.6',
    '"The fraction of the photon energy that exceeds the sum of the adiabatic double ionization energy of H2 (31.03 eV) and the KER is shared between the two electrons" ; '
    '"we restricted our investigation to fast electrons that carried more than 96% of the excess energy"', IP2_eV=31.03, frac_min=0.96)
src('grundmann_note19', 'Grundmann et al. 2020 (arXiv:2010.08298v1)', 'https://arxiv.org/pdf/2010.08298', 'reference note 19',
    '"This possible additional delay of the double ionization process is symmetric with respect to inversion of the two centers of the molecule. Therefore, our interferometric technique is blind to this absolute delay."')
src('zenodo', 'Grundmann et al. data, Zenodo 10.5281/zenodo.3899920 (record JSON saved as data/stage32/zen.json)', 'https://zenodo.org/api/records/3899920',
    'files Fig2_data.txt, Fig3_data.txt', 'title "Zeptosecond Birth Time Delay in Molecular Photoionization"; md5 Fig2 dc8986f4cc4529650f8932b8c7583e59, Fig3 f367578243fadce987bb1b48d3c74147 '
    '(re-downloaded this session; identical to data/stage31 copies)', md5_fig2='dc8986f4cc4529650f8932b8c7583e59', md5_fig3='f367578243fadce987bb1b48d3c74147')
src('fig3c_digitized', 'Grundmann et al. 2020 Fig. 3C, re-rendered at 600 dpi this session (data/stage32/grundmann2020_fig3C_600dpi_crop.png)', 'https://arxiv.org/pdf/2010.08298',
    'p.11 Fig. 3C; axis calibrated on right-axis ticks 0/200/400 zs and frame (cos beta = +-1)',
    'DIGITIZED (pixel colour tracing, +-2 zs): blue line slope 246.9 zs (check: CITED 247, Eq.3 246.8); red line slope 271.7 zs (at cos beta = 0.966 and -0.940: 271.6, 271.9); '
    'dashed purple fit 445.8 zs; purple band 326-568 zs (half-width 121)', label='DIGITIZED', blue_zs=246.9, red_zs=271.7, paper_fit_zs=445.8, band_lo_zs=326.0, band_hi_zs=568.0)
src('codata2022', 'CODATA 2022 complete listing (NIST), saved as data/stage32/codata.txt', 'https://physics.nist.gov/cuu/Constants/Table/allascii.txt', 'lines for c, a0, E_h, atomic unit of time, e, eps0',
    '"speed of light in vacuum 299 792 458 (exact)"; "Bohr radius 5.291 772 105 44 e-11"; "Hartree energy in eV 27.211 386 245 981"; "atomic unit of time 2.418 884 326 5864 e-17"; '
    '"elementary charge 1.602 176 634 e-19 (exact)"; "vacuum electric permittivity 8.854 187 8188 e-12"',
    c=J.C, a0=J.BOHR, Eh_eV=J.HARTREE_EV, t_au=J.AU_T, e=J.E_CH, eps0=J.EPS0)
src('stage31', 'Stage 31 of this repo (DERIVED there from the Zenodo Fig3 file; recomputed here)', 'docs/hyperfield/results/stage31.json', 'grundmann.slope_zs etc.',
    'slope 504.5 +- 36.2 zs (stat, window 0.4); residual over red line +235 +- 128 zs', label='DERIVED (repo)')

# consistency checks of the fetched files
assert md5(os.path.join(D31, 'grundmann2020_Fig2_data.txt')) == S['zenodo']['values']['md5_fig2']
assert md5(os.path.join(D31, 'grundmann2020_Fig3_data.txt')) == S['zenodo']['values']['md5_fig3']
_txt = open(os.path.join(D32, 'grundmann2020_arXiv2010.08298v1_text.txt'), encoding='utf-8', errors='ignore').read()
for frag in ['300.000 ion-ion-electron coincidences', '0.04 atomic units', 'subset S', '31.03 eV', '(735', 'blind to this absolute delay']:
    assert frag in _txt, frag
_cod = open(os.path.join(D32, 'codata.txt')).read()
for frag in ['299 792 458', '5.291 772 105 44 e-11', '27.211 386 245 981', '2.418 884 326 5864 e-17', '8.854 187 8188 e-12', '2022 CODATA']:
    assert frag in _cod, frag

# ------------------------------------------------------------------------------------------------ data
def load2(path, colx, coly):
    d = np.loadtxt(path, skiprows=2, delimiter=',')
    xs = np.unique(d[:, colx]); ys = np.unique(d[:, coly]); H = np.zeros((len(xs), len(ys)))
    for r in d: H[np.searchsorted(xs, r[colx]), np.searchsorted(ys, r[coly])] = r[4]
    return xs, ys, H
ca3, cb3, H3 = load2(os.path.join(D31, 'grundmann2020_Fig3_data.txt'), 2, 3)     # H3[cos a, cos b], subset S
ca2, K2, H2 = load2(os.path.join(D31, 'grundmann2020_Fig2_data.txt'), 2, 3)      # H2[cos a, KER], beta-averaged, all R
R2 = J.KE_EVA / K2; N2k = H2.sum(0); NS3 = H3.sum()
DCOS3 = float(np.mean(np.diff(ca3)))

# ------------------------------------------------------------------------------------------------ (1) reproduce Stage 31 on subset S, Fisher vs bootstrap
KER_S = J.KE_EVA / J.R_S
Ee_S = float(J.Ee_of_KER(KER_S)); k_S = float(J.zs_per_cos(J.R_S, J.EE0))
rows_S = []
for j in range(len(cb3)):
    mu, er, par = J.fit_row(ca3, H3[:, j], J.W_S); rows_S.append((mu, er, par))
muS = np.array([r[0] for r in rows_S]); erS = np.array([r[1] for r in rows_S])
bS, dbS, _, _ = J.wline(cb3, muS * k_S, erS * k_S)
st31 = json.load(open(os.path.join(RES, 'stage31.json')))['grundmann']
boot31 = np.array(st31['windows']['0.4']['dtau_zs'])
fisher_over_boot = float(np.median(erS * k_S / boot31))
tauRc_S = float(J.tau_Rc(J.R_S)); dg = S['fig3c_digitized']['values']
d1_S31 = bS - dg['red_zs']                                   # residual over the paper's red line (DERIVED)
reproduce = dict(slope_zs=bS, slope_stat_fisher_zs=dbS, stage31_slope_zs=st31['slope_zs'], stage31_slope_stat_bootstrap_zs=st31['slope_stat_zs'],
                 fisher_over_bootstrap_rowerr_median=fisher_over_boot, residual_over_red_zs=d1_S31, tauRc_S_zs=tauRc_S, nd_S_zs=dg['red_zs'] - dg['blue_zs'],
                 counts_subset_S_fig3=float(NS3))

# ------------------------------------------------------------------------------------------------ (2) R-dependence of the fringe from Fig. 2
def fig2_group(lo, hi):
    m = (R2 >= lo) & (R2 < hi); y = H2[:, m].sum(1); N = y.sum(); Rm = float((N2k[m] * R2[m]).sum() / N); Km = float((N2k[m] * K2[m]).sum() / N)
    mu, er, par = J.fit_row(ca2, y, J.W_S * J.R_S / Rm, 0.107 * J.R_S / Rm)
    A, _, sg, B = par
    return dict(N=float(N), R=Rm, KER=Km, A_per_N=A / N, B_per_N=B / N, sigma=sg, mu=mu, mu_err=er)
G_S = fig2_group(0.72, 0.76)                                      # Fig. 2 slice matching subset S
fine_edges = [0.60, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.83, 0.87, 0.96]
fine = [fig2_group(lo, hi) for lo, hi in zip(fine_edges[:-1], fine_edges[1:])]
fR = np.array([g['R'] for g in fine])
# smooth (linear in R) fits of the shape ratios vs subset S -- ASSUMED smooth; used for every binning
def lfit(y): c = np.polyfit(fR - 0.75, y, 1, w=np.array([sqrt(g['N']) for g in fine])); return lambda R: np.polyval(c, np.asarray(R) - 0.75)
rA = lfit([g['A_per_N'] / G_S['A_per_N'] for g in fine]); rB = lfit([g['B_per_N'] / G_S['B_per_N'] for g in fine])
rSR = lfit([g['sigma'] * g['R'] / (G_S['sigma'] * G_S['R']) for g in fine])     # sigma*R relative to S (1 = pure 1/R scaling)
def r_sigma(R, pe_ratio=1.0): return rSR(R) * (G_S['R'] / np.asarray(R)) / pe_ratio

# subset-S row template from Fig. 3 (per cos beta row): A_i, sigma_i, B_i
TPL = [(par[0], par[2], par[3]) for (_, _, par) in rows_S]

def pe_of(Ee): return np.sqrt(2 * Ee / J.HARTREE_EV)

def lam_rows(R, N, Eph=800.0, s_zs=None, Ee=None):
    """expected counts lam[ib, ia] for an R bin with N events (Fig. 2 normalisation), centred fringe shift from slope s_zs (zs)."""
    if Ee is None: Ee = float(J.Ee_of_KER(J.KE_EVA / R, Eph=Eph))
    pr = pe_of(Ee) / pe_of(Ee_S)
    sc = float(r_sigma(R, pr)); k = float(J.zs_per_cos(R, Ee))
    out = np.zeros((len(cb3), len(ca3)))
    for i, (A, sg, B) in enumerate(TPL):
        mu = 0.0 if s_zs is None else s_zs * cb3[i] / k
        out[i] = (N / G_S['N']) * (float(rA(R)) * A * np.exp(-0.5 * ((ca3 - mu) / (sg * sc)) ** 2) + float(rB(R)) * B)
    return out, sc, k

def fisher_slope_err(R, N, Eph=800.0, Ee=None):
    lam, sc, k = lam_rows(R, N, Eph, None, Ee)
    win = J.W_S * sc; m = np.abs(ca3) <= win; e = []
    for i, (A, sg, B) in enumerate(TPL):
        e.append(J.fisher_mu_err(float(rA(R)) * A * N / G_S['N'], 0.0, sg * sc, float(rB(R)) * B * N / G_S['N'], ca3[m]) * k)
    _, db, _, _ = J.wline(cb3, np.zeros_like(cb3), np.array(e))
    return db

def groups(edges):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (R2 >= lo) & (R2 < hi); N = N2k[m].sum()
        out.append(dict(lo=lo, hi=hi, N=float(N), R=float((N2k[m] * R2[m]).sum() / N), KER=float((N2k[m] * K2[m]).sum() / N)))
    return out

def equal_count_edges(n, lo=0.60, hi=0.96):
    o = np.argsort(R2); cw = np.cumsum(N2k[o]) / N2k.sum()
    return [lo] + [float(R2[o][np.searchsorted(cw, q)]) for q in np.arange(1, n) / n] + [hi]

def forecast(edges, Eph=800.0, nscale=1.0, Ee_fixed=None):
    gs = groups(edges)
    for g in gs:
        Ee = Ee_fixed if Ee_fixed is not None else float(J.Ee_of_KER(g['KER'], Eph=Eph))
        g['Ee'] = Ee; g['err'] = float(fisher_slope_err(g['R'], g['N'] * nscale, Eph, Ee))
    return gs

def sig_d1(R, e, mode='B', sig_delta=None):
    R = np.asarray(R); e = np.asarray(e)
    if mode == 'B': X = np.vstack([R, np.ones_like(R)]).T; idx = 1; prior = None
    else: X = np.vstack([np.ones_like(R), R]).T; idx = 0; prior = [(1, (J.SIG_DELTA_DEFAULT if sig_delta is None else sig_delta) / J.R_S)]
    P = X.T @ np.diag(1 / e ** 2) @ X
    if prior: P[prior[0][0], prior[0][0]] += 1 / prior[0][1] ** 2
    return float(sqrt(np.linalg.inv(P)[idx, idx]))

# calibration: forecast for subset S alone must reproduce the measured stat error of the slope
cal = forecast([0.72, 0.76]); cal_ratio = (cal[0]['err'] * sqrt(cal[0]['N'] / NS3)) / dbS   # Fig3 has NS3 events (vs Fig2 slice count)
NSCALE = NS3 / G_S['N']      # DERIVED: Fig. 3 subset S has fewer events than the Fig. 2 slice (extra cuts); ASSUMED same efficiency at all R

# ------------------------------------------------------------------------------------------------ (3) headline forecast
DEF_N = 4
E_def = equal_count_edges(DEF_N); g_def = forecast(E_def, nscale=NSCALE)
Rd = np.array([g['R'] for g in g_def]); ed = np.array([g['err'] for g in g_def])
sB = sig_d1(Rd, ed, 'B'); sA = sig_d1(Rd, ed, 'A'); sA_band = sig_d1(Rd, ed, 'A', sig_delta=0.5 * (dg['band_hi_zs'] - dg['band_lo_zs']))
all_one = forecast([0.60, 0.96], nscale=NSCALE)[0]
head = dict(nbins=DEF_N, edges_A=E_def, R_A=Rd.tolist(), N=[g['N'] * NSCALE for g in g_def], slope_err_zs=ed.tolist(),
            sigma_d1_B_zs=sB, sigma_d1_A_zs=sA, sigma_d1_A_authorsband_zs=sA_band,
            dmin_B_2sig_zs=2 * sB, dmin_B_5sig_zs=5 * sB, dmin_B_2sig_95power_zs=(2 + 1.645) * sB, dmin_A_2sig_zs=2 * sA, dmin_A_5sig_zs=5 * sA,
            Z_B_if_stage31_real=d1_S31 / sB, Z_A_if_stage31_real=d1_S31 / sA, Z_A_band_if_stage31_real=d1_S31 / sA_band,
            Z_B_at_stage31_95bound_485=485.0 / sB, all_R_single_bin_slope_err_zs=all_one['err'], subsetS_slope_err_zs=dbS)

# ------------------------------------------------------------------------------------------------ (4) binning scan
bscan = []
for n in range(2, 11):
    gs = forecast(equal_count_edges(n), nscale=NSCALE); R = np.array([g['R'] for g in gs]); e = np.array([g['err'] for g in gs])
    bscan.append(dict(nbins=n, sigma_d1_B_zs=sig_d1(R, e, 'B'), sigma_d1_A_zs=sig_d1(R, e, 'A')))
for name, ed_ in [('paper-like 0.04 A slices', [0.60, 0.64, 0.68, 0.72, 0.76, 0.80, 0.84, 0.88, 0.96]),
                  ('outer thirds only', None)]:
    if ed_ is None:
        e3 = equal_count_edges(3); gs = forecast(e3, nscale=NSCALE); gs = [gs[0], gs[2]]
    else: gs = forecast(ed_, nscale=NSCALE)
    R = np.array([g['R'] for g in gs]); e = np.array([g['err'] for g in gs])
    bscan.append(dict(nbins=name, sigma_d1_B_zs=sig_d1(R, e, 'B'), sigma_d1_A_zs=sig_d1(R, e, 'A')))
best = min([b for b in bscan], key=lambda b: b['sigma_d1_B_zs'])

# ------------------------------------------------------------------------------------------------ (5) second photon energy (ASSUMED: same event count, same visibility/envelope, fringe width ~ 1/p_e)
def combo(E2_list, common_dp=False):
    sets = [(800.0, forecast(E_def, nscale=NSCALE))] + [(E, forecast(E_def, Eph=E, nscale=NSCALE)) for E in E2_list]
    rows, ys = [], []
    nE = len(sets)
    for iE, (E, gs) in enumerate(sets):
        for g in gs:
            if not common_dp:   # (B): free a_E per energy, common d1
                x = np.zeros(nE + 1); x[iE] = g['R']; x[-1] = 1.0
            else:               # (A2): travel time fixed, common momentum offset dp (free, no prior): tau = dp*R/E_e, common d1
                x = np.array([1.0, g['R'] / g['Ee']])
            rows.append(x); ys.append(g['err'])
    X = np.array(rows); e = np.array(ys); P = X.T @ np.diag(1 / e ** 2) @ X
    idx = -1 if not common_dp else 0
    return float(sqrt(np.linalg.inv(P)[idx, idx]))
E2 = {}
for lab, lst in [('800 repeat (x2 events)', [800.0]), ('400', [400.0]), ('1600', [1600.0]), ('400+1600', [400.0, 1600.0])]:
    E2[lab] = dict(sigma_d1_B_zs=combo(lst), sigma_d1_A2_commondp_zs=(combo(lst, True) if lab != '800 repeat (x2 events)' else None))
E2['800 only'] = dict(sigma_d1_B_zs=sB, sigma_d1_A2_commondp_zs=None)
# energy dependence of the per-energy slope error (single all-R bin) for the record
E2_detail = {str(E): forecast([0.60, 0.96], Eph=E, nscale=NSCALE)[0] for E in (400.0, 800.0, 1600.0)}

# ------------------------------------------------------------------------------------------------ (6) systematics of the shape test (B)
def nd_bias(true_s, fit_s):
    R = Rd; e = ed; s = J.tau_Rc(R) + J.nd_term(R, true_s)
    return J.hypothesis_tests(R, s, e, nd_scaling=fit_s)['B_d1_zs']
ndb = {f'true{t:+.0f}_fit{f:+.0f}': nd_bias(t, f) for t in (-1.0, 0.0, 1.0) for f in (-1.0, 0.0, 1.0)}
# R resolution (ASSUMED sigma_R = 0.03 A, not given in the paper): mean true R per measured bin regresses to the mean
def smeared_R(edges, sR):
    o = np.argsort(R2); Rt = R2[o]; w = N2k[o]
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        from math import erf
        p = np.array([0.5 * (erf((hi - r) / (sR * sqrt(2))) - erf((lo - r) / (sR * sqrt(2)))) for r in Rt]) * w
        out.append(float((p * Rt).sum() / p.sum()))
    return np.array(out)
Rsm = smeared_R(E_def, 0.03)
s_h1 = J.tau_Rc(Rsm) + J.nd_term(Rsm) + d1_S31
t_sm = J.hypothesis_tests(Rd, s_h1, ed)        # analyst uses measured-bin R, truth used smeared R
Rres = dict(sigma_R_A=0.03, R_true_means=Rsm.tolist(), R_measured_means=Rd.tolist(), d1_recovered_zs=t_sm['B_d1_zs'], d1_true_zs=d1_S31,
            lever_shrink=float(np.std(Rsm) / np.std(Rd)))

A_FREE = (bS - J.nd_term(J.R_S, -1.0)) / J.R_S     # DERIVED: "everything scales with R", calibrated to the subset-S slope
s_h0f = A_FREE * Rsm + J.nd_term(Rsm)
Rres['d1_recovered_if_H0free_truth_zs'] = J.hypothesis_tests(Rd, s_h0f, ed)['B_d1_zs']
s_h0 = J.tau_Rc(Rsm) + J.nd_term(Rsm)
Rres['d1_recovered_if_H0_truth_zs'] = J.hypothesis_tests(Rd, s_h0, ed)['B_d1_zs']

# ------------------------------------------------------------------------------------------------ (7) Monte Carlo validation of the Fisher forecast (full pipeline: fits of simulated histograms)
def simulate_counts(edges, d1, rng, nscale=NSCALE, nd_s=-1.0, truth=None):
    """Poisson histogram, built at fine KER resolution (80 Fig. 2 KER bins), each with its own R, shift and fringe width; summed into R bins."""
    nR = len(edges) - 1; Hc = np.zeros((nR, len(cb3), len(ca3))); Rs = np.zeros(nR); Ks = np.zeros(nR); W = np.zeros(nR)
    for k, R in enumerate(R2):
        if not (edges[0] <= R < edges[-1]) or N2k[k] <= 0: continue
        j = int(np.searchsorted(edges, R, side='right') - 1)
        s = float(J.tau_Rc(R) + J.nd_term(R, nd_s) + d1) if truth is None else float(truth(R))
        lam, _, _ = lam_rows(R, N2k[k] * nscale, s_zs=s)
        Hc[j] += rng.poisson(lam); Rs[j] += N2k[k] * R; Ks[j] += N2k[k] * K2[k]; W[j] += N2k[k]
    return Hc, Rs / W, Ks / W

def one_mc(args):
    seed, d1 = args; rng = np.random.default_rng(seed)
    if d1 == 'H0free':
        Hc, Rc, Kc = simulate_counts(E_def, 0.0, rng, truth=lambda R: A_FREE * R + J.nd_term(R, -1.0))
    else:
        Hc, Rc, Kc = simulate_counts(E_def, d1, rng)
    res = J.analyse_hist(Hc, ca3, cb3, Rc, Kc, sigR_scale=[float(r_sigma(R)) for R in Rc])
    t = J.hypothesis_tests(Rc, res['slope_zs'], res['slope_err_zs'])
    return [t['B_d1_zs'], t['B_d1_err'], t['A_d1_zs'], t['A_d1_err']] + res['slope_zs']

NMC = int(os.environ.get('S32_NMC', '200'))
if __name__ == '__main__':
    with Pool(min(8, os.cpu_count() or 2)) as pool:
        mc0 = np.array(pool.map(one_mc, [(1000 + i, 0.0) for i in range(NMC)]))
        mc1 = np.array(pool.map(one_mc, [(5000 + i, d1_S31) for i in range(NMC)]))
        mcf = np.array(pool.map(one_mc, [(9000 + i, 'H0free') for i in range(NMC)]))
    def summ(m, d1):
        return dict(d1_true_zs=d1, B_mean_zs=float(m[:, 0].mean()), B_std_zs=float(m[:, 0].std()), B_mean_reported_err_zs=float(m[:, 1].mean()),
                    A_mean_zs=float(m[:, 2].mean()), A_std_zs=float(m[:, 2].std()), slope_means_zs=m[:, 4:].mean(0).tolist(), slope_stds_zs=m[:, 4:].std(0).tolist(),
                    frac_B_Z_gt_2=float(np.mean(m[:, 0] / m[:, 1] > 2)), frac_B_Z_gt_1=float(np.mean(m[:, 0] / m[:, 1] > 1)))
    MC = dict(n_per_hypothesis=NMC, H0=summ(mc0, 0.0), H1_stage31=summ(mc1, d1_S31), H0free_cal=summ(mcf, 0.0), fisher_sigma_d1_B_zs=sB, fisher_slope_err_zs=ed.tolist())

    # -------------------------------------------------------------------------------------------- (8) predicted curves + d1 scan
    Rg = np.round(np.arange(0.60, 0.961, 0.01), 3)
    h0 = J.tau_Rc(Rg) + J.nd_term(Rg, -1.0)
    a_free = (bS - J.nd_term(J.R_S, -1.0)) / J.R_S                       # H0free calibrated to the Stage 31 / subset S slope
    curves = []
    for R, v0 in zip(Rg, h0):
        k = float(J.zs_per_cos(R, float(J.Ee_of_KER(J.KE_EVA / R))))
        row = dict(R_A=R, KER_eV=J.KE_EVA / R, tauRc_zs=float(J.tau_Rc(R)), H0_zs=float(v0), H0_nd0_zs=float(J.tau_Rc(R) + J.nd_term(R, 0)),
                   H0_nd1_zs=float(J.tau_Rc(R) + J.nd_term(R, 1)), H0free_cal_zs=float(a_free * R + J.nd_term(R, -1.0)),
                   H1_stage31_zs=float(v0 + d1_S31), H1_bound485_zs=float(v0 + 485.0), zs_per_cos=k)
        for key in ['H0_zs', 'H0free_cal_zs', 'H1_stage31_zs', 'H1_bound485_zs']:
            row[key.replace('_zs', '_cosalpha_slope')] = row[key] / k
        curves.append(row)
    d1grid = np.arange(0, 1001, 25.0)
    scan = [dict(d1_zs=float(d), Z_B=float(d / sB), Z_A=float(d / sA), Z_A_authorsband=float(d / sA_band),
                 Z_B_with_400_1600=float(d / E2['400+1600']['sigma_d1_B_zs']), Z_A2_with_400_1600=float(d / E2['400+1600']['sigma_d1_A2_commondp_zs'])) for d in d1grid]

    out = dict(stage=32, title='Forecast: joint angle x bond-length split of the Grundmann 2020 H2 birth-delay data (H0 scales with R vs H1 do-si-do R-independent)',
               sources=S, reproduce_stage31=reproduce, calibration=dict(subsetS_forecast_over_measured_stat=cal_ratio, NSCALE_fig3_over_fig2slice=NSCALE,
               fig2_subsetS_slice=G_S, fig2_fine_groups=fine), headline=head, binning_scan=bscan, best_binning=best, second_energy=E2,
               second_energy_single_bin=E2_detail, nd_scaling_bias_on_B_d1_zs=ndb, R_resolution_variant=Rres, monte_carlo=MC,
               H1_stage31_d1_zs=d1_S31, H0free_cal_a_zs_per_A=a_free, H0free_cal_a_over_inv_c=float(a_free / J.tau_Rc(1.0)))
    os.makedirs(RES, exist_ok=True)
    json.dump(out, open(os.path.join(RES, 'stage32.json'), 'w'), indent=1, default=float)
    with open(os.path.join(RES, 'stage32_sources.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['key', 'label', 'ref', 'url', 'location', 'quote', 'values'])
        for k, v in S.items(): w.writerow([k, v['label'], v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
    with open(os.path.join(RES, 'stage32_curves.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(curves[0].keys())); w.writeheader(); [w.writerow(r) for r in curves]
    with open(os.path.join(RES, 'stage32_forecast.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['quantity', 'value', 'unit', 'label', 'note'])
        rows = [('slope_subsetS_refit', bS, 'zs', 'DERIVED', 'Stage 31 method, window 0.4, Fisher stat error below'),
                ('slope_subsetS_stat_fisher', dbS, 'zs', 'DERIVED', f'Stage 31 bootstrap {st31["slope_stat_zs"]:.1f}'),
                ('residual_over_red_line_d1', d1_S31, 'zs', 'DERIVED', 'slope - 271.7 (red line DIGITIZED this session)'),
                ('all_R_single_bin_slope_err', all_one['err'], 'zs', 'DERIVED', 'all ~267k Fig.2 events x Fig3 cut efficiency'),
                ('sigma_d1_B_4bins', sB, 'zs', 'DERIVED', 'shape test, a*R free: immune to momentum-zero and R-scale systematics'),
                ('dmin_B_2sigma', 2 * sB, 'zs', 'DERIVED', 'expected Z=2 (Asimov)'), ('dmin_B_5sigma', 5 * sB, 'zs', 'DERIVED', 'expected Z=5'),
                ('dmin_B_2sigma_95pct_power', (2 + 1.645) * sB, 'zs', 'DERIVED', ''),
                ('Z_B_if_stage31_residual_real', d1_S31 / sB, 'sigma', 'DERIVED', ''),
                ('sigma_d1_A_4bins', sA, 'zs', 'DERIVED', f'travel time fixed; delta*R prior {J.SIG_DELTA_DEFAULT:.1f} zs at R_S (momentum zero + R scale)'),
                ('Z_A_if_stage31_residual_real', d1_S31 / sA, 'sigma', 'DERIVED', ''),
                ('sigma_d1_A_authors_band_prior', sA_band, 'zs', 'DERIVED', 'prior = authors purple band half-width 121 zs (DIGITIZED)'),
                ('Z_A_band_if_stage31_real', d1_S31 / sA_band, 'sigma', 'DERIVED', ''),
                ('best_binning', best['nbins'], '', 'DERIVED', f"sigma_d1_B = {best['sigma_d1_B_zs']:.1f} zs"),
                ('sigma_d1_B_with_2nd_energy_400', E2['400']['sigma_d1_B_zs'], 'zs', 'DERIVED+ASSUMED', 'same events, visibility; width ~1/p_e'),
                ('sigma_d1_B_with_2nd_energy_1600', E2['1600']['sigma_d1_B_zs'], 'zs', 'DERIVED+ASSUMED', ''),
                ('sigma_d1_B_with_800_repeat', E2['800 repeat (x2 events)']['sigma_d1_B_zs'], 'zs', 'DERIVED', ''),
                ('sigma_d1_A2_commondp_400', E2['400']['sigma_d1_A2_commondp_zs'], 'zs', 'DERIVED+ASSUMED', 'travel fixed, same absolute dp0 at both energies'),
                ('sigma_d1_A2_commondp_1600', E2['1600']['sigma_d1_A2_commondp_zs'], 'zs', 'DERIVED+ASSUMED', ''),
                ('MC_B_std_H1', MC['H1_stage31']['B_std_zs'], 'zs', 'DERIVED', f'{NMC} full-pipeline simulations'),
                ('MC_B_mean_H1', MC['H1_stage31']['B_mean_zs'], 'zs', 'DERIVED', f'true {d1_S31:.1f}'),
                ('MC_B_mean_H0', MC['H0']['B_mean_zs'], 'zs', 'DERIVED', 'true 0'),
                ('MC_frac_Z_gt_2_if_H1', MC['H1_stage31']['frac_B_Z_gt_2'], '', 'DERIVED', ''),
                ('H0free_cal_coefficient', a_free / J.tau_Rc(1.0), 'x 1/c', 'DERIVED', 'if the whole excess scales with R')]
        for r in rows: w.writerow(r)
    print(json.dumps(dict(reproduce=reproduce, cal_ratio=cal_ratio, NSCALE=NSCALE, head=head, best=best, E2=E2, ndb=ndb, Rres=Rres,
                          MC={k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if 'slope' not in kk}) for k, v in MC.items()}), indent=1, default=float))
    print('bscan', json.dumps(bscan, default=float))

def write_selftest_file(path='/tmp/stage32_selftest_joint.csv', d1=None, seed=7):
    """simulated joint file in the format the analysis script accepts (KER, cos(alpha), cos(beta), count)."""
    d1 = d1_S31 if d1 is None else d1; rng = np.random.default_rng(seed)
    with open(path, 'w') as f:
        f.write('Variables = \nKinetic Energy Release (eV), cos(alpha), cos(beta), count\n')
        for k, R in enumerate(R2):
            if N2k[k] <= 0: continue
            lam, _, _ = lam_rows(R, N2k[k] * NSCALE, s_zs=float(J.tau_Rc(R) + J.nd_term(R) + d1))
            Hk = rng.poisson(lam)
            for i in range(len(cb3)):
                for a_ in range(len(ca3)): f.write(f'{K2[k]:.6f}, {ca3[a_]:.6f}, {cb3[i]:.6f}, {Hk[i, a_]:.0f}\n')
    return path
