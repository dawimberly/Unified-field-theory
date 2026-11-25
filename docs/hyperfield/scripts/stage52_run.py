"""Stage 52: audit of the Stage 31 Grundmann et al. (Science 370, 339 (2020), arXiv:2010.08298v1) slope of 504.5 zs.
All inputs fetched 2026-11-25 (this session): arXiv PDF, Zenodo 10.5281/zenodo.3899920 Fig3_data.txt, NIST CODATA 2022 allascii.txt.
Steps: (1) recompute Eq.2 conversion; (2) digitize the paper's Fig.3C points + dashed fit from a 600 dpi render of p.11;
(3) refit the public histogram with several fit models; (4) compare per-row centres with the paper's points."""
import os, json, csv
import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.optimize import curve_fit
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
DD = os.path.join(ROOT, 'data', 'stage52')
# CITED: NIST CODATA 2022 allascii.txt (fetched this session)
HARTREE_EV = 27.211386245981; BOHR = 5.29177210544e-11; AU_T = 2.4188843265864e-17; C = 299792458.0; ALPHA_INV = 137.035999177; MEC2_EV = 0.51099895069e6
# CITED: arXiv:2010.08298v1 Fig.2/3 captions "Ee = (735+-15) eV", "R = (0.74 +- 0.02) A"; p.4 "800 eV"; SM "+-0.04 atomic units"
Ee, dEe, R_A, dR_A, Eph, dp0 = 735.0, 15.0, 0.74, 0.02, 800.0, 0.04
TAU_PAPER = 247.0     # CITED p.7 "the expected birth time delay was +-247 zs for R = 0.74 A"
RED_ZS, BAND_HW = 270.0, 120.0   # Stage 31 DIGITIZED (red model line, purple band half-width), carried over

out = {}
# ---- (1) conversion audit (DERIVED)
pe = np.sqrt(2 * Ee / HARTREE_EV)                  # nonrelativistic electron momentum, a.u.
pe_rel = np.sqrt(Ee**2 + 2 * Ee * MEC2_EV) / (HARTREE_EV * ALPHA_INV)   # relativistic pc / (Eh * alpha^-1) -> a.u.
vph = (Ee / HARTREE_EV) / pe                        # Eq.2 definition v_ph = Ee/pe
R_au = R_A * 1e-10 / BOHR
ZS = R_au / vph * AU_T * 1e21                       # zs per unit cos(alpha0), Eq.2
ZS_vg = R_au / pe * AU_T * 1e21                     # if group velocity v_g = pe/m were used
ZS_halfR = ZS / 2                                   # if R/2 were used
tauRc = R_A * 1e-10 / C * 1e21
# Eq.1-consistent light-travel slope in cos(alpha0) units: electron phase Ee*tau'/hbar (paper) vs photon phase Eph*tau'/hbar
s_alpha_Ee = (Ee / HARTREE_EV) / pe / ALPHA_INV     # = v_ph/c
s_alpha_Eph = (Eph / HARTREE_EV / ALPHA_INV) / pe   # = p_gamma/p_e
out['conversion'] = dict(pe_au=pe, pe_rel_au=pe_rel, rel_corr=pe_rel / pe - 1, vph_au=vph, vg_au=pe, zs_per_cos_alpha0=ZS,
                         zs_per_cos_alpha0_if_vg=ZS_vg, zs_per_cos_alpha0_if_halfR=ZS_halfR, dZS_from_Ee_frac=0.5 * dEe / Ee,
                         R_over_c_zs=tauRc, light_slope_cos_alpha0_Ee=s_alpha_Ee, light_slope_cos_alpha0_Eph=s_alpha_Eph,
                         light_slope_zs_equiv_Eph=s_alpha_Eph * ZS)

# ---- (2) digitize Fig.3C (600 dpi render of arXiv p.11, crop saved in data/stage52)
a = np.array(Image.open(os.path.join(DD, 'fig3c_600.png')).convert('RGB')).astype(int)
r, g, b = a[..., 0], a[..., 1], a[..., 2]
purp = (r > 100) & (r < 150) & (g < 30) & (b > 100) & (b < 150)
dark = a.sum(2) < 200
ticks = np.where(dark[:, 101:125].sum(1) >= 20)[0]
tick_rows = [int(np.mean(grp)) for grp in np.split(ticks, np.where(np.diff(ticks) > 1)[0] + 1)][1:-1]   # drop frame rows
lab, n = ndimage.label(purp)
sq = []
for i in range(1, n + 1):
    ys, xs = np.where(lab == i)
    if 22 <= np.ptp(xs) <= 30 and 22 <= np.ptp(ys) <= 30: sq.append((xs.mean(), ys.mean()))
sq = sorted(sq); assert len(sq) == 15
X = np.array([s[0] for s in sq]); Y = np.array([s[1] for s in sq])
x0 = (X[0] + X[-1]) / 2; sx = (X[-1] - X[0]) / (2 * 0.933333)
Y0 = tick_rows[2]; pz = (tick_rows[-1] - tick_rows[0]) / 4 / 200.0     # right axis: 200 zs per tick
cb_d = (X - x0) / sx; tau_d = (Y0 - Y) / pz
ys, xs = np.where(purp); keep = np.ones(len(xs), bool)
for x, y in zip(X, Y): keep &= ~(np.abs(xs - x) <= 4); keep &= ~((np.abs(xs - x) <= 16) & (np.abs(ys - y) <= 16))
xs, ys = xs[keep], ys[keep]; keep = (xs > 105) & (xs < 1233) & (ys > 55) & (ys < 1228); xs, ys = xs[keep], ys[keep]
cbp = (xs - x0) / sx; tp = (Y0 - ys) / pz; m = np.ones(len(xs), bool)
for _ in range(5):
    p = np.polyfit(cbp[m], tp[m], 1); res = tp - np.polyval(p, cbp); m = (np.abs(res) < 3 * np.std(res[m])) & (np.abs(res) < 40)
dash_slope, dash_icpt = p
uw_origin_digitized = float(cb_d @ tau_d / (cb_d @ cb_d))
out['digitized'] = dict(tick_rows_px=tick_rows, zs_per_px=1 / pz, cos_beta=cb_d.tolist(), tau_zs=tau_d.tolist(),
                        dashed_slope_zs=dash_slope, dashed_intercept_zs=dash_icpt, unweighted_origin_fit_of_points_zs=uw_origin_digitized)

# ---- (3) refit public histogram
d = np.loadtxt(os.path.join(DD, 'grundmann2020_Fig3_data.txt'), skiprows=2, delimiter=',')
ca = np.unique(d[:, 2]); cb = np.unique(d[:, 3]); H = np.zeros((len(ca), len(cb)))
for row in d: H[np.searchsorted(ca, row[2]), np.searchsorted(cb, row[3])] = row[4]
def G(x, A, mu, s): return A * np.exp(-0.5 * ((x - mu) / s) ** 2)
def GC(x, A, mu, s, B): return G(x, A, mu, s) + B
def centres(model, win):
    mu, e = [], []
    for j in range(len(cb)):
        msk = np.abs(ca) <= win; x = ca[msk]; y = H[msk, j]; s = np.sqrt(np.maximum(y, 1))
        if model == 'GC': pp, cc = curve_fit(GC, x, y, p0=[y.max() - y.min(), 0, .2, y.min()], sigma=s, absolute_sigma=True, maxfev=20000)
        else: pp, cc = curve_fit(G, x, y, p0=[y.max(), 0, .2], sigma=s, absolute_sigma=True, maxfev=20000)
        mu.append(pp[1]); e.append(np.sqrt(cc[1, 1]))
    return np.array(mu) * ZS, np.array(e) * ZS
def line(t, s, origin=False):
    A = np.c_[cb] if origin else np.c_[np.ones_like(cb), cb]; W = 1 / s ** 2
    Cv = np.linalg.inv(A.T @ (A * W[:, None])); pp = Cv @ A.T @ (W * t)
    return float(pp[-1]), float(np.sqrt(Cv[-1, -1])), float(((t - A @ pp) ** 2 * W).sum())
variants = []
for model in ['GC', 'G']:
    for win in [0.25, 0.3, 0.35, 0.4, 0.45]:
        t, s = centres(model, win)
        bw, dbw, chi = line(t, s); bu, _, _ = line(t, np.ones_like(t), True)
        k = float(t @ tau_d / (tau_d @ tau_d)); rms = float(np.sqrt(np.mean((t - tau_d) ** 2)))
        variants.append(dict(model=model, window=win, slope_w=bw, dslope_w=dbw, chi2=chi, slope_unweighted_origin=bu, scale_vs_paper_points=k, rms_vs_paper_points_zs=rms))
out['variants'] = variants
best = min(variants, key=lambda v: v['rms_vs_paper_points_zs'] + 100 * abs(v['scale_vs_paper_points'] - 1))
S31 = dict(slope=504.5, dslope=36.2)     # Stage 31 headline, reproduced this session by rerunning stage31_run.py
fit_sys = abs(S31['slope'] - best['slope_w']) / 2
def resid(sl, dsl):
    err = np.sqrt(dsl ** 2 + BAND_HW ** 2 + fit_sys ** 2 + 10.0 ** 2)
    return sl - RED_ZS, err, (sl - RED_ZS) / err, abs(sl - RED_ZS) + 1.96 * err
out['summary'] = dict(paper_matched_variant=best, stage31_slope=S31, fit_model_sys_zs=fit_sys,
                      ratio_s31_to_dashed=S31['slope'] / dash_slope, ratio_s31_to_Rc=S31['slope'] / tauRc,
                      ratio_matched_to_dashed=best['slope_w'] / dash_slope, ratio_matched_to_Rc=best['slope_w'] / tauRc,
                      ratio_dashed_to_Rc=dash_slope / tauRc,
                      resid_red_matched=resid(best['slope_w'], best['dslope_w']), resid_red_s31=resid(S31['slope'], S31['dslope']))
json.dump(out, open(os.path.join(ROOT, 'results', 'stage52.json'), 'w'), indent=1, default=float)

src = [
 ('grundmann_eq12', 'Grundmann et al., arXiv:2010.08298v1 (author version of Science 370, 339)', 'https://arxiv.org/pdf/2010.08298', 'p.5, Eqs. (1)-(2)',
  '"An electron wave is emitted and propagates with the phase velocity vph = Ee/pe for the time tau_b before a second electron wave is born at a distance R. The zeroth-order interference maximum occurs under the emission angle alpha0 for which the path length difference between both electron waves vanishes: 2pi(cos(alpha0)R/lambda - tau_b vph/lambda) = 0 (1) ... tau_b = cos(alpha0) R/vph (2)"'),
 ('grundmann_eq3', 'same', 'https://arxiv.org/pdf/2010.08298', 'p.6-7, Eq. (3)',
  '"tau_b\' = cos(beta) R/c (3) ... cos(beta) = +-1 resembles the maximum possible travel time of the photon from one molecular center to the other. For this case the expected birth time delay was +-247 zs for R = 0.74 A. Between these extreme cases the birth time delay showed a linear dependence on cos(beta)."'),
 ('grundmann_fit', 'same', 'https://arxiv.org/pdf/2010.08298', 'p.6',
  '"we determined cos(alpha0), i.e., the angular position of the central maximum, for each row of the histogram in Fig. 3 B via a Gaussian fit (red curve in Fig. 3 A). Figure 3 C shows the results of these fits and the corresponding birth time delays (according to Eq. 2) on the right vertical scale."'),
 ('grundmann_fig3cap', 'same', 'https://arxiv.org/pdf/2010.08298', 'p.11, Fig. 3 caption',
  '"(C) Location of the maxima of the zeroth-order interference fringe as function of cos(beta). ... Left axis: cos(alpha0), right axis: birth time delay calculated using Eq. 2. The blue line resembles a birth time delay given by the travel time of light across the molecule [Eq. 3]."'),
 ('grundmann_fig1', 'same', 'https://arxiv.org/pdf/2010.08298', 'p.9, Fig. 1 caption', '"the angle alpha is enclosed by the electron momentum vector and the molecular axis"'),
 ('grundmann_sm', 'same', 'https://arxiv.org/pdf/2010.08298', 'SM Materials and Methods (p.12-13)',
  '"The accuracy to which we could specify the zero point of the momentum distribution of the fast electron was approximately +-0.04 atomic units ... we determined the maximum and minimum slope of the linear fit function in Fig. 3 C" ; "the molecule direction ... points from the proton detected first to the one detected second"'),
 ('grundmann_ref20', 'same', 'https://arxiv.org/pdf/2010.08298', 'p.7', '"Note that the model agreed with the prediction from Ref. (20) [Eq. 12], if one neglects the ionization potential."'),
 ('zenodo', 'Zenodo 10.5281/zenodo.3899920', 'https://zenodo.org/api/records/3899920/files/Fig3_data.txt/content', 'Fig3_data.txt header',
  '"Bin#x, Bin#y, cos(alpha), cos(beta), count" (45 x 15 bins, 47413 counts); byte-identical to the Stage 31 copy'),
 ('codata', 'NIST CODATA 2022', 'https://physics.nist.gov/cuu/Constants/Table/allascii.txt', 'table rows',
  'Hartree energy in eV 27.211386245981; Bohr radius 5.29177210544e-11 m; atomic unit of time 2.4188843265864e-17 s; c 299792458; inverse fine-structure 137.035999177; electron mass energy 0.51099895069 MeV'),
]
with open(os.path.join(ROOT, 'results', 'stage52_sources.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['key', 'ref', 'url', 'location', 'quote']); [w.writerow(s) for s in src]

cv = out['conversion']
print(f"Eq.2: pe={pe:.4f} au (rel corr {cv['rel_corr']*100:.3f}%), vph={vph:.4f} au -> {ZS:.1f} zs/unit cos(a0); with vg {ZS_vg:.1f}; with R/2 {ZS_halfR:.1f}; Ee+-15 eV -> +-{cv['dZS_from_Ee_frac']*100:.1f}%")
print(f"R/c = {tauRc:.1f} zs; light slope in cos(a0): {s_alpha_Ee:.5f} (Ee phase) = {s_alpha_Ee*ZS:.1f} zs ; {s_alpha_Eph:.5f} (photon phase) = {s_alpha_Eph*ZS:.1f} zs-equiv")
print(f"Fig.3C ticks rows {tick_rows}: {1/pz:.3f} zs/px; 200 zs tick = {200/ZS:.5f} in cos(a0) (label 0.022)")
print(f"digitized dashed fit: slope {dash_slope:.1f} zs, intercept {dash_icpt:.1f}; unweighted through-origin fit of digitized points {uw_origin_digitized:.1f}")
for v in variants: print(f"  {v['model']:2s} |cos a|<={v['window']}: weighted {v['slope_w']:6.1f} +- {v['dslope_w']:4.1f} (chi2 {v['chi2']:4.1f}/13), unweighted-origin {v['slope_unweighted_origin']:6.1f}, scale vs paper pts {v['scale_vs_paper_points']:.3f}, rms {v['rms_vs_paper_points_zs']:6.1f}")
sm = out['summary']
print(f"paper-matched: {best['model']} {best['window']}: {best['slope_w']:.1f} +- {best['dslope_w']:.1f}; fit-model sys +-{fit_sys:.1f}")
print(f"ratios: S31/dashed {sm['ratio_s31_to_dashed']:.3f}, S31/Rc {sm['ratio_s31_to_Rc']:.3f}, matched/dashed {sm['ratio_matched_to_dashed']:.3f}, matched/Rc {sm['ratio_matched_to_Rc']:.3f}, dashed/Rc {sm['ratio_dashed_to_Rc']:.3f}")
for k in ['resid_red_matched', 'resid_red_s31']: print(k, 'extra %.1f +- %.1f (%.2f sigma), 95%% bound %.0f zs' % sm[k])
