#!/usr/bin/env python3
"""Stage 31 - Dan's "do-si-do" (hooked-arm) hypothesis: a minimal phenomenological model and its bounds.

Picture (Dan): particles are locked, localized field states. Locking/unlocking works like square dancers' hooked arms:
the hook releases gradually (finite release time tau, possible plastic stage, rate dependence), hooks differ in openness
(a spread of lock strengths), they have carabiner gates (a threshold to open, spring return), and release depends on
frequency, amplitude and - above all - ANGLE OF INCIDENCE. Lock and unlock need not be mirror images (hysteresis).

Model (ASSUMED structure, one parameter per element). A probe (photon / field / particle) meets a lock with incidence
angle theta (relative to the lock axis), angular frequency w and amplitude A. Two readings are tested:
  R1 "replacement": the hooks ARE the mechanism of ordinary particle events (absorption, detection, emission).
     Then each element must reproduce standard QM inside the precision of the best experiment; any surplus is the
     hook's EXTRA effect.
  R2 "add-on": the hooks are an extra localization channel on top of QM, with rate per nucleon
        Lambda = Lambda0 * g(theta) * G(A) * L(w) * h ,
     g = 1 + a1 cos(theta) + a2 P2(cos theta)          (1) angle
     G = 1/(1+exp(-(A-A_th)/w_A))                      (2) gate threshold A_th (spring return re-closes the gate)
     L = Lorentzian(w; w0, Q)                          (3) resonant release
     amplitude: Weinberg-type energy term eps*|c|^2    (4) amplitude
     h ~ distribution of openness, rel. width sigma_h  (5) spread
     release kernel exp(-t/tau)/tau, jitter sigma_tau  (6) release time
     Gamma_lock/Gamma_unlock = exp(-dE/kT)*(1+eta)      (7) hysteresis / asymmetry
Every experimental input below was read from a source fetched in THIS session (Stage 31 standing rule); each carries
`source` (document + location) and `quote`. Labels: DERIVED (computed here), CITED (fetched source), ASSUMED (our choice).
"""
import json, os, csv, time
import numpy as np
from math import pi, sqrt, log, exp
from scipy.optimize import minimize
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
t0 = time.time()
H_PL = 6.62607015e-34; HBAR = H_PL / (2 * pi); C = 299792458.0; EV = 1.602176634e-19
AU_T = 2.4188843265857e-17; BOHR = 5.29177210903e-11; HARTREE_EV = 27.211386245988   # CODATA exact/defined-unit conversions

# ------------------------------------------------------------------------------------------------ verified sources
S = {}
def src(key, ref, url, loc, quote, **vals):
    S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)

src('papaliolios', 'Papaliolios, PRL 18, 622 (1967), as reported by Valentini & Varma, arXiv:2505.07510v3 (2025)',
    'https://arxiv.org/html/2505.07510v3', 'Sec. 1 Introduction, paragraph on ref. [48]',
    '"Papaliolios verified that the transmission probability through the third polariser varied with Theta according to (3) '
    'within an accuracy of about 1%, thereby confirming the approximate validity of the Born rule over timescales of order 10^-13 s"',
    malus_dev=0.01, t_after=1e-13)
src('valentini_f', 'Valentini & Varma, arXiv:2505.07510v3', 'https://arxiv.org/html/2505.07510v3', 'Sec. 2.3, Eqs. (30)-(31)',
    '"f(Theta) = -cos(2Theta) sin^2(2Theta) (where f(Theta)=0 at Theta=0, pi/4, pi/2). The usual Malus law (3) ... receives a correction"')
src('marx', 'Marx et al., PRL 110, 254801 (2013) (APS abstract via INSPIRE record 1280263)',
    'https://inspirehep.net/api/literature/1280263?format=json', 'abstract',
    '"The polarization purity of 6.457- and 12.914-keV x rays has been improved to the level of 2.4x10-10 and 5.7x10-10. '
    'The polarizers are channel-cut silicon crystals using six 90 deg reflections."', purity=2.4e-10, E_keV=6.457)
src('grundmann', 'Grundmann et al., Science 370, 339 (2020), arXiv:2010.08298; data Zenodo 10.5281/zenodo.3899920 (Fig3_data.txt)',
    'https://arxiv.org/pdf/2010.08298', 'main text Eqs. (2),(3); Fig. 2/3 captions; SM Materials and Methods',
    '"tau_b = cos(alpha_0) R/v_ph" ; "v_ph = E_e/p_e" ; "tau_b\' = cos(beta) R/c ... the expected birth time delay was +-247 zs for R=0.74 A" ; '
    '"fast electrons [E_e = (735+-15) eV] ... R = (0.74+-0.02) A" ; "The accuracy to which we could specify the zero point of the momentum '
    'distribution of the fast electron was approximately +-0.04 atomic units in all three directions" ; "This model is in reasonable agreement '
    'with the experimental results, but more theoretical work including electron-electron correlation is needed to clarify the deviation." ; '
    '[19] "our interferometric technique is blind to this absolute delay"', R_A=0.74, dR_A=0.02, Ee_eV=735.0, dEe_eV=15.0, dp0_au=0.04, tau_model_zs=247)
src('grundmann_fig3c', 'Grundmann et al. 2020, published-version PDF (Frankfurt), Fig. 3C, digitized from a 400 dpi render',
    'https://www.atom.uni-frankfurt.de/publications/files/Grundmann_Science_2020.pdf', 'p. 11, Fig. 3C (crop saved as data/stage31/grundmann2020_fig3C_crop.png)',
    'DIGITIZED (reading +-10 zs): blue R/c line reaches 249 zs at cos(beta)=1 (check vs 247 CITED); red line (nondipole + travel time) ~270 zs; '
    'dashed purple linear fit to the data ~440 zs at |cos(beta)|=1; purple systematic band ~335-575 zs there (half-width ~120 zs). Caption: "The error bars '
    'include statistical and systematic errors and the purple-shaded error range indicates the systematical error"',
    blue_check_zs=249, red_slope_zs=270, paper_fit_slope_zs=440, paper_sys_halfwidth_zs=120)
src('grundmann_fig2', 'Grundmann et al. 2020, Zenodo 10.5281/zenodo.3899920 Fig2_data.txt (cos(alpha) vs KER, beta-averaged)',
    'https://zenodo.org/api/records/3899920/files/Fig2_data.txt/content', 'file columns "cos(alpha), Kinetic Energy Release (eV), count"; paper text',
    '"In the reflection approximation, the internuclear distance R at the moment of photoabsorption is related to the kinetic energy release via '
    'KER = e^2/(4 pi eps0 R)"')
src('ossiander', 'Ossiander et al., Nature Phys. 13, 280 (2017) (author PDF, TUM)',
    'https://www.ph.nat.tum.de/fileadmin/w00bya/e11/Publications/Ossiander_Nature_Physics_2016.pdf', 'Table 1; p.2 text',
    'Table 1 "Time shift dtau between direct and shake-up photoemission": 93.9 eV -12.6 (s.e. 0.99) theory -12.7; 97.2 eV -10.6 (0.85) -11.0; '
    '108.2 eV -5.0 (1.01) -5.9; 113.0 eV -4.9 (1.60) -5.8 [1e-18 s]. "The remarkable agreement within the sub-attosecond standard error excludes '
    'systematic errors"', E=[93.9, 97.2, 108.2, 113.0], exp=[-12.6, -10.6, -5.0, -4.9], se=[0.99, 0.85, 1.01, 1.60], th=[-12.7, -11.0, -5.9, -5.8])
src('korzh', 'Korzh et al., Nature Photon. 14, 250 (2020), arXiv:1804.06839', 'https://arxiv.org/abs/1804.06839', 'abstract',
    '"We reduce this value to 2.7+-0.2 ps at 400 nm and 4.6+-0.2 ps at 1550 nm"', jitter_ps=2.7)
src('minev', 'Minev et al., Nature 570, 200 (2019), arXiv:1803.00545', 'https://arxiv.org/abs/1803.00545', 'abstract; main text Fig. 3/4; SM Sec. II',
    '"the jump evolution when completed is continuous, coherent, and deterministic ... we catch and reverse a quantum jump mid-flight" ; '
    '"between Delta t_catch = 2 us and 10 us, the state of the atom evolves continuously and coherently" ; "succeeds in reversing the jump to |G> '
    'with 83.1%+-0.3% fidelity, while ... completes the jump to |D>, with 82.0%+-0.3%" ; "The measured mid-flight time, 3.95 us, is predicted, in a '
    'first approximation, by Eq. (14) ... 4.3 us"', t_mid_us=3.95, t_mid_th_us=4.3, F_rev=0.831, F_comp=0.820)
src('dixit', 'Dixit et al., PRL 126, 141302 (2021), arXiv:2008.12231', 'https://arxiv.org/abs/2008.12231', 'abstract; main text (detector characterization)',
    '"6.011 GHz (24.86 ueV)" ; "the cavity state collapses to a Fock state ... upon measurement, rather than being absorbed and destroyed" ; '
    '"We fit this relationship with the function n_meas = eta n_inj + delta ... eta = 0.409+-0.055 and ... delta = (4.3+-1.1)x10^-4"',
    E_ueV=24.86, f_GHz=6.011, eta=0.409, d_eta=0.055, delta=4.3e-4)
src('heinzen', 'Heinzen, Bollinger, Itano, Gilbert & Wineland, Coherence & Quantum Optics VI (Plenum 1990) p.479 = Bollinger et al. PRL 63, 1031 (1989)',
    'https://tf.nist.gov/general/pdf/898.pdf', 'p.480, Eq. (5)',
    '"the angle of the first RF pulse alternates between theta_A = 0.325pi and theta_B = 0.675pi" ; "|eps| < 2.4 x [exponent lost in PDF text extraction] eV (5.8 uHz) on a possible '
    'nonlinear correction ... less than 4 parts in 10^27 of the binding energy per nucleon" ; "transition frequency of v = 303 MHz" ; PRL abstract '
    '(doi.org/10.1103/physrevlett.63.1031): "limit of 4x10^-27 on the fraction of binding energy per nucleon"', eps_uHz=5.8, nu_MHz=303.0, frac_BE=4e-27)
src('pruttivarasin', 'Pruttivarasin et al., Nature 517, 592 (2015), arXiv:1412.2194', 'https://arxiv.org/abs/1412.2194', 'abstract',
    '"After a 23 hour measurement, we limit the energy variations to h x 11 mHz ... verifying that Lorentz symmetry is preserved at the level of 1x10^-18"',
    dE_mHz=11.0, frac=1e-18)
src('vip2', 'Piscicchia et al. (VIP), Acta Phys. Pol. A 142, 361 (2022) [reviewing Napolitano et al., Symmetry 14, 893 (2022)]',
    'https://doi.org/10.12693/aphyspola.142.361', 'Eq. (2)',
    '"beta^2/2 <= 8.6 x 10^-31 (Bayesian), beta^2/2 <= 8.9 x 10^-31 (CL)" (90% CL, simple diffusion model)', b2=8.6e-31)
src('helou', 'Helou et al., PRD 95, 084054 (2017), arXiv:1606.03637', 'https://arxiv.org/abs/1606.03637', 'abstract; text',
    '"we bound the CSL collapse rate to be at most (2.96+-0.12)x10^-8 s^-1 ... new frequency regime, 0.7 mHz to 20 mHz" ; "r_CSL ... 100 nm ... we will '
    'follow this convention" ; "lower bound 10^(-8+-2) s^-1 proposed by Adler"', lam=2.96e-8, fmin=0.7e-3, fmax=20e-3)
src('vinante', 'Vinante et al., PRL 119, 110401 (2017), arXiv:1611.09776', 'https://arxiv.org/abs/1611.09776', 'abstract; text near Fig. 3; SM',
    '"High accuracy measurements of the cantilever thermal fluctuations reveal a nonthermal force noise of unknown origin" ; "For the standard choice '
    'r_C=10^-7 m this would imply lambda=10^-7.7 s^-1" (if all excess noise were CSL) ; "resonance frequency of 8174.0 Hz" ; "GRW ... lambda ~ 10^-16 s^-1"',
    lam=10 ** -7.7, f=8174.0, lam_grw=1e-16)
src('majorana', 'Arnquist et al. (Majorana), PRL 129, 080401 (2022), arXiv:2202.01343', 'https://arxiv.org/abs/2202.01343', 'abstract; Eqs. after (5)-(6)',
    '"no evidence of such radiation in the 19-100 keV range" ; "n-m-p CSL ... lambda/r_C^2=(5.15+-0.16)x10^-6 s^-1 m^-2" ; "only the 30 quasi-free '
    'electrons ... lambda/r_C^2=(17.4+-0.5) s^-1 m^-2" ; "coherent emission from nuclei ... (4.94+-0.15)x10^-1 s^-1 m^-2" ; "Extensions of the CSL models '
    'such as the colored noise models are required to evade such stringent constraints"', nmp=5.15e-6, mp_e=17.4, mp_n=0.494, Emin_keV=19, Emax_keV=100)
src('firas', 'Fixsen et al., ApJ 473, 576 (1996), arXiv:astro-ph/9605054', 'https://arxiv.org/abs/astro-ph/9605054', 'abstract; Sec. 1',
    '"The RMS deviations are less than 50 parts per million of the peak of the CMBR ... |y| < 15x10^-6 and |mu| < 9x10^-5 (95% CL)" ; "it could be '
    'distorted by energy release after a redshift z ~ 3x10^6"', mu=9e-5, rms_ppm=50)

# ------------------------------------------------------------------------------------------------ (1)+(6) angle-of-incidence timing: Grundmann reanalysis
g = S['grundmann']['values']
d = np.loadtxt(os.path.join(ROOT, 'data', 'stage31', 'grundmann2020_Fig3_data.txt'), skiprows=2, delimiter=',')
ca = np.unique(d[:, 2]); cb = np.unique(d[:, 3])
Hh = np.zeros((len(ca), len(cb)))
for r in d: Hh[np.searchsorted(ca, r[2]), np.searchsorted(cb, r[3])] = r[4]
pe_au = sqrt(2 * g['Ee_eV'] / HARTREE_EV)                     # DERIVED nonrelativistic (735 eV / 511 keV: 0.14% correction, neglected)
vph_au = (g['Ee_eV'] / HARTREE_EV) / pe_au                   # CITED definition v_ph = E_e/p_e
R_au = g['R_A'] * 1e-10 / BOHR
zs_per_cos = R_au / vph_au * AU_T * 1e21                      # tau_b [zs] = cos(alpha0) * this
tauRc_zs = g['R_A'] * 1e-10 / C * 1e21                        # R/c model slope

rng = np.random.default_rng(31)
def _model(q, x):
    lam = np.exp(q[0]) * np.exp(-0.5 * ((x - q[1]) / np.exp(q[2])) ** 2) + np.exp(q[3])
    if len(q) == 5: lam = lam * (1 + q[4] * x)                    # tilted (linear) envelope variant
    return np.clip(lam, 1e-9, None)
def _fit(x, y, tilt=False):
    def nll(q):
        lam = _model(q, x)
        return np.sum(lam - y * np.log(lam))
    q0 = [np.log(max(y.max() - y.min(), 1.0)), 0.0, np.log(0.2), np.log(max(y.min(), 1.0))] + ([0.0] if tilt else [])
    r = minimize(nll, q0, method='Nelder-Mead', options=dict(xatol=1e-8, fatol=1e-10, maxiter=40000))
    return r.x

def fit_peak(y, x, win, tilt=False, nboot=300):
    m = np.abs(x) <= win; x = x[m]; y = y[m]
    q = _fit(x, y, tilt)
    lam = _model(q, x)
    boots = [_fit(x, rng.poisson(lam).astype(float), tilt)[1] for _ in range(nboot)]   # parametric Poisson bootstrap for sigma(mu)
    return q[1], float(np.std(boots))

def analyse(win, tilt=False):
    mus, ers = [], []
    for j in range(len(cb)):
        m_, e_ = fit_peak(Hh[:, j], ca, win, tilt); mus.append(m_); ers.append(e_)
    mus = np.array(mus); ers = np.array(ers)
    tau = mus * zs_per_cos; dtau = ers * zs_per_cos
    X = np.vstack([np.ones_like(cb), cb, 0.5 * (3 * cb ** 2 - 1)]).T
    W = np.diag(1 / dtau ** 2)
    cov = np.linalg.inv(X.T @ W @ X); beta = cov @ X.T @ W @ tau
    chi2 = float(((tau - X @ beta) / dtau) @ ((tau - X @ beta) / dtau))
    X2 = X[:, :2]; cov2 = np.linalg.inv(X2.T @ W @ X2); b2 = cov2 @ X2.T @ W @ tau
    chi2_2 = float(((tau - X2 @ b2) / dtau) @ ((tau - X2 @ b2) / dtau))
    return dict(window=win, tilt=tilt, cos_beta=cb.tolist(), tau_zs=tau.tolist(), dtau_zs=dtau.tolist(),
                lin=dict(a=b2[0], da=sqrt(cov2[0, 0]), b=b2[1], db=sqrt(cov2[1, 1]), chi2=chi2_2, dof=len(cb) - 2),
                quad=dict(a=beta[0], b=beta[1], db=sqrt(cov[1, 1]), c=beta[2], dc=sqrt(cov[2, 2]), chi2=chi2, dof=len(cb) - 3))

WINS = [0.3, 0.4, 0.5]                                           # ASSUMED fit windows in cos(alpha); all reported, 0.4 headline
GR = {w: analyse(w) for w in WINS}
GT = analyse(0.4, tilt=True)                                     # ASSUMED variant: tilted envelope inside window
h = GR[0.4]
sys_slope = g['dp0_au'] / pe_au * zs_per_cos                    # DERIVED: momentum-zero offset along photon axis projects as cos(beta)
sys_R = h['lin']['b'] * g['dR_A'] / g['R_A']
b_spread = max(abs(GR[w]['lin']['b'] - h['lin']['b']) for w in WINS)
b_tilt = GT['lin']['b']
db_tot = sqrt(h['lin']['db'] ** 2 + sys_slope ** 2 + sys_R ** 2 + b_spread ** 2)
extra_b = h['lin']['b'] - tauRc_zs
lim_b95 = abs(extra_b) + 1.96 * db_tot
lim_c95 = abs(h['quad']['c']) + 1.96 * sqrt(h['quad']['dc'] ** 2 + sys_R ** 2)
gd = S['grundmann_fig3c']['values']
extra_red = h['lin']['b'] - gd['red_slope_zs']
err_red = sqrt(h['lin']['db'] ** 2 + gd['paper_sys_halfwidth_zs'] ** 2 + b_spread ** 2 + 10.0 ** 2)   # stat + paper sys band + window + digitizing
lim_red95 = abs(extra_red) + 1.96 * err_red
extra_red_paper = gd['paper_fit_slope_zs'] - gd['red_slope_zs']
d2 = np.loadtxt(os.path.join(ROOT, 'data', 'stage31', 'grundmann2020_Fig2_data.txt'), skiprows=2, delimiter=',')
Ku = np.unique(d2[:, 3]); wK = np.array([d2[d2[:, 3] == k][:, 4].sum() for k in Ku]); cK = np.cumsum(wK) / wK.sum()
K05, K95 = Ku[np.searchsorted(cK, 0.05)], Ku[np.searchsorted(cK, 0.95)]
KE_COUL_eVA = 1.602176634e-19 / (4 * pi * 8.8541878128e-12) / 1e-10            # e^2/(4 pi eps0) in eV*Angstrom (DERIVED, CODATA)
R_range_A = [KE_COUL_eVA / K95, KE_COUL_eVA / K05]
grundmann = dict(KER_5_95_eV=[K05, K95], R_range_A=R_range_A, Rc_range_zs=[r * 1e-10 / C * 1e21 for r in R_range_A], pe_au=pe_au, vph_au=vph_au, zs_per_unit_cos_alpha0=zs_per_cos, model_slope_R_over_c_zs=tauRc_zs,
                 windows={str(w): GR[w] for w in WINS}, tilted_envelope_win04=GT, headline_window=0.4,
                 slope_zs=h['lin']['b'], slope_stat_zs=h['lin']['db'], slope_sys_p0_zs=sys_slope, slope_sys_R_zs=sys_R,
                 slope_window_spread_zs=b_spread, slope_tilted_envelope_zs=b_tilt, slope_tilted_envelope_stat_zs=GT['lin']['db'], slope_total_err_zs=db_tot,
                 extra_slope_zs=extra_b, extra_slope_sigma=extra_b / db_tot,
                 bound_extra_cos_term_zs_95=lim_b95, bound_extra_P2_term_zs_95=lim_c95,
                 bound_extra_cos_term_fraction_of_Rc=lim_b95 / tauRc_zs,
                 extra_vs_red_model_zs=extra_red, extra_vs_red_err_zs=err_red, extra_vs_red_sigma=extra_red / err_red,
                 bound_extra_vs_red_95_zs=lim_red95, paper_slope_minus_red_zs=extra_red_paper,
                 paper_slope_minus_red_sigma=extra_red_paper / gd['paper_sys_halfwidth_zs'],
                 label='DERIVED (reanalysis of CITED Zenodo histogram; fit window ASSUMED, 3 variants reported)',
                 note='Symmetric (angle-independent) release delays are invisible to this interferometer (CITED note [19]).')

# ------------------------------------------------------------------------------------------------ (6)/(5) channel-dependent release time: Ossiander
o = S['ossiander']['values']
diff = np.array(o['exp']) - np.array(o['th']); se = np.array(o['se'])
w_ = 1 / se ** 2; md = float(np.sum(w_ * diff) / np.sum(w_)); dmd = float(1 / sqrt(np.sum(w_)))
chi2_o = float(np.sum(((diff) / se) ** 2))
ossiander = dict(diff_as=diff.tolist(), weighted_mean_as=md, err_as=dmd, chi2_vs_theory=chi2_o, dof=4, p=float(stats.chi2.sf(chi2_o, 4)),
                 bound_abs_extra_channel_delay_as_95=abs(md) + 1.96 * dmd, label='DERIVED from CITED Table 1 (theory errors ASSUMED 0)')

# ------------------------------------------------------------------------------------------------ (1) angle laws
malus = dict(generic_dev_max=S['papaliolios']['values']['malus_dev'], timescale_s=1e-13,
             crossed_extra_leak_max=S['marx']['values']['purity'],
             note='Crossed-polarizer purity bounds any extra term at Theta=90 deg (it also contains all imperfections, so it is an upper bound). '
                  'A shape f(Theta) vanishing at 0, 45, 90 deg (e.g. Valentini-Varma Eq. 30) evades it; then only the ~1% Papaliolios bound applies.',
             label='CITED')
pr = S['pruttivarasin']['values']
aniso = dict(dE_eV=pr['dE_mHz'] * 1e-3 * H_PL / EV, frac=pr['frac'], label='CITED',
             note='Bounds orientation dependence (lab frame vs fixed stars) of electron energy; for an add-on lock RATE the anisotropic part is '
                  'bounded by the total rate bound (element 3 bands), since 0 <= Lambda(theta) <= 2 x mean for |a1|,|a2| <= 1.')

# ------------------------------------------------------------------------------------------------ (2) threshold
dx = S['dixit']['values']
thresh = dict(E_quantum_eV=dx['E_ueV'] * 1e-6, E_quantum_J=dx['E_ueV'] * 1e-6 * EV, f_GHz=dx['f_GHz'],
              linear_eta=dx['eta'], linear_deta=dx['d_eta'], false_pos=dx['delta'],
              note='R1: single 24.86 ueV quanta are individually registered (locked into a Fock state) and the response is linear in injected '
                   'n = alpha^2 << 1, i.e. no amplitude gate down to sub-single-photon amplitudes; the QND readout locks the photon number without '
                   'absorbing it, so an energy-transfer threshold > 0 is not required for a lock. R2: a gate simply hides the add-on below A_th, '
                   'so A_th alone is unconstrained (degenerate with Lambda0).', label='CITED + DERIVED')

# ------------------------------------------------------------------------------------------------ (3) frequency selectivity: band bounds (CSL-equivalent rate per nucleon, r_C = 1e-7 m)
mj = S['majorana']['values']; rc = 1e-7
bands = [dict(name='LISA Pathfinder (Helou 2017)', fmin=S['helou']['values']['fmin'], fmax=S['helou']['values']['fmax'], lam=S['helou']['values']['lam'], kind='force noise'),
         dict(name='Cantilever 8174 Hz (Vinante 2017; all excess attributed)', fmin=8174.0 * (1 - 1e-3), fmax=8174.0 * (1 + 1e-3), lam=S['vinante']['values']['lam'], kind='force noise'),
         dict(name='Majorana X-ray, m-p electrons only', fmin=mj['Emin_keV'] * 1e3 * EV / H_PL, fmax=mj['Emax_keV'] * 1e3 * EV / H_PL, lam=mj['mp_e'] * rc ** 2, kind='radiation'),
         dict(name='Majorana X-ray, m-p coherent nuclei', fmin=mj['Emin_keV'] * 1e3 * EV / H_PL, fmax=mj['Emax_keV'] * 1e3 * EV / H_PL, lam=mj['mp_n'] * rc ** 2, kind='radiation'),
         dict(name='Majorana X-ray, n-m-p', fmin=mj['Emin_keV'] * 1e3 * EV / H_PL, fmax=mj['Emax_keV'] * 1e3 * EV / H_PL, lam=mj['nmp'] * rc ** 2, kind='radiation')]
gap = dict(fmin=20e-3, fmax=bands[2]['fmin'], note='Between 20 mHz and 4.6e18 Hz only the single 8174 Hz cantilever line was verified this stage; '
           'other bounds exist in the literature (cold atoms, bulk heating, other cantilevers) but were not fetched, so they are not used.')
freq = dict(bands=bands, gap=gap, rC_m=rc, ref_rates=dict(GRW=S['vinante']['values']['lam_grw'], Adler='1e-8 (x100^+-1)'),
            decades_unverified=float(np.log10(gap['fmax'] / gap['fmin']) - np.log10(bands[1]['fmax'] / bands[1]['fmin'])),
            label='CITED (lambda) + DERIVED (keV -> Hz)')

# ------------------------------------------------------------------------------------------------ (4) amplitude
hz = S['heinzen']['values']
eps_eV = hz['eps_uHz'] * 1e-6 * H_PL / EV
amp = dict(eps_eV=eps_eV, eps_over_transition=hz['eps_uHz'] * 1e-6 / (hz['nu_MHz'] * 1e6), frac_binding=hz['frac_BE'],
           label='CITED + DERIVED', note='Drive pulse area (0.325pi vs 0.675pi) and hence state amplitudes were alternated: no frequency change.')

# ------------------------------------------------------------------------------------------------ (5) spread of openness
spread = dict(intrinsic_fraction_max=S['vip2']['values']['b2'], mapping='ASSUMED O(1): a fraction f of electrons with a distinguishable '
              'openness label would escape the exclusion principle in K-shell capture, so f <~ beta^2/2',
              extrinsic_channel_delay_as_95=ossiander['bound_abs_extra_channel_delay_as_95'],
              note='If openness is set by context (extrinsic), particles stay identical and only the add-on strength bounds apply.')

# ------------------------------------------------------------------------------------------------ (6) release time summary
mv = S['minev']['values']
tau = dict(angle_dependent_zs_95=lim_red95, channel_dependent_as_95=ossiander['bound_abs_extra_channel_delay_as_95'],
           detection_jitter_ps=S['korzh']['values']['jitter_ps'], post_lock_plastic_stage=dict(born_dev_max=0.01, within_s=1e-13),
           universal_constant='not bounded by relative-timing data (cable-delay degenerate)',
           qm_already=dict(minev_t_mid_us=mv['t_mid_us'], minev_t_mid_theory_us=mv['t_mid_th_us'], reverse=mv['F_rev'], complete=mv['F_comp']),
           label='DERIVED/CITED')

# ------------------------------------------------------------------------------------------------ (7) hysteresis: detailed balance -> mu distortion
mu = S['firas']['values']['mu']
hyst = dict(eta_max=exp(mu) - 1, mu_max=mu, derivation='DERIVED: absorption ~ A n, emission ~ E (1+n); steady state n/(1+n) = E/A. '
            'QM: E/A = exp(-x). With E/A = exp(-x)/(1+eta): n = 1/((1+eta) e^x - 1) = Bose-Einstein with mu = ln(1+eta).',
            assumption='ASSUMED: the asymmetry acts in the photon-number-changing processes that thermalised the CMB at z >~ 3e6 (CITED epoch), '
                       'and is frequency independent; a frequency-dependent eta(x) is bounded by the 50 ppm-of-peak rms instead.',
            label='DERIVED + CITED')

# ------------------------------------------------------------------------------------------------ verdict table
verdict = [
 dict(element='1 angle of incidence', model='P = cos^2 Theta + delta f(Theta); birth delay = R cos(beta)/c + d1 cos(beta) + d2 P2; Lambda(theta)=Lambda0 g(theta)',
      qm_has='Malus cos^2, Fresnel, dipole patterns, photon-momentum retardation R cos(beta)/c (Grundmann measured it)',
      bound=f"generic |delta| < 1e-2 (1e-13 s after a lock); extra leak at 90 deg < 2.4e-10; extra angle-dependent birth delay beyond light-travel+nondipole |d1| < {lim_red95:.0f} zs per unit cos(beta) (95%; measured residual {extra_red:+.0f} +- {err_red:.0f} zs), |d2| < {lim_c95:.0f} zs; lab-orientation energy anisotropy < 1e-18",
      verdict=f'constrained but with real room in timing (|d1| < {lim_red95 / tauRc_zs:.1f} x the light-travel term, residual {extra_red / err_red:.1f} sigma); Malus 1e-2 generic, 2.4e-10 at 90 deg'),
 dict(element='2 gate threshold', model='G(A) = logistic((A - A_th)/w_A)', qm_has='photoelectric frequency threshold (work function); none in amplitude',
      bound='R1: E_th <= 24.86 ueV per quantum; no amplitude gate down to n<<1 (linear eta=0.41); R2: unconstrained alone (hides effect)',
      verdict='R1 excluded above 2.5e-5 eV / any amplitude gate; R2 open but untestable alone'),
 dict(element='3 frequency selectivity', model='Lorentzian L(w; w0, Q)', qm_has='resonant absorption lines',
      bound='Lambda(w0) <= 3.0e-8 s^-1 (0.7-20 mHz), <= 2.0e-8 (8174 Hz), <= 5.2e-20 (n-m-p) .. 1.7e-13 (m-p electrons) (4.6e18-2.4e19 Hz) at r_C=1e-7 m',
      verdict='constrained in 3 bands; open across ~20 unverified decades (20 mHz-4.6e18 Hz except 8174 Hz)'),
 dict(element='4 amplitude', model='Weinberg eps |c|^2 energy term', qm_has='none (QM linear)',
      bound=f'|eps| < {eps_eV:.1e} eV (h x 5.8 uHz) = {amp["eps_over_transition"]:.1e} of the 303 MHz transition = 4e-27 of binding energy per nucleon', verdict='excluded above 1.9e-14 of the transition energy'),
 dict(element='5 spread of openness', model='h ~ p(h), width sigma_h', qm_has='none for identical particles',
      bound=f"intrinsic label: fraction < 8.6e-31 (VIP-2, ASSUMED O(1) map); extrinsic: channel delay < {ossiander['bound_abs_extra_channel_delay_as_95']:.1f} as",
      verdict='intrinsic excluded; extrinsic open (degenerate with strength)'),
 dict(element='6 release time tau', model='kernel e^{-t/tau}/tau, jitter sigma_tau', qm_has='Wigner/EWS delays; continuous quantum-jump flights (Minev)',
      bound=f"angle-dependent < {lim_red95:.0f} zs; channel-dependent < {ossiander['bound_abs_extra_channel_delay_as_95']:.1f} as; jitter < 2.7 ps; plastic stage: Born holds to 1% at 1e-13 s",
      verdict='constrained; a universal constant tau is unobservable (open)'),
 dict(element='7 hysteresis', model='Gamma_lock/Gamma_unlock = e^{-dE/kT}(1+eta)', qm_has='detailed balance (Einstein A/B)',
      bound=f'|eta| < {hyst["eta_max"]:.1e} (FIRAS mu)', verdict='constrained < 9e-5'),
]

# ------------------------------------------------------------------------------------------------ formal equivalences (literature, from fetched reference lists)
equiv = {
 'add-on lock rate Lambda0 (white)': 'CSL: Ghirardi, Pearle & Rimini, PRA 42, 78 (1990) [Helou ref 20]; GRW',
 'frequency selectivity': 'colored CSL: Adler & Bassi, J. Phys. A 40, 15083 (2007); 41, 395308 (2008); Pearle PRA 48, 913 (1993) [Majorana refs 21-23]',
 'hysteresis / asymmetry (lock vs unlock, energy flow)': 'dissipative CSL: Smirne & Bassi, Sci. Rep. 5, 12518 (2015) [Helou ref 15]',
 'release time set by the lock energy': 'Diosi-Penrose / Penrose OR tau = hbar/E_G: Diosi PRA 40, 1165 (1989); Penrose GRG 28, 581 (1996) [Majorana refs 8-9]',
 'amplitude dependence': 'Weinberg nonlinear QM, PRL 62, 485 (1989) (fetched APS abstract)',
 'plastic stage after locking, Born deviations shortly after a measurement': 'Bohm & Bub, RMP 38, 453 (1966) [Valentini-Varma ref 49], tested by Papaliolios 1967',
 'spread of hidden openness relaxing to Born statistics': 'pilot-wave quantum non-equilibrium: Valentini, PLA 156, 5 (1991) [Valentini-Varma ref 18]',
 'gradual, reversible locking (no extra physics)': 'quantum trajectory theory: Carmichael (1993); Minev et al. 2019 [Minev ref 5]',
}

out = dict(stage=31, title="Dan's do-si-do (hooked-arm) hypothesis: minimal model and bounds",
           sources=S, grundmann=grundmann, ossiander=ossiander, malus=malus, anisotropy=aniso, threshold=thresh, frequency=freq,
           amplitude=amp, spread=spread, tau=tau, hysteresis=hyst, verdict=verdict, equivalences=equiv, runtime_s=round(time.time() - t0, 2))
os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
json.dump(out, open(os.path.join(ROOT, 'results', 'stage31.json'), 'w'), indent=1, default=float)
with open(os.path.join(ROOT, 'results', 'stage31_sources.csv'), 'w', newline='') as f:
    wr = csv.writer(f); wr.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for k, v in S.items(): wr.writerow([k, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])

print(f"Grundmann: zs per unit cos(alpha0) = {zs_per_cos:.0f}; R/c = {tauRc_zs:.1f} zs")
for w in WINS:
    L = GR[w]['lin']; Q = GR[w]['quad']
    print(f"  win {w}: slope {L['b']:.1f} +- {L['db']:.1f} zs (chi2 {L['chi2']:.1f}/{L['dof']}), offset {L['a']:.1f}+-{L['da']:.1f}; P2 coef {Q['c']:.1f} +- {Q['dc']:.1f}")
print(f"  tilted envelope (win 0.4): slope {GT['lin']['b']:.1f} +- {GT['lin']['db']:.1f} zs, chi2 {GT['lin']['chi2']:.1f}/13")
print(f"  sys p0 {sys_slope:.1f}, sys R {sys_R:.1f}, window {b_spread:.1f} -> total {db_tot:.1f}; extra {extra_b:.1f} zs ({extra_b / db_tot:.2f} sigma); 95% |d1| < {lim_b95:.1f} zs ({lim_b95 / tauRc_zs:.2f} R/c); |d2| < {lim_c95:.1f}")
print(f"  vs red (nondipole+travel) model 270 zs: extra {extra_red:.1f} +- {err_red:.1f} zs ({extra_red / err_red:.2f} sigma), 95% < {lim_red95:.0f} zs; paper slope-red {extra_red_paper} ({extra_red_paper / gd['paper_sys_halfwidth_zs']:.2f} sigma)")
print(f"Ossiander: exp-th = {md:.2f} +- {dmd:.2f} as, chi2 {chi2_o:.2f}/4; 95% bound {ossiander['bound_abs_extra_channel_delay_as_95']:.2f} as")
print(f"Weinberg eps = {eps_eV:.2e} eV, /transition = {amp['eps_over_transition']:.2e}")
print(f"anisotropy dE = {aniso['dE_eV']:.2e} eV")
for b in bands: print(f"  band {b['name']}: {b['fmin']:.3g}-{b['fmax']:.3g} Hz lam <= {b['lam']:.3g}")
print(f"eta_max {hyst['eta_max']:.2e}; unverified decades {freq['decades_unverified']:.1f}")
print('runtime', out['runtime_s'])
