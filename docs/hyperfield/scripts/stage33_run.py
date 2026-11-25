#!/usr/bin/env python3
"""Stage 33 - closing the frequency gaps for resonant (frequency-selective) release in Dan's do-si-do hypothesis.

Parametrization (kept EXACTLY from Stage 31, element 3): an add-on lock/release channel with rate per nucleon
    Lambda(w) = Lambda0 * L(w; w0, Q)
expressed as a CSL-equivalent (mass-proportional) collapse rate per nucleon lambda at r_C = 1e-7 m. An experiment that
samples the noise in a frequency band [fmin, fmax] bounds lambda(f) for f in that band. Narrow resonators are drawn as
+-0.1 % lines (Stage 31 convention, ASSUMED). This is formally colored CSL (Adler-Bassi 2007/2008).

Which frequency does an experiment sample?  (ASSUMED physical rule, applied uniformly, no tuning)
  * mechanical oscillators / interferometers: the measured force-noise frequency (resonance or measured band);
  * trapped-ion and phonon-mode heating: the mode frequency (excitation needs noise at w_mode);
  * bulk heating of a solid: longitudinal phonons w = v_s q weighted by q^4 exp(-q^2 r_C^2) (Eq. 17 of
    Carlesso-Ferialdi-Bassi 2018; half-maximum band, DERIVED);
  * free-falling cold atoms: the free-evolution response kernel |int_0^T (T-s) e^{iws} ds|^2 (half-power, DERIVED);
  * spontaneous radiation: photon frequency E/h.
Every experimental number is read from a source fetched in THIS session (results/stage33_sources.csv: url, location,
verbatim quote). Labels: CITED (number from fetched source), DERIVED (computed here from CITED inputs),
ASSUMED (our stated choice).
"""
import json, os, csv, time
import numpy as np
from math import pi, sqrt, exp, log10
from scipy.optimize import brentq
from scipy.integrate import quad

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
t0 = time.time()
S = {}
def src(key, ref, url, loc, quote, **vals):
    S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)

# ------------------------------------------------------------------ constants (CITED: NIST CODATA pages fetched)
src('nist', 'NIST CODATA (physics.nist.gov/cuu), pages Value?hbar, Value?h, Value?ukg, Value?evj',
    'https://physics.nist.gov/cgi-bin/cuu/Value?hbar', 'Numerical value fields',
    '"reduced Planck constant Numerical value 1.054 571 817... x 10 -34 J s" ; "Planck constant Numerical value 6.626 070 15 x 10 -34 J Hz -1" ; '
    '"unified atomic mass unit Numerical value 1.660 539 068 92 x 10 -27 kg" ; "electron volt Numerical value 1.602 176 634 x 10 -19 J" ; "electron mass in u Numerical value 5.485 799 090 441 x 10 -4 u" (Value?meu)',
    hbar=1.054571817e-34, h=6.62607015e-34, u=1.66053906892e-27, eV=1.602176634e-19, me_u=5.485799090441e-4)
HBAR, H, U, EV = (S['nist']['values'][k] for k in ('hbar', 'h', 'u', 'eV'))
M0 = U          # ASSUMED: CSL reference mass m0 = 1 u (sources say "nucleon"; u vs m_p changes bounds by 0.7 %)
RC = 1e-7       # CITED convention (Stage 31; Helou 2017 "r_CSL ... 100 nm ... we will follow this convention")
keV2Hz = 1e3 * EV / H

# ------------------------------------------------------------------ sources
src('ccb18', 'Carlesso, Ferialdi & Bassi, Eur. Phys. J. D 72, 159 (2018), arXiv:1805.10100 "Colored collapse models from the non-interferometric perspective"',
    'https://arxiv.org/pdf/1805.10100', 'Eq. (13),(17),(19); Sec. IV A',
    '"S_CSL = hbar^2 lambda r_C^3/(pi^3/2 m0^2) Int dk |mu(k)|^2 k_x^2 e^{-k^2 r_C^2}, which is independent of omega, being the noise spectrum flat" ; '
    '"LISA Pathfinder ... two identical cubic masses (length 4.6 cm, average distance 37.6 cm, mass 1.928 kg)" ; '
    '"LIGO ... its noise is monitored in the 10 - 10^3 Hz band" ; "characterized by a frequency omega_m/2pi = 8174.01 Hz" ; '
    '"lambda_eff = 2 lambda r_C^5/(3 pi^3/2) Int dq e^{-q^2 r_C^2} q^2 f(omega_L(q)) ... In the white noise case lambda_eff = lambda" ; '
    '"frequency of the X-rays (omega ~ 10^19 s^-1), at which the collapse noise is sampled, exceeds the cut-off frequency"',
    L_lpf=0.046, m_lpf=1.928)
src('review22', 'Carlesso, Donadi, Ferialdi, Paternostro, Ulbricht & Bassi, Nature Phys. 18, 243 (2022), arXiv:2203.04231 "Present status and future challenges of non-interferometric tests of collapse models"',
    'https://arxiv.org/pdf/2203.04231', 'Sec. 2.1-2.6, Eq. (3), Table 1',
    '"This in turn sets the bound lambda < 3.3 x 10^-11 s^-1 at rC = 10^-7 m" ; "Application of Eq. (4) to experimental data leads to lambda < 5.1 x 10^-8 s^-1 for the reference value of rC = 10^-7 m" ; '
    '"These tests constrained the CSL model to around lambda < 1.9 x10^-8 s^-1 at rC = 10^-7 m ... yielding an upper bound lambda < 2.0 x 10^-10 s^-1" ; '
    '"LIGO sets lambda < 1.0 x 10^-5 s^-1, AURIGA gives lambda < 4.6x10^-2 s^-1, while LISA Pathfinder provides lambda < 3.8x10^-9 s^-1" ; '
    '"yielding lambda < 4.1 x 10^-5 s^-1 and lambda < 6.7 x 10^-7 s^-1 for the CSL model at rC = 10^-7 m" [Paul trap ref 49 Pontin 2020; magneto-gravitational ref 50 Zheng 2020] ; '
    '"lambda < 5.2 x 10^-13 s^-1 at rC = 10^-7 m" (ref 51) ; "upper bound on the CSL rate of lambda < 10^-3 s^-1" (SQUID) ; '
    'Table 1: "Non-dissociation of hydrogen <1 ; CMB distorsion (COBE/FIRAS) <10^-1 ; Protons heating contribution to the CMB <10^-5 ; Heating in neutron stars <9.4x10^-7 ; Heating of the Intergalactic Medium <10^-8 ; Heating of Neptune <6.6x10^-11" ; '
    '"State-of-art interferometric experiments ... upper bound of lambda < 10^-7 s^-1 at rC = 10^-7 m"',
    bulk=3.3e-11, atoms_x=5.1e-8, cant=1.9e-8, cant_ml=2.0e-10, ligo=1.0e-5, auriga=4.6e-2, lpf=3.8e-9, paul=4.1e-5, mag=6.7e-7,
    ge=5.2e-13, squid=1e-3, h_diss=1.0, firas=1e-1, cmb_p=1e-5, ns=9.4e-7, igm=1e-8, neptune=6.6e-11, mwi=1e-7)
src('helou', 'Helou, Slagmolen, McClelland & Chen, PRD 95, 084054 (2017), arXiv:1606.03637', 'https://arxiv.org/pdf/1606.03637', 'abstract; p.1 col.2; Eq. (9)',
    '"we bound the CSL collapse rate to be at most (2.96 +- 0.12) x 10^-8 s^-1" ; "LISA pathfinder measures the relative acceleration noise between two free-falling test masses at a record accuracy of Sa = 5.2 +- 0.1 fm s^-2/sqrt(Hz) for frequencies between 0.7 mHz and 20 mHz"',
    lam=2.96e-8, Sa_fm=5.2, fmin=0.7e-3, fmax=20e-3)
src('armano18', 'Armano et al. (LPF), PRL 120, 061101 (2018) "Beyond the required LISA free-fall performance: new LISA Pathfinder results down to 20 uHz" (APS open-access PDF)',
    'https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.120.061101', 'abstract; Introduction; Discussion',
    '"the measured differential acceleration noise figure, which is at (1.74 +- 0.05) fm s^-2/sqrt(Hz) above 2 mHz and (6 +- 1) x 10 fm s^-2/sqrt(Hz) at 20 uHz" ; '
    '"limited by Brownian noise at (5.2 +- 0.1) fm s^-2/sqrt(Hz), for frequencies 1 mHz <~ f <~ 30 mHz; rising above the Brownian noise floor for frequencies f <~ 1 mHz, increasing to <~12 fm s^-2/sqrt(Hz) at f = 0.1 mHz" [2016 data, Ref. 3] ; '
    '"At frequencies above 10 mHz, the noise is dominated by displacement measurement noise, at a level of roughly 100 fm/sqrt(Hz)"',
    a_floor=1.74, a_20uHz=60.0, a_01mHz_2016=12.0, x_noise_fm=100.0)
src('cbfv16', 'Carlesso, Bassi, Falferi & Vinante, PRD 94, 124036 (2016), arXiv:1606.04581', 'https://arxiv.org/pdf/1606.04581', 'Sec. III',
    '"the effective force noise reaches a minimum SF(w)~95 fN Hz^-1/2 at w/2pi ~ 30-35 Hz" (LIGO) ; "Shh(w) = 1.6 x 10[..] Hz^-1/2 at w/2pi = 931 Hz" (AURIGA) ; "the bound from LISA Pathfinder lambda < 3 x 10^-8 s^-1"',
    ligo_f=(30.0, 35.0), auriga_f=931.0)
src('vinante16', 'Vinante et al., PRL 116, 090402 (2016), arXiv:1510.05791', 'https://arxiv.org/pdf/1510.05791', 'text p.2, p.3',
    '"The resonant frequency of the cantilever is f0 = 3084 Hz" ; "upper bound coming from matter-wave interferometry, lambda <= 5.0 x 10^-6 s^-1 for rC = 10^-7 m" ; "It also improves the bound coming from matter wave interferometry by 2 orders of magnitude"',
    f0=3084.0, mwi=5.0e-6)
src('vinante17', 'Vinante et al., PRL 119, 110401 (2017), arXiv:1611.09776', 'https://arxiv.org/pdf/1611.09776', 'abstract; text',
    '"thermal excess noise of unknown origin" ; "the cantilever is f0 = 8174.01 Hz" ; "rC = 10^-7 m this would imply lambda = 10^-7.7"', f0=8174.01, lam=10 ** -7.7)
src('vinante20', 'Vinante et al., PRL 125, 100404 (2020), arXiv:2002.09782', 'https://arxiv.org/pdf/2002.09782', 'p.2; p.4',
    '"The resonance frequency is f0 = 3532.7 Hz, while it was measured as f0 = 8174 Hz before attaching the multilayer test mass" ; '
    '"the absolute value of the excess force noise in the previous experiment, featuring the same cantilever, magnet, and SQUID, is compatible with the error bar of the new experiment. The improved bound on the CSL parameters arises entirely from the largest mass load and the specific multilayer structure."',
    f0=3532.7)
src('zheng20', 'Zheng et al., PRR 2, 013057 (2020), arXiv:1907.06896', 'https://arxiv.org/pdf/1907.06896', 'Introduction; Fig. caption',
    '"for rC = 10^-7 m, we estimate the upper bound lambda = 10^-6.4 s^-1 on the collapse rate at 95% confidence level" ; "resonance frequencies w2/2pi ~ 9.3 Hz and w1/2pi ~ 12.9 Hz"',
    f=(9.3, 12.9), lam_paper=10 ** -6.4)
src('bilardello16', 'Bilardello, Donadi, Vinante & Bassi, Physica A 462, 764 (2016), arXiv:1605.01891', 'https://arxiv.org/pdf/1605.01891', 'Sec. II',
    '"Step 1: The harmonic trap is removed and the gas evolves freely for a relatively long time, Dt1 = 1.1 s" ; "the delta-kick ... lasts for a short time, dt2 ~ 35 ms" ; "time, Dt3 = 1.8 s. The position variance of the gas is"',
    t1=1.1, t2=0.035, t3=1.8)
src('adler_vinante18', 'Adler & Vinante, PRA 97, 052119 (2018), arXiv:1801.06857 "Bulk heating effects as tests for collapse models"', 'https://arxiv.org/pdf/1801.06857', 'p.2-3',
    '"the residual heating from unknown sources is limited to roughly 10^-11 W/kg. This means that the effective lambda for bulk heating is at most 3.1 x 10^-11 s^-1" ; '
    '"heating takes place only if the noise has frequency components of at least wL(|q| = rc^-1) ~ vs |q| ~ 0.4 x 10^11 s^-1 ... vs the speed of sound (which for this estimate we have taken as 4000 m/s characteristic of copper at low temperature)"',
    lam=3.1e-11, vs=4000.0)
src('schrinski23_conv', 'Schrinski et al. PRL 130, 133604 (2023) (same PDF), CSL mapping and diffusion rate', 'https://www.klaus-hornberger.de/pdf/prl19.pdf', 'p.133604-2 col.2; p.133604-3',
    '"usually parametrized in terms of a rate per atomic mass unit lambda_CSL = (1 u/me)^2/tau_e, and a localization length r_CSL = hbar/(sqrt2 sigma_q)" ; "we find from Eq. (4) a maximal diffusion rate of Gamma_max = 3.5 x 10^13/tau_e at a critical length scale hbar/sigma_q ~ 0.5 um"',
    G_tau=3.5e13, hbar_over_sq=0.5e-6)
src('schrinski23', 'Schrinski, Yang, von Luepke, Bild, Chu, Hornberger, Nimmrichter & Fadel, PRL 130, 133604 (2023) "Macroscopic quantum test with bulk acoustic wave resonators"',
    'https://www.klaus-hornberger.de/pdf/prl19.pdf', 'p.133604-3 "Experiment" paragraph',
    '"It oscillates over the bulk length L = 435 um at angular frequency w = 2pi . 5.961 GHz with wavelength lambda = 1.8 um" ; "Assuming a homogenous mass density 3.98 g/cm3, and given the transverse waist w0 = 27 um, we obtain the effective oscillator mass meff = pi w0^2 L rho/4 = 1.0 ug" ; '
    '"a small one-phonon steady-state population of 1.6 +- 0.2 %, by the relaxation time T1 = 85.8 +- 1.5 us"',
    f=5.961e9, lam_ac=1.8e-6, rho=3980.0, n=0.016, dn=0.002, T1=85.8e-6, w0=27e-6, L=435e-6)
src('hbar_next', 'Ultracold mechanical quantum sensor for tests of new physics, arXiv:2507.02653v3 (Table 2)', 'https://arxiv.org/html/2507.02653v3', 'Table 2',
    '"Frequency 5 GHz | 3-10 GHz | 15-800 MHz" (current / next generation / MHz device) ; "Estimated CSL rate 5.7x10^-8 s^-1 | 10^-10 s^-1 | 10^-13 s^-1" ; "Estimated CSL radius 3.0x10^-7 m" (all columns)',
    note='used only as a PROPOSAL for 4.6 MHz-5.6 GHz; bound itself dropped (r_C = 3e-7 m, not 1e-7 m)')
src('kaufmann17', 'H. Kaufmann, PhD thesis, Univ. Mainz (2017)', 'https://www.quantenbit.physik.uni-mainz.de/files/2019/03/PhD_HKaufmann.pdf', 'Abstract',
    '"Zeeman sublevels ... of 40Ca+ ions" ; "We achieve a low heating rate of three motional quanta per second at a radial mode frequency of 2pi x 4.6 MHz"',
    A=40, f=4.6e6, ndot=3.0)
src('bruzewicz15', 'Bruzewicz, Sage & Chiaverini, PRA 91, 041402 (2015), arXiv:1412.5119', 'https://arxiv.org/pdf/1412.5119', 'abstract; Table 1',
    '"a 88Sr+ ion trapped in a surface-electrode ion trap" ; "trap frequencies and temperatures between ~0.6 and 1.5 MHz and ~4 and 295 K" ; '
    'Table 1 (fit n(T) = n0(1+(T/T0)^beta)): "0.60 21(2) 14(2) 1.6(1); 0.80 18(2) 21(4) 1.6(1); 1.0 9(2) 13(4) 1.5(1); 1.3 7(1) 19(2) 1.6(1); 1.5 8(1) 24(4) 1.5(1)"',
    A=88, rows=[(0.60e6, 21, 2, 14, 1.6), (0.80e6, 18, 2, 21, 1.6), (1.0e6, 9, 2, 13, 1.5), (1.3e6, 7, 1, 19, 1.6), (1.5e6, 8, 1, 24, 1.5)], Tmin=4.0)
src('igex17', 'Piscicchia et al., Entropy 19, 319 (2017), arXiv:1710.01973 (IGEX data)', 'https://arxiv.org/pdf/1710.01973', 'Sec. 3, Eq. (18)',
    '"We restricted our analysis to the energy range DE = (14.5 / 48.5) keV" ; "setting rC to the value 10^-7 m, corresponding to a probability level of 95% lambda <= 6.8 . 10^-12 s^-1 mass prop."',
    Emin=14.5, Emax=48.5, lam=6.8e-12)
src('majorana22', 'Arnquist et al. (Majorana), PRL 129, 080401 (2022), arXiv:2202.01343', 'https://arxiv.org/pdf/2202.01343', 'abstract; Eqs. (5)-(6) and text',
    '"no evidence of such radiation in the 19-100 keV range" ; "lambda/rC^2 = (5.15 +- 0.16) x 10^-6 s^-1 m^-2" (n-m-p) ; "When only the 30 quasi-free electrons are considered ... lambda/rC^2 = (17.4 +- 0.5) s^-1 m^-2" ; "(4.94 +- 0.15) x 10^-1" (m-p coherent nuclei) ; '
    '"Extensions of the CSL models such as the colored noise models are required to evade such stringent constraints"',
    Emin=19.0, Emax=100.0, nmp=5.15e-6, mp_e=17.4, mp_n=0.494)
src('donadi21', 'Donadi et al., Eur. Phys. J. C 81, 773 (2021), arXiv:2107.11237 "Novel CSL bounds from the noise-induced radiation emission from atoms"', 'https://arxiv.org/pdf/2107.11237', 'Fig. 1 caption; Eq. (15)',
    '"Measured X-rays spectrum ... in the energy range DE = (1000 / 3800) keV" ; "We obtain the upper limit on lambda: lambda < 5.2 . 10^-13 s^-1"', Emin=1000.0, Emax=3800.0, lam=5.2e-13)

# ------------------------------------------------------------------ physics helpers (DERIVED)
def S_slab(lam, rho, A):
    """Two-sided CSL force PSD on a rigid body with two flat faces of area A normal to x, thickness >> r_C (DERIVED
    from Eq. 13 of ccb18: Int|k_x mu|^2 e^{-k^2 rc^2} -> rho^2 (2pi)^2 A * 2 sqrt(pi)/rc)."""
    return 8 * pi * HBAR ** 2 * lam * rho ** 2 * A * RC ** 2 / M0 ** 2

# Calibration check: reproduce Helou's LPF number from the same formula (one-sided Sa = 2 x two-sided; relative accel of
# two independent masses doubles the PSD).
L_lpf, m_lpf = S['ccb18']['values']['L_lpf'], S['ccb18']['values']['m_lpf']
rho_lpf, A_lpf = m_lpf / L_lpf ** 3, L_lpf ** 2
def lam_lpf_from_asd(asd_fm):
    Sa1 = (asd_fm * 1e-15) ** 2
    return Sa1 * m_lpf ** 2 / (2 * 2 * S_slab(1.0, rho_lpf, A_lpf))
lam_check = lam_lpf_from_asd(S['helou']['values']['Sa_fm'])
calib = dict(formula_lambda_from_5p2fm=lam_check, helou_cited=S['helou']['values']['lam'], ratio=lam_check / S['helou']['values']['lam'],
             label='DERIVED check: Eq.13 slab limit + LPF geometry reproduces Helou 2017 to <1 %, so lambda scales as ASD^2')

# LPF 2018 curve (DERIVED: lambda = lambda_Helou * (ASD/5.2)^2; ASSUMED log-log interpolation between quoted points;
# above 10 mHz ASSUMED quadrature sum of 1.74 fm s^-2 floor and 100 fm/rtHz displacement noise x (2 pi f)^2)
a = S['armano18']['values']
def lpf18_asd(f):
    pts_f = np.array([20e-6, 0.1e-3, 2e-3]); pts_a = np.array([a['a_20uHz'], a['a_01mHz_2016'], a['a_floor']])
    if f <= 2e-3:
        return float(10 ** np.interp(np.log10(f), np.log10(pts_f), np.log10(pts_a)))
    xa = a['x_noise_fm'] * (2 * pi * f) ** 2 if f > 10e-3 else 0.0
    return sqrt(a['a_floor'] ** 2 + xa ** 2)
f_lpf = np.logspace(np.log10(20e-6), np.log10(30e-3), 200)
lam_lpf18 = np.array([S['helou']['values']['lam'] * (lpf18_asd(f) / S['helou']['values']['Sa_fm']) ** 2 for f in f_lpf])

# Cold atoms: half-power frequency of the free-evolution kernel (DERIVED), T = Dt1 + Dt3 (ASSUMED: kick neglected)
b = S['bilardello16']['values']; T_free = b['t1'] + b['t2'] + b['t3']
def kern(w):
    re = quad(lambda s: (T_free - s) * np.cos(w * s), 0, T_free, limit=200)[0]
    im = quad(lambda s: (T_free - s) * np.sin(w * s), 0, T_free, limit=200)[0]
    return (re ** 2 + im ** 2) / (T_free ** 4 / 4)
w_half = brentq(lambda w: kern(w) - 0.5, 1e-3, 10 / T_free)
f_atoms = w_half / (2 * pi)

# Bulk heating band: weight q^4 exp(-q^2 rc^2) dq (3D), w = vs q; half-maximum points (DERIVED)
wfun = lambda x: x ** 4 * np.exp(-x ** 2); wmax = wfun(sqrt(2))
x_lo = brentq(lambda x: wfun(x) - wmax / 2, 0.1, sqrt(2)); x_hi = brentq(lambda x: wfun(x) - wmax / 2, sqrt(2), 5)
vs = S['adler_vinante18']['values']['vs']
f_bulk = (x_lo * vs / (2 * pi * RC), x_hi * vs / (2 * pi * RC)); f_bulk_peak = sqrt(2) * vs / (2 * pi * RC)

# HBAR single phonon mode (DERIVED): excitation rate per unit lambda for a mode u = z sin(k0 z) exp(-r^2/w0^2)
#   Gamma = D/(2 m_eff hbar w),  D = hbar^2 lam rc^3/(pi^1.5 m0^2) Int d^3k |k.(mu u)~|^2 e^{-k^2 rc^2}
#   -> Gamma/lam = 4 pi^1.5 hbar rho rc^3 k0^2 exp(-k0^2 rc^2)/(m0^2 w)   (independent of w0, L)
hb = S['schrinski23']['values']; k0 = 2 * pi / hb['lam_ac']; w_hb = 2 * pi * hb['f']
G_per_lam = 4 * pi ** 1.5 * HBAR * hb['rho'] * RC ** 3 * k0 ** 2 * exp(-(k0 * RC) ** 2) / (M0 ** 2 * w_hb)
Gamma_max = (hb['n'] + 2 * hb['dn']) / hb['T1']          # ASSUMED conservative: all population attributed, +2 sigma
lam_hbar = Gamma_max / G_per_lam
# cross-check vs the paper's own convention: Gamma_max = 3.5e13/tau_e at hbar/sigma_q ~ 0.5 um, lambda_CSL=(u/me)^2/tau_e
#   our formula at r = 0.5um/sqrt2 for comparison (DERIVED; flags a convention factor)
sc = S['schrinski23_conv']['values']; r_alt = sc['hbar_over_sq'] / sqrt(2)
G_per_lam_alt_ours = 4 * pi ** 1.5 * HBAR * hb['rho'] * r_alt ** 3 * k0 ** 2 * exp(-(k0 * r_alt) ** 2) / (M0 ** 2 * w_hb)
G_per_lam_alt_paper = sc['G_tau'] * S['nist']['values']['me_u'] ** 2      # Gamma/lambda = 3.5e13 * (me/u)^2
hbar_check = dict(r=r_alt, ours=G_per_lam_alt_ours, paper=G_per_lam_alt_paper, ratio=G_per_lam_alt_ours / G_per_lam_alt_paper,
                  label='DERIVED cross-check of the single-mode formula against the paper at its own optimum length')

# Trapped ions (DERIVED): point mass m = A u, per-axis CSL heating (1/3)(3 lam hbar^2 m/(4 m0^2 rc^2)) [review Eq. 3],
#   ndot = that/(hbar w)  =>  lam <= 4 m0^2 rc^2 w ndot/(hbar m)
def lam_ion(A, f, ndot):
    return 4 * M0 ** 2 * RC ** 2 * (2 * pi * f) * ndot / (HBAR * A * U)
ions = []
k = S['kaufmann17']['values']; ions.append(dict(name=f"40Ca+ heating 4.6 MHz (Kaufmann 2017)", f=k['f'], lam=lam_ion(k['A'], k['f'], k['ndot'])))
bz = S['bruzewicz15']['values']
for (f, n0, dn0, T0, beta) in bz['rows']:
    nd = (n0 + dn0) * (1 + (bz['Tmin'] / T0) ** beta)    # DERIVED: heating at ~4 K from the fit, +1 sigma on n0
    ions.append(dict(name=f"88Sr+ heating {f/1e6:.2f} MHz (Bruzewicz 2015)", f=f, lam=lam_ion(bz['A'], f, nd), ndot=nd))

# ------------------------------------------------------------------ the band list
line = lambda f: (f * (1 - 1e-3), f * (1 + 1e-3))     # ASSUMED Stage 31 narrow-line convention
R = S['review22']['values']
B = []
def band(name, fmin, fmax, lam, label, keys, kind, mp=True, note='', curve=None, stage='33'):
    B.append(dict(name=name, fmin=float(fmin), fmax=float(fmax), lam=float(lam), label=label, sources=keys, kind=kind,
                  mass_prop=mp, note=note, curve=curve, stage=stage))
band('LISA Pathfinder 2016 (Helou 2017)', 0.7e-3, 20e-3, S['helou']['values']['lam'], 'CITED', ['helou'], 'free-fall force noise', stage='31')
band('LISA Pathfinder 2018 (Armano 2018 ASD -> lambda, DERIVED)', 20e-6, 30e-3, lam_lpf18.min(), 'DERIVED (ASD^2 scaling of Helou); interpolation ASSUMED',
     ['armano18', 'helou', 'ccb18'], 'free-fall force noise', note=f'min {lam_lpf18.min():.2e} at 2-10 mHz; {lam_lpf18[0]:.2e} at 20 uHz; review quotes 3.8e-9',
     curve=dict(f=f_lpf.tolist(), lam=lam_lpf18.tolist()))
band('Cold atoms, 87Rb pK free expansion (Kovachy data, Bilardello 2016)', 1e-6, f_atoms, R['atoms_x'], 'CITED lambda; band DERIVED (kernel half-power, T=2.94 s)',
     ['review22', 'bilardello16'], 'position diffusion', note='responds to all f below the half-power point; lower edge = plot edge')
for f0 in S['zheng20']['values']['f']:
    band(f'Diamagnetic levitated micro-sphere {f0} Hz (Zheng 2020)', *line(f0), R['mag'], 'CITED (review value; paper says 10^-6.4)', ['review22', 'zheng20'], 'oscillator')
band('LIGO test-mass force noise 30-35 Hz', *S['cbfv16']['values']['ligo_f'], R['ligo'], 'CITED', ['review22', 'cbfv16'], 'interferometer force noise')
band('AURIGA bar 931 Hz', *line(S['cbfv16']['values']['auriga_f']), R['auriga'], 'CITED', ['review22', 'cbfv16'], 'resonant bar')
band('Nanocantilever 3084 Hz (Vinante 2016)', *line(S['vinante16']['values']['f0']), S['vinante16']['values']['mwi'] / 100, 'DERIVED from text ("improves 5.0e-6 by 2 orders")', ['vinante16'], 'oscillator')
band('Multilayer cantilever 3532.7 Hz (Vinante 2020)', *line(S['vinante20']['values']['f0']), R['cant_ml'], 'CITED', ['review22', 'vinante20'], 'oscillator')
band('Cantilever 8174 Hz (Vinante 2017)', *line(S['vinante17']['values']['f0']), R['cant'], 'CITED', ['review22', 'vinante17'], 'oscillator', stage='31')
for io in ions:
    band(io['name'], *line(io['f']), io['lam'], 'DERIVED (review Eq. 3, point mass, per axis)', ['kaufmann17' if 'Ca' in io['name'] else 'bruzewicz15', 'review22', 'nist'], 'trapped-ion heating')
band('Bulk heating of solids (10 pW/kg residual; Adler-Vinante 2018)', *f_bulk, S['adler_vinante18']['values']['lam'], 'CITED lambda; band DERIVED (q^4 e^{-q^2 rc^2}, v_s=4000 m/s)',
     ['adler_vinante18', 'ccb18', 'review22'], 'phonon heating', note=f'peak {f_bulk_peak:.2e} Hz; review quotes 3.3e-11')
band('HBAR phonon mode 5.961 GHz (Schrinski 2023)', *line(hb['f']), lam_hbar, 'DERIVED (single-mode CSL excitation, all population attributed)', ['schrinski23', 'ccb18', 'nist'], 'phonon-mode heating')
band('IGEX Ge X-rays 14.5-48.5 keV (Piscicchia 2017)', S['igex17']['values']['Emin'] * keV2Hz, S['igex17']['values']['Emax'] * keV2Hz, S['igex17']['values']['lam'], 'CITED', ['igex17'], 'radiation')
mj = S['majorana22']['values']
band('Majorana 19-100 keV, m-p electrons only', mj['Emin'] * keV2Hz, mj['Emax'] * keV2Hz, mj['mp_e'] * RC ** 2, 'CITED', ['majorana22'], 'radiation', stage='31')
band('Majorana 19-100 keV, m-p coherent nuclei', mj['Emin'] * keV2Hz, mj['Emax'] * keV2Hz, mj['mp_n'] * RC ** 2, 'CITED', ['majorana22'], 'radiation', stage='31')
band('Majorana 19-100 keV, NON-mass-proportional', mj['Emin'] * keV2Hz, mj['Emax'] * keV2Hz, mj['nmp'] * RC ** 2, 'CITED', ['majorana22'], 'radiation', mp=False, stage='31',
     note='different model (n-m-p); not comparable with mechanical bounds; excluded from coverage stats')
band('Ge underground 1000-3800 keV (Donadi 2021)', S['donadi21']['values']['Emin'] * keV2Hz, S['donadi21']['values']['Emax'] * keV2Hz, S['donadi21']['values']['lam'], 'CITED', ['donadi21', 'review22'], 'radiation')

# ------------------------------------------------------------------ coverage statistics (mass-proportional only)
def union(iv):
    iv = sorted(iv); out = []
    for lo, hi in iv:
        if out and lo <= out[-1][1]: out[-1][1] = max(out[-1][1], hi)
        else: out.append([lo, hi])
    return out
mpB = [b_ for b_ in B if b_['mass_prop']]
U_all = union([(b_['fmin'], b_['fmax']) for b_ in mpB])
U_31 = union([(b_['fmin'], b_['fmax']) for b_ in mpB if b_['stage'] == '31'])
GAP = (0.02, mj['Emin'] * keV2Hz)            # Stage 31 gap: 20 mHz .. 19 keV
def dec_in(Uv, lo, hi):
    return sum(max(0.0, log10(min(h_, hi)) - log10(max(l_, lo))) for l_, h_ in Uv if h_ > lo and l_ < hi)
gap_dec = log10(GAP[1] / GAP[0])
cov_gap = dec_in(U_all, *GAP); cov_gap31 = dec_in(U_31, *GAP)
SPAN = (20e-6, S['donadi21']['values']['Emax'] * keV2Hz)   # lowest measured LPF point .. highest Ge energy
open_windows = []
prev = SPAN[0]
for lo, hi in U_all:
    if lo > prev * 1.0000001 and lo > SPAN[0]:
        open_windows.append((prev, min(lo, SPAN[1])))
    prev = max(prev, hi)
if prev < SPAN[1]: open_windows.append((prev, SPAN[1]))
open_windows = [(lo, hi) for lo, hi in open_windows if log10(hi / lo) > 1e-6]

# per-decade best bound in the gap (envelope)
def best_at(f):
    c = [b_['lam'] if not b_['curve'] else float(np.interp(np.log10(f), np.log10(b_['curve']['f']), b_['curve']['lam']))
         for b_ in mpB if b_['fmin'] <= f <= b_['fmax']]
    return min(c) if c else None
tight = min(mpB, key=lambda b_: b_['lam']); weak = max(mpB, key=lambda b_: b_['lam'])

unassigned = [
    dict(item='Heating of Neptune', lam=R['neptune'], why='broadband white-noise heating; no fetched sound speed / frequency mapping', src='review22'),
    dict(item='Heating of the intergalactic medium', lam=R['igm'], why='free-proton heating; frequency = recoil/Doppler transfer, not given', src='review22'),
    dict(item='Heating in neutron stars (PSR J1840-1419)', lam=R['ns'], why='no frequency mapping in fetched source', src='review22'),
    dict(item='Proton heating contribution to the CMB', lam=R['cmb_p'], why='emission redshifted from unknown z; band not assignable', src='review22'),
    dict(item='CMB distortion (COBE/FIRAS)', lam=R['firas'], why='same: emission epoch spread, band not assignable', src='review22'),
    dict(item='Non-dissociation of cosmic hydrogen', lam=R['h_diss'], why='energy transfer >= 13.6 eV implied but the source gives no spectral treatment', src='review22'),
    dict(item='SQUID supercurrent decay', lam=R['squid'], why='pair breaking needs noise above 2*Delta/h, gap value not fetched; source warns bound could be weaker', src='review22'),
    dict(item='Matter-wave interferometry', lam=R['mwi'], why='decoherence integrates over the flight time and recoil spectrum; no band', src='review22'),
    dict(item='Paul-trap levitated nanoparticle (Pontin 2020)', lam=R['paul'], why='secular frequency not stated in fetched text (only "sample rate 500-800 Hz")', src='review22'),
]
dropped = [
    dict(item='Fu 1997 / TWIN Ge 11 keV (lambda <= 2.2e-10, corrected 8e-10)', why='source says data "underestimated by a factor of about 50 at 10 keV"; energy band of the corrected analysis not stated'),
    dict(item='HBAR next-generation paper (arXiv:2507.02653) lambda < 5.7e-8', why='given only at r_CSL = 3.0e-7 m, not the Stage 31 r_C = 1e-7 m'),
    dict(item='Earth heat balance (~1e-11, Adler-Vinante)', why='stated as approximate ("~"), no frequency mapping'),
    dict(item='LIGO outside 30-35 Hz (10-1000 Hz)', why='force-noise spectrum not fetched; only the 30-35 Hz minimum is quoted'),
    dict(item='Torsion pendulums, atomic clocks, NMR / comagnetometers, axion haloscopes, SQUID noise', why='no CSL-parametrised bound found in a source fetched this session; CSL couples to mass density, so clocks/NMR see it only through motional heating (covered by the ion entries)'),
    dict(item='Cold-atom energy-based bound lambda < 1e-7+-1', why='uncertain by +-1 decade in source; position-diffusion bound used instead'),
]

# ------------------------------------------------------------------ outputs
os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
with open(os.path.join(ROOT, 'results', 'stage33_sources.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for kk, v in S.items(): w.writerow([kk, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
with open(os.path.join(ROOT, 'results', 'stage33_bounds.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['name', 'kind', 'fmin_Hz', 'fmax_Hz', 'lambda_max_per_s_rC_1e-7m', 'label', 'mass_proportional', 'stage', 'sources', 'note'])
    for b_ in sorted(B, key=lambda x: x['fmin']):
        w.writerow([b_['name'], b_['kind'], f"{b_['fmin']:.4g}", f"{b_['fmax']:.4g}", f"{b_['lam']:.3g}", b_['label'], b_['mass_prop'], b_['stage'], ';'.join(b_['sources']), b_['note']])
with open(os.path.join(ROOT, 'results', 'stage33_lpf18_curve.csv'), 'w') as fh:
    fh.write('f_Hz,lambda_max_per_s\n'); [fh.write(f'{f:.6e},{l:.6e}\n') for f, l in zip(f_lpf, lam_lpf18)]
out = dict(stage=33, title='Resonant-release frequency coverage (colored-CSL parametrization, r_C = 1e-7 m)',
           constants=dict(hbar=HBAR, h=H, u=U, eV=EV, m0='= 1 u (ASSUMED)', rC=RC), calibration=calib,
           lpf18=dict(lam_min=float(lam_lpf18.min()), lam_20uHz=float(lam_lpf18[0]), review=R['lpf'], ratio_to_review=float(lam_lpf18.min() / R['lpf'])),
           atoms=dict(T=T_free, f_half=f_atoms), bulk=dict(f_lo=f_bulk[0], f_hi=f_bulk[1], f_peak=f_bulk_peak, x_lo=x_lo, x_hi=x_hi),
           hbar_mode=dict(G_per_lambda=G_per_lam, Gamma_max=Gamma_max, lam=lam_hbar, k0rc=k0 * RC, crosscheck=hbar_check),
           ions=ions, bands=[{k_: v_ for k_, v_ in b_.items() if k_ != 'curve'} for b_ in B],
           coverage=dict(gap=GAP, gap_decades=gap_dec, covered_decades_stage31=cov_gap31, covered_decades_now=cov_gap,
                         union=U_all, open_windows=[(lo, hi, log10(hi / lo)) for lo, hi in open_windows],
                         total_span=SPAN, total_covered_decades=dec_in(U_all, *SPAN)),
           tightest=dict(name=tight['name'], lam=tight['lam'], f=(tight['fmin'], tight['fmax'])),
           weakest=dict(name=weak['name'], lam=weak['lam'], f=(weak['fmin'], weak['fmax'])),
           unassigned=unassigned, dropped=dropped, sources=S, runtime_s=time.time() - t0)
json.dump(out, open(os.path.join(ROOT, 'results', 'stage33.json'), 'w'), indent=1, default=float)

print(f"calibration: formula {lam_check:.3e} vs Helou {S['helou']['values']['lam']:.3e} (ratio {calib['ratio']:.4f})")
print(f"LPF18: min {lam_lpf18.min():.3e} (review 3.8e-9), at 20 uHz {lam_lpf18[0]:.3e}, at 30 mHz {lam_lpf18[-1]:.3e}")
print(f"cold atoms: T={T_free:.3f} s, f_half={f_atoms:.4f} Hz")
print(f"bulk band {f_bulk[0]:.3e}-{f_bulk[1]:.3e} Hz (peak {f_bulk_peak:.3e}); x {x_lo:.3f}-{x_hi:.3f}")
print(f"HBAR: k0 rc={k0*RC:.3f}, Gamma/lambda={G_per_lam:.3e} s^-1 per s^-1, Gamma_max={Gamma_max:.1f}/s -> lambda<={lam_hbar:.3e}")
print(f"HBAR cross-check at r={r_alt:.3e}: ours {G_per_lam_alt_ours:.3e} vs paper {G_per_lam_alt_paper:.3e} ratio {hbar_check['ratio']:.3f}")
for io in ions: print(f"  ion {io['name']}: lambda <= {io['lam']:.3e}")
for b_ in sorted(B, key=lambda x: x['fmin']): print(f"  {b_['fmin']:.3e}-{b_['fmax']:.3e} Hz  {b_['lam']:.3e}  [{b_['stage']}] {b_['name']}")
print(f"gap {GAP[0]:.3g}-{GAP[1]:.3e} Hz = {gap_dec:.2f} decades; covered stage31 {cov_gap31:.4f}, now {cov_gap:.3f} decades")
print("open windows (>1e-3 dec):")
for lo, hi in open_windows:
    if log10(hi / lo) > 1e-3: print(f"   {lo:.4g} - {hi:.4g} Hz  ({log10(hi/lo):.2f} dec)")
print(f"tightest {tight['name']} {tight['lam']:.3e}; weakest {weak['name']} {weak['lam']:.3e}")
