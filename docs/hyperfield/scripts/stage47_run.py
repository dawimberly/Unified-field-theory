"""Stage 47: hysteresis (lock/unlock asymmetry) vs detailed balance in lab thermometry.

Stage 31 definition (re-used unchanged): emission/absorption rate ratio per quantum E/A = e^{-x}/(1+eta),
x = h nu / kT  ->  steady-state occupation n = 1/((1+eta) e^x - 1) = 1/(e^{x+ln(1+eta)} - 1).
FIRAS |mu| < 9e-5 (95%)  ->  |eta| < 9.0e-5 (Stage 31).

Lab readings (DERIVED from that definition):
 * any mode at x << 1:  n ~ 1/(x+eta), so a thermometer reading T_meas = T/(1+eta/x).  Fractional deviation eps -> |eta| < eps*x.
 * two-level Boltzmann ratio at any x:  ln(N_l/N_u) = x + ln(1+eta), so a population-ratio test at fractional precision eps gives |eta| < eps.
"""
import json, csv, numpy as np
h = 6.62607015e-34; kB = 1.380649e-23; c2 = 1.438776877  # cm K (second radiation constant hc/k, derived from exact SI h,c,k)
OUT = '/workspace/uft-stage33/docs/hyperfield/'
S = []
def src(key, ref, url, loc, quote, **vals):
    S.append(dict(key=key, ref=ref, url=url, location=loc, quote=quote, values=json.dumps(vals)))
src('firas', 'Fixsen et al., ApJ 473, 576 (1996), arXiv:astro-ph/9605054', 'https://arxiv.org/pdf/astro-ph/9605054',
    'abstract; Sec. 1 (lines 76, 84 of text)', '"|y| <15 x 10^-6 and |mu| <9 x 10^-5 (95% CL)" ; "this frequency range (2 to 21 cm-1)" ; "2.728 +- 0.004 K (95% CL)"',
    mu95=9e-5, nu_lo_cm=2, nu_hi_cm=21, T=2.728)
src('stage31', 'This repo, HYPERFIELD_RESULTS.md Stage 31 section (7)', 'docs/hyperfield/HYPERFIELD_RESULTS.md', 'Stage 31, element (7)',
    '"If E/A = e^(-x)/(1+eta), then n = 1/((1+eta)e^x - 1). That is a Bose-Einstein spectrum with mu = ln(1+eta)." ; "|eta| < 9.0x10^-5"', eta=9e-5)
src('qu2017', 'Qu, Benz, Coakley et al., Metrologia 54, 549 (2017), NIST pub 923355 / PMC5621608', 'https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=923355',
    'abstract; Sec. 4 (fit, Fig. 5-6)', '"sensing resistor immersed in a triple-point-of-water cell is compared to the noise power of a quantum-accurate pseudo-random noise waveform" ; '
    '"k = 1.380 649 7(37) x 10-23 J K-1 with a relative standard uncertainty of 2.7 x 10-6 and a relative offset of 0.89 x 10-6 from the CODATA 2014 recommended value" ; '
    '"= a0 + a2 f2 + a4 f4 ..." ; "bandwidths of 10 kHz to different upper cut-off frequencies fmax" ; "d = 2 model when fmax = 368.75 kHz"',
    u_rel=2.7e-6, offset=0.89e-6, f_lo=10e3, f_max=368.75e3, T=273.16)
src('purdy2015', 'Purdy, Yu, Kampel et al., PRA 92, 031802(R) (2015), arXiv:1406.7247', 'https://arxiv.org/pdf/1406.7247',
    'p.1-3 text; Fig. 4 discussion', '"Rsa = e^(hbar wm/kb T) = (n+1)/n" ; "temperature, T0, that is set in the range of 4.8 to 50 K" ; "(3,2) mode ... wm/2pi = 2.637 MHz" ; '
    '"extrapolated temperatures are compared to the temperature as measured by a silicon diode thermometer attached to the cryostat. The results agree within the statistical uncertainty with an average deviation of less than 10%"',
    f=2.637e6, T_lo=4.8, T_hi=50, dev=0.10)
src('li2025', 'Li, Wang, ... Tennyson, Li, Hu, Sci. Adv. 11, eadz6560 (2025)', 'https://www.science.org/doi/full/10.1126/sciadv.adz6560',
    'Results (Fig. 3 discussion, LRT section); abstract',
    '"Line intensities in the (3-0) vibrational band of 12C16O at 296 K" ; "UCL2025 (cyan line) reduces these deviations to <0.02%" ; "Discrepancies of up to 0.02% in transition probability ratios" ; '
    '"For the R(0)/R(20) pair, the lower-state energies are 0 and 806.3828 cm-1" ; "we derived an empirical model of the Einstein-A ratios of transitions through a fitting" ; '
    '"The blue line indicates the temperature (T) calibrated from the readings of platinum sensors"', T=296.0, dev_abinitio=2e-4, E_hi_cm=806.3828)
R = {}
# --- FIRAS (Stage 31) ---
R['firas'] = dict(eta95=9e-5, x_range=[c2*2/2.728, c2*21/2.728])
# --- JNT: ratio(f) = 1/(1+eta/x(f)) ~ 1 - (eta kT/h)/f ; fitted with a0 + a2 f^2 over 10..368.75 kHz ---
T = 273.16; f = np.arange(10e3, 368.75e3 + 1, 6.25e3 / 4)
A = np.vstack([np.ones_like(f), f**2]).T
coef = np.linalg.lstsq(A, 1.0 / f, rcond=None)[0]     # response of a0 to a pure 1/f term (uniform weights, ASSUMED)
a0_per_c = coef[0]                                  # a0 = a0_per_c * c, c = eta kT/h
resp = a0_per_c * 368.75e3                          # in units of 1/f_max
a0_lim = 0.89e-6 + 1.96 * 2.7e-6                    # 95% two-sided on |a0 - 1| (offset + 1.96 sigma)
c_lim = a0_lim / a0_per_c
eta_jnt = c_lim * h / (kB * T)
resid_lo = (1/10e3 - (coef[0] + coef[1]*(10e3)**2)) * c_lim   # residual the fit would leave at 10 kHz at this eta
R['jnt'] = dict(T=T, x_fmax=h*368.75e3/(kB*T), a0_response_in_units_of_c_over_fmax=float(resp), a0_lim95=a0_lim,
                eta95=float(eta_jnt), fit_residual_at_10kHz=float(resid_lo),
                note='a0 of a quadratic fit to c/f over 10-368.75 kHz equals %.2f c/f_max; |a0-1| < %.2e -> |eta| < %.1e. The fit residual at 10 kHz would be only %.1e, so the a0 bound is the conservative one.' % (resp, a0_lim, eta_jnt, resid_lo),
                eta_firas_effect='at eta = 9e-5 the predicted noise ratio at 368.75 kHz would be %.2e of nominal' % (1/(1 + 9e-5/(h*368.75e3/(kB*T)))))
# --- Purdy Raman ratio: T_meas = T0/(1+eta/x); conservative: |eta| <= dev95 * x(T_lo) ---
x_lo = h*2.637e6/(kB*4.8); x_hi = h*2.637e6/(kB*50)
R['purdy'] = dict(x_range=[x_hi, x_lo], dev95=0.20, eta95=0.20*x_lo,
                  note='ASSUMED 95% = 2 x stated 10% average deviation; average over unknown T points bounded by the largest-x point (4.8 K), so conservative')
# --- CO line-ratio vs ab initio (Boltzmann factor at x up to 3.9; same-lower-state R/P pairs are eta-free) ---
R['co'] = dict(x_range=[0.0, c2*806.3828/296.0], eta95=2e-4,
               note='ASSUMED: eta multiplies every population ratio of two different rotational levels once (single-step reading); measured ratios agree with UCL2025 ab initio to <0.02%, so |eta| < 2e-4. The 6 ppm LRT-vs-sensor agreement is NOT used: the Einstein-A ratios were fitted to the same lab data at 296 K with the sensor temperature, so a constant eta is absorbed (circular).')
# --- eta proportional to x (pure temperature rescale) ---
R['eta_prop_x'] = dict(note='if eta = eps*x the effect is a rescaled temperature; only cross-thermometer comparisons see it: JNT vs CODATA2014 (acoustic-gas dominated) |eps| < %.1e; Purdy |eps| < 0.2' % a0_lim, eps95_jnt=a0_lim)
R['verdict'] = dict(
    freq_independent_eta=dict(best_lab='JNT', eta95=float(eta_jnt), factor_vs_firas=float(9e-5/eta_jnt)),
    x_ge_1_only=dict(best_lab='CO line ratios vs ab initio', eta95=2e-4, factor_vs_firas=2e-4/9e-5, beats_firas=False))
json.dump(R, open(OUT+'results/stage47.json', 'w'), indent=1)
with open(OUT+'results/stage47_sources.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(S[0].keys())); w.writeheader(); [w.writerow(s) for s in S]
print(json.dumps(R, indent=1))

# ---------------- figure ----------------
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(14, 5.4))
a = ax[0]
def box(xr, e, col, lab):
    a.fill_between(xr, [e, e], [1, 1], color=col, alpha=0.25); a.plot(xr, [e, e], color=col, lw=3, label=lab)
box(R['firas']['x_range'], 9e-5, 'k', 'FIRAS CMB (Stage 31): 9.0e-5')
box([R['jnt']['x_fmax']/36.875, R['jnt']['x_fmax']], eta_jnt, 'tab:blue', 'Johnson noise thermometry (Qu 2017): %.1e' % eta_jnt)
box(R['purdy']['x_range'], R['purdy']['eta95'], 'tab:green', 'optomech. Raman ratio (Purdy 2015): %.1e' % R['purdy']['eta95'])
box([0.05, R['co']['x_range'][1]], 2e-4, 'tab:red', 'CO line ratios vs ab initio (Li 2025): 2e-4')
a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(1e-9, 30); a.set_ylim(1e-14, 1e-2)
a.set_xlabel('x = h nu / kT probed'); a.set_ylabel('|eta| excluded above line (95%)')
a.set_title('hysteresis eta (Stage 31 definition): shaded = excluded')
a.legend(loc='lower right', fontsize=8); a.grid(alpha=0.3, which='both')
a.annotate('lab wins only if eta is\nfrequency-independent', (1e-6, 3e-4), fontsize=8)
a.annotate('x >= 1: FIRAS still best', (0.02, 1.5e-5), fontsize=8)
b = ax[1]
ff = np.linspace(10e3, 368.75e3, 400); xx = h*ff/(kB*T)
for e, col in [(9e-5, 'k'), (1e-11, 'tab:orange'), (1e-12, 'tab:purple'), (eta_jnt, 'tab:blue')]:
    b.plot(ff/1e3, 1/(1 + e/xx) - 1, color=col, label='eta = %.1e' % e)
b.axhspan(-a0_lim, a0_lim, color='tab:blue', alpha=0.2, label='JNT |a0-1| < %.1e (95%%)' % a0_lim)
b.set_yscale('symlog', linthresh=1e-6); b.set_xlabel('frequency (kHz), JNT fit band'); b.set_ylabel('Johnson-noise ratio - 1')
b.set_title('Johnson noise at 273.16 K: eta suppresses low-f noise as 1/(1+eta kT/hf)')
b.text(15, 3e-6, 'bound uses the fitted a0 (quadratic fit absorbs part of the 1/f shape)', fontsize=7)
b.legend(fontsize=8, loc='lower right'); b.grid(alpha=0.3)
fig.suptitle('Stage 47 - hysteresis vs detailed balance: freq-independent eta < %.0e (JNT, %.1e x tighter than FIRAS); at x >= 1 lab gives 2e-4, FIRAS (9e-5) still wins' % (eta_jnt, 9e-5/eta_jnt), fontsize=11)
fig.tight_layout(); fig.savefig(OUT+'figures/stage47_hysteresis.png', dpi=110)
