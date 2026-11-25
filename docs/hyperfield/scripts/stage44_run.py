#!/usr/bin/env python3
"""Stage 44 - do-si-do release-time spread and gate threshold from long-baseline decay-law tests and quantum-jump statistics.

MODEL (ASSUMED; QM values in brackets). A 'hook' (unstable nucleus or shelved ion) releases after a waiting time t.
 D1 shape:   t ~ Weibull(tau, beta): hazard h(t) = (beta/tau)(t/tau)^(beta-1); spread sigma_tau/tau = CV(beta)   [beta = 1, CV = 1]
 D2 static spread: each hook has its own fixed rate G = Gbar(1 + s*eps), Var(eps) = 1 (gamma-distributed);
             population hazard Gbar/(1 + s^2 Gbar t)  ->  apparent half-life T(t) = T0 (1 + s^2 Gbar t)            [s = 0]
 D3 slow openness (one system, memory): rate spread s, constant over >= 2 consecutive jumps, re-drawn within < 1000 jumps
             -> consecutive intervals correlated                                                                  [s = 0]
 D4 gate:    no release possible for t < t_g after the hook is formed (dead time), exponential afterwards        [t_g = 0]
Data: Norman et al. 1988 (LBL-24562 preprint of PRL 60, 2246), Berkeland et al. 2004 (arXiv:physics/0304013).
"""
import json, os, csv
import numpy as np
from math import log, sqrt, gamma as G
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
S = {}
def src(key, ref, url, loc, quote, **vals): S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)
src('norman88', 'Norman, Gazes, Crane & Bennett, "Tests of the Exponential Decay Law at Short and Long Times", LBL-24562 (Dec 1987) = PRL 60, 2246 (1988)',
    'https://escholarship.org/content/qt5568p03n/qt5568p03n.pdf', 'Table II; text p.4-6 (OCR of scanned preprint)',
    'Table II (1811 keV, half-life free): sample 1 "0.3 t1/2 ~ 13 t1/2" "2.583 +- 0.008"; sample 2 "11 t1/2 ~ 22 t1/2" "2.585 +- 0.005"; sample 3 "20 t1/2 ~ 33 t1/2" "2.576 +- 0.005"; '
    'sample 4 "34 t1/2 ~ 45 t1/2" "2.573 +- 0.019" [hours] ; "In all four cases the half-life determined from our data agreed with the known value, thus indicating no dependence on the age of the source" ; '
    '"56Mn ... 2.5785 hours" ; 60Co: "irradiated ... in a burst with a full width at half maximum of 12.5 milliseconds" ; "After allowing this activity to decay for 3.5 hours, the sample contained 15 micro-Curies of 60Co" ; '
    '"Counting was begun within 8x10^-5 t1/2(60Co) after the activation" (exponent OCR-garbled; consistent with 3.5 h / 5.271 y = 7.6e-5) ; "60Co ... 5.271 years" ; "Half of the 60Co activity that we observed was produced initially as the t1/2 = 10.5 minute 60Com"',
    T=[2.583, 2.585, 2.576, 2.573], dT=[0.008, 0.005, 0.005, 0.019], lo=[0.3, 11, 20, 34], hi=[13, 22, 33, 45], Tknown=2.5785, co_wait_h=3.5, co_T_y=5.271)
src('berkeland04', 'Berkeland, Raymondson & Tassin, "Tests for nonrandomness in quantum jumps", PRA 69, 052103 (2004), arXiv:physics/0304013v2', 'https://arxiv.org/pdf/physics/0304013v2',
    'abstract; p.2-3, Eqs. (1)-(5)', '"we have observed 228 000 quantum jumps of a single trapped and laser cooled 88Sr+ ion" ; "We measure the fractional reduction in the entropy of information to be < 6.5x10^-4 when the value of any interval between quantum jumps is known" ; '
    '"divide our data sets into continuous strings of 1000 events and normalize each interval time by the average interval time of its string" ; "choose tm(n) at intervals that result in roughly uniform values of Nm,n" ; '
    '"We set the number of bins in each row and column (mmax = nmax) to be 13" ; "U(x|y) = (H(x) - H(x|y))/H(x)" ; "intervals are uncorrelated when they are separated by 1 to 20 intervals, which corresponds to approximately 2 ms to several seconds"',
    U=6.5e-4, nb=13)
src('stage31', 'HYPERFIELD_RESULTS.md Stage 31 (this repo)', 'docs/hyperfield/HYPERFIELD_RESULTS.md', '31.2 (5),(6)', 'openness spread: VIP-2 beta^2/2 <= 8.6e-31 -> intrinsic distinguishable label <~ 1e-30; release time: jitter < 2.7 ps, universal tau unbounded')
out = {}
# ---------- D2 / D1 from Norman Table II ----------
n = S['norman88']['values']; T = np.array(n['T']); dT = np.array(n['dT']); lo = np.array(n['lo']); hi = np.array(n['hi'])
L2 = log(2)
def wls(x, y, w):
    W = np.sum(w); xm = np.sum(w * x) / W; ym = np.sum(w * y) / W
    b = np.sum(w * (x - xm) * (y - ym)) / np.sum(w * (x - xm) ** 2); a = ym - b * xm
    sb = 1 / sqrt(np.sum(w * (x - xm) ** 2)); chi2 = float(np.sum(w * (y - a - b * x) ** 2))
    return a, b, sb, chi2, xm
d2 = {}; d1 = {}
for tag, age in (('midpoint', (lo + hi) / 2), ('start', lo)):          # ASSUMED effective age of each fit
    t_life = age * L2                                                     # age in mean lives
    a, b, sb, chi2, _ = wls(t_life, T, 1 / dT ** 2)
    s2 = b / a; ss2 = sb / a; s2_up = max(s2, 0) + 1.645 * ss2           # one-sided 95% (s^2 >= 0)
    d2[tag] = dict(T0=a, slope_h_per_life=b, s2=s2, s2_err=ss2, s2_95=s2_up, s_95=sqrt(s2_up), chi2=chi2, dof=2)
    x = np.log(t_life); y = np.log(T); w = (T / dT) ** 2
    c, m, sm, chi2w, _ = wls(x, y, w)                                     # ln T = c + (1-beta) ln t
    beta = 1 - m; d1[tag] = dict(beta=beta, beta_err=sm, beta_lo95=beta - 1.96 * sm, beta_hi95=beta + 1.96 * sm, chi2=chi2w)
def cv_weib(b): return sqrt(G(1 + 2 / b) / G(1 + 1 / b) ** 2 - 1)
for tag in d1:
    d1[tag]['CV_range95'] = [cv_weib(d1[tag]['beta_hi95']), cv_weib(d1[tag]['beta_lo95'])]
out['D2_static_spread'] = d2; out['D1_weibull_shape'] = d1
# ---------- D3 from Berkeland mutual information ----------
b = S['berkeland04']['values']; nb = b['nb']; Imax = b['U'] * log(nb)     # all measured U attributed (it is ~ finite-sample bias)
edges = -np.log(1 - np.arange(nb + 1) / nb); edges[-1] = np.inf          # equiprobable bins for unit-mean exponential
gh_x, gh_w = np.polynomial.hermite_e.hermegauss(60); gh_w = gh_w / gh_w.sum()
def MI(s):                                                               # lognormal rate, mean-normalised (interval rescaling, as in the paper)
    sig = sqrt(log(1 + s * s)); r = np.exp(sig * gh_x - sig ** 2 / 2)
    P = np.exp(-np.outer(r, edges[:-1])) - np.exp(-np.outer(r, edges[1:]))   # P[rate, bin]
    pm = gh_w @ P; pmn = np.einsum('k,km,kn->mn', gh_w, P, P)
    return float(np.sum(pmn * np.log(pmn / np.outer(pm, pm))))
ss = np.logspace(-2, 0.3, 200); mis = np.array([MI(s) for s in ss])
s_lim = float(np.interp(Imax, mis, ss))
out['D3_slow_openness'] = dict(I_max_nats=Imax, s_95=s_lim, note='U < 6.5e-4 is itself about the finite-sample bias (nb-1)^2/(2N ln nb) for N ~ 4e4 pairs; all attributed. Rates varying slower than ~1000 jumps are removed by the paper\'s normalisation, faster than one interval average out.')
# ---------- D4 gate from Norman 60Co ----------
tau_co = n['co_T_y'] * 365.25 * 86400 / L2; tg = n['co_wait_h'] * 3600
out['D4_gate'] = dict(t_g_max_s=tg, t_g_over_tau=tg / tau_co, logic='all 60Co nuclei are <= 3.5 h old when counting starts (12.5 ms burst; 60mCo feeding finished after 20 of its half-lives); counts with flat exponential time spectra were recorded, so a dead time t_g > ~3.5 h is excluded')
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
for tag, mk in (('midpoint', 'o'),):
    age = (lo + hi) / 2 * L2
    ax[0].errorbar(age, T, yerr=dT, xerr=[(age - lo * L2), (hi * L2 - age)], fmt=mk, color='k', capsize=3, label='Norman 1988, Table II (fit per sample)')
tt = np.linspace(0, 32, 50); dd = d2['midpoint']
ax[0].plot(tt, dd['T0'] + dd['slope_h_per_life'] * tt, 'b-', label=f"linear fit: s² = {dd['s2']:.1e} ± {dd['s2_err']:.1e}")
ax[0].plot(tt, dd['T0'] * (1 + dd['s2_95'] * tt), 'r--', label=f"D2 95% limit s < {dd['s_95']:.3f}")
ax[0].axhline(n['Tknown'], color='g', lw=0.8, ls=':', label='known 2.5785 h')
ax[0].set_xlabel('age of the 56Mn population (mean lives)'); ax[0].set_ylabel('fitted half-life (h)'); ax[0].legend(fontsize=7); ax[0].set_title('D1/D2: does the decay rate change with age?', fontsize=9)
ax[1].loglog(ss, mis, 'b-', label='D3 model: mutual information of consecutive intervals'); ax[1].axhline(Imax, color='r', ls='--', label=f'Berkeland 2004: U<6.5e-4 → I<{Imax:.2e} nats')
ax[1].axvline(s_lim, color='r', lw=0.8); ax[1].set_xlabel('relative rate spread s (persisting ≥ 2 jumps)'); ax[1].set_ylabel('I (nats), 13×13 equiprobable bins'); ax[1].legend(fontsize=7)
ax[1].set_title(f'D3: memory in quantum jumps → s < {s_lim:.2f}', fontsize=9)
names = ['D2 static rate spread s\n(56Mn, 45 t½)', 'D1 Weibull |β−1|\n(56Mn)', 'D3 slow openness s\n(Sr⁺ jumps)', 'D4 gate t_g/τ\n(60Co)']
vals = [dd['s_95'], max(abs(d1['midpoint']['beta_lo95'] - 1), abs(d1['midpoint']['beta_hi95'] - 1)), s_lim, tg / tau_co]
ax[2].barh(names, vals, color=['tab:blue', 'tab:cyan', 'tab:purple', 'tab:orange']); ax[2].set_xscale('log')
for i, v in enumerate(vals): ax[2].text(v, i, f' {v:.1e}', va='center', fontsize=8)
ax[2].set_xlim(1e-6, 10); ax[2].set_title('95% upper limits (QM: all zero)', fontsize=9)
fig.suptitle('Stage 44 — do-si-do release-time spread & gate: no deviation from exponential/memoryless release', fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, 'figures', 'stage44_release_time.png'), dpi=110)
with open(os.path.join(ROOT, 'results', 'stage44_sources.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for kk, v in S.items(): wr.writerow([kk, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
json.dump(dict(stage=44, results=out, sources=S), open(os.path.join(ROOT, 'results', 'stage44.json'), 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
