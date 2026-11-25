#!/usr/bin/env python3
"""Stage 32 figure: figures/stage32_split_forecast.png (reads results/stage32.json and results/stage32_curves.csv)."""
import os, json, csv
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
d = json.load(open(os.path.join(ROOT, 'results', 'stage32.json')))
cv = list(csv.DictReader(open(os.path.join(ROOT, 'results', 'stage32_curves.csv'))))
g = lambda k: np.array([float(r[k]) for r in cv])
R = g('R_A'); h = d['headline']; rp = d['reproduce_stage31']; d1 = d['H1_stage31_d1_zs']
fig, ax = plt.subplots(2, 2, figsize=(13, 10))
a = ax[0, 0]
a.fill_between(R, g('H0_nd1_zs'), g('H0_nd0_zs'), color='0.85', label='H0 nondipole-scaling variants (R$^0$, R$^{+1}$)')
a.plot(R, g('tauRc_zs'), 'b-', label='light travel R/c (Eq. 3)')
a.plot(R, g('H0_zs'), 'r-', label='H0: R/c + nondipole (red line, ~1/R ASSUMED)')
a.plot(R, g('H0free_cal_zs'), 'm--', label=f"H0free: all of excess scales with R (a = {d['H0free_cal_a_over_inv_c']:.2f}/c)")
a.plot(R, g('H1_stage31_zs'), 'g-', lw=2, label=f'H1: H0 + d$_1$ = {d1:.0f} zs (Stage 31 residual)')
a.plot(R, g('H1_bound485_zs'), 'g:', label='H1 at Stage 31 95% bound (485 zs)')
Rb = np.array(h['R_A']); eb = np.array(h['slope_err_zs'])
from_curve = lambda key: np.interp(Rb, R, g(key))
a.errorbar(Rb - 0.004, from_curve('H1_stage31_zs'), eb, fmt='go', ms=5, capsize=3, label='forecast stat errors, 4 equal-count R bins (on H1)')
a.errorbar(Rb + 0.004, from_curve('H0free_cal_zs'), eb, fmt='ms', ms=5, capsize=3, label='same errors on H0free')
a.errorbar([0.74], [rp['slope_zs']], [rp['slope_stat_fisher_zs']], fmt='k*', ms=12, capsize=3, label='subset S refit (Zenodo Fig. 3), stat')
a.set_xlabel('internuclear distance R (Å)'); a.set_ylabel('birth-delay slope dτ/dcosβ (zs)'); a.set_title('(a) predicted slope vs R: the curve to compare the joint data with')
a.legend(fontsize=7, loc='upper left'); a.set_ylim(0, 1000)
a = ax[0, 1]
for key, st, lab in [('H0_cosalpha_slope', 'r-', 'H0'), ('H0free_cal_cosalpha_slope', 'm--', 'H0free (flat: R-scaling = R-independent fringe shift)'),
                     ('H1_stage31_cosalpha_slope', 'g-', 'H1 (d$_1$ = %.0f zs): fringe shift ∝ 1/R' % d1), ('H1_bound485_cosalpha_slope', 'g:', 'H1 (485 zs)')]:
    a.plot(R, g(key), st, label=lab)
k = np.interp(Rb, R, g('zs_per_cos'))
a.errorbar(Rb, from_curve('H1_stage31_cosalpha_slope'), eb / k, fmt='go', capsize=3, label='forecast errors')
a.set_xlabel('R (Å)'); a.set_ylabel('fringe-shift slope dcos(α$_0$)/dcosβ'); a.set_title('(b) same, in the raw observable (cos α$_0$ per unit cos β)'); a.legend(fontsize=8)
a = ax[1, 0]
dd = np.linspace(0, 1000, 201)
sB = h['sigma_d1_B_zs']; E2 = d['second_energy']
a.plot(dd, dd / sB, 'k-', lw=2, label=f'(B) shape test, 800 eV, ~267k events: σ = {sB:.0f} zs')
a.plot(dd, dd / E2['800 repeat (x2 events)']['sigma_d1_B_zs'], 'k:', label=f"(B) + repeat at 800 eV (x2 events): σ = {E2['800 repeat (x2 events)']['sigma_d1_B_zs']:.0f}")
a.plot(dd, dd / E2['1600']['sigma_d1_B_zs'], 'b-', label=f"(B) + 1600 eV (same events, ASSUMED): σ = {E2['1600']['sigma_d1_B_zs']:.0f}")
a.plot(dd, dd / E2['400']['sigma_d1_B_zs'], 'c-', label=f"(B) + 400 eV (same events, ASSUMED): σ = {E2['400']['sigma_d1_B_zs']:.0f}")
a.plot(dd, dd / h['sigma_d1_A_zs'], 'r--', label=f"(A) excess over fixed R/c (NOT a scaling test): σ = {h['sigma_d1_A_zs']:.0f}")
a.axhline(2, color='0.5', lw=0.8); a.axhline(5, color='0.5', lw=0.8)
a.axvline(d1, color='g', lw=1); a.text(d1 + 5, 0.4, 'Stage 31\nresidual', color='g', fontsize=8)
a.axvline(485, color='g', ls=':'); a.text(490, 0.4, '95% bound', color='g', fontsize=8)
a.set_xlabel('R-independent extra release term d$_1$ (zs per unit cos β)'); a.set_ylabel('expected significance Z (Asimov)'); a.set_ylim(0, 12)
a.set_title('(c) can the split tell H1 from H0?'); a.legend(fontsize=7, loc='upper right')
a = ax[1, 1]
bs = [b for b in d['binning_scan'] if isinstance(b['nbins'], int)]
a.plot([b['nbins'] for b in bs], [b['sigma_d1_B_zs'] for b in bs], 'ko-', label='equal-count R bins')
for b in d['binning_scan']:
    if not isinstance(b['nbins'], int): a.axhline(b['sigma_d1_B_zs'], ls='--', lw=0.8, color='C1' if 'paper' in b['nbins'] else 'C2', label=b['nbins'])
mc = d['monte_carlo']
a.errorbar([4.15], [mc['H1_stage31']['B_std_zs']], fmt='r^', label=f"full-pipeline MC spread ({mc['n_per_hypothesis']} sims, H1)")
a.errorbar([3.85], [mc['H0']['B_std_zs']], fmt='bv', label='MC spread (H0)')
a.set_xlabel('number of R bins'); a.set_ylabel('σ(d$_1$), test (B) (zs)'); a.set_title('(d) beyond ~4 bins the bin choice matters little'); a.legend(fontsize=8)
fig.suptitle('Stage 32 - forecast for a joint angle × bond-length split of Grundmann et al. 2020 (H$_2$, 800 eV)', fontsize=12)
fig.tight_layout(); out = os.path.join(ROOT, 'figures', 'stage32_split_forecast.png'); fig.savefig(out, dpi=130); print('wrote', out)
