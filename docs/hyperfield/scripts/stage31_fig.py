#!/usr/bin/env python3
"""Stage 31 figure: do-si-do hypothesis bounds (reads results/stage31.json)."""
import json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
d = json.load(open(os.path.join(ROOT, 'results', 'stage31.json')))
fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
# (a) angle-of-incidence birth delay
g = d['grundmann']; w = g['windows']['0.4']
cb = np.array(w['cos_beta']); t = np.array(w['tau_zs']); e = np.array(w['dtau_zs'])
a = ax[0, 0]; x = np.linspace(-1, 1, 50)
a.errorbar(cb, t, e, fmt='s', color='purple', label='Stage 31 reanalysis (Zenodo Fig3 data, stat. only)')
a.plot(x, w['lin']['a'] + w['lin']['b'] * x, '--', color='purple', label=f"fit slope {w['lin']['b']:.0f} zs")
a.plot(x, g['model_slope_R_over_c_zs'] * x, 'b-', label=f"light travel R cos(beta)/c ({g['model_slope_R_over_c_zs']:.0f} zs)")
red = d['sources']['grundmann_fig3c']['values']['red_slope_zs']
a.plot(x, red * x, 'r-', label=f'+ nondipole (paper red line, digitized ~{red} zs)')
a.fill_between(x, (red - g['bound_extra_vs_red_95_zs']) * x, (red + g['bound_extra_vs_red_95_zs']) * x, color='orange', alpha=0.15,
               label=f"room for extra hook term (95%): |d1| < {g['bound_extra_vs_red_95_zs']:.0f} zs")
a.set_xlabel('cos(beta)  (angle of incidence: photon direction vs molecular axis)'); a.set_ylabel('birth-time delay (zs)')
a.set_title('(1)+(6) Angle of incidence vs release time: H$_2$, 800 eV (Grundmann 2020)'); a.legend(fontsize=7.5, loc='upper left'); a.grid(alpha=.3)
# (b) frequency map
a = ax[0, 1]
cols = ['C0', 'C1', 'C2', 'C3', 'C4']
for b, c in zip(d['frequency']['bands'], cols):
    a.plot([b['fmin'], b['fmax']], [b['lam'], b['lam']], lw=6 if b['fmax'] / b['fmin'] > 1.01 else 0, color=c, solid_capstyle='butt', label=f"{b['name']}: {b['lam']:.1e}")
    if b['fmax'] / b['fmin'] <= 1.01: a.plot([b['fmin']], [b['lam']], 'v', ms=12, color=c, zorder=5)
gp = d['frequency']['gap']; a.axvspan(gp['fmin'], gp['fmax'], color='grey', alpha=0.12, label='no bound verified this stage (except 8174 Hz)')
a.axhline(1e-16, color='k', ls=':', lw=1, label='GRW reference 1e-16 s$^{-1}$')
a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(1e-4, 1e20); a.set_ylim(1e-21, 1e-5)
a.set_xlabel('release (resonance) frequency f0 (Hz)'); a.set_ylabel('max extra lock rate per nucleon (s$^{-1}$, r$_C$=100 nm)')
a.set_title('(3) Frequency-selective release: where a resonance could hide'); a.legend(fontsize=7, loc='lower left', framealpha=0.9); a.grid(alpha=.3)
# (c) channel-dependent release time
a = ax[1, 0]; o = d['ossiander']; s = d['sources']['ossiander']['values']
a.errorbar(s['E'], o['diff_as'], s['se'], fmt='o', color='k', label='experiment - ab initio theory (Table 1)')
a.axhline(0, color='grey'); a.axhspan(o['weighted_mean_as'] - 1.96 * o['err_as'], o['weighted_mean_as'] + 1.96 * o['err_as'], color='C2', alpha=.2,
                                      label=f"mean {o['weighted_mean_as']:.2f} +- {o['err_as']:.2f} as (95% band)")
a.set_xlabel('photon energy (eV)'); a.set_ylabel('extra release-time difference (as)')
a.set_title('(5)+(6) Shake-up vs direct photoemission in He (Ossiander 2017)'); a.legend(fontsize=8); a.grid(alpha=.3)
# (d) summary of dimensionless bounds
a = ax[1, 1]
items = [('1 Malus, generic shape', 1e-2), ('1 Malus, extra leak at 90 deg', d['malus']['crossed_extra_leak_max']),
         ('1 lab-orientation anisotropy', d['anisotropy']['frac']), ('1+6 angle-dep. delay / (R/c)', d['grundmann']['bound_extra_vs_red_95_zs'] / d['grundmann']['model_slope_R_over_c_zs']),
         ('2 gate threshold E_th / 1 eV (R1)', d['threshold']['E_quantum_eV']), ('4 amplitude eps / transition', d['amplitude']['eps_over_transition']),
         ('5 intrinsic openness fraction', d['spread']['intrinsic_fraction_max']), ('6 channel delay / 12.6 as', d['ossiander']['bound_abs_extra_channel_delay_as_95'] / 12.6),
         ('7 hysteresis eta', d['hysteresis']['eta_max'])]
y = np.arange(len(items)); v = [np.log10(i[1]) for i in items]
a.barh(y, v, color=['C1' if x > -1 else 'C0' for x in v]); a.set_yticks(y); a.set_yticklabels([i[0] for i in items], fontsize=8)
for yi, vi, it in zip(y, v, items): a.text(vi + (0.3 if vi < 0 else 0.3), yi, f'{it[1]:.1e}', va='center', fontsize=7)
a.axvline(0, color='k'); a.set_xlabel('log10(max allowed extra effect, dimensionless)'); a.set_xlim(-32, 3)
a.set_title('Summary: orange = order-unity room left'); a.grid(alpha=.3, axis='x')
fig.suptitle("Stage 31: Dan's do-si-do (hooked-arm) hypothesis vs verified data", fontsize=13)
fig.tight_layout(); os.makedirs(os.path.join(ROOT, 'figures'), exist_ok=True)
fig.savefig(os.path.join(ROOT, 'figures', 'stage31_dosido.png'), dpi=130)
print('saved figures/stage31_dosido.png')
