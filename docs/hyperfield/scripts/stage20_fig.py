#!/usr/bin/env python3
"""Stage 20 figure: (left) required impingement flux per design vs demonstrated/assumed vacuum, in flux units;
(right) equilibrium flux from saturated vapour / adsorbed layers vs temperature."""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import stage20_run as s20  # noqa: E402
d = json.load(open(os.path.join(ROOT, 'results/stage20.json')))
fig, (a, b) = plt.subplots(1, 2, figsize=(15, 6.5), gridspec_kw=dict(width_ratios=[1.2, 1]))
L = d['limits']
cols = {'upper limit (lifetime, 68% CL)': 'tab:green', 'upper limit (lifetime)': 'tab:green', 'upper limit (collision rate)': 'tab:green',
        'upper limit (conservative)': 'tab:green', 'measured (gauge-limited)': 'tab:blue', 'proposed': 'tab:orange',
        'study estimate': 'tab:orange', 'assumed': 'tab:red', 'Stage 12/16 convention': 'tab:gray'}
y = np.arange(len(L) + 1)
for i, l in enumerate(L):
    a.barh(i, np.log10(l['J']), color=cols[l['kind']], alpha=0.8)
    a.text(np.log10(l['J']) + 0.05, i, f"{l['label']}  (n={l['n']:.1e} m$^{{-3}}$, {l['species']}, {l['T']:g} K)", va='center', fontsize=6.5)
a.barh(len(L), np.log10(d['solar_wind']['J']), color='k', alpha=0.5)
a.text(np.log10(d['solar_wind']['J']) + 0.05, len(L), 'solar wind protons at 1 AU, unshielded (n v)', va='center', fontsize=6.5)
for (k, dd), c in zip(d['designs'].items(), ['tab:purple', 'tab:brown', 'navy']):
    a.axvline(np.log10(dd['J_req']), color=c, ls='--', lw=1.5, label=f"required: {k} (J={dd['J_req']:.2g})")
a.set_xlim(7, 16.5); a.set_yticks([]); a.set_xlabel(r'log$_{10}$ impingement flux $J=n\bar v/4$ (m$^{-2}$ s$^{-1}$)')
a.set_title('Vacuum as a species/temperature-independent flux\n(green: demonstrated upper limits; blue: gauge; orange/red: proposed/assumed)', fontsize=9)
a.legend(fontsize=7, loc='lower right')
T = np.geomspace(0.01, 6, 400)
def J(Tv, kind):
    out = []
    for t in Tv:
        if kind == 'He3': p, m = s20.svp(t, 'He3')[0], s20.M['He3']
        elif kind == 'He4': p, m = s20.svp(t, 'He4')[0], s20.M['He4']
        elif kind == 'H2': p, m = s20.svp_H2(t), s20.M['H2']
        else:
            out.append(s20.n_ads_eq(t, kind) * s20.vbar(t, s20.M['He4']) / 4); continue
        out.append(p / (s20.K_B * t) * s20.vbar(t, m) / 4)
    return np.array(out)
for kind, lab, c in [('He3', '$^3$He bulk SVP', 'tab:red'), ('He4', '$^4$He bulk SVP', 'tab:orange'), ('H2', 'solid H$_2$ SVP (Souers, extrap. <4 K)', 'tab:blue'),
                     (1.0, 'He on graphite, 1 monolayer', 'tab:green'), (1e-6, 'He on graphite, 1e-6 monolayer', 'olive')]:
    j = J(T, kind)
    b.loglog(T, np.where(j > 1e-30, j, np.nan), color=c, label=lab)
for (k, dd), c in zip(d['designs'].items(), ['tab:purple', 'tab:brown', 'navy']):
    b.axhline(dd['J_req'], color=c, ls='--', lw=1)
b.axhline(L[0]['J'], color='tab:green', ls=':', lw=1.5, label='BASE 2017 upper limit (He)')
for Tm, lab in [(0.015, '15 mK'), (0.05, '50 mK'), (1.0, '1 K'), (4.0, '4 K')]:
    b.axvline(Tm, color='gray', lw=0.5); b.text(Tm, 1e20, lab, rotation=90, fontsize=7)
b.set_ylim(1e4, 1e32); b.set_xlabel('surface temperature (K)'); b.set_ylabel(r'equilibrium flux $J$ (m$^{-2}$ s$^{-1}$)')
b.set_title('Equilibrium sources: requirement (dashed) is met below\n$^3$He 61 mK, $^4$He 0.17 K, H$_2$ 1.9 K, adsorbed He ~3 K', fontsize=9)
b.legend(fontsize=7, loc='lower right')
plt.tight_layout(); plt.savefig(os.path.join(ROOT, 'figures/stage20_vacuum.png'), dpi=150); print('ok')
