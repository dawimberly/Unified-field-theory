#!/usr/bin/env python3
"""Stage 16 figure: heavier trapped test masses vs the full budget (reads results/stage16.json)."""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
D = json.load(open(os.path.join(ROOT, 'results/stage16.json')))
B = D['baseline']
sty = {'diamond_spin': ('C0', 'o'), 'lead_meissner': ('C3', 's'), 'niobium_meissner': ('C1', 'D'),
       'silica_trap': ('C2', '^'), 'osmium_space': ('C4', 'v'), 'osmium_space_cold': ('C5', 'P')}


def series(key, field):
    rows = [r for r in D['candidates'][key]['rows'] if r.get('feasible')]
    return np.array([r['m'] for r in rows]), np.array([float(r[field]) for r in rows])


fig, ax = plt.subplots(2, 3, figsize=(17.5, 10.2))
A = ax.ravel()
for key, c in D['candidates'].items():
    col, mk = sty[key]
    lab = c['label']
    m, dx = series(key, 'dx')
    if len(m) == 0:
        continue
    A[0].loglog(m, dx, color=col, marker=mk, ms=3.5, lw=1.4, label=lab)
    A[1].semilogx(m, series(key, 'joint_gap')[1], color=col, marker=mk, ms=3.5, lw=1.4, label=lab)
    A[2].loglog(m, series(key, 'sqrt_Sa_req')[1] * 1e15, color=col, marker=mk, ms=3.5, lw=1.4)
    A[3].loglog(m, series(key, 'heat_budget_W')[1], color=col, marker=mk, ms=3.5, lw=1.4)
    A[4].loglog(m, series(key, 'P_req')[1], color=col, marker=mk, ms=3.5, lw=1.4)
    if c['split'] == 'deloc':
        A[5].loglog(m, series(key, 'deloc_factor')[1], color=col, marker=mk, ms=3.5, lw=1.4)

# (a) dx vs m with SOA
for s in D['soa']:
    A[0].loglog(s['m'], s['dx'], 'k*' if s['com'] else 'kx', ms=11)
    A[0].annotate(s['label'] + ('' if s['com'] else ' (not c.o.m.)'), (s['m'], s['dx']), fontsize=7,
                  xytext=(4, 4), textcoords='offset points')
A[0].loglog(B['m'], B['dx'], marker='*', color='gold', mec='k', ms=15, ls='none', label='Stage 12 cheapest point')
A[0].set_xlabel('test mass  m  [kg]'); A[0].set_ylabel(r'smallest workable superposition  $\Delta x$  [m]')
A[0].set_title('(a) required delocalisation vs mass (best hold per mass)')
A[0].legend(handles=[A[0].lines[-1]], fontsize=7, loc='upper right')
A[0].set_xlim(1e-23, 1e-4)

# (b) joint gap
A[1].axhline(D['baseline']['joint_gap'][0], color='k', ls='--', lw=0.9)
A[1].text(2e-16, D['baseline']['joint_gap'][0] + 0.15, 'Stage 12 cheap point: %.2f' % D['baseline']['joint_gap'][0], fontsize=8)
mm = np.geomspace(1e-16, 1e-5, 50)
A[1].semilogx(mm, np.log10(mm / D['soa'][0]['m']), color='0.5', ls=':', lw=1.2, label='mass-only gap vs 172 kDa (Stage 12 metric)')
A[1].set_xlabel('m [kg]'); A[1].set_ylabel('orders of magnitude')
A[1].set_title(r'(b) joint gap: $\Sigma\,\log_{10}$(req/SOA) over $m,\ \Delta x,\ t$')
A[1].legend(fontsize=6.5, loc='upper left'); A[1].set_ylim(6, 19)

# (c) acceleration noise
A[2].axhline(1.74, color='k', ls='--', lw=0.9)
A[2].text(2e-16, 2.1, 'LISA Pathfinder 1.74 fm s$^{-2}$/$\\sqrt{Hz}$ (space, >2 mHz)', fontsize=7.5)
f = D['cross_checks']
A[2].loglog(D['fuchs']['m'], f['fuchs_thermal_acc_fm'], 'k^', ms=8)
A[2].annotate('Fuchs 2024 thermal (0.43 mg)', (D['fuchs']['m'], f['fuchs_thermal_acc_fm']), fontsize=7, xytext=(-120, -12), textcoords='offset points')
A[2].loglog(D['hofer']['m'], f['hofer_thermal_acc_fm'], 'kv', ms=8)
A[2].annotate('Hofer 2023 thermal (5.6 ug)', (D['hofer']['m'], f['hofer_thermal_acc_fm']), fontsize=7, xytext=(-110, 6), textcoords='offset points')
A[2].loglog(B['m'], B['sqrt_Sa_req'] * 1e15, marker='*', color='gold', mec='k', ms=14, ls='none')
A[2].set_xlabel('m [kg]'); A[2].set_ylabel(r'allowed uncorrelated $\sqrt{S_a}$  [fm s$^{-2}$ Hz$^{-1/2}$]')
A[2].set_title('(c) allowed acceleration noise (10% of budget)')

# (d) heat budget
A[3].axhline(1.8161706e-18, color='k', ls='--', lw=0.9)
A[3].text(2e-16, 2.6e-18, 'Stage 14: 1.82 aW (1 K to 4 K in 3 s)', fontsize=7.5)
A[3].loglog(B['m'], B['heat_budget_to_ceiling_aW'] * 1e-18, marker='*', color='gold', mec='k', ms=14, ls='none')
A[3].set_xlabel('m [kg]'); A[3].set_ylabel('allowed internal heating  [W]')
A[3].set_title(r'(d) internal-heating budget (to witness or T$_c$/2 ceiling)')

# (e) pressure
A[4].axhline(1e-15, color='k', ls='--', lw=0.9)
A[4].text(2e-16, 1.4e-15, 'QGEM-design XHV 1e-15 Pa', fontsize=7.5)
A[4].loglog(B['m'], B['P_req'], marker='*', color='gold', mec='k', ms=14, ls='none')
A[4].set_xlabel('m [kg]'); A[4].set_ylabel('required residual-He pressure [Pa]')
A[4].set_title('(e) vacuum (gas allowed 40% of what is left)')

# (f) expansion factor
A[5].axhline(3.5, color='k', ls='--', lw=0.9)
A[5].text(2e-16, 4.5, 'Rossi 2025: 3.5 x initial coherence length', fontsize=7.5)
for k, v in D['scaling'].items():
    if k in ('lead', 'diamond', 'osmium') and isinstance(v, dict):
        A[5].axhline(v['deloc_factor_large_R'], color='0.6', lw=0.6, ls=':')
A[5].text(1e-15, D['scaling']['lead']['deloc_factor_large_R'] * 0.45, 'dotted: analytic large-R limits (Os, Pb, diamond)', fontsize=7, color='0.4')
A[5].set_xlabel('m [kg]'); A[5].set_ylabel(r'$\Delta x / (2 z_{\rm zpf})$ at 100 Hz')
A[5].set_title('(f) expansion factor needed (no-spin routes)')

for a in A:
    a.grid(alpha=0.25, which='both', lw=0.4)
fig.suptitle('Stage 16: heavier test masses trade the mass gap for delocalisation, vacuum and noise; they do not close it', fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.965))
out = os.path.join(ROOT, 'figures', 'stage16_mass.png')
fig.savefig(out, dpi=140)
print('wrote', out)
