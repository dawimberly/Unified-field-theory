#!/usr/bin/env python3
"""Stage 12 figures: budget comparison, Pareto front, mass-time / delocalisation trade."""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
os.chdir(ROOT)
D = json.load(open('results/stage12.json'))

# ---------- fig 1: literature + cheapest, key numbers ----------
fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.1))

# left: phase, decoherence budget, shots for named points + cheapest
names, Phi, gt, C, logN, labels = [], [], [], [], [], []
for p in D['literature']['points']:
    names.append(p['name'].replace('_', '\n'))
    Phi.append(abs(p['Phi']))
    gt.append(p['gamma_t'])
    C.append(p['concurrence'] or 0.0)
    n = p['N_shots_5sig_witness']
    logN.append(np.log10(n) if n and np.isfinite(n) else np.nan)
# cheapest
ch = D['cheapest']
names.append('cheapest\n(this work)')
Phi.append(abs(ch['Phi']))
gt.append(ch['gamma_t'])
C.append(ch['concurrence'])
logN.append(np.log10(ch['N_shots']))

ax = axes[0]
x = np.arange(len(names))
w = 0.36
ax.bar(x - w/2, Phi, w, color='C0', label=r'$|\Phi|=|\phi_1+\phi_2|$ (rad)')
ax.bar(x + w/2, gt, w, color='C3', label=r'decoherence $\gamma t$')
ax.axhline(0.05, color='0.4', ls='--', lw=0.8, label='Schut 5σ phase at 10k shots (0.05 rad)')
ax.set_xticks(x)
ax.set_xticklabels(names, fontsize=7)
ax.set_ylabel('radians  (dimensionless)')
ax.set_title('Entangling phase vs decoherence budget')
ax.set_ylim(0, 0.72)
ax.legend(fontsize=7, loc='upper left')
# concurrence as text
for i, c in enumerate(C):
    ax.text(i, max(Phi[i], gt[i]) + 0.03, f'C={c:.3f}', ha='center', fontsize=7, color='0.2')

# right: orders-of-magnitude gaps at the cheapest point
ax = axes[1]
terms = ch['difficulty_terms']
# display names, skip negative
order = [
    ('mass vs 172 kDa interference', terms['mass']),
    ('mass vs $10^{-18}$ kg cooled NP', 3.75),
    (r'$\Delta x$ vs 133 nm grating', terms['dx']),
    ('duration vs 10 ms flight', terms['time']),
    ('vacuum vs $10^{-10}$ Pa UHV', terms['vacuum_vs_UHV']),
    ('vacuum vs $10^{-15}$ Pa (QGEM design)', terms['vacuum_vs_XHV']),
    (r'$\partial_x B$ vs $10^4$ T/m', terms['gradient']),
    ('Humpty-Dumpty vs $10^{-6}$ rel. precision', terms['humpty_gs']),
    ('shot count vs 10 000', terms['shots']),
]
labs = [o[0] for o in order]
vals = [o[1] if o[1] is not None else 0 for o in order]
colors = ['C3' if v > 0.05 else 'C2' for v in vals]
ax.barh(range(len(labs))[::-1], vals, color=colors)
ax.axvline(0, color='k', lw=0.6)
ax.set_yticks(range(len(labs))[::-1])
ax.set_yticklabels(labs, fontsize=8)
ax.set_xlabel('orders of magnitude harder than SOA  (negative = easier)')
ax.set_title('Gaps at the cheapest workable point'
             '\nparallel + plate + trap; 5.6e-15 kg, 3.9 μm, 2.5 s')
fig.tight_layout()
fig.savefig('figures/stage12_budget.png', dpi=120)
print('wrote figures/stage12_budget.png')

# ---------- fig 2: Pareto in (m, dx) and (m, t) ----------
fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.1))
pts = D['workable_thin']
front = D['pareto']
# unique front by (m,dx,t,d,geom)
seen=set(); front_u=[]
for p in front:
    key=(p['geometry'], p['trapped'], round(p['m'],30), round(p['dx'],20), round(p['t_int'],5), round(p['d'],12))
    if key in seen: continue
    seen.add(key); front_u.append(p)

def colour(p):
    if p['trapped']:
        return 'C0'
    if p.get('shielded') and p['geometry']=='parallel':
        return 'C1'
    if p.get('shielded'):
        return 'C2'
    return '0.5'

ax = axes[0]
for p in pts:
    ax.scatter(p['m'], p['dx'], s=8, c=colour(p), alpha=0.25, linewidths=0)
for p in front_u:
    ax.scatter(p['m'], p['dx'], s=28, c=colour(p), edgecolors='k', linewidths=0.4, zorder=3)
ax.scatter([ch['m']], [ch['dx']], marker='*', s=180, c='C3', zorder=5, label='cheapest')
# literature markers
ax.scatter([1e-14], [250e-6], marker='x', s=60, c='k', zorder=4, label='Bose 2017')
ax.scatter([1e-15], [23e-6], marker='+', s=80, c='k', zorder=4, label='VDK 2020')
ax.scatter([1e-14], [1.7e-6], marker='D', s=30, c='k', zorder=4, label='Schut 2024 Table I (d=21 μm)')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('mass  [kg]'); ax.set_ylabel(r'superposition size $\Delta x$  [m]')
ax.set_title('Workable points and Pareto front')
handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='C0', markersize=8, label='parallel + plate + trap'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='C1', markersize=8, label='parallel + plate, free-fall'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='C2', markersize=8, label='linear + plate, free-fall'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='0.5', markersize=8, label='unshielded'),
    Line2D([0], [0], marker='*', color='C3', markersize=12, linestyle='None', label='cheapest (this work)'),
]
ax.legend(handles=handles, fontsize=7, loc='upper right')

ax = axes[1]
for p in pts:
    ax.scatter(p['m'], p['t_int'] + 2*(p['tau_acc'] if 'tau_acc' in p else 0),
               s=8, c=colour(p), alpha=0.25, linewidths=0)
# t_tot may be missing on some thin records — use t_int+2 tau if present
for p in front_u:
    ttot = p.get('t_tot', p['t_int'] + 0.5)
    ax.scatter(p['m'], ttot, s=28, c=colour(p), edgecolors='k', linewidths=0.4, zorder=3)
ax.scatter([ch['m']], [ch['t_tot']], marker='*', s=180, c='C3', zorder=5)
# mass-time trade curves at Phi=0.05
trade = D['mass_time_trade']
from collections import defaultdict
g=defaultdict(list)
for r in trade:
    g[r['label']].append((r['m'], r['t_for_Phi_0p05']))
styles = {'linear_noshield_Bose_like': ('k', '--', 'Bose-like, Φ=0.05 (gravity only)'),
          'linear_shield_VDK_like': ('C2', '--', 'VDK-like, Φ=0.05 (gravity only)'),
          'parallel_trap_d21': ('C0', '-', 'trap d=21 μm, Φ=0.05 (gravity only)'),
          'parallel_trap_d61': ('C0', ':', 'trap d=61 μm, Φ=0.05 (gravity only)')}
for lab, pts_ in g.items():
    pts_ = sorted(pts_)
    c, ls, name = styles.get(lab, ('0.5', '-', lab))
    ax.plot([p[0] for p in pts_], [p[1] for p in pts_], color=c, ls=ls, lw=1.2, label=name)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('mass  [kg]'); ax.set_ylabel('time  [s]')
ax.set_ylim(1e-3, 1e4)
ax.set_title('Mass–time trade at fixed phase, and workable points')
ax.legend(fontsize=6.5, loc='upper right')
fig.tight_layout()
fig.savefig('figures/stage12_pareto.png', dpi=120)
print('wrote figures/stage12_pareto.png')

# ---------- fig 3: Pedernales / smaller-dx, Phi contours ----------
fig, ax = plt.subplots(figsize=(6.4, 5.0))
grid = D['pedernales']['grid']
# build arrays
from collections import defaultdict
byd=defaultdict(list)
for r in grid:
    byd[r['d']].append((r['dx'], abs(r['Phi'])))
for d, rows in sorted(byd.items()):
    rows=sorted(rows)
    ax.loglog([r[0] for r in rows], [r[1] for r in rows],
              label=f'd = {d*1e6:.0f} μm', lw=1.8)
ax.axhline(0.05, color='0.4', ls='--', lw=0.9, label='0.05 rad (10k-shot phase 5σ)')
ax.axhline(0.3, color='0.4', ls=':', lw=0.9, label='0.3 rad (Bose-like, comfortable witness)')
# mark cheapest
ax.scatter([ch['dx']], [abs(ch['Phi']) * (1.0/ch['t_int'])],  # no, Phi is already for t=2.5
           marker='*', s=160, c='C3', zorder=5)
# Phi at cheapest is for t=2.5, contours in pedernales grid are for t=1, m=1e-14
# scale cheapest to t=1 for visual comparison: Phi/t * 1
ax.scatter([ch['dx']], [abs(ch['Phi'])/ch['t_int']], marker='*', s=160, c='C3',
           zorder=5, label=r'cheapest (rescaled to 1 s, $\Phi/t$)')
ax.set_xlabel(r'superposition size $\Delta x$  [m]')
ax.set_ylabel(r'$|\Phi|$ in 1 s at $m=10^{-14}$ kg  [rad]')
ax.set_title('Smaller cats only help if $d$ shrinks too\n(Pedernales–Plenio: rate is set by delocalisation)')
ax.legend(fontsize=7, loc='upper left')
fig.tight_layout()
fig.savefig('figures/stage12_trade.png', dpi=120)
print('wrote figures/stage12_trade.png')
