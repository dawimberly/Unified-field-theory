#!/usr/bin/env python3
"""Stage 24 figure: sensitivity of the non-local window and the Stage 22 closing mass to the smearing length and conventions."""
import json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
R = json.load(open(os.path.join(ROOT, 'results/stage24.json')))
fig, ax = plt.subplots(1, 3, figsize=(19, 6.5))
marks = [(1.5e-22, 'cluster\nCompton'), (9.2e-18, 'Na\nCompton'), (2.1e-16, 'nucleon\nCompton'), (1e-15, 'r_N\n(Opp.)'),
         (1e-12, 'R_lam\n(Opp.)'), (0.54e-10, 'Donadi\nDP'), (8.2e-9, 'cluster\nradius'), (1e-7, 'r_C')]
sty = {'d': '-', 'd_half': '--', 'xi_d': ':'}
for k, g in R['grid'].items():
    mode, dm = k.split('|'); c = 'C0' if mode == 'floor-only' else 'C3'
    ax[0].plot(g['l'], g['window'], sty[dm], color=c, label=f'{mode}, dx={dm}')
    mc = np.array(g['closing_MDa'], dtype=float); ax[1].plot(g['l'], mc, sty[dm], color=c, label=f'{mode}, dx={dm}')
for a in ax[:2]:
    a.set_xscale('log')
    for x, t in marks:
        a.axvline(x, color='grey', lw=.6, alpha=.6)
    a.grid(alpha=.3)
ymax = 40
for x, t in marks:
    ax[0].text(x, ymax * 0.93, t, fontsize=6, ha='center', va='top')
ax[0].axhline(0, color='k'); ax[0].axhline(3.66, color='g', lw=.8, ls='--')
ax[0].set_ylim(-25, ymax); ax[0].set_xlabel('smearing length l [m]'); ax[0].set_ylabel('window width [orders] (<0 = already closed)')
ax[0].set_title('(a) non-local window vs smearing length\nblue: floor only; red: same l in floor and ceiling')
ax[0].legend(fontsize=7, loc='lower right')
ax[1].set_yscale('log'); ax[1].axhspan(1e-7, R['thresholds']['M_reach_MDa'], color='C2', alpha=.12, label='reachable in Stage 22 geometry (sag <= 31.7 mm)')
ax[1].axhline(0.172, color='g', ls=':', lw=1, label='Pedalino 172 kDa (already demonstrated)')
ax[1].set_ylim(1e-3, 1e8); ax[1].set_xlabel('smearing length l [m]'); ax[1].set_ylabel('closing mass at cheapest MUSCLE point (xi=0.8) [MDa]')
ax[1].set_title('(b) Stage 22 closing mass vs l\nStage 22 survives only for l <= 1.7e-14 m (floor only) / 1.65e-15 m (consistent)')
ax[1].legend(fontsize=7, loc='upper left')
rows = R['rows']
lab = [r['label'][:70] for r in rows]; w = [r['window_orders'] for r in rows]
col = ['C2' if r.get('reachable_in_stage22_constraints') and r['window_orders'] > 0 else ('C7' if r['window_orders'] <= 0 else 'C3') for r in rows]
y = np.arange(len(rows)); ax[2].barh(y, w, color=col); ax[2].axvline(0, color='k'); ax[2].axvline(3.66, color='g', ls='--', lw=.8)
ax[2].set_yticks(y); ax[2].set_yticklabels(lab, fontsize=5.5); ax[2].invert_yaxis()
ax[2].set_xlabel('window [orders]'); ax[2].set_title('(c) all table rows\ngreen: Stage 22 closes it; red: unreachable; grey: closed', fontsize=10)
ax[2].grid(alpha=.3, axis='x')
fig.suptitle('Stage 24 - how much of the Stage 22 result is convention: smearing length, superposition size, constituent counting', fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, 'figures/stage24_smearing.png'), dpi=130); print('ok')
