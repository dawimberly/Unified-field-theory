#!/usr/bin/env python3
"""Stage 14 figures: energy budget, cooldown, gas cooling, absorbed-power ceiling."""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
D = json.load(open(os.path.join(ROOT, 'results/stage14.json')))

# ---------- fig 1: four-panel budget ----------
fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.6))

# (a) T vs deposited energy
ax = axes[0, 0]
Q = np.array(D['T_vs_energy']['Q_J'])
ax.loglog(Q * 1e18, D['T_vs_energy']['T_from_1K'], 'C0', lw=2.0, label='starting at 1 K')
ax.loglog(Q * 1e18, D['T_vs_energy']['T_from_4K'], 'C1', lw=2.0, label='starting at 4 K')
ax.axhline(4, color='k', ls='--', lw=0.7)
ax.axvline(D['energy_budget']['Q_1K_to_4K_aJ'], color='C0', ls=':', lw=0.8)
ax.axhline(12, color='0.5', ls=':', lw=0.7)
ax.text(8, 4.4, r'4 K witness-safe', fontsize=8)
ax.text(D['energy_budget']['Q_1K_to_4K_aJ']*1.15, 1.3, r'5.45 aJ', fontsize=8, color='C0')
ax.set_xlim(1e-3, 1e9)
ax.set_ylim(0.8, 400)
ax.set_xlabel('deposited energy  [aJ]')
ax.set_ylabel(r'$T_{\mathrm{int}}$  [K]')
ax.set_title(r'(a)  Debye $T^3$ heat capacity: 5.45 aJ takes 1 K $\to$ 4 K')
ax.legend(fontsize=8, loc='upper left')

# (b) radiative cooldown to a 4 K bath
ax = axes[0, 1]
col = {'window': 'C0', 'millar22': 'C1', 'dirty': 'C3', 'stage12': '0.45'}
lab = {'window': r'window-grade, $\tan\delta=10^{-5}$',
       'millar22': r'cleaner ND / Millar 22 GHz, $10^{-3}$',
       'dirty': r'dirty surface, $10^{-2}$',
       'stage12': r'Stage 12 Im $\chi=0.01$'}
ls = {'window': '-', 'millar22': '-', 'dirty': '-', 'stage12': '--'}
for k in ('window', 'millar22', 'dirty', 'stage12'):
    c = D['cooldown_curves'][k]
    ax.loglog(c['t_s'], c['T'], color=col[k], ls=ls[k], lw=1.9, label=lab[k])
ax.axhline(4, color='k', ls='--', lw=0.7)
ax.axvline(3, color='0.5', ls=':', lw=0.8)
ax.text(4, 6, '3 s run', fontsize=7, color='0.3')
ax.set_xlim(1e-1, 1e8)
ax.set_ylim(3.5, 350)
ax.set_xlabel('time  [s]')
ax.set_ylabel(r'$T_{\mathrm{int}}$  [K]')
ax.set_title(r'(b)  Radiative cooldown into a 4 K bath  (lasers off)')
ax.legend(fontsize=7, loc='upper right')

# (c) gas-cooling timescale vs pressure
ax = axes[1, 0]
pts = D['gas']['points']
P = np.array([p['P_Pa'] for p in pts])
ax.loglog(P, [p['tau_He_4K_s'] for p in pts], 'C0', lw=2.0, label=r'He, 4 K, $\alpha_E=0.5$')
ax.loglog(P, [p['tau_He_4K_alpha01'] for p in pts], 'C0', ls='--', lw=1.4, label=r'He, 4 K, $\alpha_E=0.1$')
ax.loglog(P, [p['tau_N2_4K_s'] for p in pts], 'C2', lw=1.6, label=r'N$_2$, 4 K, $\alpha_E=0.5$')
ax.loglog(P, [p['tau_He_300K_s'] for p in pts], 'C1', lw=1.6, label=r'He, 300 K, $\alpha_E=0.5$')
ax.axhline(3, color='k', ls='--', lw=0.7)
ax.axvline(D['gas']['p_He_4K_tau_3s_Pa'], color='C0', ls=':', lw=0.8)
ax.axvline(1e-15, color='0.4', ls=':', lw=0.8)
ax.text(2e-15, 2e3, r'$10^{-15}$ Pa', fontsize=7, color='0.3', rotation=90, va='bottom')
ax.set_xlim(1e-16, 1e4)
ax.set_ylim(1e-6, 1e10)
ax.set_xlabel('gas pressure  [Pa]')
ax.set_ylabel('thermalisation time $C/G$  [s]')
ax.set_title(r'(c)  Buffer-gas cooling  (Knudsen formula, monatomic $2k$)')
ax.legend(fontsize=7, loc='upper right')

# (d) T_ss vs absorbed extra power, T_env = 1 K
# invert P_em(T) - P_em(1K) = P_extra
ax = axes[1, 1]
Ts = np.linspace(1.05, 40, 300)
from numpy import array
# reuse radiative P_em ∝ T^5: P(T) = P4 * (T/4)^5
for k, color, label in (
    ('window', 'C0', r'window-grade, $T_{\mathrm{env}}=1$ K'),
    ('dirty', 'C3', r'dirty surface, $T_{\mathrm{env}}=1$ K'),
):
    P4 = D['radiative'][k]['P_em_4K']
    P1 = D['radiative'][k]['P_em_1K']
    Pextra = P4 * (Ts / 4.0) ** 5 - P1
    ax.loglog(Pextra, Ts, color=color, lw=2.0, label=label)
# impulse-budget equivalent power over 3 s, starting at 1 K
# (ignore radiation during the run for the clean case; it is a tiny correction)
ax.axvline(D['energy_budget']['P_budget_impulse_1to4_W'], color='k', ls='--', lw=0.8)
ax.text(D['energy_budget']['P_budget_impulse_1to4_W']*1.4, 1.4,
        r'3 s impulse budget' '\n' r'(1 K$\to$4 K, no cooling)', fontsize=7)
ax.axhline(4, color='k', ls=':', lw=0.7)
ax.set_xlim(1e-26, 1e-12)
ax.set_ylim(1.05, 40)
ax.set_xlabel('steady absorbed power  [W]')
ax.set_ylabel(r'steady $T_{\mathrm{int}}$  [K]')
ax.set_title(r'(d)  Steady $T_{\mathrm{int}}$ versus leftover absorbed power')
ax.legend(fontsize=7, loc='upper left')

fig.suptitle(r'Stage 14: keeping $T_{\mathrm{int}}\leq 4$ K for the $0.73\,\mu$m, $5.6\times10^{-15}$ kg diamond',
             fontsize=11, y=1.01)
fig.tight_layout()
out1 = os.path.join(ROOT, 'figures/stage14_thermal.png')
fig.savefig(out1, dpi=130, bbox_inches='tight')
print('wrote', out1)
