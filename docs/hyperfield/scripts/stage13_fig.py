#!/usr/bin/env python3
"""Stage 13 figure: blackbody budget vs measured diamond Im χ."""
import json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
D = json.load(open(os.path.join(ROOT, 'results/stage13.json')))

# colours / styles
style = {
    'stage12_assumed':        dict(c='0.35', ls='--', lw=1.6, z=3),
    'silica_romero':          dict(c='0.55', ls=':',  lw=1.2, z=2),
    'bulk_best_cvd':          dict(c='C0',  ls='-',  lw=2.0, z=4),
    'bulk_typical_window':    dict(c='C0',  ls='-.', lw=1.6, z=4),
    'bulk_scd_subkelvin':     dict(c='C0',  ls=':',  lw=1.4, z=3),
    'bulk_millar_thz_300K':   dict(c='C2',  ls='-',  lw=1.5, z=3),
    'bulk_millar_22ghz_cryo': dict(c='C1',  ls='-',  lw=1.8, z=4),
    'nanodiamond_cleaner':    dict(c='C3',  ls='--', lw=1.8, z=4),
    'nanodiamond_surface':    dict(c='C3',  ls='-',  lw=2.0, z=5),
}
short = {
    'stage12_assumed': 'Stage 12 assumption, Im χ=0.01',
    'silica_romero': 'silica (Romero-Isart), Im χ=0.1',
    'bulk_best_cvd': r'best window CVD, tanδ=$4.5\times10^{-6}$',
    'bulk_typical_window': r'typical window CVD, tanδ=$10^{-5}$',
    'bulk_scd_subkelvin': r'single-crystal, tanδ=$2\times10^{-6}$',
    'bulk_millar_thz_300K': r'Millar 1 THz, tanδ=$6.3\times10^{-5}$',
    'bulk_millar_22ghz_cryo': r'Millar 22 GHz / 35–100 K, tanδ=$9.5\times10^{-4}$',
    'nanodiamond_cleaner': r'cleaner nanodiamond, tanδ=$10^{-3}$ (assumed)',
    'nanodiamond_surface': r'dirty nanodiamond / sp² shell, tanδ=$10^{-2}$',
}

fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.55))

# --- panel 0: C vs T_env = T_int
ax = axes[0]
for cid, st in style.items():
    cur = D['curves'][cid]
    ax.plot(cur['T'], cur['C_eq'], color=st['c'], ls=st['ls'], lw=st['lw'],
            zorder=st['z'], label=short[cid])
ax.axhline(0.1, color='k', lw=0.7, ls='--')
ax.axhline(0.0, color='k', lw=0.5)
ax.text(0.28, 0.108, 'C = 0.1', fontsize=7, color='k')
ax.set_xlim(0.5, 22)
ax.set_ylim(-0.01, 0.16)
ax.set_xlabel(r'$T_{\mathrm{env}}=T_{\mathrm{int}}$  [K]')
ax.set_ylabel('concurrence $C$')
ax.set_title('Equal temperatures\n(air + GGN fixed at the cheap-point budget)')
ax.legend(fontsize=5.6, loc='upper right', framealpha=0.92)

# --- panel 1: C vs T_int at T_env = 1 K
ax = axes[1]
for cid, st in style.items():
    cur = D['curves'][cid]
    ax.plot(cur['T'], cur['C_int_Tenv1'], color=st['c'], ls=st['ls'], lw=st['lw'],
            zorder=st['z'])
ax.axhline(0.1, color='k', lw=0.7, ls='--')
ax.axhline(0.0, color='k', lw=0.5)
ax.set_xlim(0.5, 30)
ax.set_ylim(-0.01, 0.16)
ax.set_xlabel(r'$T_{\mathrm{int}}$  [K]   ($T_{\mathrm{env}}=1$ K)')
ax.set_ylabel('concurrence $C$')
ax.set_title('Internal heating\n(environment held at 1 K)')

# --- panel 2: T_max vs Im χ
ax = axes[2]
sw = D['im_chi_sweep']
ims = [s['im_chi'] for s in sw]
ax.plot(ims, [s['Tmax_eq_Cgt0'] for s in sw], color='C0', lw=1.8, label=r'$T_{\mathrm{env}}=T_{\mathrm{int}}$, $C>0$')
ax.plot(ims, [s['Tmax_eq_Cgt0p1'] for s in sw], color='C0', lw=1.4, ls='--', label=r'$T_{\mathrm{env}}=T_{\mathrm{int}}$, $C>0.1$')
ax.plot(ims, [s['Tmax_int_Cgt0'] for s in sw], color='C3', lw=1.8, label=r'$T_{\mathrm{int}}$ at $T_{\mathrm{env}}=1$ K, $C>0$')
ax.plot(ims, [s['Tmax_int_Cgt0p1'] for s in sw], color='C3', lw=1.4, ls='--', label=r'$T_{\mathrm{int}}$ at $T_{\mathrm{env}}=1$ K, $C>0.1$')
# scatter-only ceiling
sc = D['scatter_only_ceiling']
ax.axhline(sc['Tmax_eq_Cgt0'], color='C0', lw=0.7, ls=':')
ax.text(1.2e-8, sc['Tmax_eq_Cgt0']+0.25, r'scattering ceiling (Im $\chi=0$)', fontsize=6.5, color='C0')

# vertical markers for cases
markers = [
    ('bulk_typical_window', 'C0'),
    ('bulk_millar_22ghz_cryo', 'C1'),
    ('nanodiamond_surface', 'C3'),
    ('stage12_assumed', '0.35'),
]
for cid, col in markers:
    rec = next(x for x in D['cases'] if x['id']==cid)
    ax.axvline(rec['im_chi'], color=col, lw=0.7, ls=':', alpha=0.8)

ax.set_xscale('log')
ax.set_xlim(1e-8, 1.2e-1)
ax.set_ylim(0, 28)
ax.set_xlabel(r'Im$[(\varepsilon-1)/(\varepsilon+2)]$')
ax.set_ylabel(r'maximum temperature  [K]')
ax.set_title('Temperature ceiling versus loss')
ax.legend(fontsize=6, loc='upper right')

fig.suptitle('Stage 13: blackbody budget at the cheapest QGEM point  '
             r'($m=5.6\times10^{-15}$ kg, $\Delta x=3.9$ μm, 3 s)',
             fontsize=10, y=1.02)
fig.tight_layout()
out = os.path.join(ROOT, 'figures/stage13_bbr.png')
fig.savefig(out, dpi=130, bbox_inches='tight')
print('wrote', out)
