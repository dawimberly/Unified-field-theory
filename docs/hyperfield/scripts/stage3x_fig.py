#!/usr/bin/env python3
"""Generic coverage figure for Stages 35+ : bound on the resonant-release rate lambda(f0) (CSL-equivalent, per nucleon, r_C = 1e-7 m) vs frequency.
usage: stage3x_fig.py <stage>   (reads results/stage<N>.json, writes figures/stage<N>_frequency_coverage.png)"""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
N = sys.argv[1]
d = json.load(open(os.path.join(ROOT, 'results', f'stage{N}.json')))
fig, ax = plt.subplots(figsize=(13, 7))
cols = {'free-fall force noise': 'tab:green', 'position diffusion': 'tab:olive', 'oscillator': 'tab:purple', 'interferometer force noise': 'tab:blue',
        'resonant bar': 'tab:cyan', 'trapped-ion heating': 'tab:brown', 'phonon heating': 'tab:red', 'phonon-mode heating': 'darkred', 'radiation': 'tab:orange',
        'cosmic background': 'tab:pink', 'membrane / resonator heating': 'magenta', 'torsion pendulum': 'teal', 'levitated': 'olive'}
for lo, hi, dec in d['coverage']['open_windows']:
    if dec > 1e-3: ax.axvspan(lo, hi, color='0.85', alpha=0.6, lw=0, zorder=0)
ax.axvspan(1, 1, color='0.85', alpha=0.6, label='open window (no verified bound)')
ax.axhspan(1e-10, 1e-6, color='gold', alpha=0.08, zorder=0, label='Adler $10^{-8\\pm2}$ s$^{-1}$')
ax.axhline(1e-16, color='k', ls=':', lw=1, label='GRW $10^{-16}$ s$^{-1}$')
seen = set()
for b in d['bands']:
    c = cols.get(b['kind'], 'k'); lab = b['kind'] if b['kind'] not in seen else None
    if b['name'].startswith('LISA Pathfinder 2018'):
        continue
    if b.get('curve'):
        arr = np.loadtxt(os.path.join(ROOT, 'results', b['curve']), delimiter=',', skiprows=1)
        ax.fill_between(arr[:, 0], arr[:, 2], arr[:, 1], color=c, alpha=0.35, lw=0)
        ax.plot(arr[:, 0], arr[:, 1], color=c, lw=1.4, label=lab); seen.add(b['kind']); continue
    seen.add(b['kind'])
    if not b['mass_prop']:
        ax.hlines(b['lam'], b['fmin'], b['fmax'], colors='0.5', linestyles='--', lw=2, label='Majorana non-mass-proportional (other model)')
        continue
    if b['fmax'] / b['fmin'] < 1.3:
        f0 = np.sqrt(b['fmin'] * b['fmax'])
        ax.plot([f0, f0], [b['lam'], 3 * b['lam']], color=c, lw=1, alpha=0.6)
        ax.plot(f0, b['lam'], 'v', color=c, ms=8, label=lab, zorder=5)
    else:
        ax.fill_between([b['fmin'], b['fmax']], b['lam'], b['lam'] * 30, color=c, alpha=0.15, lw=0)
        ax.hlines(b['lam'], b['fmin'], b['fmax'], colors=c, lw=2.5, label=lab, zorder=4)
spec_path = os.path.join(ROOT, 'results', 'stage33_lpf18_curve.csv')
if os.path.exists(spec_path):
    arr = np.loadtxt(spec_path, delimiter=',', skiprows=1)
    ax.plot(arr[:, 0], arr[:, 1], color='darkgreen', lw=2, label='LPF 2018, ASD$^2$-scaled')
ax.axvline(12.4e3 * 2.418e14, color='tab:orange', ls='--', lw=1)
ax.text(12.4e3 * 2.418e14, 1e3, ' 12.4 keV: below this, neutral\n atoms emit no CSL photons', fontsize=7.5, color='tab:orange', va='bottom', ha='right')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e-6, 1e21); ax.set_ylim(1e-21, 1e9)
ax.set_xlabel('release (resonance) frequency $f_0$ (Hz)')
ax.set_ylabel('max extra release rate per nucleon $\\lambda(f_0)$ (s$^{-1}$, CSL-equiv., $r_C$=100 nm)')
cv = d['coverage']
ax.set_title(d.get('fig_title', f"Stage {N}"), fontsize=10.5)
ax.grid(alpha=0.3, which='major'); ax.legend(fontsize=7.5, loc='lower left', ncol=2, framealpha=0.95)
fig.text(0.5, 0.005, 'Narrow resonators drawn as ▼. Curves are frequency-resolved. Grey: no verified bound.', ha='center', fontsize=7.5)
plt.tight_layout(rect=(0, 0.02, 1, 1))
os.makedirs(os.path.join(ROOT, 'figures'), exist_ok=True)
plt.savefig(os.path.join(ROOT, 'figures', f'stage{N}_frequency_coverage.png'), dpi=150)
print('saved')
