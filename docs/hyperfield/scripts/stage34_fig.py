#!/usr/bin/env python3
"""Stage 34 figure (Stage 33 + LIGO 10 Hz-5 kHz curve): bound on the resonant-release rate lambda(f0) (CSL-equivalent, per nucleon, r_C = 1e-7 m) vs frequency."""
import json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
d = json.load(open(os.path.join(ROOT, 'results', 'stage34.json')))
import csv
fig, ax = plt.subplots(figsize=(13, 7))
cols = {'free-fall force noise': 'tab:green', 'position diffusion': 'tab:olive', 'oscillator': 'tab:purple', 'interferometer force noise': 'tab:blue',
        'resonant bar': 'tab:cyan', 'trapped-ion heating': 'tab:brown', 'phonon heating': 'tab:red', 'phonon-mode heating': 'darkred', 'radiation': 'tab:orange'}
for lo, hi, dec in d['coverage']['open_windows']:
    if dec > 1e-3: ax.axvspan(lo, hi, color='0.85', alpha=0.6, lw=0, zorder=0)
ax.axvspan(1, 1, color='0.85', alpha=0.6, label='open window (no verified bound)')
ax.axhspan(1e-10, 1e-6, color='gold', alpha=0.08, zorder=0, label='Adler $10^{-8\\pm2}$ s$^{-1}$ (CITED, ccb18)')
ax.axhline(1e-16, color='k', ls=':', lw=1, label='GRW $10^{-16}$ s$^{-1}$ (CITED)')
seen = set()
# LPF 2018 curve
import importlib.util
for b in d['bands']:
    c = cols[b['kind']]; lab = b['kind'] if b['kind'] not in seen else None; seen.add(b['kind'])
    narrow = b['fmax'] / b['fmin'] < 1.3
    if b['name'].startswith('LISA Pathfinder 2018') or b['stage'] == '34':
        continue
    if not b['mass_prop']:
        ax.hlines(b['lam'], b['fmin'], b['fmax'], colors='0.5', linestyles='--', lw=2, label='Majorana non-mass-proportional (other model)')
        continue
    if narrow:
        f0 = np.sqrt(b['fmin'] * b['fmax'])
        ax.plot([f0, f0], [b['lam'], 3 * b['lam']], color=c, lw=1, alpha=0.6)
        ax.plot(f0, b['lam'], 'v', color=c, ms=8, label=lab, zorder=5)
    else:
        ax.fill_between([b['fmin'], b['fmax']], b['lam'], b['lam'] * 30, color=c, alpha=0.15, lw=0)
        ax.hlines(b['lam'], b['fmin'], b['fmax'], colors=c, lw=2.5, label=lab, zorder=4)
# overlay the LPF 2018 frequency-resolved curve
spec_path = os.path.join(ROOT, 'results', 'stage33_lpf18_curve.csv')
lg = np.loadtxt(os.path.join(ROOT, 'results', 'stage34_ligo_curve.csv'), delimiter=',', skiprows=1)
ax.fill_between(lg[:, 0], lg[:, 2], lg[:, 1], color='navy', alpha=0.35, lw=0)
ax.plot(lg[:, 0], lg[:, 1], color='navy', lw=1.2, label='LIGO O1/O3 force noise, 10 Hz–5 kHz (DERIVED, Stage 34)')
ax.axvspan(2000, 5000, ymin=0.0, ymax=1.0, color='navy', alpha=0.05, hatch='//', lw=0)
if os.path.exists(spec_path):
    arr = np.loadtxt(spec_path, delimiter=',', skiprows=1)
    ax.plot(arr[:, 0], arr[:, 1], color='darkgreen', lw=2, label='LPF 2018, ASD$^2$-scaled (DERIVED)')
ann = [(3e-4, 3e-11, 'LPF'), (0.02, 3e-9, 'cold atoms'), (11, 1e-8, 'levitated'), (300, 3e-1, 'LIGO (34)'), (931, 3e-2, 'AURIGA'), (5e3, 3e-12, 'cantilevers'), (1.5e6, 2e1, 'ions'),
       (6e9, 3e-13, 'bulk heating'), (6e9, 3e-3, 'HBAR'), (1e19, 3e-11, 'X-ray'), (5e20, 3e-14, 'γ (Ge)')]
for x, y, t in ann: ax.text(x, y, t, ha='center', va='bottom', fontsize=8)
ax.text(3e14, 1e-4, 'OPEN: 13 GHz – 14.5 keV\n(THz, IR, optical, UV, soft X)\n8.4 decades', ha='center', fontsize=10, color='0.3')
ax.text(5e7, 1e-4, 'OPEN:\n4.6 MHz–5.6 GHz', ha='center', fontsize=9, color='0.3')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e-6, 1e21); ax.set_ylim(1e-21, 1e6)
ax.set_xlabel('release (resonance) frequency $f_0$ (Hz)')
ax.set_ylabel('max extra release rate per nucleon $\\lambda(f_0)$ (s$^{-1}$, CSL-equiv., $r_C$=100 nm)')
cv = d['coverage']
ax.set_title(f"Stage 34 — Stage-31 gap 20 mHz–19 keV ({cv['gap_decades']:.1f} dec): "
             f"Stage 33 {cv['covered33']:.2f} → Stage 34 {cv['covered34']:.2f} dec (hatched: 2–5 kHz, calibration unquantified)", fontsize=10.5)
ax.grid(alpha=0.3, which='major'); ax.legend(fontsize=7.5, loc='lower left', ncol=2, framealpha=0.95)
fig.text(0.5, 0.005, 'Bounds are on the noise level averaged over each experiment\'s band; a release line narrower than that band is bounded more weakly by (band width)/(line width). '
         'Narrow resonators drawn as ▼ (±0.1 %, Stage 31 convention).', ha='center', fontsize=7.5)
plt.tight_layout(rect=(0, 0.02, 1, 1))
os.makedirs(os.path.join(ROOT, 'figures'), exist_ok=True)
plt.savefig(os.path.join(ROOT, 'figures', 'stage34_frequency_coverage.png'), dpi=150)
print('saved')
