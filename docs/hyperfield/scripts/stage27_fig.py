#!/usr/bin/env python3
"""Stage 27 figure: signature (correlation vs separation, scaling), tuning gains, discrete-kernel source modulation."""
import json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
r = json.load(open(os.path.join(ROOT, 'results/stage27.json')))
fig, ax = plt.subplots(2, 2, figsize=(13, 10))
# (1) correlation coefficient vs separation
a = ax[0, 0]
s = r['spheres']; p = r['points']
a.loglog(s['d_over_R'], np.abs(s['rho_long']), 'C0-', label='spheres R >> l: |rho| longitudinal (anti-corr.) vs d/R')
a.loglog(s['d_over_R'], np.abs(s['rho_trans']), 'C0--', label='spheres: rho transverse (+) vs d/R')
pl = np.array(p['rho_long']); dl = np.array(p['d_over_l'])
a.loglog(dl[pl > 0], pl[pl > 0], 'C3-', label='point-like R << l: rho_long > 0 (common mode) vs d/l')
a.loglog(dl[pl < 0], -pl[pl < 0], 'C3:', label='point-like: |rho_long| < 0 branch vs d/l')
for g in r['geometries'][:4]:
    R = (3 * g['V'] / (4 * np.pi)) ** (1 / 3); x = g['d'] / R
    a.plot(x, abs(g['rho_long']), 'ko', ms=4); a.annotate(g['geometry'].split(',')[0][:28], (x, abs(g['rho_long'])), fontsize=7)
a.set_xlabel('separation d / R  (or d / l for point-like bodies)'); a.set_ylabel('|correlation coefficient|')
a.set_ylim(1e-6, 1.5); a.set_title('Signature 1: two-body correlation (independent of kappa)\ncross = -2Għ/d³ along d (Kafri-Taylor-Milburn)'); a.legend(fontsize=7)
# (2) scaling exponents
a = ax[0, 1]
rows = r['signature']['scaling_table']; keys = ['R', 'rho', 'T', 'P', 'f']
M = np.array([[v if isinstance(v, (int, float)) else np.nan for v in (row[k] for k in keys)] for row in rows], float)
im = a.imshow(M, cmap='coolwarm', vmin=-5, vmax=5, aspect='auto')
for i in range(len(rows)):
    for j, k in enumerate(keys):
        v = rows[i][k]; a.text(j, i, str(v) if not isinstance(v, str) else v[:9], ha='center', va='center', fontsize=7)
a.set_xticks(range(len(keys))); a.set_xticklabels(['radius R', 'density', 'temperature', 'pressure', 'frequency'])
a.set_yticks(range(len(rows))); a.set_yticklabels([row['source'][:34] for row in rows], fontsize=7)
a.set_title('Signature 2: exponent of acceleration PSD vs knob\n(CQ: only volume matters; flat in T, P, f, composition)')
# (3) tuning
a = ax[1, 0]
sm = r['summary']; names = list(sm)
b0 = [sm[n]['stage26'] for n in names]; b1 = [sm[n]['best'] for n in names]
x = np.arange(len(names))
a.bar(x - 0.2, b0, 0.4, label='Stage 26 shortfall (log10 u)'); a.bar(x + 0.2, b1, 0.4, label='best Stage 27 tuning')
for i in range(len(names)): a.text(x[i] + 0.2, b1[i] + 0.3, '%.1f' % b1[i], ha='center', fontsize=8); a.text(x[i] - 0.2, b0[i] + 0.3, '%.1f' % b0[i], ha='center', fontsize=8)
a.axhline(0, color='k'); a.set_xticks(x); a.set_xticklabels([n + '\n' + sm[n]['route'] for n in names], fontsize=8)
a.set_ylabel('orders short (log10 u); < 0 = reaches'); a.set_title('Tuning: gains are 1-2.3 orders; only MAQRO crosses zero'); a.legend(fontsize=8)
# (4) discrete kernel window
a = ax[1, 1]
d = r['tuning']['b_modulation']['discrete']; sc = d['scan']
E = np.array(sc['ell']); lo = np.array(sc['Dhat_min']); hi = np.array(sc['Dhat_max'])
a.loglog(E, lo, 'C1-', label='Dhat >= (XENONnT heating)'); a.loglog(E, hi, 'C2-', label='Dhat <= (LPF self-sourced noise; LXe env.)')
ok = np.isfinite(np.array(sc['window'], float))
a.fill_between(E[ok], lo[ok], hi[ok], color='C2', alpha=0.25, label='discrete-kernel window (l >= %.1f um)' % (d['l_min'] * 1e6))
a.axhline(d['Dhat_detectable_LPFclass_kg'], color='C3', ls='--', label='source-mass lock-in reach: LPF-class sensor, 10 kg at 5 cm')
a.axhline(d['Dhat_detectable_torsion_1kg_1cm_kg'], color='C4', ls=':', label='same with Westphal torsion, 1 kg at 1 cm')
a.axvline(52e-6, color='gray', lw=0.8)
a.set_xlim(1e-7, 1e-4); a.set_ylim(1e-20, 1e-4)
a.set_xlabel('l (m)'); a.set_ylabel('Dhat (kg)'); a.set_title('Modulation works only for state-dependent kernels'); a.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, 'figures/stage27_signature.png'), dpi=130)
print('ok')
