#!/usr/bin/env python3
"""Stage 28 figure: LPF white floor vs time/temperature, zero-pressure extrapolation, profile likelihood, CQ bound."""
import json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
r = json.load(open(os.path.join(ROOT, 'results/stage28.json')))
runs = r['runs']; c = r['cold_run']; TV = 62.4; T0 = 293.0
t = np.array([q['t_minus_tv_table'] for q in runs]); T = np.array([q['T'] for q in runs])
S = np.array([q['S'] for q in runs]); sg = np.sqrt(np.array([q['sigma'] for q in runs]) ** 2 + (r['digitization_rel_uncertainty_assumed'] * S) ** 2)
A = r['variants']['A: all 13 runs (baseline)']
gas = lambda p, t, T: p['a'] * (TV / t) ** p['gamma'] * np.exp(-p['T_a'] * (1 / T - 1 / T0))
fig, ax = plt.subplots(2, 2, figsize=(13, 10))
a = ax[0, 0]
cold = T < 290
a.errorbar(t[~cold], S[~cold], sg[~cold], fmt='o', color='navy', label='warm runs (~295-297 K)')
a.errorbar(t[cold], S[cold], sg[cold], fmt='s', color='C0', label='runs 10, 11 (~284 K)')
a.errorbar([c['t_minus_tv']], [c['S']], [c['sigma']], fmt='D', color='c', label='cold run 274.8 K (not in baseline fit)')
tt = np.logspace(np.log10(40), np.log10(560), 200)
for TT, ls in ((296, '-'), (284.5, '--')):
    a.plot(tt, gas(A, tt, TT) + A['S0_fit'], 'k' + ls, lw=1, label='fit gas(T=%.1f K) + S0' % TT)
a.axhline(A['S0_fit'], color='C3', label='fitted floor S0 = %.2f' % A['S0_fit'])
a.axhspan(*A['S0_68'], color='C3', alpha=0.15)
a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel('t - t_v (days since venting)'); a.set_ylabel('white Δg PSD  [fm² s⁻⁴ Hz⁻¹]')
a.set_title('LPF white (Brownian) floor per run\n(digitized from Armano et al. 2024 Fig. 3; T from Table 1)'); a.legend(fontsize=7)
a = ax[0, 1]
g = gas(A, t, T)
a.errorbar(g, S, sg, fmt='o', color='navy', label='13 runs')
gc = gas(A, np.array([c['t_minus_tv']]), np.array([c['T']]))[0]
a.errorbar([gc], [c['S']], [c['sigma']], fmt='D', color='c', label='cold run (model-extrapolated gas term)')
x = np.linspace(0, 40, 50); a.plot(x, x + A['S0_fit'], 'k-', lw=1, label='S = gas + S0 (slope 1)')
a.errorbar([0], [A['S0_fit']], yerr=[[A['S0_fit'] - A['S0_68'][0]], [A['S0_68'][1] - A['S0_fit']]], fmt='*', color='C3', ms=14, label='zero-pressure intercept (68%)')
a.axhline(A['S0_95_upper'], color='C3', ls=':', label='95%% upper limit %.2f' % A['S0_95_upper'])
a.set_xlim(-1, 40); a.set_ylim(0, 40); a.set_xlabel('modelled gas term  a (t_v/t)^γ exp[-T_a(1/T-1/T0)]  (∝ pressure)'); a.set_ylabel('measured white PSD [fm² s⁻⁴ Hz⁻¹]')
ins = a.inset_axes([0.55, 0.08, 0.42, 0.35]); ins.errorbar(g, S, sg, fmt='o', color='navy', ms=3); ins.errorbar([gc], [c['S']], [c['sigma']], fmt='D', color='c', ms=3)
ins.plot(x, x + A['S0_fit'], 'k-', lw=1); ins.errorbar([0], [A['S0_fit']], yerr=[[A['S0_fit'] - A['S0_68'][0]], [A['S0_68'][1] - A['S0_fit']]], fmt='*', color='C3')
ins.set_xlim(-0.3, 8); ins.set_ylim(0, 9); ins.tick_params(labelsize=7)
a.set_title('Extrapolation to zero pressure'); a.legend(fontsize=7, loc='upper left')
a = ax[1, 0]
for k, v in r['variants'].items():
    gr = np.array(v['profile']['grid']); ch = np.array(v['profile']['chi2'])
    a.plot(gr, ch - ch.min(), label='%s: S0 = %.2f, 95%% < %.2f' % (k.split(':')[0] + ':' + k.split(':')[1][:26], v['S0_best'], v['S0_95_upper']))
a.axhline(2.71, color='gray', ls='--', lw=0.8); a.axhline(1, color='gray', ls=':', lw=0.8)
a.set_ylim(0, 12); a.set_xlabel('pressure-free floor S0 [fm² s⁻⁴ Hz⁻¹]'); a.set_ylabel('Δχ² (profile)'); a.set_title('Profile likelihood of the floor (a, γ, T_a floated)'); a.legend(fontsize=7)
a = ax[1, 1]
B = r['bounds']; names = list(B); vals = [B[n]['log10u'] for n in names]
y = np.arange(len(names))
a.barh(y, vals, color=['C0'] * 4 + ['C2', 'C7'])
for i, v in enumerate(vals): a.text(v + 0.01, i, '%.2f' % v, va='center', fontsize=8)
a.axvline(9.40, color='k', ls='-', label='Stage 26 (run 10 total ASD 1.74 fm s⁻²/√Hz): 9.40')
a.axvline(8.40, color='C3', ls='--', label='Stage 27 projection: 8.4')
a.set_yticks(y); a.set_yticklabels([n[:46] for n in names], fontsize=7); a.set_xlim(8.0, 9.6)
a.set_xlabel('log10 u = log10[Γ_bound / Γ_DP(κ=1)] at the plateau (orders short)')
a.set_title('CQ volume-law bound: Γ·V ≤ %.2g m⁵ s⁻³ (baseline)\nκ ∈ [%.1e, %.1e]' % (B[names[0]]['GammaV_bound'], B[names[0]]['kappa_min'], B[names[0]]['kappa_max']))
a.legend(fontsize=7, loc='lower left')
fig.tight_layout(); fig.savefig(os.path.join(ROOT, 'figures/stage28_lpf.png'), dpi=130); print('ok')
