#!/usr/bin/env python3
"""Stage 19 figure: required delocalisation F = dx/2z_zpf(100 Hz) vs mass (Stage 16 budget, Pb Meissner)
against what each spin-free route provides (demonstrated vs proposed)."""
import json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
d = json.load(open(os.path.join(ROOT, 'results/stage19.json')))
P = d['transient_penalty']
fig, (a, b) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw=dict(width_ratios=[1.35, 1]))
cur = d['requirements']['lead_curve']
for t, c in zip([0.5, 1.0, 2.5, 5.0], ['C0', 'C1', 'C2', 'C3']):
    pts = sorted([(float(x['m']), float(x['F'])) for x in cur if float(x['t']) == t])
    m, F = np.array(pts).T
    a.loglog(m, F, '-', color=c, label=f'required, Pb Meissner, hold {t:g} s')
    a.loglog(m, F * P, ':', color=c, lw=0.8)
lp = d['requirements']['stage16_lead_1e14kg_1s']; cp = d['requirements']['stage12_cheap_point']
a.loglog(lp['m'], lp['F'], 'k*', ms=14, label='Stage 16 Pb lead point (1e-14 kg, 1 s)')
a.loglog(cp['m'], cp['F'], 'kD', ms=8, label='Stage 12 cheap diamond point (2.5 s)')
R = {r['key']: r for r in d['routes']}
lines = [('rossi_2025', 'Rossi 2025 (demonstrated) 3.5', 'k', '-'),
         ('weiss_2021', 'Weiss 2021 loop 1e3', 'm', '--'),
         ('roda_llordes_2024_XL', 'Roda-Llordes XL 7.1e3', 'purple', '--'),
         ('roda_llordes_2024_XXL', 'Roda-Llordes XXL 7.1e4', 'r', '--'),
         ('hexapole_2026', 'Hexapole double well ~50', 'brown', '--')]
for k, lab, c, ls in lines:
    a.axhline(R[k]['F'], color=c, ls=ls, lw=1.3, label=lab + (' [demonstrated]' if ls == '-' else ' [proposed]'))
rp = R['pino_2018']; ri = R['romero_isart_2012']
a.loglog(rp['m'], rp['F'], 'o', mfc='none', mec='navy', ms=11, mew=2, label='Pino 2018 xi=607 nm (proposed, Nb)')
a.loglog(ri['m'], ri['F'], 's', mfc='none', mec='teal', ms=10, mew=2, label='Romero-Isart 2012 flux-qubit tag (proposed)')
a.set_ylim(3e-3, 3e6); a.set_xlabel('mass (kg)'); a.set_ylabel(r'$F=\Delta x/2z_{\rm zpf}$ (100 Hz)')
a.set_title('Spin-free delocalisation: needed (Stage 16 budget) vs available\n(dotted: x sqrt(3) transient penalty)')
a.legend(fontsize=7, loc='lower right'); a.grid(True, which='major', alpha=0.3)
# right: gap bars
items = [('Pino 2018\n(dx)', rp['gap_dx_orders'], 'proposed'),
         ('Pino\nvacuum', rp['gap_vacuum_orders'], 'proposed'),
         ('Pino\nnoise vs LPF', rp['gap_noise_vs_LPF_orders'], 'proposed'),
         ('RL XXL\n@Pb lead', R['roda_llordes_2024_XXL']['gap_dx_orders'], 'proposed'),
         ('RL XXL\n@cheap', R['roda_llordes_2024_XXL']['gap_dx_stage12_cheap'], 'proposed'),
         ('tag+XXL', d['combinations']['tag_then_inflate']['gap_orders'], 'proposed'),
         ('Hexapole\nQ', R['hexapole_2026']['gap_Q_orders'], 'proposed'),
         ('LPP\nT/Q', R['lami_pedernales_plenio_2024']['gap_TQ_vs_hofer_orders'], 'proposed'),
         ('Demonstrated\n@Pb lead', d['demonstrated_gap']['lead_1e14'], 'demonstrated'),
         ('RI 2012\ntag dx', ri['gap_dx_orders'], 'proposed')]
x = np.arange(len(items))
b.bar(x, [i[1] for i in items], color=['tab:green' if i[2] == 'demonstrated' else 'tab:orange' for i in items])
for xi, i in zip(x, items):
    b.text(xi, i[1] + 0.1, f'{i[1]:.2f}', ha='center', fontsize=8)
b.set_xticks(x); b.set_xticklabels([i[0] for i in items], fontsize=7, rotation=45, ha='right')
b.set_ylabel('gap (orders of magnitude)'); b.set_title('Remaining gaps (orange: vs proposals, green: vs demonstrated)')
b.grid(True, axis='y', alpha=0.3)
plt.tight_layout(); plt.savefig(os.path.join(ROOT, 'figures/stage19_tagging.png'), dpi=150)
print('ok')
