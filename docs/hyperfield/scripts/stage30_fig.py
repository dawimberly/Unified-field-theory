#!/usr/bin/env python3
"""Stage 30 figure (verification pass: all geometries from fetched sources; open square = V from mass x handbook density): (a) G residual vs test-body volume; (b) |z| vs volume (the actual CQ prediction); (c) CQ-predicted ppm vs 500 ppm."""
import json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
d = json.load(open(os.path.join(ROOT, 'results/stage30.json')))
T = [t for t in d['table'] if t['method'] != 'AI' and t['V_eff'] == t['V_eff']]
col = {'TB-dyn': 'C0', 'TB-comp': 'C1', 'SUSP': 'C2', 'STAT': 'C3'}
fig, ax = plt.subplots(1, 3, figsize=(17, 5))
for t in T:
    mk = 'o' if t['dims_V'] else 's'
    ax[0].errorbar(t['V_eff'], t['ppm'], yerr=t['u_ppm'], fmt=mk, color=col[t['method']], mfc='white' if mk == 's' else col[t['method']], capsize=2)
    ax[0].annotate(t['id'], (t['V_eff'], t['ppm']), fontsize=7, xytext=(3, 3), textcoords='offset points')
    ax[1].plot(t['V_eff'], abs(t['ppm'] / t['u_ppm']), mk, color=col[t['method']], mfc='white' if mk == 's' else col[t['method']])
for m, c in col.items(): ax[0].plot([], [], 'o', color=c, label=m)
ax[0].plot([], [], 'ks', mfc='white', label='V = cited mass / handbook density'); ax[0].plot([], [], 'ko', label='V from cited dimensions')
st = d['statistics']['macro_subset']['log10_Veff']
ax[0].axhline(0, color='gray', lw=0.8); ax[0].set_xscale('log')
ax[0].set_xlabel('test-body volume V_eff [m^3]'); ax[0].set_ylabel('(G - CODATA)/CODATA [ppm]')
sv = d['statistics']['sensitivity_verified_volumes']['drop_3_old_dumbbells_NIST_TRD_LANL']
ax[0].set_title('signed: Spearman %.2f (p=%.2f), WLS perm p=%.3f, LEE %.2f;\nwithout NIST-82/TR&D-96/LANL-97: perm p=%.2f' % (st['spearman_rho_ppm'], st['spearman_p_ppm'], st['wls_perm_p'],
                d['statistics']['look_elsewhere']['bonferroni'], sv['wls_perm_p']), fontsize=9)
ax[0].legend(fontsize=7)
a = d['statistics_extra']['abs_residual_vs_logV']
ax[1].set_xscale('log'); ax[1].set_xlabel('V_eff [m^3]'); ax[1].set_ylabel('|z| = |G - CODATA|/u')
ax[1].set_title('CQ predicts more scatter at small V: Spearman %.2f (p=%.2f) - null' % (a['spearman_rho_absz'], a['p_absz']), fontsize=9)
P = d['physics']['per_experiment']
ids = [p['id'] for p in P]; xx = np.arange(len(P))
ax[2].semilogy(xx, [p['scatter_ppm_at_max'] for p in P], 'o', label='CQ scatter at max allowed (Stage 28 / X-ray for atoms)')
ax[2].semilogy(xx, [p['scatter_ppm_DP_floor'] for p in P], 'v', label='CQ scatter at DP floor')
ax[2].semilogy(xx, [p['matter_sourced_ppm'] for p in P], 'x', label='matter-sourced kernel (Stage 29 ceiling)')
ax[2].semilogy(xx, [p['rectification_bias_ppm'] for p in P], 'd', label='rectification bias bound (time-of-swing)')
ax[2].axhline(500, color='r', label='500 ppm spread'); ax[2].axhline(d['scatter']['tau_ppm'], color='r', ls=':', label='excess scatter tau = %.0f ppm' % d['scatter']['tau_ppm'])
ax[2].set_xticks(xx); ax[2].set_xticklabels(ids, rotation=70, fontsize=7); ax[2].set_ylabel('ppm (1 day integration)')
ax[2].set_title('CQ effect is >= %.1f orders below 500 ppm' % min(p['orders_below_500ppm'] for p in P), fontsize=9); ax[2].set_ylim(1e-12, 1e5); ax[2].legend(fontsize=6, loc='lower right')
plt.tight_layout(); plt.savefig(os.path.join(ROOT, 'figures/stage30_bigG.png'), dpi=130); print('ok')
