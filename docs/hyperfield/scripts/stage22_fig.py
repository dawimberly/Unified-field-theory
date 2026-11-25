#!/usr/bin/env python3
"""Stage 22 figure: closing the non-local Oppenheim window (reads results/stage22.json)."""
import json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
R = json.load(open(os.path.join(ROOT, 'results/stage22.json')))
H = 6.62607015e-34; U = 1.66053906660e-27
fig, ax = plt.subplots(2, 2, figsize=(13, 10))

a = ax[0, 0]
for k, c in R['closing_curve'].items():
    xs = [r['xi'] for r in c['rows']]; ms = [r['M_kDa'] / 1e3 for r in c['rows']]
    a.plot(xs, ms, 'o-', label=f"{k} (d={c['d_nm']:.0f} nm)")
ch = R['designs'][R['cheapest']]
a.plot(ch['xi'], ch['M_kDa'] / 1e3, 'k*', ms=18, label='cheapest robust point')
a.plot(0.983 * H / (1e6 * U * 25 * (133e-9) ** 2), 1.0, 'rs', ms=9, label='MUSCLE plan 1 MDa @25 m/s (0.56 left)')
a.plot(0.983 * H / (172e3 * U * 160 * (133e-9) ** 2), 0.172, 'g^', ms=9, label='Pedalino 2025, 172 kDa (floor)')
a.set_xlabel(r'$\xi = L/L_T$'); a.set_ylabel('cluster mass needed to close [MDa]'); a.set_yscale('log')
a.set_title(r'(a) closing curve $2\xi M^4 d^{5/3}/(h m_a r_N^{2/3}) \geq$ Gisler ceiling (v drops out)')
a.legend(fontsize=7); a.grid(alpha=.3)

a = ax[0, 1]
v = np.logspace(0.8, 2.3, 200); L = 0.983; d = 133e-9
Rl = (1e-30 * d) ** (1 / 3); m_na = 22.98976928 * U
Mc = (R['U_ceiling'] * m_na * Rl * v / (2 * L)) ** (1 / 3) / U / 1e6
a.plot(v, Mc, 'k-', lw=2, label='closing line (L = 0.983 m, d = 133 nm)')
for xi, ls in ((0.5, ':'), (0.8, '--'), (1.5, '-.')):
    a.plot(v, L * H / (xi * v * d ** 2) / U / 1e6, ls, color='C1', label=fr'$\xi$ = {xi}')
a.fill_between(v, L * H / (1.5 * v * d ** 2) / U / 1e6, L * H / (0.5 * v * d ** 2) / U / 1e6, color='C1', alpha=.12,
               label=r'robust q/cl discrimination band (0.5 $\leq\xi\leq$ 1.5, model)')
a.plot(ch['v_m_s'], ch['M_kDa'] / 1e3, 'k*', ms=18)
a.plot(25, 1.0, 'rs', ms=9); a.plot(160, 0.172, 'g^', ms=9)
a.axvspan(120, 170, color='C2', alpha=.1, label='demonstrated neutral beam 120-170 m/s')
a.set_xscale('log'); a.set_yscale('log'); a.set_ylim(0.1, 10); a.set_xlabel('velocity [m/s]'); a.set_ylabel('mass [MDa]')
a.set_title('(b) L = 0.983 m: above black line and inside band = closes, q/cl distinct', fontsize=10)
a.legend(fontsize=7); a.grid(alpha=.3, which='both')

a = ax[1, 0]
ks = [k for k in R['designs'] if k.count('_') == 2]
xs = [R['designs'][k]['xi'] for k in ks]
dv = [R['designs'][k]['best']['dV_worst_alpha'] for k in ks]
t5 = [R['designs'][k]['best']['T5sigma_s'] for k in ks]
a.plot(xs, dv, 'o-', color='C0', label=r'worst-case $\Delta V$ (q - cl), $\alpha\pm$11%, $\sigma_{ion}\pm$20%')
a.axhline(0.05, color='C0', ls=':', lw=1)
a.set_xlabel(r'$\xi$ along the closing curve (ion route, L = 0.983 m)'); a.set_ylabel(r'$\Delta V$', color='C0')
b = a.twinx(); b.semilogy(xs, t5, 's-', color='C3', label='5-sigma statistical time [s]'); b.set_ylabel('T(5 sigma) [s]', color='C3')
a.set_title('(c) discrimination and run time along the curve'); a.legend(loc='upper left', fontsize=7); b.legend(loc='upper right', fontsize=7)
a.grid(alpha=.3)

a = ax[1, 1]
names = [r['name'][:48] for r in R['routes']] + [c['name'][:48] for c in R['ceiling']['candidates']]
vals = [-r['orders_left'] for r in R['routes']] + [c['margin_orders'] for c in R['ceiling']['candidates']]
cols = ['C2' if not r['flag'] else 'C7' for r in R['routes']] + ['C4' if 'differential' not in c['name'] else 'C7' for c in R['ceiling']['candidates']]
y = np.arange(len(names))
a.barh(y, vals, color=cols); a.axvline(0, color='k')
a.set_yticks(y); a.set_yticklabels(names, fontsize=6); a.invert_yaxis()
a.set_xlabel('orders beyond closing (>0 closes; <0 = orders still short)')
a.set_title('(d) floor routes (green), ceiling sensors (purple), grey = flagged', fontsize=10)
a.set_xlim(-5, 7); a.grid(alpha=.3, axis='x')
fig.suptitle('Stage 22 - closing the non-local Oppenheim window (3.66 orders: Pedalino 2025 floor vs Gisler 2022 ceiling)', fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, 'figures/stage22_closing.png'), dpi=130)
print('ok')
