#!/usr/bin/env python3
import json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
d = json.load(open(os.path.join(ROOT, 'results/stage25.json')))
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
cols = {'Pedalino 172 kDa': 'C0', 'Stage-22 1.415 MDa': 'C3'}
# (a) lattice sum vs regimes, dx = 133 nm
a = ax[0]
for lab, v in d['lattice'].items():
    rows = v['cases']['dx = d = 133 nm']
    l = np.array([r['ell'] for r in rows]); S = np.array([r['S_lattice'] for r in rows])
    a.loglog(l, S, 'o', color=cols[lab], label=f'{lab}: bcc lattice sum (N={v["N"]})')
    a.loglog(l, [r['S_incoherent'] for r in rows], '--', color=cols[lab], lw=1, label='incoherent  2N/(sqrt(pi) l)  [= M^2/N]')
    a.loglog(l, [r['S_coherent_sphere'] for r in rows], ':', color=cols[lab], lw=1.5, label='coherent sphere  2.4 N^2/R  [= M^2]')
    a.loglog(l, [r['S_blob_l_gg_R'] for r in rows], '-.', color=cols[lab], lw=0.8, alpha=0.6, label='l >> R blob  [M^2/l]')
    a.loglog(l, [r['S_lattice'] / r['ratio_to_NM2_form'] for r in rows], '-', color=cols[lab], lw=2.5, alpha=0.35, label='published form N*M^2/R_lam (r_N=1e-15)')
a.axvline(1e-15, color='k', lw=0.6); a.text(1.1e-15, 1e10, 'l = r_N', fontsize=8)
a.set_ylim(1e9, 1e30); a.set_xlabel('smearing length l [m]'); a.set_ylabel('S = sum_ij [2K(r) - K(r+dx) - K(r-dx)]  [1/m]')
a.set_title('(A) Eq.137 on a Na bcc cluster, dx = 133 nm'); a.legend(fontsize=6.3, loc='upper right')
# (b) effective N exponent
b = ax[1]
for dlab, mk in (('dx = d = 133 nm', 'o-'), ('dx = 1 A (< lattice spacing)', 's--')):
    p = d['A']['N_power'][dlab]
    b.semilogx([x['ell'] for x in p], [x['p_N'] for x in p], mk, label=dlab)
for y, t in ((1, 'N^1 (incoherent, M^2/N)'), (5 / 3, 'N^5/3 (coherent, M^2/R)'), (2, 'N^2 (coherent, l >> R)'), (3, 'N^3 = N*M^2 (Eq. 47 as published)')):
    b.axhline(y, color='gray', lw=0.7, ls=':'); b.text(1.2e-16, y + 0.04, t, fontsize=7.5)
b.set_ylim(0.8, 3.3); b.set_xlabel('smearing length l [m]'); b.set_ylabel('p in S ~ N^p (fixed l, dx; 7,482 -> 61,600 atoms)')
b.set_title('(A) N-scaling of the decoherence rate: never N*M^2'); b.legend(fontsize=8, loc='center right')
# (c) floors vs ceilings
c = ax[2]; sc = d['B']['scan']; l = np.array([r['ell'] for r in sc])
c.loglog(l, [r['floor_interf_published'] for r in sc], 'C0-', label='interferometric floor, Eq.47 as published (Pedalino)')
c.loglog(l, [r['floor_interf_derived'] for r in sc], 'C0--', label='interferometric floor, derived counting (Pedalino)')
c.loglog(l, [r['floor_xray_Donadi'] for r in sc], 'C1-', lw=2, label='X-ray floor, Donadi 2021 Ge (Diosi conv. R0>0.18 A)')
c.loglog(l, [r['floor_xray_XENONnT'] for r in sc], 'C3-', lw=2, label='X-ray floor, XENONnT 2026 (R0>4.5 A)')
c.loglog(l, [r['ceiling_floor_only'] for r in sc], 'k-', label='ceiling (Gisler), r_N fixed at 1e-15 (floor-only)')
c.loglog(l, [r['ceiling_consistent'] for r in sc], 'k--', label='ceiling, r_N = l (consistent)')
B = d['B']
for k, m in (('l_min_consistent_Donadi', 'C1'), ('l_min_consistent_XENONnT', 'C3')):
    c.axvline(B[k], color=m, ls=':', lw=1)
c.axvline(1e-15, color='gray', lw=0.6)
c.annotate('Donadi excludes l=1e-15\nby %.1f orders\n(XENONnT %.1f)' % (B['exclusion_orders_Donadi'], B['exclusion_orders_XENONnT']),
           xy=(1e-15, B['floor_xray_Donadi_l1em15']), xytext=(2e-14, 1e5), fontsize=8, arrowprops=dict(arrowstyle='->'))
c.set_ylim(1e-40, 1e15); c.set_xlabel('smearing length l [m]'); c.set_ylabel('l_P^2 D2 bound  [kg^2 s m^-1]  (Eq. 47 units)')
c.set_title('(B) radiation bounds carried over via the trade-off'); c.legend(fontsize=6.5, loc='lower left')
fig.suptitle('Stage 25: testing the Stage 24 flags - (A) Eq.47 N-counting, (B) DP spontaneous-radiation carry-over', fontsize=12)
fig.tight_layout()
os.makedirs(os.path.join(ROOT, 'figures'), exist_ok=True)
fig.savefig(os.path.join(ROOT, 'figures/stage25_tests.png'), dpi=140)
print('ok')
