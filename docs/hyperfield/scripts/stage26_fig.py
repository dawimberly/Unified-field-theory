#!/usr/bin/env python3
"""Stage 26 figure: surviving CQ parameter space in (l, kappa) and the reach of each experiment."""
import json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
d = json.load(open(os.path.join(ROOT, 'results/stage26.json')))
E = np.array(d['scan']['ell']); U = {k: np.array(v, float) for k, v in d['scan']['u'].items()}
arr = lambda v: np.array([np.nan if x is None else x for x in v], float)
def kpm(u):
    u = np.asarray(u, float); ok = u >= 1; kp = np.where(ok, (u + np.sqrt(np.where(ok, u * u - 1, 0))) / 2, np.nan)
    return 1 / (4 * kp), kp
HB2G = d['D2_map']['hbar_over_2G']
fig, ax = plt.subplots(1, 3, figsize=(19, 6.0))
sty = {'XENONnT 2026 Xe (PRL 136, 120201)': ('C3', '-', 2.4, 'XENONnT 2026 X-rays'),
       'Donadi 2021 Ge (Nat Phys 17, 74)': ('C1', '-', 1.6, 'Donadi 2021 X-rays'),
       'LISA Pathfinder 2018 (Armano)': ('C0', '-', 2.0, 'LISA Pathfinder (differential, L1)'),
       'Gisler 2022 SiN nanowire': ('k', '-', 1.2, 'Gisler 2022 nanowire (best absolute)'),
       'Pedalino 2025 172 kDa Na cluster': ('C2', '-', 1.2, 'Pedalino 2025 interferometer'),
       'MAQRO baseline 1e10 u silica, dx=R, 100 s': ('C4', '--', 1.6, 'MAQRO baseline 1e10 u (planned)'),
       'MAQRO max 1e11 u silica, dx=R, 100 s': ('C4', ':', 1.8, 'MAQRO max 1e11 u (planned)'),
       'Pino 2018 Nb 1 um, dx 0.5 um, 0.596 s': ('C5', '--', 1.8, 'Pino 2018 Nb sphere (proposal)'),
       'LISA requirement (single TM)': ('C9', '--', 1.2, 'LISA requirement (planned)'),
       'MAGIS-100 Sr atom, 1 m, 9 s': ('C8', '--', 1.0, 'MAGIS-100 (planned)'),
       'Stage-12 BMV diamond 5.6e-15 kg, 3.85 um, 2.5 s (our design)': ('C7', '-.', 1.0, 'our Stage-12 BMV diamond (design)')}
# (a) allowed region in (l, kappa)
a = ax[0]
S1 = d['scenarios']['S1_markov_full_existing']; km, kp = arr(S1['kappa_minus']), arr(S1['kappa_plus'])
a.fill_between(E, km, kp, where=np.isfinite(km), color='C2', alpha=0.25, label='SURVIVES: all existing data (Markovian)')
S4 = d['scenarios']['S4_markov_full_existing_plus_planned']; km4, kp4 = arr(S4['kappa_minus']), arr(S4['kappa_plus'])
a.fill_between(E, km4, kp4, where=np.isfinite(km4), color='C5', alpha=0.35, label='survives MAQRO + Pino + LISA too')
for n, (c, ls, lw, lab) in sty.items():
    if n.startswith('MAGIS') or n.startswith('Stage-12') or n.startswith('LISA req'): continue
    m, p = kpm(U[n]); a.loglog(E, p, color=c, ls=ls, lw=lw, label=lab); a.loglog(E, m, color=c, ls=ls, lw=lw)
a.axvline(d['inputs']['L_ISL'], color='gray', ls='-', lw=1); a.text(d['inputs']['L_ISL'] * 0.5, 1e-25, '1/r^2 tested\nto 52 um', fontsize=7, ha='right')
a.axhline(0.5, color='gray', lw=0.6, ls=':'); a.text(2e-16, 0.6, 'kappa = 1/2: least-noise point (Tilloy-Diosi)', fontsize=7)
a.axvline(1e-15, color='gray', lw=0.6); a.text(1.1e-15, 1e22, 'Oppenheim r_N', fontsize=7, rotation=90)
a.set_xlim(1e-16, 1e-4); a.set_ylim(1e-30, 1e30)
a.set_xlabel('smearing length l of the mass density [m]'); a.set_ylabel('kappa = D0 / D0(Diosi)   (decoherence strength)')
a.set_title('(a) allowed (l, kappa) region, DP kernel (most permissive)'); a.legend(fontsize=6.4, loc='lower left')
a2 = a.twinx(); a2.set_yscale('log'); a2.set_ylim(HB2G / 1e-30, HB2G / 1e30)
a2.set_ylabel('saturated l_P^2 D2 = hbar/(2 G kappa)  [kg^2 s/m]  (Stage 25 normalisation)', fontsize=8)
# (b) reach: log10 u per experiment
b = ax[1]
for n, (c, ls, lw, lab) in sty.items():
    b.semilogx(E, np.log10(U[n]), color=c, ls=ls, lw=lw, label=lab)
b.axhline(0, color='k', lw=1); b.fill_between(E, -10, 0, color='C3', alpha=0.06)
b.text(2e-13, -9.3, 'u < 1: experiment alone excludes this l for EVERY kernel (any kappa)', fontsize=7.5)
b.axvline(4.5e-10, color='C3', ls=':', lw=1); b.text(5e-10, -6, 'l_min = R0(XENONnT) = 4.5 A', fontsize=7.5, color='C3')
b.axvline(d['inputs']['L_ISL'], color='gray')
b.set_xlim(1e-16, 1e-4); b.set_ylim(-10, 30)
b.set_xlabel('smearing length l [m]'); b.set_ylabel('log10 u_E = log10(measured bound / Diosi(kappa=1) prediction)')
b.set_title('(b) orders short (u>1) or excluding (u<1), fixed published sensitivities'); b.legend(fontsize=6.4, loc='upper left')
# (c) window size in orders of kappa
c = ax[2]
for s, col, lab in (('S1_markov_full_existing', 'C2', 'existing, X-rays see D0 + field noise (primary)'),
                    ('S2_lindblad_only_xray_existing', 'C1', 'existing, X-rays see D0 term only (Stage 25 assumption)'),
                    ('S3_no_xray_existing', 'C0', 'existing, noise coloured below X-ray band (no X-ray bound)'),
                    ('S4_markov_full_existing_plus_planned', 'C5', '+ MAQRO, Pino, LISA, MAGIS (planned/proposed)'),
                    ('S5_markov_full_plus_planned_plus_our_designs', 'C7', '+ our Stage 12/16/22 designs')):
    W = arr(d['scenarios'][s]['window_orders']); c.semilogx(E, W, color=col, lw=2, label=lab)
c.axvline(d['inputs']['L_ISL'], color='gray'); c.axvline(4.5e-10, color='C3', ls=':', lw=1)
c.set_xlim(1e-16, 1e-4); c.set_ylim(0, 32)
c.set_xlabel('smearing length l [m]'); c.set_ylabel('surviving window  log10(kappa+/kappa-)  [orders]')
c.set_title('(c) size of the surviving window vs l'); c.legend(fontsize=7, loc='upper left')
fig.suptitle('Stage 26: what survives of classical-quantum Newtonian gravity - every probe sees A_E(l)(kappa + 1/4kappa); X-rays set l >= 4.5 A for any kernel', fontsize=11.5)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, 'figures/stage26_window.png'), dpi=140)
print('ok')
