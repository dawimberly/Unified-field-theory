#!/usr/bin/env python3
"""Stage 29 figure: (a) volume-law scope (shape, finite l/R); (b) cross-correlation checks; (c) discrete-kernel window, Stage 27 vs corrected."""
import json, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import erf
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
d = json.load(open(os.path.join(ROOT, 'results/stage29.json')))
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
v = d['claim1_volume_law']
sh = v['shape_spheroids']; m = [s['c_over_a'] for s in sh]
ax[0].semilogx(m, [s['GammaV_axis_over_claim'] for s in sh], 'o-', label='spheroid, along symmetry axis (3N_zz)')
ax[0].semilogx(m, [s['GammaV_transverse_over_claim'] for s in sh], 's-', label='spheroid, transverse (3N_xx)')
f = v['finite_l_over_R']
ax[0].semilogx(f['l_over_R'], f['Gamma_ratio_to_plateau'], 'k^--', label='sphere, finite l/R (x-axis = l/R)')
ax[0].axhline(1, color='gray', lw=0.8)
ax[0].set_xlabel('aspect ratio c/a   (or l/R for the dashed curve)'); ax[0].set_ylabel(r'$\Gamma V / [(4\pi/3)G\hbar(\kappa+1/4\kappa)]$')
ax[0].set_title('(1) volume law holds for N = 1/3 (sphere, cube) only'); ax[0].legend(fontsize=8)
dd = np.logspace(np.log10(2.0), 2, 200)
ax[1].loglog(dd, 3 * (4 / 3 * np.pi) / (2 * np.pi * dd ** 3), 'C0-', label=r'DP family: $|\rho_{long}| = (3/2\pi)V/d^3$ (any $\kappa$)')
fc = d['claim2_cross']['fourier_checks']
for k, c in fc.items():
    dr = float(k.split('=')[1][:-1]); ax[1].loglog(dr, -c['rho_long_numeric'], 'ko')
ax[1].loglog([], [], 'ko', label='independent Fourier integral (Stage 29)')
ax[1].axhline(1e-3, color='gray', ls=':', lw=0.8)
ax[1].text(2.1, 2e-6, 'KTM (LOCC): rho = 0\nTD local gamma: rho_long = 0', fontsize=8)
ax[1].set_xlabel('d / R (equal spheres)'); ax[1].set_ylabel(r'$-\rho_{long}$ (anti-correlation)')
ax[1].set_title('(2) -2Għ/d³ confirmed; sign pattern is kernel-specific'); ax[1].legend(fontsize=8)
s = d['claim3_discrete']['scan']; E = np.array(s['ell'])
ax[2].loglog(E, s['Dmin_isolated'], 'C3--', label='Dhat_min, X-ray, isolated nucleus (Stage 27)')
ax[2].loglog(E, s['Dmin_corrected'], 'C3-', label='Dhat_min corrected (total density; + interferometry)')
ax[2].loglog(E, s['Dmax'], 'C0-', label='Dhat_max (LPF self-noise)')
ax[2].axhline(d['claim3_discrete']['lockin_Dhat_detectable_kg_s'], color='C2', ls=':', label='source-mass lock-in reach (Stage 27)')
ax[2].axvline(d['claim3_discrete']['l_min_stage27_isolated_reproduced'], color='C3', ls='--', lw=0.7)
ax[2].axvline(d['claim3_discrete']['l_min_dense_corrected'], color='C3', lw=0.7)
ax[2].axvline(52e-6, color='gray', lw=0.7)
ax[2].set_ylim(1e-33, 1e-8); ax[2].set_xlabel('smearing length l [m]'); ax[2].set_ylabel('Dhat [kg s]')
ax[2].set_title('(3) discrete kernel: l_min 2.5 um -> 5.5 nm'); ax[2].legend(fontsize=7, loc='upper right')
plt.tight_layout(); plt.savefig(os.path.join(ROOT, 'figures/stage29_checks.png'), dpi=130)
print('ok')
