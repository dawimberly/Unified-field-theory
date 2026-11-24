import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from math import pi
MBAR, MPL = 2.435e18, 1.22089e19
xi = np.linspace(-1, 0.4, 400)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for nW, lab in [(45, 'SM (45 Weyl)'), (48, 'SM + nu_R (48 Weyl)')]:
    N = 4*6*(1/6-xi) + nW - 48
    ax[0].plot(xi, N, label=lab)
    L = np.where(N > 0, MBAR*np.sqrt(96*pi**2/np.where(N > 0, N, 1)), np.nan)
    ax[1].semilogy(xi, L, label=lab)
N2 = 4*6*(1/6-xi) + 45 + 24
ax[0].plot(xi, N2, '--', label='vectors counted as 2 scalars (Kabat-type ambiguity)')
ax[1].semilogy(xi, MBAR*np.sqrt(96*pi**2/N2), '--', label='vectors as 2 scalars')
ax[0].axhline(0, c='k', lw=.6); ax[0].axvline(0, c='gray', ls=':'); ax[0].axvline(1/6, c='gray', ls=':')
ax[0].text(0.005, 30, 'xi=0', fontsize=8); ax[0].text(0.17, 30, 'xi=1/6', fontsize=8)
ax[0].set_xlabel('Higgs non-minimal coupling xi'); ax[0].set_ylabel('N_eff (units of one minimal real scalar)')
ax[0].set_title('Induced 1/G sign: SM one loop (vectors -4 each)'); ax[0].legend(fontsize=7); ax[0].set_ylim(-30, 100)
for y, t in [(MPL, 'M_Pl'), (MBAR, 'reduced M_Pl'), (2e16, 'GUT')]:
    ax[1].axhline(y, c='gray', ls=':'); ax[1].text(-0.98, y*1.2, t, fontsize=8)
ax[1].set_xlabel('xi'); ax[1].set_ylabel('cutoff Lambda needed [GeV]'); ax[1].set_title('Cutoff needed for 1/(8 pi G) = (2.435e18 GeV)^2'); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig('figures/stage8_induced_G.png', dpi=110)
