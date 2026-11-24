import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = json.load(open('results/stage9.json')); s = d['kk_spectrum']
HBARC = 1.973269804e-16; R = d['kk_reduction']['MD_rows']['6.0']['R_m']
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
mE = [np.sqrt(n/2)*np.sqrt(2)*HBARC/R for n, c in s['E8_shells_norm_count']]; cE = np.cumsum([c for n, c in s['E8_shells_norm_count']])
mZ = [np.sqrt(n)*HBARC/R for n, c in s['Z8_shells_norm_count']]; cZ = np.cumsum([c for n, c in s['Z8_shells_norm_count']])
ax[0].step(mE, cE, where='post', label='E8 torus: 240, 2160, 6720, ...'); ax[0].step(mZ, cZ, where='post', label='Z8 torus: 16, 112, 448, ...')
for m, (n, c) in zip(mE[:3], s['E8_shells_norm_count'][:3]): ax[0].annotate(str(c), (m, 0), xytext=(m, 30), fontsize=8, color='C0')
ax[0].set_yscale('log'); ax[0].set_xlabel('KK graviton mass [GeV] (n=8, M_D=6 TeV, equal volume)'); ax[0].set_ylabel('cumulative number of modes')
ax[0].set_title('KK spectrum: E8 vs hypercubic torus'); ax[0].legend(fontsize=8)
t = np.linspace(0.2, 5, 400)
ax[1].plot(t, np.sqrt(2)/t, label='point-particle KK (momentum only): sqrt2/R, never 0')
ax[1].plot(t, np.sqrt(2)*np.abs(1/t-t), label='string momentum+winding: ~sqrt2|1/t - t|, 0 at self-dual point')
ax[1].axhline(0, c='k', lw=.6, ls='--', label='heterotic internal Gamma8 roots: massless at any R')
ax[1].set_xlabel('torus scale R / sqrt(alpha\')'); ax[1].set_ylabel('mass of the 240 root states [1/sqrt(alpha\')]'); ax[1].set_ylim(-0.3, 6)
ax[1].set_title('When do the 240 E8 roots become massless gauge bosons?'); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig('figures/stage9_e8_torus.png', dpi=110)
