import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import importlib.util, sys
spec = importlib.util.spec_from_file_location('g', 'scripts/stage10_gw.py'); g = importlib.util.module_from_spec(spec); sys.argv=['x']
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(g)
out = g.out
# extra checks: stage-5 headline (peak 0.035 Hz) -> beta/H, and g*, v_w sensitivity
bH5 = 0.035/(1.65e-5*60*(106.75/100)**(1/6))
ex = {'stage5_betaH_implied_by_0.035Hz': bH5}
for a in (0.05, 0.1, 0.3, 1.0):
    ex[f'stage5_alpha{a}'] = g.all_snr(6000, a, bH5)
for gs in (106.75, 150, 200):
    ex[f'gstar{gs}_alpha0.1_bH35'] = g.all_snr(6000, 0.1, bH5, gs=gs)['LISA']
for vw in (0.3, 0.6, 0.95):
    ex[f'vw{vw}_alpha0.1_bH35'] = g.all_snr(6000, 0.1, bH5, vw=vw)['LISA']   # kappa_v kept at v_w~1 fit (ASSUMED)
ex['with_collisions_kphi1_alpha0.1_bH35'] = g.all_snr(6000, 0.1, bH5, kappa_phi=1.0)['LISA']
ex['no_sw_suppression_alpha0.1_bH35'] = g.all_snr(6000, 0.1, bH5, suppress=False)['LISA']
# hysteresis-derived lower bound on beta/H: delta <= 1/(2 a1)
Sc = out['lock_scan'][0]['S3_over_T']
ex['betaH_min_vs_a1'] = {a1: max(Sc*(2*(1-0.5/a1)/(0.5/a1)-1), 0) for a1 in (0.6, 0.75, 1, 1.5, 2, 3)}
# T at which LISA SNR peaks for alpha=0.1 beta/H=100
out['extra'] = ex; g.save()
print(json.dumps(ex, indent=1, default=float))

fig, ax = plt.subplots(1, 3, figsize=(18, 5.5))
f = g.FG
for n, c in [('LISA', 'k'), ('DECIGO', 'tab:purple'), ('BBO', 'tab:blue'), ('Einstein Telescope (crude)', 'tab:green'), ('Cosmic Explorer (crude)', 'tab:olive')]:
    Sn, (lo, hi), _ = g.DET[n]; ff = f[(f >= lo) & (f <= hi)]; ax[0].loglog(ff, g.omega_n(Sn, ff), c, lw=1.5, label=n+' noise')
ax[0].hlines(g.LVK_h2, 20, 76.6, 'r', lw=3, label='LVK O3 bound (flat)')
for (a, b, ls, lab) in [(0.1, bH5, '-', f'α=0.1, β/H={bH5:.0f} (Stage 5 0.035 Hz)'), (0.5, 10, '--', 'α=0.5, β/H=10'),
                        (0.0277, 401, ':', 'lock natural median α=0.028, β/H=400'), (0.12, 4614, '-.', 'lock strongest α=0.12, β/H=4600')]:
    h2, _ = g.total(f, 6000, a, b); ax[0].loglog(f, h2, ls, color='tab:orange', label=lab)
ax[0].set_ylim(1e-18, 1e-6); ax[0].set_xlim(1e-5, 1e3); ax[0].set_xlabel('f [Hz]'); ax[0].set_ylabel('h² Ω'); ax[0].legend(fontsize=7); ax[0].set_title('T*=6 TeV spectra vs sensitivity')
m = out['map']; A, B = np.meshgrid(m['beta_over_H'], m['alpha'])
cs = ax[1].contourf(A, B, np.log10(np.maximum(np.array(m['LISA']), 1e-6)), levels=np.arange(-6, 7), cmap='viridis')
plt.colorbar(cs, ax=ax[1], label='log10 SNR (LISA)')
for k, c in [('LISA', 'w'), ('DECIGO', 'tab:purple'), ('BBO', 'tab:cyan')]:
    ax[1].contour(A, B, np.array(m[k]), levels=[10], colors=c, linewidths=2)
    ax[1].plot([], [], color=c, label=f'{k} SNR=10')
pts = [r for r in out['lock_scan'] if r['alpha']]
ax[1].scatter([r['beta_over_H'] for r in pts], [r['alpha'] for r in pts], c='r', marker='x', s=50, label='lock model, natural O(1) inputs')
ax[1].axvline(bH5, color='orange', ls='--', label='Stage-5 β/H≈35')
ax[1].set_xscale('log'); ax[1].set_yscale('log'); ax[1].set_xlabel('β/H*'); ax[1].set_ylabel('α'); ax[1].legend(fontsize=7, loc='lower left'); ax[1].set_title('T*=6 TeV, v_w=0.95, g*=106.75')
ts = [r for r in out['T_scan'] if r['alpha'] == 0.1]
for n, c in [('LISA', 'k'), ('DECIGO', 'tab:purple'), ('BBO', 'tab:blue'), ('Cosmic Explorer (crude)', 'tab:olive')]:
    ax[2].loglog([r['T'] for r in ts], [max(r['snr'][n], 1e-8) for r in ts], c, label=n)
ax[2].axhline(10, color='gray', ls='--'); ax[2].axvline(6000, color='orange', ls=':')
ax[2].set_ylim(1e-6, 1e7); ax[2].set_xlabel('T* [GeV]'); ax[2].set_ylabel('SNR'); ax[2].set_title('α=0.1, β/H=100: SNR vs T*'); ax[2].legend(fontsize=8)
plt.tight_layout(); plt.savefig('figures/stage10_gw.png', dpi=120); print('ok')
