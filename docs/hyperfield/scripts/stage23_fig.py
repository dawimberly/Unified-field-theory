"""Stage 23 figure: space design with measured inputs (reads results/stage23.json)."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import stage23_run as s23
d = json.load(open(os.path.join(ROOT, 'results', 'stage23.json')))
fig, ax = plt.subplots(2, 2, figsize=(14, 10.5))

# (a) measured vs required (log10 ratio)
a = ax[0, 0]
dp = d['dipole']; jt = d['jitter']; rl = d['release']
rows = [('permanent dipole, diamond (silica-scaled R^2), 1 mm', np.log10(dp['diamond_scaled_R2_eum'] / dp['required_eum']['0.001'])),
        ('permanent dipole, same, 1 cm', np.log10(dp['diamond_scaled_R2_eum'] / dp['required_eum']['0.01'])),
        ('  ... with spin averaging x120 (Afek 2021)', np.log10(dp['diamond_scaled_R2_eum'] / 120 / dp['required_eum']['0.01'])),
        ('S/C jitter 8 nm/rtHz vs required (LPF stiffness)', np.log10(8e-9 / jt['J_required'])),
        ('S/C jitter 30 nm/rtHz vs required', np.log10(30e-9 / jt['J_required'])),
        ('release velocity 10 um/s (LPF) vs 0.9 um/s', np.log10(rl['v_lpf'] / rl['v_required_1pct'])),
        ('free flight demonstrated 6 us vs 3 s', -np.log10(rl['t_free_demo'] / rl['t_free_needed'])),
        ('stray E at 1 mm: LPF-measured vs ion-trap model', d['efield']['lpf_over_iontrap_orders'])]
y = np.arange(len(rows)); v = [r[1] for r in rows]
a.barh(y, v, color=['tab:red' if x > 0 else 'tab:green' for x in v]); a.axvline(0, color='k')
a.set_yticks(y); a.set_yticklabels([r[0] for r in rows], fontsize=7); a.invert_yaxis()
for yi, x in zip(y, v): a.text(x + 0.1, yi, '%+.1f' % x, va='center', fontsize=7)
a.set_xlabel('log10(measured / required)'); a.set_title('(a) measured inputs vs what the space design needs', fontsize=10)

# (b) spectra (acceleration ASD, diamond design) vs requirement
b = ax[0, 1]
f = np.array(s23.s21.FGRID)
des = s23.design('diamond'); m = des['m']; p = s23.dipole('diamond', des['R'], 2)
for z, ls in ((1e-3, '-'), (1e-2, '--')):
    b.loglog(f, [p * s23.delta_asd(x) / z ** 2 / m for x in f], ls=ls, color='tab:red', label='permanent dipole %.0f e um, electrodes %g mm' % (p / s23.EUM, z * 1e3))
b.loglog(f, [p / 120 * s23.delta_asd(x) / 1e-2 ** 2 / m for x in f], ':', color='tab:red', label='same, 1 cm, spin-averaged x120')
for J, c in ((8e-9, 'tab:blue'), (30e-9, 'tab:cyan')):
    b.loglog(f, [s23.OMEGA2_SC * s23.jitter_asd(x, J) for x in f], color=c, label='S/C jitter %d nm x LPF stiffness 5.25e-7 s^-2' % round(J * 1e9))
b.axhline(des['req'], color='k', lw=2, label='requirement %.2g m s^-2/rtHz' % des['req'])
b.axvspan(1e-3, 0.02, color='grey', alpha=0.12, label='LPF-measured band')
b.set_xlabel('frequency (Hz)'); b.set_ylabel('acceleration ASD (m s^-2/rtHz)'); b.legend(fontsize=6.5, loc='lower left')
b.set_title('(b) measured-input spectra (D = 273 um diamond design)', fontsize=10)

# (c) scenarios
c = ax[1, 0]
sc = d['scenarios']
keys = ['diamond_measured_1mm', 'diamond_measured_1mm_R3', 'diamond_1cm', 'diamond_1cm_spin', 'diamond_1cm_spin_30nm',
        'silica_1mm', 'silica_1cm_spin', 'osmium_1mm', 'osmium_1cm', 'osmium_1cm_spin', 'diamond_1mm_one_electron']
y = np.arange(len(keys))
terms = ['permanent_dipole', 'induced_dipole', 'sc_jitter', 'charge']
cols = dict(permanent_dipole='tab:red', induced_dipole='tab:orange', sc_jitter='tab:blue', charge='tab:purple')
for i, k in enumerate(keys):
    for t in terms:
        if t in sc[k]['terms']:
            c.plot(sc[k]['terms'][t], i, 'o', color=cols[t], label=t if i == 0 or (t == 'charge' and 'electron' in k) else None)
    c.barh(i, sc[k]['total'], color='lightgrey', zorder=0)
    c.text(sc[k]['total'] + 0.15, i, '%+.1f' % sc[k]['total'], va='center', fontsize=7)
c.axvline(2.0, color='k', ls=':', label='Stage 21 (assumed inputs) +2.0')
c.set_yticks(y); c.set_yticklabels(keys, fontsize=7); c.invert_yaxis(); c.axvline(0, color='k')
c.set_xlabel('residual orders (grey = total)'); c.legend(fontsize=7, loc='lower right')
c.set_title('(c) recomputed space residual with measured inputs', fontsize=10)

# (d) D stability vs release velocity
dd = ax[1, 1]
v = np.geomspace(1e-9, 3e-5, 200); t = 3.0; D = des['D']
dd.loglog(v * 1e6, v * t / D, label='dD/D after 3 s')
dd.loglog(v * 1e6, 3 * v * t / D, label='dPhi/Phi = 3 dD/D')
dd.axvline(10, color='tab:red', ls='--', label='LPF in-flight release ~10 um/s')
dd.axvline(rl['v_required_1pct'] * 1e6, color='k', ls=':', label='needed for 1%%: %.2f um/s (T_eff <= %.1f mK)' % (rl['v_required_1pct'] * 1e6, rl['T_eff_required'] * 1e3))
dd.axhline(0.01, color='grey', lw=0.8)
dd.set_xlabel('release velocity spread (um/s)'); dd.set_ylabel('fractional change'); dd.legend(fontsize=7)
dd.set_title('(d) holding D = 273 um vs release velocity', fontsize=10)
fig.suptitle('Stage 23 - Stage 21 space design with measured inputs: +8.1 as drawn, +4.0 best (diamond), floor +2.6 from S/C jitter', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(ROOT, 'figures', 'stage23_space.png'); fig.savefig(out, dpi=140); print(out)
