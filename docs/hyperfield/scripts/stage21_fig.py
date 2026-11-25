"""Stage 21 figure: levers that loosen the BMV acceleration-noise requirement (reads results/stage21.json)."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import stage21_run as s21
d = json.load(open(os.path.join(ROOT, 'results', 'stage21.json')))
fig, ax = plt.subplots(2, 2, figsize=(14, 10.5))

a = ax[0, 0]
f = np.array(d['lever1_echo']['f'])
for r in d['lever1_echo']['same_hold'][:3]:
    a.loglog(f, r['W2'], lw=1.5, label='|W(f)|^2, %d branch swap(s), t_tot %.1f s' % (r['n_swap'], r['t_tot']))
a.set_ylim(1e-6, 30); a.set_xlabel('frequency (Hz)'); a.set_ylabel('phase filter |W(f)|^2 (s^2)')
a2 = a.twinx()
for k, ls in (('seismic_nlnm', '-'), ('seismic_iso_flat', '--'), ('tilt', ':')):
    fn = s21.SPECTRA[k][0]
    a2.loglog(f, [fn(x) for x in f], color='k', ls=ls, lw=1, label=k)
a2.set_ylabel('acceleration ASD (m s^-2/rtHz), black'); a2.legend(fontsize=7, loc='upper right')
a.legend(fontsize=7, loc='lower left')
a.set_title('(a) echo moves the response off DC onto ~1/T, i.e. onto the 0.1-0.3 Hz microseism', fontsize=10)

b = ax[0, 1]
rows = d['lever1_hold']
for n, mk in ((0, 'o-'), (2, 's--')):
    xs = [r['t_hold'] for r in rows if r['n_swap'] == n and r['req_white'] > 0]
    ys = [r['req_white'] * 1e15 for r in rows if r['n_swap'] == n and r['req_white'] > 0]
    b.semilogx(xs, ys, mk, label='N5<=1e4, %d swaps' % n)
l2 = [r for r in d['lever2_shots']['rows'] if r['req_white'] > 0]
b.semilogx([r['t_hold'] for r in l2], [r['req_white'] * 1e15 for r in l2], '^:', label='fixed wall-clock (lever 2)')
b.axhline(d['scaling']['grossardt_sqrt_2Ghbar_D3'] * 1e15, color='r', ls='-.', label='sqrt(2 G hbar/D^3) = %.2f fm (ceiling)' % (d['scaling']['grossardt_sqrt_2Ghbar_D3'] * 1e15))
b.axhline(0.442, color='k', ls=':', label='Stage 18 (N5=5049) 0.44 fm')
b.set_xlabel('hold time (s)'); b.set_ylabel('allowed white sqrt(S_a) (fm s^-2/rtHz)')
b.set_title('(b) hold-time trade (holds < 2.5 s: N5 cap missed; 20 s+swaps: gas)', fontsize=10)
b.legend(fontsize=7)

c = ax[1, 0]
base = d['baseline']['residual_orders']
e1 = d['lever1_echo']['same_hold'][1]['residual_orders']
L3 = d['lever3_plate']; ff = d['lever5_freefall']
items = [
    ('N5 cap 5049 -> 1e4 (more shots)', np.log10(d['baseline']['req_white'] / d['baseline']['req_white_stage18_criterion'])),
    ('hold 2.5 -> 5 s', np.log10([r for r in d['lever1_hold'] if r['t_hold'] == 5.0 and r['n_swap'] == 0][0]['req_white'] / d['baseline']['req_white'])),
    ('echo N=1: tilt (1/f)', base['tilt'] - e1['tilt']),
    ('echo N=1: isolated seismic (flat)', base['seismic_iso_flat'] - e1['seismic_iso_flat']),
    ('echo N=1: NLNM (microseism)', base['seismic_nlnm'] - e1['seismic_nlnm']),
    ('shorter shots (0.1-1 s), same wall-clock', -3.0),
    ('remove plate: requirement (D 11->273 um)', -L3['tighten_orders']),
    ('neutral dipole + non-metal plate vs E-field', base['efield_white'] - L3['neutral_dipole_orders']),
    ('stiffer trap', 0.0),
    ('drop tower, trapped (vs isolated ground)', base['seismic_iso_flat'] - ff['drop_tower_trapped']['orders']),
    ('ISS, trapped (vs isolated ground)', base['seismic_iso_flat'] - ff['iss_spec_trapped']['orders']),
    ('one-trap relative mode / DFS', 0.0),
]
y = np.arange(len(items))
vals = [v for _, v in items]
c.barh(y, vals, color=['tab:green' if v > 0.05 else ('tab:grey' if abs(v) <= 0.05 else 'tab:red') for v in vals])
c.set_yticks(y); c.set_yticklabels([k for k, _ in items], fontsize=7); c.invert_yaxis(); c.axvline(0, color='k')
for yi, v in zip(y, vals):
    c.text(v + (0.1 if v >= 0 else -0.1), yi, '%+.2f' % v, va='center', ha='left' if v >= 0 else 'right', fontsize=7)
c.text(-2.9, 5, 'infeasible', fontsize=7, va='center')
c.set_xlabel('orders loosened (+) or tightened (-)'); c.set_xlim(-6.5, 3.5)
c.set_title('(c) lever ledger (exact Stage 18 model, N5 <= 1e4)', fontsize=10)

dd = ax[1, 1]
cb = d['combined']
designs = [('ground baseline (same criterion)\nseismic-iso + tilt', cb['ground_baseline_same_criterion']['inertial']),
           ('ground baseline + plate surface noise', cb['ground_baseline_same_criterion']['with_surface']),
           ('ground: echo N=1, plate kept', cb['ground_plate_echo']['inertial']),
           ('ground: echo + neutral dipole + non-metal plate', cb['ground_plate_echo_neutral']['total']),
           ('ground: no plate (D=273 um), echo', cb['ground_noplate_echo']['inertial']),
           ('drop tower, trapped', ff['drop_tower_trapped']['orders']),
           ('space, trapped, LPF-grade readout', cb['space_trapped_lpf']['orders']),
           ('space, released, no plate (E-field white)', cb['space_free_floating_noplate']['orders']),
           ('space, released, no plate (E-field 1/f)', cb['space_free_floating_noplate']['orders_1overf'])]
y = np.arange(len(designs))
dd.barh(y, [v for _, v in designs], color=['grey', 'grey', 'tab:blue', 'tab:blue', 'tab:blue', 'tab:orange',
                                           'tab:purple', 'tab:purple', 'tab:purple'])
dd.set_yticks(y); dd.set_yticklabels([k for k, _ in designs], fontsize=7); dd.invert_yaxis()
for yi, (_, v) in zip(y, designs):
    dd.text(v + 0.1, yi, '%+.1f' % v, va='center', fontsize=7)
dd.set_xlabel('residual: log10(noise / requirement)'); dd.axvline(0, color='k')
dd.set_title('(d) combined: best ground +6.0, best space +2.0', fontsize=10)
fig.suptitle('Stage 21 - loosening the BMV acceleration-noise requirement: what survives', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(ROOT, 'figures', 'stage21_redesign.png'); fig.savefig(out, dpi=140); print(out)
