"""Stage 18 figure: common-mode vs differential noise for the trapped BMV point (reads results/stage18.json)."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import stage18_run as s18

d = json.load(open(os.path.join(ROOT, 'results', 'stage18.json')))
REQ = d['REQ']
fig, ax = plt.subplots(2, 2, figsize=(13.5, 10))

# (a) witness vs correlated noise
a = ax[0, 0]
w = d['witness']; sw = w['sigma_sweep']
for k, lab, st in (('uncorrelated', 'independent (c=0)', '-'), ('common', 'common-mode (c=+1)', '--'),
                   ('anti', 'anti-correlated (c=-1)', ':')):
    a.plot(sw['sigma'], sw[k], st, lw=2, label=lab)
a.axhline(w['lambda_target'], color='k', lw=0.8, ls='-.', label='Stage 16 budget (N5=5.0e3)')
a.axhline(0, color='r', lw=0.8)
a.set_xscale('log'); a.set_xlabel('which-path phase noise per mass  sigma (rad)')
a.set_ylabel('PPT  lambda_min  (<0 = entangled)')
a.set_title('(a) PPT witness has no built-in common-mode rejection\nallowed sigma^2: c=0 %.4f, c=+1 %.4f  ->  CMRR = %.3f'
            % (w['var_allowed']['uncorrelated'], w['var_allowed']['common'], w['intrinsic_cmrr']['common']), fontsize=10)
a.legend(fontsize=8, loc='upper left')

# (b) spectra
b = ax[0, 1]
f = np.geomspace(1e-3, 100, 400)
S, _ = s18.src_table()
show = ['seismic_nlnm', 'seismic_aligo_iso', 'tilt', 'rotation_diff', 'trap_current_persist', 'johnson_plate',
        'efield_white', 'newtonian_seismic']
for s in S:
    if s['key'] in show:
        b.loglog(f, [s['asd'](x) for x in f], ls='-' if s['kind'] == 'common' else '--', lw=1.6, label=s['label'])
b.loglog(f, (2 * np.pi * f) ** 2 * s18.LPF_READ, color='grey', lw=1, ls=':', label='LPF readout x(2 pi f)^2 (space)')
b.axhline(REQ, color='k', lw=2, label='requirement %.2f fm (exact, this stage)' % (REQ * 1e15))
b.axhline(0.899e-15, color='k', lw=1, ls='--', label='Stage 16 model 0.90 fm')
b.axhline(d['PINO_REQ'], color='m', lw=1, ls='--', label='Pino 2018 point 0.68 fm (Stage 19)')
W2 = np.array([s18.W2(x) for x in f])
b2 = b.twinx(); b2.loglog(f, W2, color='orange', alpha=0.35); b2.set_ylabel('|W(f)|^2 (s^2), split profile', color='orange'); b2.set_ylim(1e-6, 1e2)
b.set_xlabel('frequency (Hz)'); b.set_ylabel('acceleration ASD per mass along split axis (m s^-2 /rtHz)')
b.set_title('(b) per-mass noise along the split axis (solid common, dashed differential)', fontsize=10)
b.legend(fontsize=6.5, loc='lower left'); b.set_ylim(1e-18, 1e-3)

# (c) ledger in orders
c = ax[1, 0]
rows = [r for r in d['sources'] if r['key'] != 'efield_1overf']
labels = [r['label'] for r in rows]; vals = [r['orders_over_req'] for r in rows]
cols = ['tab:blue' if r['kind'] == 'common' else 'tab:red' for r in rows]
g = d['ggn_person']['rows'][0]
labels.append('one person moving 0.3 m at 2 m (absolute g, common)'); vals.append(g['orders_common']); cols.append('tab:blue')
for k, v in d['reference_floors'].items():
    labels.append('reference floor: ' + k); vals.append(v['orders_over_req']); cols.append('tab:green')
for k, v in d['scenarios'].items():
    labels.append('TOTAL ' + k); vals.append(v['orders_over_req']); cols.append('black')
y = np.arange(len(vals))
c.barh(y, vals, color=cols)
c.set_yticks(y); c.set_yticklabels(labels, fontsize=6.5); c.invert_yaxis()
c.axvline(0, color='k', lw=1)
c.set_xlabel('log10( white-equivalent ASD / %.2f fm )' % (REQ * 1e15))
c.set_title('(c) ledger: blue common, red differential, green best sensors, black totals', fontsize=10)
for yi, v in zip(y, vals):
    c.text(v + 0.1, yi, '%+.1f' % v, va='center', fontsize=6)

# (d) allowed region
dd = ax[1, 1]
reg = w['region']
dd.loglog(np.array(reg['sqrt_Sc']) * 1e15, np.maximum(np.array(reg['sqrt_Su_allowed']), 1e-20) * 1e15, lw=2,
          label='PPT witness (exact)')
x = np.geomspace(1e-2, 10, 100)
ideal = np.sqrt(np.clip(REQ ** 2 * 1e30 - 0 * x, 0, None)) * np.ones_like(x)
dd.loglog(x, ideal, 'k:', label='if common mode were free (ideal gradiometer)')
dd.axvline(REQ * 1e15, color='grey', ls='--')
dd.set_xlabel('common-mode ASD (fm s^-2 /rtHz)'); dd.set_ylabel('allowed independent ASD (fm s^-2 /rtHz)')
dd.set_xlim(1e-2, 10); dd.set_ylim(1e-2, 2)
dd.set_title('(d) allowed region: common + independent add in quadrature\n(needs common <= %.2f fm, not a CMRR-relaxed limit)'
             % (REQ * 1e15), fontsize=10)
dd.legend(fontsize=8)
fig.suptitle('Stage 18 - trapped BMV point (diamond 5.6e-15 kg, dx 3.85 um, D 11 um, 3 s): common-mode noise is not rejected',
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(ROOT, 'figures', 'stage18_noise.png')
fig.savefig(out, dpi=140)
print(out)
