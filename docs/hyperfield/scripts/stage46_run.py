#!/usr/bin/env python3
"""Stage 46 - release/clock gate from 57Fe nuclear forward scattering (NFS) locked to the synchrotron pulse.

MODEL (ASSUMED), excitation at the prompt (t = 0, non-resonant photons visible in the same spectrum):
 clock gate (a): the whole coherent response is delayed by t_g -> every beat minimum moves by +t_g relative to a calculation with t = 0 at excitation
 release gate (b): no emission for t < t_g, then the normal response with free normalisation -> zero nuclear counts in [prompt tail, t_g]
Standard QM/NFS: coherent forward emission starts at t = 0+ with the *fastest* (speed-up) decay; no delayed onset; minima set by hyperfine + thickness beats.
Data: Marx-Glowna et al., J. Synchrotron Rad. 28, 120 (2021), Fig. 6 (PMC7842219, open access): 9.3 um 57Fe foil, polarimeter, fit with CONUSS from 1.5 ns.
DIGITIZED from the embedded 1385x724 jpeg -> data/stage46_marxglowna2021_fig6_digitized.csv (cols t, data median/low/high per pixel column, calculated curve).
"""
import json, os, csv
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
S = {}
def src(key, ref, url, loc, quote, **vals): S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)
src('marxglowna21', 'Marx-Glowna, Uschmann, Schulze, Marschner, Wille, Schlage, Stoehlker, Roehlsberger, Paulus, "Advanced X-ray polarimeter design for nuclear resonant scattering", J. Synchrotron Rad. 28, 120 (2021), PMC7842219',
    'https://pmc-oa-opendata.s3.amazonaws.com/PMC7842219.1/PMC7842219.1.pdf', 'Fig. 6 and text p.4',
    '"Time spectrum of a 9.3 um-thick foil of 57Fe (95% enrichment) magnetized at 45 deg to the polarization direction. The red line is a fit to the data via the program package CONUSS" ; '
    '"The use of the polarimeter makes it possible to fit the time spectrum already 1.5 ns after excitation" ; "The total data collection time was 69 min" ; axis "time after excitation [ns]", "photons/channel"; arrow "non-resonantly scattered photons" at 0')
src('esrf_tau', 'ESRF Spotlight 47, Nuclear forward scattering in pulsed magnetic fields', 'https://www.esrf.fr/news/spotlight/spotlight47/spotlight47nfs', 'para 1', '"(tau 57Fe = 141 ns)"', tau=141.0)
src('digitization', 'DIGITIZED (this stage)', 'data/stage46_marxglowna2021_fig6_digitized.csv', 'Fig. 6 raster', 'x ticks 0/50/100/150 ns at px 300/638/976/1313 (6.757 px/ns, 0.148 ns/px); y decade ticks (right frame) 1e4/1e3/1e2/1 at px 72/189/307/541.5 (117.4 px/decade); black = data, red = calculated')
d = np.genfromtxt(os.path.join(ROOT, 'data', 'stage46_marxglowna2021_fig6_digitized.csv'), delimiter=',', skip_header=1)
t, med, low, high, calc = d.T
tau = S['esrf_tau']['values']['tau']; out = {}
# ---- clock gate: beat minima, data vs calculation ----
cents = [8.3, 19.9, 31, 51, 60, 70, 77, 81.5, 99.5, 106, 116, 142]; mins = []
for c0 in cents:
    w = (t > c0 - 2.5) & (t < c0 + 2.5)
    mins.append((c0, float(t[w][np.nanargmin(calc[w])]), float(t[w][np.nanargmin(med[w])])))
mins = np.array(mins); off = mins[:, 2] - mins[:, 1]
m, s = off.mean(), off.std(ddof=1) / np.sqrt(len(off))
PROMPT_SYS = 0.30   # ns: prompt maximum digitized at -0.15 / +0.30 ns depending on envelope -> excitation instant known to +-0.3 ns (ASSUMED)
out['clock'] = dict(offsets_ns=off.tolist(), mean_ns=m, sem_ns=s, scatter_ns=float(off.std(ddof=1)), tg95_stat_ns=[m - 1.96 * s, m + 1.96 * s],
                    tg95_ns=[m - 1.96 * s - PROMPT_SYS, m + 1.96 * s + PROMPT_SYS], note='ASSUMED: the CONUSS calculation has t = 0 at excitation and no fitted time offset')
out['clock']['tg_over_tau'] = max(abs(x) for x in out['clock']['tg95_ns']) / tau
# ---- release gate: earliest clearly nuclear counts ----
w = (t > 1.3) & (t < 3.0)
dip = float(np.nanmin(med[(t > 1.0) & (t < 1.6)]))           # prompt tail + background ceiling at the dip
sig = [(float(a_), float(b_), float((b_ - dip) / np.sqrt(dip))) for a_, b_ in zip(t[w], med[w])]   # Poisson vs dip ceiling, per-column median
first = next(sig[i] for i in range(len(sig)) if all(x[2] > 5 for x in sig[i:]))   # >5 sigma here and at every later column to 3 ns
out['release'] = dict(dip_ceiling=dip, per_channel=sig, tg95_ns=first[0], tg_over_tau=first[0] / tau,
                      logic='a release gate t_g removes all nuclear counts before t_g; the data exceed the prompt-tail/background ceiling by >5 sigma in every column from t = %.2f ns to 3 ns, so t_g < that time' % first[0])
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
ax[0].semilogy(t, med, 'k-', lw=0.6, label='data (DIGITIZED, Fig. 6)'); ax[0].semilogy(t, calc, 'r-', lw=0.8, label='CONUSS calculation (paper)')
ax[0].set_xlabel('time after excitation (ns)'); ax[0].set_ylabel('photons / channel'); ax[0].legend(fontsize=7); ax[0].set_title('⁵⁷Fe NFS, 9.3 µm foil (Marx-Glowna 2021)', fontsize=9)
z = (t > -2) & (t < 7)
ax[1].semilogy(t[z], med[z], 'k.-', lw=0.6, ms=3, label='data'); ax[1].semilogy(t[z], calc[z], 'r-', label='calculation, no gate')
tg = out['release']['tg95_ns']; ax[1].axvspan(-2, 0.9, color='grey', alpha=0.2, label='prompt (non-resonant)')
ax[1].axvline(tg, color='b', ls='--', label=f'release-gate limit {tg:.2f} ns'); ax[1].axhline(dip, color='g', ls=':', label=f'tail+bg ceiling {dip:.0f}')
ax[1].set_xlabel('time (ns)'); ax[1].legend(fontsize=7); ax[1].set_title('early window: nuclear counts present by ~2 ns', fontsize=9)
ax[2].errorbar(mins[:, 1], off, yerr=out['clock']['scatter_ns'], fmt='ko', ms=4, capsize=2)
ax[2].axhline(0, color='k', lw=0.5); ax[2].axhspan(*out['clock']['tg95_ns'], color='r', alpha=0.15, label=f"clock gate 95%: [{out['clock']['tg95_ns'][0]:.2f}, {out['clock']['tg95_ns'][1]:.2f}] ns")
ax[2].set_xlabel('beat minimum (ns)'); ax[2].set_ylabel('data − calculation (ns)'); ax[2].legend(fontsize=7); ax[2].set_title('beat minima do not shift', fontsize=9)
fig.suptitle(f"Stage 46 — ⁵⁷Fe NFS gate: clock shift |t_g| < {max(abs(x) for x in out['clock']['tg95_ns']):.2f} ns ({out['clock']['tg_over_tau']:.1e} τ); release gate t_g < {tg:.1f} ns ({out['release']['tg_over_tau']:.1e} τ)", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, 'figures', 'stage46_nfs_gate.png'), dpi=110)
with open(os.path.join(ROOT, 'results', 'stage46_sources.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for kk, v in S.items(): wr.writerow([kk, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
json.dump(dict(stage=46, results=out, minima=mins.tolist(), sources=S), open(os.path.join(ROOT, 'results', 'stage46.json'), 'w'), indent=1, default=float)
print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != 'per_channel'} for k, v in out.items()}, indent=1, default=float)); print(sig)
