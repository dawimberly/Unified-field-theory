#!/usr/bin/env python3
"""Stage 40 - extend the Stage 36 membrane result into the 8 kHz - 4.6 MHz gap with a second public membrane dataset:
Yuan, Cohen & Steele, APL 107, 263501 (2015), arXiv:1510.07468 (SiN membranes at 14 mK, 121-242 kHz).
Same parametrization (CSL-equivalent lambda per nucleon, mass-proportional, r_C = 1e-7 m) and the Stage 36 thin-plate formula
  Gamma/lambda = 4 pi hbar rho r_C^2 (1 - e^{-h^2/4 r_C^2}) / (m0^2 h w).
Bound: CSL heating <= measured total heating n_mode * gamma_m = kB T_mode/(hbar Q) (all attributed; T_mode is the MEASURED mode temperature).
Every input fetched this session (results/stage40_sources.csv).
"""
import json, os, csv
from math import pi, exp, log10, sqrt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
S = {}
def src(key, ref, url, loc, quote, **vals): S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)
src('nist', 'NIST CODATA 2022 Value?hbar, ukg, k (fetched this session)', 'https://physics.nist.gov/cgi-bin/cuu/Value?hbar', 'Numerical value fields',
    '"1.054 571 817... x 10 -34 J s" ; "1.660 539 068 92 x 10 -27 kg" ; "1.380 649 x 10 -23 J K -1"', hbar=1.054571817e-34, u=1.66053906892e-27, kB=1.380649e-23)
src('yuan15', 'Yuan, Cohen & Steele, "Silicon nitride membrane resonators at millikelvin temperatures with quality factors exceeding 10^8", APL 107, 263501 (2015), arXiv:1510.07468',
    'https://arxiv.org/pdf/1510.07468', 'p.2 device text; Table I; p.4-5 Q and mode temperature',
    '"The membranes are 50 nm thick" ; "We deposit a metal layer of 20 nm of Al on the membrane" ; "taking an effective density of 3.0 g/cm3" ; Table I "wm/2pi (kHz) 242 [Device I (1,1)] 121 193 192 [Device II (1,1),(1,2),(2,1)]" ; '
    '"reaching a record of Q = 1.27x10^8 for a fundamental mode at 14 mK" ; "Quality factors of all three modes are above 5 x 10^6" ; "Mode (2,1) of Device II has a much lower Q = 1.2 x 10^5 that is independent of temperature" ; '
    '"the mode temperature saturates at 210 mK for Device I and 180 mK for Device II" ; "The high mode temperature in the experiment is likely related to mechanical vibration noise in the setup"',
    h=50e-9, h_al=20e-9, rho_eff=3000.0)
src('tsaturyan17', 'Tsaturyan et al., Nature Nanotech. 12, 776 (2017), arXiv:1608.00937 (fetched for Stage 36 this session)', 'https://arxiv.org/pdf/1608.00937', 'p.3', '"(rho = 3200 kg/m3 is the density of SiN)"', rho=3200.0)
N = S['nist']['values']; HBAR, U, KB = N['hbar'], N['u'], N['kB']; M0, RC = U, 1e-7
y = S['yuan15']['values']; rho_sin = S['tsaturyan17']['values']['rho']
def G_plate(rho, h, f): return 4 * pi * HBAR * rho * RC ** 2 * (1 - exp(-h ** 2 / (4 * RC ** 2))) / (M0 ** 2 * h * 2 * pi * f)
modes = [('Device I (1,1)', 242e3, 1.27e8, 0.210, 'Q CITED (14 mK)'),
         ('Device II (1,1)', 121e3, 5e6, 0.180, 'Q >= 5e6 CITED (lower bound used -> conservative)'),
         ('Device II (1,2)', 193e3, 5e6, 0.180, 'Q >= 5e6 CITED (lower bound used -> conservative)'),
         ('Device II (2,1)', 192e3, 1.2e5, 0.180, 'Q = 1.2e5 CITED; mode T of Device II ASSUMED to apply')]
B = []; rows = []
for nm, f, Q, Tm, note in modes:
    heat = KB * Tm / (HBAR * Q)
    G = min(G_plate(rho_sin, y['h'], f), G_plate(y['rho_eff'], y['h'] + y['h_al'], f))   # conservative: smaller prediction
    lam = heat / G
    rows.append(dict(mode=nm, f=f, Q=Q, T_mode=Tm, heat=heat, G_per_lam=G, lam=lam))
    B.append(dict(name=f'SiN membrane {f/1e3:.0f} kHz, {nm} (Yuan 2015)', fmin=f * 0.999, fmax=f * 1.001, lam=lam, label='DERIVED (Stage 36 thin-plate formula; heating kB T_mode/(hbar Q), all attributed)',
                  sources=['yuan15', 'tsaturyan17', 'nist'], kind='membrane / resonator heating', mass_prop=True, note=note, curve=None, stage='40'))
s36 = json.load(open(os.path.join(ROOT, 'results', 'stage36.json')))
Ball = [dict(b) for b in s36['bands']] + B
def union(iv):
    iv = sorted(iv); out = []
    for lo, hi in iv:
        if out and lo <= out[-1][1]: out[-1][1] = max(out[-1][1], hi)
        else: out.append([lo, hi])
    return out
U40 = union([(b['fmin'], b['fmax']) for b in Ball if b['mass_prop']]); SPAN = (2e-05, 9.188359119922688e+20)
open_w = []; prev = SPAN[0]
for lo, hi in U40:
    if lo > prev and log10(lo / prev) > 1e-3: open_w.append((prev, lo, log10(lo / prev)))
    prev = max(prev, hi)
if SPAN[1] > prev: open_w.append((prev, SPAN[1], log10(SPAN[1] / prev)))
cov = dict(s36['coverage']); cov['open_windows'] = open_w
best = min(rows, key=lambda r: r['lam'])
fig_title = f"Stage 40 — Yuan 2015 mK SiN membranes fill 121–242 kHz points: best {best['lam']:.1e} s$^{{-1}}$ at {best['f']/1e3:.0f} kHz (plus Seis 2.9e-9 at 1.49 MHz)"
with open(os.path.join(ROOT, 'results', 'stage40_sources.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for kk, v in S.items(): wr.writerow([kk, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
with open(os.path.join(ROOT, 'results', 'stage40_bounds.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['name', 'kind', 'fmin_Hz', 'fmax_Hz', 'lambda_max_per_s_rC_1e-7m', 'label', 'stage', 'sources', 'note'])
    for b in sorted(Ball, key=lambda x: x['fmin']):
        wr.writerow([b['name'], b['kind'], f"{b['fmin']:.4g}", f"{b['fmax']:.4g}", f"{b['lam']:.3g}", b['label'], b['stage'], ';'.join(b['sources']), b['note']])
json.dump(dict(stage=40, fig_title=fig_title, results=rows, bands=Ball, coverage=cov, sources=S), open(os.path.join(ROOT, 'results', 'stage40.json'), 'w'), indent=1, default=float)
for r in rows: print({k: (f'{v:.3g}' if isinstance(v, float) else v) for k, v in r.items()})
