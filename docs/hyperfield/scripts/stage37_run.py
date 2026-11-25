#!/usr/bin/env python3
"""Stage 37 - deep-dive on the Vinante 2017 8174 Hz cantilever excess (resonant-release coincidence candidate).
What is known since, whether anything else tested near 8174 Hz, and what a frequency-selective release would predict elsewhere.
Every input fetched in THIS session; quotes in results/stage37_sources.csv.
"""
import json, os, csv
import numpy as np
from math import pi, sqrt, log
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
S = {}
def src(key, ref, url, loc, quote, **vals): S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)
src('vinante17', 'Vinante, Mezzena, Falferi, Carlesso & Bassi, PRL 119, 110401 (2017), arXiv:1611.09776', 'https://arxiv.org/pdf/1611.09776', 'abstract; p.3-5; Table I; Suppl. "Additional information on fitting"',
    '"reveal a nonthermal force noise of unknown origin. This excess noise is compatible with the CSL heating predicted by Adler" ; "The resonant frequency of the fundamental mode of the cantilever is f0 = 8174.01 Hz" ; '
    'Table I: "Off 116+-1 1.87+-0.16+-0.1(sys)" ; "Off 347+-3 2.12+-0.20+-0.1(sys)" ; "On 114+-2 2.58+-0.20+-0.1(sys)" (SF0 in aN^2/Hz) ; '
    '"For the standard choice rC = 10^-7 m this would imply lambda = 10^-7.7" ; "fitted ... on the fixed interval 8100 - 8240 Hz" ; "the cantilever noise available over the whole noise bandwidth (roughly 20 Hz)" ; '
    '"the residual force noise ... is only slightly increased with respect to the measurements with the pulse tube off"',
    f0=8174.01, SF_off=(1.87, 0.16, 0.1), SF_off_hi=(2.12, 0.20, 0.1), SF_on=(2.58, 0.20, 0.1), log_lam=-7.7, fit=(8100, 8240), bw=20.0)
src('vinante20', 'Vinante, Carlesso, Bassi, Chiasera, Varas, Falferi, Margesin, Mezzena & Ulbricht, PRL 125, 100404 (2020), arXiv:2002.09782', 'https://arxiv.org/pdf/2002.09782', 'p.2-3; Suppl.',
    '"The resonance frequency is f0 = 3532.7 Hz, while it was measured as f00 = 8174 Hz before attaching the multilayer test mass" ; "SF0 = (-1.51 +- 1.77) x 10^-36 N^2/Hz" ; "SF0,CSL <= 2.07x10^-36 N^2/Hz at the 95% confidence level" ; '
    '"excess noise observed in a previous related experiment [31] is incompatible with a CSL effect for rC = 10^-7 m. On the other hand, they do not provide substantial new insight into the origin of that excess noise. Indeed, the absolute value of the excess force noise in the previous experiment, featuring the same cantilever, magnet, and SQUID, is compatible with the error bar of the new experiment" ; '
    '"[3515, 3550] Hz is chosen"', f0=3532.7, SF=(-1.51e-36, 1.77e-36), SF_ul=2.07e-36)
src('review22', 'Carlesso et al., Nature Phys. 18, 243 (2022), arXiv:2203.04231', 'https://arxiv.org/pdf/2203.04231', 'Sec. 2.2',
    '"tailored for testing the CSL model at rC = 10^-7 m, yielding an upper bound lambda < 2.0 x 10^-10 s^-1, which completely covers the values suggested by Adler"', lam20=2.0e-10)
src('adler_vinante18', 'Adler & Vinante, PRA 97, 052119 (2018), arXiv:1801.06857', 'https://arxiv.org/pdf/1801.06857', 'p.3',
    '"This means that the effective lambda for bulk heating is at most 3.1 x 10^-11 s^-1, ruling out a collapse model white noise interpretation of the excess noise reported in [3]"')
src('search', 'WebSearch "Vinante 2017 cantilever excess noise ... follow-up 2021 2022 2023" (this session)', 'https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.100404', 'search synthesis',
    'no later explanation found; 2020 result is the latest direct test located (search is not exhaustive)')
src('ligo_O3a', 'LIGO-G1902347-v1 L1 O3a displacement ASD (Stage 34 data file, re-read)', 'https://dcc.ligo.org/public/0165/G1902347/001/O3-L1-C01_CLEAN_SUB60HZ-1240573680.0_sensitivity_displacement_asd.txt',
    'last rows', '"4.9900000e+03 2.0452922e-19" (file ends at 5 kHz; uncalibrated above)', x5k=2.0e-19)
src('ligo_geom', 'Carlesso et al. PRD 94 124036 (2016) and ccb18 (Stage 34 sources): m = 40 kg, silica 2200 kg/m^3, R = 0.17 m', 'https://arxiv.org/pdf/1606.04581', 'Sec. III A',
    'as quoted in results/stage34_sources.csv', m=40.0, rho=2200.0, R=0.17)

v17, v20 = S['vinante17']['values'], S['vinante20']['values']
lam17 = 10 ** v17['log_lam']; lam20 = S['review22']['values']['lam20']
ratio_needed = lam20 / lam17                   # spectral weight at 3532.7 Hz relative to 8174 Hz allowed for a mass-coupled release (DERIVED)
df = v17['f0'] - v20['f0']
hwhm_lor_max = df * sqrt(ratio_needed / (1 - ratio_needed))
sig_gauss_max = df / sqrt(2 * log(1 / ratio_needed))
hwhm_min = v17['bw'] / 2                       # ASSUMED: excess flat across the ~20 Hz cantilever noise band (fits pass chi^2 with a constant)
# pulse-tube on vs off (environment sensitivity)
a, b = v17['SF_on'], v17['SF_off']
d_pt = a[0] - b[0]; s_pt = sqrt(a[1] ** 2 + b[1] ** 2); z_pt = d_pt / s_pt
# 2020 absolute test of an instrument-like (not mass-coupled) noise of the 2017 size at 3532.7 Hz
SF17 = b[0] * 1e-36
z20 = (SF17 - v20['SF'][0]) / v20['SF'][1]
# LIGO at 8174 Hz: prediction for lambda = 10^-7.7 (DERIVED; Stage 34 K)
HBAR, U = 1.054571817e-34, 1.66053906892e-27
g = S['ligo_geom']['values']; K = 16 * pi * HBAR ** 2 * g['rho'] ** 2 * pi * g['R'] ** 2 * 1e-14 / U ** 2
w0 = 2 * pi * v17['f0']; x_pred = sqrt(4 * K * lam17) / (g['m'] * w0 ** 2)
gap_ligo = S['ligo_O3a']['values']['x5k'] / x_pred
res = dict(lam17=lam17, lam20=lam20, ratio_needed=ratio_needed, df=df, hwhm_lorentz_max_Hz=hwhm_lor_max, sigma_gauss_max_Hz=sig_gauss_max, hwhm_min_Hz=hwhm_min,
           pulse_tube=dict(delta_aN2=d_pt, sigma=s_pt, z=z_pt), z_2020_vs_2017_absolute=z20, ligo=dict(x_pred_m_rtHz=x_pred, ligo_5kHz=S['ligo_O3a']['values']['x5k'], orders_short=np.log10(gap_ligo)),
           harmonics=dict(f_half=v17['f0'] / 2, f_2=2 * v17['f0'], f_3=3 * v17['f0']))
# figure: allowed line shapes vs existing points
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
s36 = json.load(open(os.path.join(ROOT, 'results', 'stage36.json')))
fig, ax = plt.subplots(figsize=(10, 6))
f = np.logspace(np.log10(500), np.log10(3e4), 4000)
for hw, ls in ((hwhm_min, ':'), (hwhm_lor_max, '-')):
    ax.plot(f, lam17 * hw ** 2 / ((f - v17['f0']) ** 2 + hw ** 2), ls, color='crimson', label=f'Lorentzian release line at 8174 Hz, HWHM {hw:.0f} Hz (λ = 10$^{{-7.7}}$ on line)')
ax.plot(f, lam17 * np.exp(-(f - v17['f0']) ** 2 / (2 * sig_gauss_max ** 2)), '--', color='darkorange', label=f'Gaussian line, σ = {sig_gauss_max:.0f} Hz (widest allowed)')
lg = np.loadtxt(os.path.join(ROOT, 'results', 'stage34_ligo_curve.csv'), delimiter=',', skiprows=1)
ax.plot(lg[:, 0], lg[:, 1], color='navy', lw=1, label='LIGO O1/O3 bound (Stage 34)')
for bnd in s36['bands']:
    fc = sqrt(bnd['fmin'] * bnd['fmax'])
    if 500 < fc < 3e4 and bnd['fmax'] / bnd['fmin'] < 1.3:
        ax.plot(fc, bnd['lam'], 'v', ms=9, color='purple'); ax.text(fc * 1.03, bnd['lam'] * 1.5, bnd['name'].split('(')[0].strip(), fontsize=7.5)
ax.axvspan(*v17['fit'], color='crimson', alpha=0.15, label='2017 fit window 8100–8240 Hz')
for h_, nm in ((v17['f0'] / 2, 'f₀/2'), (2 * v17['f0'], '2f₀'), (3 * v17['f0'], '3f₀')):
    ax.axvline(h_, color='0.5', ls='-.', lw=0.8); ax.text(h_, 1e-12, nm, fontsize=8, ha='center')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(1e-13, 1e5); ax.set_xlim(500, 3e4)
ax.set_xlabel('frequency (Hz)'); ax.set_ylabel('λ(f) (s$^{-1}$, CSL-equiv., r$_C$=100 nm)')
ax.set_title('Stage 37 — what a resonant release at 8174 Hz must look like to explain Vinante 2017 and survive Vinante 2020', fontsize=10)
ax.legend(fontsize=7.5, loc='upper left'); ax.grid(alpha=0.3, which='both')
plt.tight_layout(); plt.savefig(os.path.join(ROOT, 'figures', 'stage37_8174Hz.png'), dpi=150)
with open(os.path.join(ROOT, 'results', 'stage37_sources.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for kk, v in S.items(): wr.writerow([kk, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
json.dump(dict(stage=37, results=res, sources=S), open(os.path.join(ROOT, 'results', 'stage37.json'), 'w'), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
