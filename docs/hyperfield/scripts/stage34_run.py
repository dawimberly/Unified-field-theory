#!/usr/bin/env python3
"""Stage 34 - filling the two cheapest Stage-33 open windows for resonant release (do-si-do element 3).

Same parametrization as Stages 31/33: extra release rate per nucleon expressed as the CSL-equivalent, mass-proportional
collapse rate lambda at r_C = 1e-7 m; an experiment sampling the noise at frequency f bounds lambda(f).
(1) LIGO: public calibrated DARM displacement spectra (O3a L1, O3b H1, LIGO DCC) and the GWTC-1 GW150914 PSDs (O1),
    converted to per-test-mass force noise and to lambda(f) over 10 Hz - 5 kHz; validated against the published
    30-35 Hz LIGO CSL bound.
(2) Germanium 100 keV - 1 MeV: no fetchable count-rate spectrum with detection efficiency -> dropped (honest negative).
Every input was fetched in THIS session; url/location/quote in results/stage34_sources.csv.
Labels: CITED (number from fetched source), DERIVED (computed here), ASSUMED (stated choice).
"""
import json, os, csv, time
import numpy as np
from math import pi, log10

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
DATA = os.path.join(ROOT, 'data', 'stage34')
t0 = time.time()
S = {}
def src(key, ref, url, loc, quote, **vals):
    S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)

src('nist', 'NIST CODATA pages Value?hbar, Value?ukg (fetched this session)', 'https://physics.nist.gov/cgi-bin/cuu/Value?hbar', 'Numerical value fields',
    '"Numerical value 1.054 571 817... x 10 -34 J s" ; "Numerical value 1.660 539 068 92 x 10 -27 kg"', hbar=1.054571817e-34, u=1.66053906892e-27)
src('cbfv16', 'Carlesso, Bassi, Falferi & Vinante, PRD 94, 124036 (2016), arXiv:1606.04581', 'https://arxiv.org/pdf/1606.04581', 'Sec. III A, Eq. (15)',
    '"the effective force noise reaches a minimum SF ~ 95 fN Hz^-1/2 at w/2pi ~ 30-35 Hz. We have used the numerical values m = 40 kg for the test mass and a = 4 km for the arm length" ; '
    '"Shh(w) = Sxx(w)/a^2" ; "Sxx(w) = (4/m^2) SFF(w)/((w0^2-w^2)^2 + (w w0/Q)^2)" ; "resonant frequency w0/2pi below 1 Hz ... the test masses can be considered to a good approximation in the free-mass limit"',
    m=40.0, a=4000.0, SF_fN=95.0, f_ref=(30.0, 35.0))
src('ccb18', 'Carlesso, Ferialdi & Bassi, EPJD 72, 159 (2018), arXiv:1805.10100', 'https://arxiv.org/pdf/1805.10100', 'Sec. IV A',
    '"a cylindrical silica mirror (density 2200 kg/m3, radius 17 cm and length 20 cm) are suspended and oscillate at a frequency below 1 Hz, while its noise is monitored in the 10 - 10^3 Hz band" ; Eq. (13) S_CSL',
    rho=2200.0, R=0.17, L=0.20)
src('review22', 'Carlesso et al., Nature Phys. 18, 243 (2022), arXiv:2203.04231', 'https://arxiv.org/pdf/2203.04231', 'Sec. 2.3',
    '"LIGO sets lambda < 1.0 x 10^-5 s^-1, AURIGA gives lambda < 4.6x10^-2 s^-1, while LISA Pathfinder provides lambda < 3.8x10^-9 s^-1"', ligo=1.0e-5)
src('L1_O3a', 'LIGO-G1902347-v1, "L1 Calibrated Sensitivity Spectra Apr 29 2019 (Representative best of O3a -- C01_CLEAN_SUB60HZ)", E. Goetz',
    'https://dcc.ligo.org/public/0165/G1902347/001/O3-L1-C01_CLEAN_SUB60HZ-1240573680.0_sensitivity_displacement_asd.txt', 'DCC abstract; file columns f [Hz], ASD',
    '"These are calibrated strain and displacement spectra from a representative best of O3a" ; "We expect the uncertainty in calibration at this time 5% and 3 deg, valid between 10 and 2000 [Hz] (the uncertainty is larger and unquantified outside this band)" ; file rows 1.0000000e+01 ... 5.0000000e+03',
    file='O3-L1-C01_CLEAN_SUB60HZ-1240573680.0_sensitivity_displacement_asd.txt', cal_band=(10.0, 2000.0))
src('H1_O3b', 'LIGO-G2100674-v1, H1 Calibrated Sensitivity Spectra Jan 04 2020 (Representative best of O3b -- C01_CLEAN_SUB60HZ)',
    'https://dcc.ligo.org/public/0174/G2100674/001/O3-H1-C01_CLEAN_SUB60HZ-1262197260.0_sensitivity_displacement_asd.txt', 'DCC abstract; file columns f [Hz], ASD',
    '"These are calibrated strain and displacement spectra from a representative best of O3b, taken on Jan 04 2020" ; "We expect the uncertainty in calibration at this time 5% and 3 deg, valid between 10 and 2000 [Hz] (th[e uncertainty is larger and unquantified outside this band])"',
    file='O3-H1-C01_CLEAN_SUB60HZ-1262197260.0_sensitivity_displacement_asd.txt', cal_band=(10.0, 2000.0))
src('gwtc1', 'LIGO-P1900011-v1, "Power Spectral Densities (PSD) release for GWTC-1", GWTC1_GW150914_PSDs.dat + README.txt',
    'https://dcc.ligo.org/public/0158/P1900011/001/GWTC1_GW150914_PSDs.dat', 'README; file header',
    '"# Freq (Hz) LIGO_Hanford_PSD (1/Hz) LIGO_Livingston_PSD (1/Hz)" ; "the LIGO and Virgo data are only properly calibrated in the frequency range between 10 Hz and 5 kHz" ; "only data in the frequency range from 20 Hz ... and 1 kHz ... was used. Thus the PSDs in this data release only have values in that range"',
    file='GWTC1_GW150914_PSDs.dat')
src('ge_search', 'Search for tabulated low-background Ge spectra 100 keV-1 MeV (WebSearch + Donadi 2021 arXiv:2107.11237)', 'https://arxiv.org/pdf/2107.11237', 'Sec. III B; Fig. 2',
    '"Since, due to self-absorption, the efficiency function strongly depends on the geometry of the apparatus, the efficiency is obtained by means of a MC simulation" ; "Measured X-rays spectrum ... in the energy range DE = (1000/3800) keV" ; GERDA/CDEX spectra found only as plots (no tables) and only with ROI background indices near 2039 keV',
    note='no fetchable count-rate table with efficiency for 100-1000 keV')

HBAR, U = S['nist']['values']['hbar'], S['nist']['values']['u']
M0, RC = U, 1e-7                     # ASSUMED m0 = 1 u (as Stage 33); r_C convention
c = S['cbfv16']['values']; g = S['ccb18']['values']
m, rho, A = c['m'], g['rho'], pi * g['R'] ** 2
# CSL per-mass one-sided force PSD (DERIVED, Stage 33 slab formula, faces of area pi R^2 normal to the beam, L >> r_C):
#   S_F1 = 2 * 8 pi hbar^2 lam rho^2 A rc^2 / m0^2
K = 16 * pi * HBAR ** 2 * rho ** 2 * A * RC ** 2 / M0 ** 2
lam_from_SF = lambda SF1: SF1 / K
# DARM displacement -> per-mass force (DERIVED, free-mass limit, 4 independent test masses; = Carlesso Eq. 15 convention):
#   S_x = 4 S_F1/(m^2 w^4)  ->  S_F1 = m^2 w^4 S_x / 4
SF_from_Sx = lambda f, Sx: m ** 2 * (2 * pi * f) ** 4 * Sx / 4

# --- validation 1: the published 95 fN/rtHz reproduces the review's 1.0e-5
lam_95 = lam_from_SF((c['SF_fN'] * 1e-15) ** 2)
# --- validation 2: the O1 GW150914 PSD at 30-35 Hz (strain x a = displacement)
d = np.loadtxt(os.path.join(DATA, S['gwtc1']['values']['file']))
o1 = {}
for j, det in [(1, 'H1'), (2, 'L1')]:
    f = d[:, 0]; Sx = d[:, j] * c['a'] ** 2
    lam = lam_from_SF(SF_from_Sx(f, Sx))
    mref = (f >= 30) & (f <= 35)
    o1[det] = dict(f=f, lam=lam, lam_ref_min=float(lam[mref].min()), lam_ref_median=float(np.median(lam[mref])),
                   SF_ref_min_fN=float(np.sqrt(SF_from_Sx(f[mref], Sx[mref]).min()) * 1e15))
valid = dict(lam_from_95fN=lam_95, review=S['review22']['values']['ligo'], ratio_95=lam_95 / S['review22']['values']['ligo'],
             O1_H1_min_30_35=o1['H1']['lam_ref_min'], O1_L1_min_30_35=o1['L1']['lam_ref_min'],
             O1_H1_SFmin_fN=o1['H1']['SF_ref_min_fN'], O1_L1_SFmin_fN=o1['L1']['SF_ref_min_fN'],
             ratio_O1_H1_to_review=o1['H1']['lam_ref_min'] / S['review22']['values']['ligo'])

# --- O3 curves
o3 = {}
for key in ('L1_O3a', 'H1_O3b'):
    dd = np.loadtxt(os.path.join(DATA, S[key]['values']['file']))
    f, asd = dd[:, 0], dd[:, 1]
    o3[key] = dict(f=f, lam=lam_from_SF(SF_from_Sx(f, asd ** 2)))
# envelope on a common grid (min over independent detectors/runs; each bound alone is valid)
fg = o3['L1_O3a']['f']
env = np.minimum(o3['L1_O3a']['lam'], np.interp(fg, o3['H1_O3b']['f'], o3['H1_O3b']['lam']))
for det in ('H1', 'L1'):
    inside = (fg >= o1[det]['f'][0]) & (fg <= o1[det]['f'][-1])
    env[inside] = np.minimum(env[inside], np.interp(fg[inside], o1[det]['f'], o1[det]['lam']))
def stat(lo, hi):
    mm = (fg >= lo) & (fg <= hi)
    return dict(min=float(env[mm].min()), f_at_min=float(fg[mm][np.argmin(env[mm])]), median=float(np.median(env[mm])), max=float(env[mm].max()))
bands_stats = {k: stat(*v) for k, v in dict(b10_30=(10, 30), b30_35=(30, 35), b35_1k=(35, 1000), b1k_2k=(1000, 2000), b2k_3k=(2000, 3000), b3k_5k=(3000, 5000)).items()}
pts = {str(q): float(env[np.argmin(abs(fg - q))]) for q in (10, 20, 50, 100, 300, 1000, 2000, 3000, 5000)}
# lowest-noise point overall
i_best = int(np.argmin(env))
# down-sampled curve for CSV / figure (log grid, take the MAX inside each log bin = conservative display)
edges = np.logspace(1, np.log10(5000), 301)
cur_f, cur_lam_max, cur_lam_min = [], [], []
for lo, hi in zip(edges[:-1], edges[1:]):
    mm = (fg >= lo) & (fg < hi)
    if mm.any():
        cur_f.append(float(np.sqrt(lo * hi))); cur_lam_max.append(float(env[mm].max())); cur_lam_min.append(float(env[mm].min()))
os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
with open(os.path.join(ROOT, 'results', 'stage34_ligo_curve.csv'), 'w') as fh:
    fh.write('f_Hz,lambda_bin_max_per_s,lambda_bin_min_per_s\n')
    for a_, b_, c_ in zip(cur_f, cur_lam_max, cur_lam_min): fh.write(f'{a_:.5e},{b_:.4e},{c_:.4e}\n')

# --- merge with Stage 33 bands and recompute coverage
s33 = json.load(open(os.path.join(ROOT, 'results', 'stage33.json')))
B = [dict(b) for b in s33['bands']]
for b in B: b['curve'] = None
B.append(dict(name='LIGO O3a L1 / O3b H1 / O1 GW150914 force noise (DERIVED)', fmin=10.0, fmax=5000.0, lam=float(env.min()),
              label='DERIVED (DCC displacement ASD -> per-mass force -> lambda); free-mass + rigid-body ASSUMED', sources=['L1_O3a', 'H1_O3b', 'gwtc1', 'cbfv16', 'ccb18', 'nist'],
              kind='interferometer force noise', mass_prop=True, stage='34', note=f'curve; min {env.min():.2e} at {fg[i_best]:.1f} Hz; calibration quantified 10-2000 Hz only',
              curve='stage34_ligo_curve.csv'))
def union(iv):
    iv = sorted(iv); out = []
    for lo, hi in iv:
        if out and lo <= out[-1][1]: out[-1][1] = max(out[-1][1], hi)
        else: out.append([lo, hi])
    return out
def dec_in(Uv, lo, hi):
    return sum(max(0.0, log10(min(h_, hi)) - log10(max(l_, lo))) for l_, h_ in Uv if h_ > lo and l_ < hi)
mp = [b for b in B if b['mass_prop']]
U33 = union([(b['fmin'], b['fmax']) for b in mp if b['stage'] != '34'])
U34 = union([(b['fmin'], b['fmax']) for b in mp])
U34_cal = union([(b['fmin'], b['fmax']) if b['stage'] != '34' else (10.0, 2000.0) for b in mp])
GAP = tuple(s33['coverage']['gap']); SPAN = tuple(s33['coverage']['total_span'])
cov = dict(gap=GAP, gap_decades=log10(GAP[1] / GAP[0]), covered33=dec_in(U33, *GAP), covered34=dec_in(U34, *GAP), covered34_cal_only=dec_in(U34_cal, *GAP),
           span=SPAN, span_decades=log10(SPAN[1] / SPAN[0]), span_covered33=dec_in(U33, *SPAN), span_covered34=dec_in(U34, *SPAN))
cov['newly'] = cov['covered34'] - cov['covered33']; cov['newly_cal_only'] = cov['covered34_cal_only'] - cov['covered33']
cov['open_gap'] = cov['gap_decades'] - cov['covered34']
cov['window_0p19_30'] = dict(covered=dec_in(U34, 0.1885, 30.0), total=log10(30 / 0.1885))
cov['window_35_3k'] = dict(covered=dec_in(U34, 35.0, 3000.0), total=log10(3000 / 35))
open_w = []; prev = SPAN[0]
for lo, hi in U34:
    if lo > prev and log10(lo / prev) > 1e-3: open_w.append((prev, lo, log10(lo / prev)))
    prev = max(prev, hi)
if SPAN[1] > prev: open_w.append((prev, SPAN[1], log10(SPAN[1] / prev)))
cov['open_windows'] = open_w

dropped = [
    dict(item='Germanium 100 keV - 1 MeV re-analysis', why='no fetchable count-rate table: GERDA/CDEX/Majorana spectra exist only as plots; the detection efficiency for internally emitted photons comes from collaboration MC (Donadi 2021 gives fitted efficiency only for 1000-3800 keV). Digitizing a plot would still leave the efficiency unverifiable.'),
    dict(item='LIGO below 10 Hz (0.19-10 Hz)', why='spectra not calibrated below 10 Hz ("data are not calibrated below 10 Hz and above 5 kHz")'),
    dict(item='LIGO above 5 kHz', why='not calibrated; files end at 5 kHz'),
]
flags = [
    f"O1 GW150914 PSD reproduces the published LIGO bound: H1 min over 30-35 Hz gives {o1['H1']['lam_ref_min']:.2e} vs 1.0e-5 (force {o1['H1']['SF_ref_min_fN']:.0f} fN/rtHz vs 95)",
    'Bound grows as f^4 x S_x: above ~1 kHz lambda > 1 s^-1, i.e. formally covered but physically weak (no coincidence; it is the free-mass response).',
    '2-5 kHz: calibration uncertainty "larger and unquantified" (DCC abstract) - coverage there flagged; counted separately.',
    'Spectral lines (violin modes, calibration lines, mains harmonics) make lambda(f) locally weaker; raw ASD used, no line removal (no tuning).',
    'Rigid-body response ASSUMED below 5 kHz (internal test-mass modes not modelled; no fetched mode frequency).',
]

with open(os.path.join(ROOT, 'results', 'stage34_sources.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for k, v in S.items(): w.writerow([k, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
with open(os.path.join(ROOT, 'results', 'stage34_bounds.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['name', 'kind', 'fmin_Hz', 'fmax_Hz', 'lambda_max_per_s_rC_1e-7m', 'label', 'mass_proportional', 'stage', 'sources', 'note'])
    for b in sorted(B, key=lambda x: x['fmin']):
        w.writerow([b['name'], b['kind'], f"{b['fmin']:.4g}", f"{b['fmax']:.4g}", f"{b['lam']:.3g}", b['label'], b['mass_prop'], b['stage'], ';'.join(b['sources']), b['note']])
out = dict(stage=34, title='LIGO force-noise conversion 10 Hz-5 kHz; Ge 100 keV-1 MeV dropped', validation=valid, ligo_bands=bands_stats, ligo_points=pts,
           ligo_best=dict(f=float(fg[i_best]), lam=float(env[i_best])), bands=B, coverage=cov, dropped=dropped, flags=flags, sources=S, runtime_s=time.time() - t0)
json.dump(out, open(os.path.join(ROOT, 'results', 'stage34.json'), 'w'), indent=1, default=float)

print(json.dumps(valid, indent=1))
print('LIGO band stats:'); [print(' ', k, {kk: f'{vv:.3g}' for kk, vv in v.items()}) for k, v in bands_stats.items()]
print('points:', {k: f'{v:.2e}' for k, v in pts.items()})
print(f"best {env[i_best]:.3e} at {fg[i_best]:.2f} Hz")
print(json.dumps({k: v for k, v in cov.items() if k != 'open_windows'}, indent=1, default=float))
for w_ in open_w: print(f"  open {w_[0]:.4g}-{w_[1]:.4g} Hz ({w_[2]:.2f} dec)")
