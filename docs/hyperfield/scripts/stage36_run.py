#!/usr/bin/env python3
"""Stage 36 - the 4.6 MHz-5.6 GHz and 0.19-9.3 Hz open windows (plus the 8 kHz-4.6 MHz gap) for resonant release.
Same parametrization as Stages 31-35 (CSL-equivalent lambda per nucleon, mass-proportional, r_C = 1e-7 m).
Every input fetched in THIS session; url/location/quote in results/stage36_sources.csv.
"""
import json, os, csv, time
import numpy as np
from math import pi, log10, sqrt, exp

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
t0 = time.time(); S = {}
def src(key, ref, url, loc, quote, **vals): S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)

src('nist', 'NIST CODATA 2022 Value?hbar, Value?ukg, Value?k (fetched this session)', 'https://physics.nist.gov/cgi-bin/cuu/Value?hbar', 'Numerical value fields',
    '"1.054 571 817... x 10 -34 J s" ; "1.660 539 068 92 x 10 -27 kg" ; "1.380 649 x 10 -23 J K -1"', hbar=1.054571817e-34, u=1.66053906892e-27, kB=1.380649e-23)
src('omahen25', 'Omahen, Storz, Bild, Scheiwiller, Fadel & Chu, "Ultracold Mechanical Quantum Sensor for Tests of New Physics", arXiv:2507.02653', 'https://arxiv.org/pdf/2507.02653',
    'abstract; p.3 Fig. 2 text; Table I; App. F',
    '"we use the upper bound, Pp404,max = 6.7x10^-5" ; "Pp404 = (1.2+-5.5)x10^-5" ; Table I: "Phonon mode frequency wp/2pi 5048.630(1) MHz", "Phonon mode number n 404(13)", "Phonon mode 404 lifetime T1,p 112(9) us", '
    '"sapphire density rho 3980 kg m-3", "HBAR height L 435 um" ; "lambda_CSL = (1 u/me)^2/tau_e, giving lambda_max_CSL = 5.7 x 10^-8 s^-1"',
    f=5048.630e6, n=404, T1=112e-6, dT1=9e-6, P=1.2e-5, dP=5.5e-5, rho=3980.0, L=435e-6, lam_paper=5.7e-8)
src('schrinski23', 'Schrinski et al., PRL 130, 133604 (2023) (formula check, Stage 33)', 'https://www.klaus-hornberger.de/pdf/prl19.pdf', 'p.3',
    'Stage 33 single-phonon-mode formula Gamma/lambda = 4 pi^1.5 hbar rho rc^3 k0^2 e^{-k0^2 rc^2}/(m0^2 w), validated to 1.009 vs this paper')
src('rossi18', 'Rossi, Mason, Chen, Tsaturyan & Schliesser, Nature 563, 53 (2018), arXiv:1805.05087', 'https://arxiv.org/pdf/1805.05087', 'p.2-3',
    '"millimetre-sized but 20-nm-thin, Si3N4 membrane" ; "The frequency Omega_m/2pi = 1.14 MHz of the defect mode" ; "we find Q = 1.03 x 10^9 in ringdown measurements" ; '
    '"mechanical coherence time gamma^-1 ~ (n_th Gamma_m)^-1 = hbar Q/(kB T) on the order of 1 ms ... (T ~ 10 K, n_th ~ O(10^5)) at which all reported experiments are conducted"',
    h=20e-9, f=1.14e6, Q=1.03e9, T=10.0)
src('seis22', 'Seis, Capelle, Langman, Saarinen, Planz & Schliesser, Nat. Commun. 13, 1507 (2022), arXiv:2107.05552', 'https://arxiv.org/pdf/2107.05552', 'Sec. 2.1-2.3, Eq. (1)',
    '"It consists of a 63-nm thick soft-clamped membrane made of silicon nitride" ; "A square portion of its central defect (an area of approximately 60 x 60 um^2) is covered with a 50-nm thick layer of aluminum" ; '
    '"a mechanical mode at Omega_m/2pi = 1.486 MHz" ; "tcoh = hbar Q/(kB Tbath) = 1/(n_th Gamma_m)" (Eq. 1) ; "Following eq. (1), we extract tcoh ~ 140 ms"',
    h=63e-9, f=1.486e6, tcoh=0.140, h_al=50e-9)
src('tsaturyan17', 'Tsaturyan, Barg, Polzik & Schliesser, Nature Nanotech. 12, 776 (2017), arXiv:1608.00937', 'https://arxiv.org/pdf/1608.00937', 'p.3',
    '"(rho = 3200 kg/m3 is the density of SiN)"', rho=3200.0)
src('komori20', 'Komori, Enomoto, Ooi, Miyazaki, Matsumoto, Sudhir, Michimura & Ando, PRA 101, 011802(R) (2020) (APS accepted manuscript)', 'https://link.aps.org/accepted/10.1103/PhysRevA.101.011802',
    'abstract; p.2; Eq. (4)',
    '"a 10 mg bar mirror (silica substrate of dimensions 15 mm x 1.5 mm x 0.2 mm)" ; "torque sensitivity of 20 aNm/sqrt(Hz)" ; "S_tau^CSL(w) = lambda_CSL 8 pi hbar^2 r_CSL^2 rho I/(m0^2 d)" ; '
    '"vertical seismic noise ... limits the sensitivity at low frequencies (< 50 Hz)" ; "excess in the 50-150 Hz interval"', m=10e-6, dims=(15e-3, 1.5e-3, 0.2e-3), tau=20e-18)
src('armano18', 'Armano et al. (LPF), PRL 120, 061101 (2018)', 'https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.120.061101', 'Introduction; p.3',
    '"limited, for f >~ 30 mHz, by the interferometer readout noise" ; "At frequencies above 10 mHz, the noise is dominated by displacement measurement noise" (no upper frequency of the published ASD stated in the text)')
src('fuchs24', 'Fuchs et al., Sci. Adv. 10, eadk2949 (2024), arXiv:2303.03545 (via search abstract/PubMed text)', 'https://arxiv.org/html/2303.03545', 'Results',
    '"The resulting force noise of this z-mode is approximately 0.5 fN/sqrt(Hz)" at "26.7 Hz" ; "six different modes corresponding to the three translational modes, and three rotational modes" (frequencies of the other modes not given in text)')

N = S['nist']['values']; HBAR, U, KB = N['hbar'], N['u'], N['kB']
M0, RC = U, 1e-7
B = []
def band(name, fmin, fmax, lam, label, keys, kind, note=''):
    B.append(dict(name=name, fmin=float(fmin), fmax=float(fmax), lam=float(lam), label=label, sources=keys, kind=kind, mass_prop=True, note=note, curve=None, stage='36'))
nar = lambda f: (f * 0.999, f * 1.001)          # Stage 31 convention for narrow resonators

# (A) HBAR mode 404 (DERIVED with the Stage-33 formula, validated vs Schrinski 2023 at 1.009)
o = S['omahen25']['values']
lam_ac = 2 * o['L'] / o['n']; k0 = 2 * pi / lam_ac; w = 2 * pi * o['f']
G_per_lam_hbar = 4 * pi ** 1.5 * HBAR * o['rho'] * RC ** 3 * k0 ** 2 * exp(-(k0 * RC) ** 2) / (M0 ** 2 * w)
Gmax_hbar = (o['P'] + 2 * o['dP']) / (o['T1'] - 2 * o['dT1'])     # ASSUMED conservative: +2 sigma population, -2 sigma T1, all attributed
lam_hbar = Gmax_hbar / G_per_lam_hbar
Gpub_hbar = 6.7e-5 / o['T1']; lam_hbar_1s = Gpub_hbar / G_per_lam_hbar
band('HBAR mode 404, 5.0486 GHz (Omahen 2025)', *nar(o['f']), lam_hbar, 'DERIVED (Stage 33 phonon-mode formula; +2 sigma)', ['omahen25', 'schrinski23', 'nist'], 'phonon-mode heating',
     note=f'paper own convention: 5.7e-8; with paper 1-sigma P_max: {lam_hbar_1s:.2e}')

# (B) soft-clamped SiN membranes: thin plate, mode wavelength >> r_C (DERIVED here):
#   D = 8 pi hbar^2 lam rho r^2 m_eff (1 - e^{-h^2/4r^2})/(m0^2 h);  Gamma = D/(2 m_eff hbar w)
#   Gamma/lam = 4 pi hbar rho r^2 (1 - e^{-h^2/4r^2})/(m0^2 h w)   (independent of mode shape and m_eff)
#   thick limit -> per-area slab formula of Stage 33/34 (validated vs Helou LPF and LIGO); same D convention as the HBAR formula.
def G_plate(rho, h, f):
    return 4 * pi * HBAR * rho * RC ** 2 * (1 - exp(-h ** 2 / (4 * RC ** 2))) / (M0 ** 2 * h * 2 * pi * f)
rho_sin = S['tsaturyan17']['values']['rho']
r = S['rossi18']['values']; G_r = G_plate(rho_sin, r['h'], r['f'])
heat_r = KB * r['T'] / (HBAR * r['Q'])                             # n_th Gamma_m (CITED formula, T ~ 10 K ASSUMED exact)
lam_rossi = heat_r / G_r
band('SiN membrane 1.14 MHz, 20 nm (Rossi 2018), heating n_th Gamma_m at T~10 K', *nar(r['f']), lam_rossi, 'DERIVED (thin-plate formula; all thermal heating attributed)', ['rossi18', 'tsaturyan17', 'nist'], 'membrane / resonator heating',
     note='T only "~10 K" in text: factor-2 uncertainty')
s = S['seis22']['values']; heat_s = 1 / s['tcoh']
G_s_sin = G_plate(rho_sin, s['h'], s['f'])
rho_bi = (rho_sin * s['h'] + 2700.0 * s['h_al']) / (s['h'] + s['h_al'])   # Al density 2700 ASSUMED (not fetched) - only used for a conservative alternative
G_s_bi = G_plate(rho_bi, s['h'] + s['h_al'], s['f'])
G_s = min(G_s_sin, G_s_bi)                                         # conservative: take the smaller predicted rate
lam_seis = heat_s / G_s
band('SiN membrane 1.486 MHz, 63 nm + Al pad (Seis 2022), 1/t_coh = n_th Gamma_m', *nar(s['f']), lam_seis, 'DERIVED (thin-plate formula; all heating attributed; min over SiN-only / bilayer)', ['seis22', 'tsaturyan17', 'nist'], 'membrane / resonator heating')

# (C) torsion pendulum (Komori 2020): Eq. (4) CITED; I = m L^2/12, rho = m/V (DERIVED)
k = S['komori20']['values']; Lb, Hb, db = k['dims']; rho_k = k['m'] / (Lb * Hb * db); I_k = k['m'] * Lb ** 2 / 12
lam_komori = k['tau'] ** 2 / (8 * pi * HBAR ** 2 * RC ** 2 * rho_k * I_k / (M0 ** 2 * db))
# frequency of the 20 aNm/rtHz minimum not stated in text (between the <50 Hz seismic and 50-150 Hz excess regions) -> not added as a band

# coverage
s35 = json.load(open(os.path.join(ROOT, 'results', 'stage35.json')))
Ball = [dict(b) for b in s35['bands']] + B
def union(iv):
    iv = sorted(iv); out = []
    for lo, hi in iv:
        if out and lo <= out[-1][1]: out[-1][1] = max(out[-1][1], hi)
        else: out.append([lo, hi])
    return out
def dec_in(Uv, lo, hi): return sum(max(0.0, log10(min(h_, hi)) - log10(max(l_, lo))) for l_, h_ in Uv if h_ > lo and l_ < hi)
U35 = union([(b['fmin'], b['fmax']) for b in Ball if b['mass_prop'] and b['stage'] != '36'])
U36 = union([(b['fmin'], b['fmax']) for b in Ball if b['mass_prop']])
GAP = tuple(s35['coverage']['gap']); SPAN = tuple(s35['coverage']['span'])
cov = dict(gap=GAP, gap_decades=log10(GAP[1] / GAP[0]), covered35=dec_in(U35, *GAP), covered36=dec_in(U36, *GAP), span_covered36=dec_in(U36, *SPAN), span_decades=log10(SPAN[1] / SPAN[0]),
           w_MHz_GHz=dict(total=log10(5.554e9 / 4.605e6), covered=dec_in(U36, 4.605e6, 5.554e9)), w_low=dict(total=log10(9.2907 / 0.18853), covered=dec_in(U36, 0.18853, 9.2907)))
cov['newly'] = cov['covered36'] - cov['covered35']
open_w = []; prev = SPAN[0]
for lo, hi in U36:
    if lo > prev and log10(lo / prev) > 1e-3: open_w.append((prev, lo, log10(lo / prev)))
    prev = max(prev, hi)
if SPAN[1] > prev: open_w.append((prev, SPAN[1], log10(SPAN[1] / prev)))
cov['open_windows'] = open_w
res = dict(hbar=dict(lam_ac_um=lam_ac * 1e6, G_per_lam=G_per_lam_hbar, Gmax=Gmax_hbar, lam=lam_hbar, lam_pub_Pmax=lam_hbar_1s, lam_paper=o['lam_paper']),
           rossi=dict(G_per_lam=G_r, heat=heat_r, lam=lam_rossi), seis=dict(G_per_lam_sin=G_s_sin, G_per_lam_bilayer=G_s_bi, heat=heat_s, lam=lam_seis),
           komori=dict(rho=rho_k, I=I_k, lam=lam_komori))
dropped = [
    dict(item='0.19-9.3 Hz: torsion pendula / levitated', why='no published bound with a stated frequency inside the window found. Komori 2020 torsion pendulum gives ~%.1e (DERIVED, Eq. 4) at its best torque (~50-150 Hz region, frequency not stated) - weaker than LIGO there. Fuchs 2024 levitated 26.7 Hz (outside). Altamura 2024 (Eot-Wash rotational) samples mHz. LPF 2018 ASD: published range above 30 mHz not stated in text -> not extended.' % lam_komori),
    dict(item='Wollack 2022 LiNbO3 phononic-crystal resonators (2.05/2.34 GHz)', why='mechanical thermal population not given in text; sub-micron mode geometry ~ r_C needs the full mode shape.'),
    dict(item='Quartz BAW (MHz) / Si optomechanical crystals (5 GHz)', why='no occupation/heating numbers with geometry fetched; shear modes (k.u = 0) couple weakly; not attempted this stage.'),
]
flags = [f'HBAR 5.0486 GHz: {lam_hbar:.2e} at r_C=100 nm vs paper {o["lam_paper"]:.1e} at its own convention (optimum hbar/sigma_q, not r_C = 100 nm) - same factor as Stage 33 Schrinski.',
         'Membrane bounds attribute all measured thermal heating to the effect (upper limits only, no excess search).',
         'Seis 2022 bath runs at ~80 mK while the cryostat is at 30 mK (attributed to vibration-isolation limits by the authors) - an excess heating, but not a candidate: bound uses all of it.']
with open(os.path.join(ROOT, 'results', 'stage36_sources.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for kk, v in S.items(): wr.writerow([kk, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
with open(os.path.join(ROOT, 'results', 'stage36_bounds.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['name', 'kind', 'fmin_Hz', 'fmax_Hz', 'lambda_max_per_s_rC_1e-7m', 'label', 'mass_proportional', 'stage', 'sources', 'note'])
    for b in sorted(Ball, key=lambda x: x['fmin']):
        wr.writerow([b['name'], b['kind'], f"{b['fmin']:.4g}", f"{b['fmax']:.4g}", f"{b['lam']:.3g}", b['label'], b['mass_prop'], b['stage'], ';'.join(b['sources']), b['note']])
fig_title = (f"Stage 36 — MHz membranes + 5.05 GHz HBAR: Stage-31 gap {cov['covered35']:.2f} → {cov['covered36']:.2f} dec (narrow lines); "
             f"4.6 MHz–5.6 GHz {cov['w_MHz_GHz']['covered']:.3f}/{cov['w_MHz_GHz']['total']:.2f} dec, 0.19–9.3 Hz {cov['w_low']['covered']:.2f}/{cov['w_low']['total']:.2f}")
out = dict(stage=36, title='MHz-GHz and sub-10 Hz windows', fig_title=fig_title, results=res, bands=Ball, coverage=cov, dropped=dropped, flags=flags, sources=S, runtime_s=time.time() - t0)
json.dump(out, open(os.path.join(ROOT, 'results', 'stage36.json'), 'w'), indent=1, default=float)
print(json.dumps(res, indent=1, default=float)); print(json.dumps({kk: v for kk, v in cov.items() if kk != 'open_windows'}, indent=1, default=float))
