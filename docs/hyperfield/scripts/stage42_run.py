#!/usr/bin/env python3
"""Stage 42 - narrow lines inside the 4.6 MHz - 5.6 GHz gap from published heating / occupancy numbers. No account registration.
 (a) thin-plate formula (Stage 36) for Al drums: Gamma/lambda = 4 pi hbar rho r_C^2 (1 - e^{-h^2/4r_C^2})/(m0^2 h w)
 (b) bulk single-phonon-mode formula (Stage 33, validated vs Schrinski 2023 at 1.009) for quartz BAW overtones and sapphire HBAR:
     Gamma/lambda = 4 pi^1.5 hbar rho r_C^3 k0^2 e^{-k0^2 r_C^2}/(m0^2 w), k0 = w / v.
     Bulk-only: free-surface terms are left out, which ASSUMES no destructive bulk/surface interference (conservative if they add).
Bound: CSL heating <= measured (or bath-implied) heating, all of it attributed. Inputs in results/stage42_sources.csv.
"""
import json, os, csv
from math import pi, exp
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
S = {}
def src(key, ref, url, loc, quote, **vals): S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)
src('nist', 'NIST CODATA 2022 Value?hbar, ukg (fetched this session)', 'https://physics.nist.gov/cgi-bin/cuu/Value?hbar', 'Numerical value fields',
    '"1.054 571 817... x 10 -34 J s" ; "1.660 539 068 92 x 10 -27 kg"', hbar=1.054571817e-34, u=1.66053906892e-27)
src('teufel11', 'Teufel et al., Nature 475, 359 (2011), arXiv:1103.2144 (Stage 41)', 'https://arxiv.org/pdf/1103.2144', 'p.1-2',
    '"100 nm thick aluminum membrane with a diameter of 15 um" ; "The oscillator mass m = 48 pg" -> rho_Al = 2716 kg/m3 (DERIVED)', rho_al=48e-15 / (pi * 7.5e-6 ** 2 * 100e-9))
src('cattiaux21', 'Cattiaux, Golokolenov, Kumar, Sillanpaa, Mercier de Lepinay, Gazizulin, Zhou, Armour, Bourgeois, Fefferman & Collin, Nat. Commun. 12, 6182 (2021), arXiv:2104.09541',
    'https://arxiv.org/pdf/2104.09541', 'abstract; p.2-3; Suppl. "Chip properties"',
    '"an aluminium drum-head ... of diameter 15 um, 100 nm thickness" ; "resonating at omega_m = 2pi x 15.1 MHz in its first flexure" ; '
    '"reducing the average number of quanta in the fundamental vibrational mode at 15 MHz to just 0.3" ; "timescale tau_m = 2/Gamma_m ~ 1 ms, the mechanical relaxation time"',
    f=15.1e6, h=100e-9, n=0.3, tau=1e-3)
src('goryachev12', 'Goryachev, Creedon, Ivanov, Galliou, Bourquin & Tobar, Appl. Phys. Lett. 100, 243504 (2012), arXiv:1202.4556', 'https://arxiv.org/pdf/1202.4556', 'abstract; Table I',
    '"When cooled to 20 mK ... (with Q-factors of order billions) at frequencies of 15.6 and 65.4 MHz" ; Table I: "Quartz BAW - 5th overtone 1.56x10^7 [Hz] 2x10^9 [Q] 0.018 [K] 23.5 [nTH]" ; '
    '"Quartz BAW - 21st overtone 6.54x10^7 1.2x10^9 0.021 6.2"', f5=1.56e7, Q5=2e9, n5=23.5, f21=6.54e7, Q21=1.2e9, n21=6.2)
src('galliou13', 'Galliou, Goryachev, Bourquin, Abbe, Aubry & Tobar, Sci. Rep. 3, 2132 (2013), arXiv:1309.4832', 'https://arxiv.org/pdf/1309.4832', 'p.4',
    '"the quasi-longitudinal (V_L = 6757 m/s)" (same BVA SC-cut family; ASSUMED to apply to the Goryachev resonator) ; "1.08 mm thick ... designed to sustain shear vibration of 5 MHz" '
    '(check: 2 x 1.08 mm x 15.6 MHz/5 = 6740 m/s, DERIVED, agrees to 0.3%)', vL=6757.0)
src('quartz_rho', 'NCBI Bookshelf, ATSDR Toxicological Profile for Silica, Table 4-2 (citing Haynes et al. 2014)', 'https://www.ncbi.nlm.nih.gov/books/NBK592829/table/ch4.tab2/', 'Density row',
    '"Density (g/cm3) at 20 C: 2.648 (alpha-quartz)"', rho=2648.0)
src('chu17', 'Chu, Kharel, Renninger, Burkhart, Frunzio, Rakich & Schoelkopf, Science 358, 199 (2017), arXiv:1703.00342', 'https://arxiv.org/pdf/1703.00342', 'p.5-10',
    '"has a frequency of ~6.65 GHz" ; "longitudinal sound velocity vl = 1.11x10^4 m/s" ; "the qubit\'s ground state population increases to 98%. This value ... represents a lower bound on the phonon ground state population" ; '
    '"T1 = 17 +- 1 us"', f=6.65e9, vl=1.11e4, Pexc=0.02, T1=17e-6, dT1=1e-6)
src('omahen25', 'Omahen et al., arXiv:2507.02653, Table I (fetched in Stage 36 of this overnight run)', 'https://arxiv.org/pdf/2507.02653', 'Table I', '"sapphire density rho 3980 kg m-3"', rho=3980.0)
src('maccabe20', 'MacCabe, Ren, Luo, Cohen, Zhou, Sipahigil, Mirhosseini & Painter, Science 370, 840 (2020), arXiv:1901.04129', 'https://arxiv.org/pdf/1901.04129', 'abstract; p.3; Methods',
    '"fundamental 5 GHz acoustic mode" ; "increasing up to tau_ph,0 ~ 1.5 seconds" ; "220 nm device layer" -- NOT converted: no validated CSL formula for a localized optomechanical-crystal breathing mode '
    '(k ~ lattice scale, k r_C ~ 1) and no dark-state occupancy/heating number quoted')
src('lpf_archive', 'ESA COSMOS LISA Pathfinder Archive (Stage 41)', 'https://www.cosmos.esa.int/web/lisa-pathfinder-archive', 'portal menu', '"SIGN IN" / "Self Registration" only -- not used, no account created')
N = S['nist']['values']; HBAR, M0, RC = N['hbar'], N['u'], 1e-7
def G_plate(rho, h, f): return 4 * pi * HBAR * rho * RC ** 2 * (1 - exp(-h ** 2 / (4 * RC ** 2))) / (M0 ** 2 * h * 2 * pi * f)
def G_bulk(rho, f, v):
    w = 2 * pi * f; k0 = w / v
    return 4 * pi ** 1.5 * HBAR * rho * RC ** 3 * k0 ** 2 * exp(-(k0 * RC) ** 2) / (M0 ** 2 * w)
rows = []
c = S['cattiaux21']['values']; rho_al = S['teufel11']['values']['rho_al']
heat = c['n'] * 2 / c['tau']; G = G_plate(rho_al, c['h'], c['f'])
rows.append(dict(name='Al drum 15.1 MHz, 500 uK (Cattiaux 2021)', f=c['f'], heat=heat, G_per_lam=G, lam=heat / G, formula='thin plate',
                 label='DERIVED (thin-plate; measured n = 0.3, Gamma_m = 2/tau_m)', keys=['cattiaux21', 'teufel11', 'nist'], note='measured occupancy; rho_Al from Teufel'))
g = S['goryachev12']['values']; q = S['quartz_rho']['values']['rho']; vL = S['galliou13']['values']['vL']
for tag in ('5', '21'):
    f, Q, n = g['f' + tag], g['Q' + tag], g['n' + tag]
    heat = n * 2 * pi * f / Q; G = G_bulk(q, f, vL)
    rows.append(dict(name=f'quartz BAW OT{tag} {f/1e6:.1f} MHz, 20 mK (Goryachev 2012)', f=f, heat=heat, G_per_lam=G, lam=heat / G, formula='bulk phonon mode',
                     label='DERIVED (bulk formula); ASSUMED mode thermalised at bath T (nTH from T, not measured) and surface terms add', keys=['goryachev12', 'galliou13', 'quartz_rho', 'nist'],
                     note='bath-implied heating n_TH*w/Q; no occupancy measurement'))
h = S['chu17']['values']; heat = h['Pexc'] / (h['T1'] - 2 * h['dT1']); G = G_bulk(S['omahen25']['values']['rho'], h['f'], h['vl'])
rows.append(dict(name='HBAR ~6.65 GHz (Chu 2017)', f=h['f'], heat=heat, G_per_lam=G, lam=heat / G, formula='bulk phonon mode',
                 label='DERIVED (Stage 33 formula; P_exc <= 2%, T1 - 2 sigma)', keys=['chu17', 'omahen25', 'nist'], note='above Omahen 5.05 GHz; weak'))
B = [dict(name=r['name'], fmin=r['f'] * 0.999, fmax=r['f'] * 1.001, lam=r['lam'], label=r['label'], sources=r['keys'], kind='membrane / resonator heating' if r['formula'] == 'thin plate' else 'phonon-mode heating',
          mass_prop=True, note=r['note'], curve=None, stage='42') for r in rows]
s41 = json.load(open(os.path.join(ROOT, 'results', 'stage41.json')))
Ball = [dict(b) for b in s41['bands']] + B
fig_title = 'Stage 42 — 4.6 MHz–5.6 GHz gap: Al drum 15.1 MHz {:.1e}; quartz BAW 15.6/65.4 MHz {:.1e}/{:.1e}; HBAR 6.65 GHz {:.1e} s$^{{-1}}$'.format(*[r['lam'] for r in rows])
with open(os.path.join(ROOT, 'results', 'stage42_sources.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for kk, v in S.items(): wr.writerow([kk, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
with open(os.path.join(ROOT, 'results', 'stage42_bounds.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['name', 'kind', 'fmin_Hz', 'fmax_Hz', 'lambda_max_per_s_rC_1e-7m', 'label', 'stage', 'sources', 'note'])
    for b in sorted(Ball, key=lambda x: x['fmin']):
        wr.writerow([b['name'], b['kind'], f"{b['fmin']:.4g}", f"{b['fmax']:.4g}", f"{b['lam']:.3g}", b['label'], b['stage'], ';'.join(b['sources']), b['note']])
json.dump(dict(stage=42, fig_title=fig_title, results=rows, bands=Ball, coverage=s41['coverage'], sources=S), open(os.path.join(ROOT, 'results', 'stage42.json'), 'w'), indent=1, default=float)
for r in rows: print(r['name'], f"heat={r['heat']:.3g}/s G/lam={r['G_per_lam']:.3g} lam<={r['lam']:.3g}")
