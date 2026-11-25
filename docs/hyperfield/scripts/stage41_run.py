#!/usr/bin/env python3
"""Stage 41 - fallback narrow-line scan (the ESA LISA Pathfinder archive requires self-registration/login -> not used).
Two published superconducting drum resonators with measured heating rates, via the Stage 36 thin-plate formula
  Gamma/lambda = 4 pi hbar rho r_C^2 (1 - e^{-h^2/4 r_C^2}) / (m0^2 h w)   (mode wavelength >> r_C).
Bound: CSL heating <= measured total heating (all attributed). Every input fetched this session (results/stage41_sources.csv).
"""
import json, os, csv
from math import pi, exp, log10
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
S = {}
def src(key, ref, url, loc, quote, **vals): S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)
src('nist', 'NIST CODATA 2022 Value?hbar, ukg (fetched this session)', 'https://physics.nist.gov/cgi-bin/cuu/Value?hbar', 'Numerical value fields',
    '"1.054 571 817... x 10 -34 J s" ; "1.660 539 068 92 x 10 -27 kg"', hbar=1.054571817e-34, u=1.66053906892e-27)
src('lpf_archive', 'ESA COSMOS "LISA Pathfinder Archive" portal (checked this session)', 'https://www.cosmos.esa.int/web/lisa-pathfinder-archive', 'portal menu',
    'menu shows "SIGN IN" and "Self Registration" only; no anonymous data link -> login required, item switched to fallback per instructions')
src('teufel11', 'Teufel, Donner, Li, Harlow, Allman, Cicak, Sirois, Whittaker, Lehnert & Simmonds, Nature 475, 359 (2011), arXiv:1103.2144', 'https://arxiv.org/pdf/1103.2144', 'p.1-2',
    '"The oscillator is a 100 nm thick aluminum membrane with a diameter of 15 um" ; "Omega_m = 2pi x 10.56 MHz, and we find an intrinsic damping rate of Gamma_m = 2pi x 32 Hz" ; "The oscillator mass m = 48 pg" ; '
    '"for the lowest achievable temperature of 15 mK. This temperature corresponds to a thermal occupancy nm = 30"',
    h=100e-9, d=15e-6, f=10.56e6, Gm=2 * pi * 32.0, m=48e-15, n=30.0)
src('youssefi23', 'Youssefi, Kono, Chegnizadeh & Kippenberg, "A squeezed mechanical oscillator with milli-second quantum decoherence", Nat. Phys. 19, 1697 (2023), arXiv:2208.13082', 'https://arxiv.org/pdf/2208.13082',
    'abstract; p.2; Suppl. Fig. S2', '"Omega_m/2pi = 1.8 MHz, corresponding to Qm = 40 x 10^6" ; "a thermal decoherence rate of 20.5 Hz (corresponding to 130 quanta/sec motional heating rate)" ; '
    '"The radius of the drumhead used in the main work is R = 75um" ; "considering t = 180 nm and sigma_Al ~ 350 MPa" (Al top plate thickness, from the FEM caption)',
    f=1.8e6, heat=130.0, h=180e-9)
N = S['nist']['values']; HBAR, U = N['hbar'], N['u']; M0, RC = U, 1e-7
def G_plate(rho, h, f): return 4 * pi * HBAR * rho * RC ** 2 * (1 - exp(-h ** 2 / (4 * RC ** 2))) / (M0 ** 2 * h * 2 * pi * f)
t = S['teufel11']['values']; rho_al = t['m'] / (pi * (t['d'] / 2) ** 2 * t['h'])      # DERIVED Al density from quoted mass and geometry
heat_t = t['n'] * t['Gm']                                                             # n_th Gamma_m at 15 mK (all attributed)
G_t = G_plate(rho_al, t['h'], t['f']); lam_t = heat_t / G_t
y = S['youssefi23']['values']; G_y = G_plate(rho_al, y['h'], y['f']); lam_y = y['heat'] / G_y   # rho_Al from Teufel (DERIVED); h = 180 nm (FEM caption, ASSUMED = device)
rows = [dict(name='Al drum 10.56 MHz (Teufel 2011)', f=t['f'], h=t['h'], rho=rho_al, heat=heat_t, G_per_lam=G_t, lam=lam_t, note='first membrane point in the 4.6 MHz-5.6 GHz window'),
        dict(name='Al vacuum-gap drum 1.8 MHz (Youssefi 2023)', f=y['f'], h=y['h'], rho=rho_al, heat=y['heat'], G_per_lam=G_y, lam=lam_y, note='thickness from FEM caption (ASSUMED equal to device)')]
B = [dict(name=r['name'], fmin=r['f'] * 0.999, fmax=r['f'] * 1.001, lam=r['lam'], label='DERIVED (Stage 36 thin-plate formula; all measured heating attributed)', sources=['teufel11' if 'Teufel' in r['name'] else 'youssefi23', 'teufel11', 'nist'],
          kind='membrane / resonator heating', mass_prop=True, note=r['note'], curve=None, stage='41') for r in rows]
s40 = json.load(open(os.path.join(ROOT, 'results', 'stage40.json')))
Ball = [dict(b) for b in s40['bands']] + B
cov = dict(s40['coverage'])
fig_title = f"Stage 41 — superconducting drums: {lam_y:.1e} s$^{{-1}}$ at 1.8 MHz, {lam_t:.1e} at 10.56 MHz (first membrane point in 4.6 MHz–5.6 GHz); LPF archive needs login"
with open(os.path.join(ROOT, 'results', 'stage41_sources.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for kk, v in S.items(): wr.writerow([kk, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
with open(os.path.join(ROOT, 'results', 'stage41_bounds.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['name', 'kind', 'fmin_Hz', 'fmax_Hz', 'lambda_max_per_s_rC_1e-7m', 'label', 'stage', 'sources', 'note'])
    for b in sorted(Ball, key=lambda x: x['fmin']):
        wr.writerow([b['name'], b['kind'], f"{b['fmin']:.4g}", f"{b['fmax']:.4g}", f"{b['lam']:.3g}", b['label'], b['stage'], ';'.join(b['sources']), b['note']])
json.dump(dict(stage=41, fig_title=fig_title, results=rows, bands=Ball, coverage=cov, sources=S), open(os.path.join(ROOT, 'results', 'stage41.json'), 'w'), indent=1, default=float)
for r in rows: print({k: (f'{v:.3g}' if isinstance(v, float) else v) for k, v in r.items()})
