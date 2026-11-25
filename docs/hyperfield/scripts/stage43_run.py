#!/usr/bin/env python3
"""Stage 43 - full (bulk + free-surface) CSL heating of a thickness-type mode, validated before use.

Model (mass-proportional CSL, linearised in the displacement u = q*phi(x)):
  momentum diffusion  eta = lam (4 pi rC^2)^{3/2} / (m0^2 (2pi)^3) * Int d^3k e^{-k^2 rC^2} |F(k)|^2,  F(k) = i Int rho(x) k.phi(x) e^{ik.x} d^3x
  phonon heating      Gamma = hbar eta / (2 M_eff w),  M_eff = Int rho |phi|^2
Slab 0<z<L with free faces, laterally wide (lateral scale >> rC), phi = z_hat cos(n pi z/L):
  Gamma/lam = hbar (4pi)^{3/2} rC^3 rho / (2 pi m0^2 L w) * Int_{-inf}^{inf} dk k^2 |I(k)|^2 e^{-k^2 rC^2},  I(k) = Int_0^L cos(k0 z) e^{ikz} dz
Limits: bulk peaks -> Stage 33 formula 4 pi^1.5 hbar rho rC^3 k0^2 e^{-k0^2 rC^2}/(m0^2 w); far tails -> surface term 8 pi hbar rho rC^2/(m0^2 L w).
Validation A (user criterion, within ~10%): Omahen 2025 published lambda_max = 5.7e-8 s^-1 at r_CSL = 3.0e-7 m.
Validation B (surface machinery): same integral for a rigid cube vs the closed form of Carlesso et al. 2016 Eq. (12b).
"""
import json, os, csv
import numpy as np
from math import pi, exp, sqrt, erf
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
S = {}
def src(key, ref, url, loc, quote, **vals): S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)
src('nist', 'NIST CODATA 2022 (fetched this session)', 'https://physics.nist.gov/cgi-bin/cuu/Value?hbar', 'Numerical value fields',
    '"1.054 571 817... x 10 -34 J s" ; "1.660 539 068 92 x 10 -27 kg"', hbar=1.054571817e-34, u=1.66053906892e-27)
src('omahen25', 'Omahen, Storz, Bild, Scheiwiller, Fadel & Chu, arXiv:2507.02653', 'https://arxiv.org/pdf/2507.02653', 'Table I; Table II; App. F',
    '"Phonon mode frequency 5048.630(1) MHz", "Phonon mode number n 404(13)", "Phonon mode 404 lifetime T1,p 112(9) us", "sapphire density rho 3980 kg m-3", "HBAR height L 435 um" ; '
    'Table II: "Phonon T1 100 us", "Excited state population 6.7 x 10-5", "Estimated CSL rate 5.7 x 10-8 s-1", "Estimated CSL radius 3.0 x 10-7 m" ; App. F: "lambda_max_CSL = 5.7 x 10^-8 s^-1" ; "rCSL = hbar/sqrt(2) sigma_q = 3.0 x 10-7 m"',
    f=5048.630e6, n=404, T1=112e-6, dT1=9e-6, T1_tab=100e-6, P=1.2e-5, dP=5.5e-5, Pmax=6.7e-5, rho=3980.0, L=435e-6, lam_pub=5.7e-8, r_pub=3.0e-7)
src('carlesso16', 'Carlesso, Bassi, Falferi & Vinante, "Experimental bounds on collapse models from gravitational wave detectors", PRD 94, 124036 (2016), arXiv:1606.04581',
    'https://arxiv.org/pdf/1606.04581', 'Eqs. (1), (11), (12b)',
    'Eq. (11): "S_FF(w) = hbar^2 lam rC^3/(2 pi^{3/2} m0^2) Int dk |mu(k)|^2 |1 - e^{i a kx}|^2 kx^2 e^{-rC^2 k^2}" ; Eq. (12b): "S^LISA_FF(w) = 16 hbar^2 lam m^2/(L^2 m0^2) (rC/L)^4 (1 - e^{-L^2/4rC^2} + fcorr) (1 - e^{-L^2/4rC^2} - sqrt(pi) L/(2rC) erf(L/2rC))^2" '
    '(relative force F_rel = (F1-F2)/2 of two cubes; single cube = 2x this with fcorr = 0)')
src('goryachev12', 'Goryachev et al., APL 100, 243504 (2012), arXiv:1202.4556 (Stage 42)', 'https://arxiv.org/pdf/1202.4556', 'Table I',
    '"Quartz BAW - 5th overtone 1.56x10^7 2x10^9 0.018 23.5" ; "Quartz BAW - 21st overtone 6.54x10^7 1.2x10^9 0.021 6.2"', f5=1.56e7, Q5=2e9, n5=23.5, f21=6.54e7, Q21=1.2e9, n21=6.2)
src('galliou13', 'Galliou et al., Sci. Rep. 3, 2132 (2013), arXiv:1309.4832 (Stage 42)', 'https://arxiv.org/pdf/1309.4832', 'p.4',
    '"quasi-longitudinal (V_L = 6757 m/s)" ; "the excitation electrodes are separated from the crystal plate by a tiny vacuum gap" (free faces)', vL=6757.0)
src('quartz_rho', 'NCBI Bookshelf ATSDR Silica Table 4-2 (Stage 42)', 'https://www.ncbi.nlm.nih.gov/books/NBK592829/table/ch4.tab2/', 'Density row', '"2.648 (alpha-quartz)"', rho=2648.0)
src('chu17', 'Chu et al., Science 358, 199 (2017), arXiv:1703.00342 (Stage 42)', 'https://arxiv.org/pdf/1703.00342', 'p.5-10',
    '"~6.65 GHz" ; "vl = 1.11x10^4 m/s" ; "For a measured substrate thickness of 420 um" ; "ground state population increases to 98% ... a lower bound" ; "T1 = 17 +- 1 us"',
    f=6.65e9, vl=1.11e4, L=420e-6, Pexc=0.02, T1=17e-6, dT1=1e-6)
N = S['nist']['values']; HBAR, M0 = N['hbar'], N['u']

def kgrid(L, rC, per=24):
    kmax = 8.0 / rC; dk = min(pi / L / per, 0.02 / rC)
    return np.arange(dk / 2, kmax, dk), dk
def ex(x, L):            # (e^{ixL}-1)/(ix) = L e^{ixL/2} sinc(xL/2)
    return L * np.exp(1j * x * L / 2) * np.sinc(x * L / (2 * pi))
def G_full(rho, L, n, f, rC):
    k0 = n * pi / L; w = 2 * pi * f; k, dk = kgrid(L, rC)
    I = 0.5 * (ex(k + k0, L) + ex(k - k0, L))
    integ = 2 * np.sum(k ** 2 * np.abs(I) ** 2 * np.exp(-(k * rC) ** 2)) * dk          # even integrand
    return HBAR * (4 * pi) ** 1.5 * rC ** 3 * rho / (2 * pi * M0 ** 2 * L * w) * integ
def G_bulk(rho, L, n, f, rC):
    k0 = n * pi / L; w = 2 * pi * f
    return 4 * pi ** 1.5 * HBAR * rho * rC ** 3 * k0 ** 2 * exp(-(k0 * rC) ** 2) / (M0 ** 2 * w)
def G_surf(rho, L, f, rC): return 8 * pi * HBAR * rho * rC ** 2 / (M0 ** 2 * L * 2 * pi * f)

out = {}
# ---- Validation B: rigid cube vs Carlesso Eq. 12b (units: hbar = lam = m0 = rho = 1, rC = 1) ----
vb = []
for Lr in (0.5, 1, 2, 5, 20, 100):
    k, dk = kgrid(Lr, 1.0, per=40)
    sinc2 = np.abs(ex(k, Lr)) ** 2
    lon = 2 * np.sum(k ** 2 * sinc2 * np.exp(-k ** 2)) * dk
    tr = 2 * np.sum(sinc2 * np.exp(-k ** 2)) * dk
    mine = (1.0 / pi ** 1.5) * lon * tr ** 2                     # hbar^2 lam rC^3/(pi^1.5 m0^2) * Int d^3k k_x^2 |mu~|^2 e^{-k^2 rC^2}
    m = Lr ** 3; e = exp(-Lr ** 2 / 4)
    carl = 2 * 16 * m ** 2 / Lr ** 2 * (1 / Lr) ** 4 * (1 - e) * (1 - e - sqrt(pi) * Lr / 2 * erf(Lr / 2)) ** 2
    vb.append(dict(L_over_rC=Lr, numeric=mine, carlesso=carl, ratio=mine / carl))
out['validation_B_rigid_cube'] = vb
# ---- internal: thickness mode = bulk + surface in the separated limit ----
o = S['omahen25']['values']
out['internal_check'] = [dict(rC=r, full=G_full(o['rho'], o['L'], o['n'], o['f'], r), bulk=G_bulk(o['rho'], o['L'], o['n'], o['f'], r), surf=G_surf(o['rho'], o['L'], o['f'], r)) for r in (1e-7, 3e-7)]
for d in out['internal_check']: d['ratio_full_to_sum'] = d['full'] / (d['bulk'] + d['surf']); d['surf_fraction'] = d['surf'] / d['full']
# ---- Validation A: Omahen published number ----
Gf = G_full(o['rho'], o['L'], o['n'], o['f'], o['r_pub'])
va = dict(G_per_lam_full_r3e7=Gf, lam_T1_112us=o['Pmax'] / o['T1'] / Gf, lam_T1_100us=o['Pmax'] / o['T1_tab'] / Gf, lam_pub=o['lam_pub'])
va['ratio_112'] = va['lam_T1_112us'] / o['lam_pub']; va['ratio_100'] = va['lam_T1_100us'] / o['lam_pub']
va['surface_fraction_at_r3e7'] = out['internal_check'][1]['surf_fraction']
out['validation_A_omahen'] = va
passA = abs(va['ratio_112'] - 1) <= 0.10 or abs(va['ratio_100'] - 1) <= 0.10
passB = all(abs(v['ratio'] - 1) <= 0.01 for v in vb)
out['validated'] = bool(passA and passB)
# ---- apply (only if validated) at rC = 1e-7 ----
RC = 1e-7; upd = []
if out['validated']:
    Gm = G_full(o['rho'], o['L'], o['n'], o['f'], RC); Gmax = (o['P'] + 2 * o['dP']) / (o['T1'] - 2 * o['dT1'])
    upd.append(dict(match='Omahen 2025', f=o['f'], L=o['L'], n=o['n'], G_full=Gm, G_bulk=G_bulk(o['rho'], o['L'], o['n'], o['f'], RC), lam=Gmax / Gm))
    g = S['goryachev12']['values']; vL = S['galliou13']['values']['vL']; rq = S['quartz_rho']['values']['rho']
    for tag, nn in (('5', 5), ('21', 21)):
        f = g['f' + tag]; L = nn * vL / (2 * f)                                   # DERIVED thickness from overtone number
        heat = g['n' + tag] * 2 * pi * f / g['Q' + tag]; Gq = G_full(rq, L, nn, f, RC)
        upd.append(dict(match=f'quartz BAW OT{tag}', f=f, L=L, n=nn, G_full=Gq, G_bulk=G_bulk(rq, L, nn, f, RC), lam=heat / Gq))
    c = S['chu17']['values']; nn = round(2 * c['L'] * c['f'] / c['vl']); Gc = G_full(S['omahen25']['values']['rho'], c['L'], nn, c['f'], RC)
    upd.append(dict(match='Chu 2017', f=c['f'], L=c['L'], n=nn, G_full=Gc, G_bulk=G_bulk(3980.0, c['L'], nn, c['f'], RC), lam=c['Pexc'] / (c['T1'] - 2 * c['dT1']) / Gc))
out['updates'] = upd
s42 = json.load(open(os.path.join(ROOT, 'results', 'stage42.json')))
Ball = []
for b in s42['bands']:
    b = dict(b)
    for u in upd:
        if u['match'] in b['name']:
            b['lam_old_bulk_only'] = b['lam']; b['lam'] = u['lam']; b['stage'] = '43'
            b['label'] = 'DERIVED (Stage 43 full bulk+surface thickness-mode formula, validated A: Omahen, B: Carlesso Eq. 12b)'
            b['name'] = b['name'] + ' [bulk+surface]'
    Ball.append(b)
q = [u for u in upd if 'quartz' in u['match']]
fig_title = ('Stage 43 — bulk+surface thickness-mode formula validated (Omahen {:.2f}, Carlesso cube to {:.0e}); quartz BAW → {:.1e}/{:.1e} s$^{{-1}}$'
             .format(va['ratio_112'], max(abs(v['ratio'] - 1) for v in vb), q[0]['lam'], q[1]['lam']) if out['validated'] else 'Stage 43 — surface formula NOT validated; bounds unchanged')
with open(os.path.join(ROOT, 'results', 'stage43_sources.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for kk, v in S.items(): wr.writerow([kk, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
with open(os.path.join(ROOT, 'results', 'stage43_bounds.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['name', 'kind', 'fmin_Hz', 'fmax_Hz', 'lambda_max_per_s_rC_1e-7m', 'lambda_old_bulk_only', 'label', 'stage', 'sources', 'note'])
    for b in sorted(Ball, key=lambda x: x['fmin']):
        wr.writerow([b['name'], b['kind'], f"{b['fmin']:.4g}", f"{b['fmax']:.4g}", f"{b['lam']:.3g}", f"{b.get('lam_old_bulk_only', float('nan')):.3g}", b['label'], b['stage'], ';'.join(b['sources']), b['note']])
json.dump(dict(stage=43, fig_title=fig_title, results=out, bands=Ball, coverage=s42['coverage'], sources=S), open(os.path.join(ROOT, 'results', 'stage43.json'), 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=lambda x: float(x)))
