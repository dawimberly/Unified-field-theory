#!/usr/bin/env python3
"""Stage 35 - the largest open window for resonant release, 13 GHz - 14.5 keV (do-si-do element 3).

Same parametrization as Stages 31-34: extra release rate per nucleon expressed as the CSL-equivalent, mass-proportional
collapse rate lambda at r_C = 1e-7 m; an experiment sampling the noise at frequency f bounds lambda(f).
Items (Dan's list): qubit quasiparticles, low-threshold cryogenic detectors (and their low-energy excess, LEE),
CMB/FIRAS, optical/IR/X-ray spontaneous emission.
Every input fetched in THIS session; url/location/quote in results/stage35_sources.csv.
Labels: CITED (number from fetched source), DERIVED (computed here), ASSUMED (stated choice).
"""
import json, os, csv, time
import numpy as np
from math import pi, log10, sqrt
from scipy.special import erfcx

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
t0 = time.time()
S = {}
def src(key, ref, url, loc, quote, **vals):
    S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)

src('nist', 'NIST CODATA 2022 pages Value?hbar,e,ep0,c,bg,ukg,h,k (fetched this session)', 'https://physics.nist.gov/cgi-bin/cuu/Value?hbar', 'Numerical value fields',
    '"1.054 571 817... x 10 -34 J s" ; "1.602 176 634 x 10 -19 C" ; "8.854 187 8188 x 10 -12 F m -1" ; "299 792 458 m s -1" ; "6.674 30 x 10 -11 m 3 kg -1 s -2" ; '
    '"1.660 539 068 92 x 10 -27 kg" ; "6.626 070 15 x 10 -34 J Hz -1" ; "1.380 649 x 10 -23 J K -1"',
    hbar=1.054571817e-34, e=1.602176634e-19, ep0=8.8541878188e-12, c=299792458.0, G=6.67430e-11, u=1.66053906892e-27, h=6.62607015e-34, kB=1.380649e-23)
src('adler_vinante18', 'Adler & Vinante, PRA 97, 052119 (2018), arXiv:1801.06857', 'https://arxiv.org/pdf/1801.06857', 'p.2-3, after Eq. (3)',
    '"the residual heating from unknown sources is limited to roughly 10^-11 W/kg. This means that the effective lambda for bulk heating is at most 3.1 x 10^-11 s^-1" ; '
    '"heating takes place only if the noise has frequency components of at least wL(|q| = rc^-1) ~ vs |q| ~ 0.4 x 10^11 s^-1 ... (which for this estimate we have taken as 4000 m/s characteristic of copper at low temperature)"',
    lam=3.1e-11, vs=4000.0, P_resid=1e-11)
src('ccb18', 'Carlesso, Ferialdi & Bassi, EPJD 72, 159 (2018), arXiv:1805.10100 (as fetched for Stage 33/34; Eq. 19)', 'https://arxiv.org/pdf/1805.10100', 'Sec. IV, Eq. (19)',
    '"lambda_eff = 2 lambda r_C^5/(3 pi^3/2) Int dq e^{-q^2 r_C^2} q^2 f(omega_L(q)) ... In the white noise case lambda_eff = lambda"')
src('majorana22', 'Majorana Collaboration (Arnquist et al.), PRL 129, 080401 (2022), arXiv:2202.01343', 'https://arxiv.org/pdf/2202.01343', 'p.2, Eqs. (1)-(2)',
    '"dGamma(E)/dE = A_f x hbar lambda/(4 pi^2 eps0 m0^2 c^3 rC^2 E)" (Eq. 2, m-p) ; "In the m-p CSL, the X-ray emission rate is independent of m_X^2" ; '
    '"the wavelength of the emitted photon lambda_k is smaller than the distance between electrons and the nucleus (0.1 nm) but much larger than the typical size of the nucleus (10^-5 nm). This corresponds to the energy range of 10 - 10^5 keV"')
src('donadi21', 'Donadi, Piscicchia, Del Grande et al., "Novel CSL bounds from the noise-induced radiation emission from atoms", EPJC 81, 773 (2021), arXiv:2107.11237', 'https://arxiv.org/pdf/2107.11237',
    'Sec. II p.5-6, Eq. (5); Appendix Eqs. (44)-(49)',
    '"one needs to find a region where the wavelength of the emitted radiation lambda_k is small enough to guarantee that the electrons and the nucleus emit independently from each other, so that their contribution does not add coherently to zero, since they have opposite charges. Therefore one has to consider wavelengths lambda_k << 0.1 nm" ; '
    '"when the distance between two particles is larger than lambda_k the emission is always incoherent. On the other hand, particles at a distance smaller than lambda_k and rC always emit coherently" ; '
    'Eq. (44) "dGamma_t/dE = hbar lambda/(4 pi^2 eps0 c^3 m0^2 rC^2 E) Sum_i q_i^2" (incoherent)', d_atom=1e-10)
src('review22', 'Carlesso et al., Nature Phys. 18, 243 (2022), arXiv:2203.04231', 'https://arxiv.org/pdf/2203.04231', 'Sec. 2.5, Eq. (6)',
    '"Eq. (6) is valid for E in [10, 10^5] keV, which corresponds to photon wavelengths larger than the size of a nucleus but smaller than that of an atom"')
src('excess22', 'Adari et al., "EXCESS workshop: Descriptions of rising low-energy spectra", SciPost Phys. Proc. 9, 001 (2022), arXiv:2202.05097', 'https://arxiv.org/pdf/2202.05097',
    'Sec. 1; CRESST-III section; EDELWEISS RED20/RED30 section; Sec. on scaling',
    '"many experiments are observing previously undetected excesses of low-energy events of unknown origin ... Typically starting at energies below 1 keV, the obtained energy spectra rise sharply towards the detector thresholds" ; '
    '"Due to variations in the shape and rate of the excess signal among experiments, detectors and measurements, a dark matter (DM) explanation seems unlikely" ; '
    '"the rate and spectral shape of the excess contributions is not compatible between these different detectors which disfavours a common origin" ; '
    '"The event rate at 200 eV is 10^5 count/keV/kg/day, decreasing to 10^4 count/keV/kg/day at 1 keV" (EDELWEISS RED20) ; '
    '"the shape of the spectra does not vary significantly in time, while the absolute rate decreases slowly over time. Sudden increases were observed at times after warming-up the detectors above 10 K"',
    R200=1e5, R1000=1e4)
src('serniak18', 'Serniak et al., PRL 121, 157701 (2018), arXiv:1803.00476 "Hot nonequilibrium quasiparticles in transmon qubits"', 'https://arxiv.org/pdf/1803.00476', 'p.3',
    '"Here Delta = 205 ueV, consistent with DC measurements"', Delta_eV=205e-6)
src('planck18', 'Planck Collaboration 2018 VI, A&A 641, A6 (2020), arXiv:1807.06209, abstract', 'https://arxiv.org/abs/1807.06209', 'abstract',
    '"baryon density Omega_b h^2 = 0.0224+-0.0001" ; "Hubble constant H_0 = (67.4+-0.5) km/s/Mpc"', Obh2=0.0224, H0=67.4)
src('parsec', 'Wikipedia "Parsec" (IAU 2015 Resolution B2 definition)', 'https://en.wikipedia.org/wiki/Parsec', 'infobox; Sec. History',
    '"1 pc ... 3.0857x10^16 m" ; "the parsec is defined as exactly (648 000/pi) au"', pc=3.0857e16)
src('hickox06', 'Hickox & Markevitch, ApJ 645, 95 (2006), arXiv:astro-ph/0512542', 'https://arxiv.org/pdf/astro-ph/0512542', 'abstract',
    '"we obtain a total intensity of the extragalactic CXB of (4.6 +- 0.3) x 10^-12 ergs cm^-2 s^-1 deg^-2 for 1-2 keV and (1.7 +- 0.2)x 10^-11 ergs cm^-2 s^-1 deg^-2 for 2-8 keV"',
    I12=(4.6e-12, 0.3e-12), I28=(1.7e-11, 0.2e-11))
src('fixsen96', 'Fixsen et al., ApJ 473, 576 (1996), arXiv:astro-ph/9605054 (COBE FIRAS)', 'https://arxiv.org/pdf/astro-ph/9605054', 'abstract; Sec. 5',
    '"temperature, 2.728 +- 0.004 K" ; "The RMS deviations are less than 50 parts per million of the peak of the CMBR" ; "the weighted rms deviation in the frequency range 2 to 21 cm-1 is 50 ppm (parts per million) of the peak brightness"',
    T=2.728, ppm=50e-6, nu_cm=(2.0, 21.0))

N = S['nist']['values']
HBAR, E_, EP0, C, G, U, H, KB = N['hbar'], N['e'], N['ep0'], N['c'], N['G'], N['u'], N['h'], N['kB']
M0, RC = U, 1e-7                       # ASSUMED m0 = 1 u, r_C = 100 nm (Stages 31-34)
keV = 1e3 * E_
Hz_of_eV = lambda eV: eV * E_ / H
B = []
def band(name, fmin, fmax, lam, label, keys, kind, mp=True, note='', curve=None):
    B.append(dict(name=name, fmin=float(fmin), fmax=float(fmax), lam=float(lam), label=label, sources=keys, kind=kind, mass_prop=mp, note=note, curve=curve, stage='35'))

# ---------------- (A) bulk-heating tail above 13 GHz (DERIVED) ----------------
# ccb18 Eq. 19 with a band-limited spectrum f(w) = 1 inside [w1,w2], 0 outside, and w = v_s q:
#   lambda_eff/lambda = Int_band x^4 e^{-x^2} dx / (3 sqrt(pi)/8),  x = q r_C
# so a release confined to the band is bounded by lambda_band <= lambda_bulk * (3 sqrt(pi)/8) / Int_band.
av = S['adler_vinante18']['values']; vs, lam_bulk = av['vs'], av['lam']
def logG(x):   # log of Int_x^inf t^4 e^{-t^2} dt = e^{-x^2}[(3 sqrt(pi)/8) erfcx(x) + x^3/2 + 3x/4]
    return -x * x + np.log(3 * sqrt(pi) / 8 * erfcx(x) + x ** 3 / 2 + 3 * x / 4)
NORM = 3 * sqrt(pi) / 8
xof = lambda f: 2 * pi * f * RC / vs
def lam_band(f1, f2):
    a, b = xof(f1), xof(f2)
    lg = logG(a) + np.log1p(-np.exp(logG(b) - logG(a)))
    return lam_bulk * NORM * np.exp(-lg)
# Stage-33 half-max band, read band-limited (consistency check of the Stage 33 convention):
f33 = (5554449066.265411, 12978178989.862732)
lam33_bandlimited = lam_band(*f33)
edges = 12978178989.862732 * 10 ** (np.arange(0, 1.01, 0.05))
tail = []
for f1, f2 in zip(edges[:-1], edges[1:]):
    tail.append((sqrt(f1 * f2), f1, f2, float(lam_band(f1, f2))))
tail = np.array(tail)
LAM_CAP = 1.0                          # ASSUMED reporting cap: count coverage only where lambda <= 1 s^-1
keep = tail[:, 3] <= LAM_CAP
f_tail_hi = tail[keep, 2].max()
f_cross = {}
for thr in (1e-8, 1e-6, 1e-3, 1.0):
    i = np.argmax(tail[:, 3] > thr); f_cross[f'{thr:g}'] = float(tail[i, 1]) if tail[i, 3] > thr else None
with open(os.path.join(ROOT, 'results', 'stage35_bulk_tail.csv'), 'w') as fh:
    fh.write('f_center_Hz,f_lo_Hz,f_hi_Hz,lambda_band_max_per_s\n')
    for r in tail: fh.write(f'{r[0]:.5e},{r[1]:.5e},{r[2]:.5e},{r[3]:.4e}\n')
with open(os.path.join(ROOT, 'results', 'stage35_bulk_tail_curve.csv'), 'w') as fh:
    fh.write('f_Hz,lambda_max,lambda_min\n')
    for r in tail[keep]:
        fh.write(f'{r[1]:.5e},{r[3]:.4e},{r[3]:.4e}\n{r[2]:.5e},{r[3]:.4e},{r[3]:.4e}\n')
band('Bulk-heating tail above 13 GHz, band-limited (0.05-dec bins, lambda<=1 shown)', 12978178989.862732, f_tail_hi, tail[0, 3],
     'DERIVED (ccb18 Eq. 19 band-limited; Adler-Vinante 10 pW/kg, v_s = 4000 m/s ASSUMED as Stage 33)', ['adler_vinante18', 'ccb18'], 'phonon heating',
     note=f'curve; {tail[0,3]:.2e} at 13 GHz rising as exp(x^2); 1 s^-1 reached at {f_cross["1"]:.3g} Hz', curve='stage35_bulk_tail_curve.csv')

# ---------------- (B) neutral matter: no CSL photons below ~12.4 keV (DERIVED) ----------------
E_atom = H * C / S['donadi21']['values']['d_atom'] / keV        # photon energy whose wavelength equals 0.1 nm
# Donadi Eq. (45)-(49): for |r_i - r_j| << lambda_k and << r_C the emission is coherent, amplitude ~ Sum q_i (m-p: masses cancel)
# -> neutral atom: Sum q_i = 0 -> leading-order rate 0 for lambda_k >> 0.1 nm, i.e. E << 12.4 keV.
# IGEX Stage-33 band starts at 14.5 keV: lambda_k = hc/E
lamk_145 = H * C / (14.5 * keV) * 1e9

# ---------------- (C) qubits: pair-breaking needs phonons >= 2 Delta ----------------
f_pb = 2 * S['serniak18']['values']['Delta_eV'] * E_ / H
x_pb = xof(f_pb)
supp_pb_log10 = (-x_pb ** 2 + 4 * np.log(x_pb)) / np.log(10)     # kernel x^4 e^{-x^2} at 2 Delta, relative to O(1)

# ---------------- (D) free charges in the diffuse ionized cosmic gas -> CXB and FIRAS (DERIVED) ----------------
p = S['planck18']['values']; pc = S['parsec']['values']['pc']
H100 = 100e3 / (1e6 * pc)                                       # 100 km/s/Mpc in s^-1
rho_b = p['Obh2'] * 3 * H100 ** 2 / (8 * pi * G)                # kg/m^3 today
n_b = rho_b / U                                                 # baryons per m^3 (DERIVED)
q2_per_baryon = 1.5                                             # e^2 units; ASSUMED conservative: all-helium (H gives 2.0)
f_ion = 0.5                                                     # ASSUMED: >= half of baryons in diffuse ionized gas (not fetched -> conservative choice)
H0 = p['H0'] * 1e3 / (1e6 * pc)
dt = 0.09 / H0                                                  # ASSUMED: only emission since z ~ 0.1 counted (< ln 1.1 / H0), redshift neglected
Cdim = HBAR * E_ ** 2 / (4 * pi ** 2 * EP0 * C ** 3 * M0 ** 2 * RC ** 2)   # E dGamma/dE per unit lambda per e^2 (dimensionless) - Majorana Eq. 2 / Donadi Eq. 44
I_E_per_lam = C / (4 * pi) * n_b * f_ion * q2_per_baryon * Cdim * dt      # W m^-2 sr^-1 per J of photon energy, per (s^-1 of lambda)
cgs_deg2_to_SI = 1e-3 * (180 / pi) ** 2
hk = S['hickox06']['values']; cxb = {}
for nm, (Ilo, Ihi), (I, dI) in (('1-2 keV', (1, 2), hk['I12']), ('2-8 keV', (2, 8), hk['I28'])):
    Iobs = (I + 2 * dI) * cgs_deg2_to_SI                         # all of the measured CXB attributed (+2 sigma)
    lam_max = Iobs / (I_E_per_lam * (Ihi - Ilo) * keV)
    cxb[nm] = dict(I_obs_SI=Iobs, lam=lam_max)
    band(f'CXB total {nm} vs free-charge emission of diffuse ionized gas (Hickox 2006)', Hz_of_eV(Ilo * 1e3), Hz_of_eV(Ihi * 1e3), lam_max,
         'DERIVED (Majorana Eq. 2 / Donadi Eq. 44 per free charge; f_ion = 0.5, z<0.1 ASSUMED)', ['hickox06', 'planck18', 'parsec', 'majorana22', 'donadi21', 'nist'], 'cosmic background',
         note='all measured CXB attributed (+2 sigma); free charges spaced >> lambda_k -> incoherent')
fx = S['fixsen96']['values']; T = fx['T']
nus = np.linspace(fx['nu_cm'][0], fx['nu_cm'][1], 2000) * 100 * C
Bnu = 2 * H * nus ** 3 / C ** 2 / np.expm1(H * nus / (KB * T)); Bpk = Bnu.max()
I_nu_per_lam = I_E_per_lam * H                                  # W m^-2 sr^-1 Hz^-1 per lambda (flat in nu)
lam_firas = 2 * fx['ppm'] * Bpk / I_nu_per_lam                  # ASSUMED: 2 x the 50 ppm weighted rms as the limit
band('COBE/FIRAS 2-21 cm^-1 residual vs free-charge emission (Fixsen 1996)', nus[0], nus[-1], lam_firas,
     'DERIVED (2 x 50 ppm of peak; same emission model)', ['fixsen96', 'planck18', 'parsec', 'majorana22', 'donadi21', 'nist'], 'cosmic background',
     note='weak: CMB is bright')

# ---------------- (E) LEE coincidence check ----------------
ex = S['excess22']['values']
alpha = log10(ex['R200'] / ex['R1000']) / log10(1000 / 200)     # local power law between the two quoted points (DERIVED)
Egrid = np.linspace(0.2, 1.0, 2001)                             # keV
R = ex['R200'] * (Egrid / 0.2) ** (-alpha)                      # counts/keV/kg/day
P_lee = np.trapezoid(Egrid * R, Egrid) * keV / 86400            # W/kg (0.2-1 keV)
phonon_quantum_eV = H * xof(13e9) / (2 * pi * RC) * vs / E_     # energy of a 13 GHz phonon (= h f)
phonon_quantum_eV = H * 13e9 / E_
n_quanta_30eV = 30.0 / (H * f_tail_hi / E_)

# ---------------- coverage (merge with Stage 34) ----------------
s34 = json.load(open(os.path.join(ROOT, 'results', 'stage34.json')))
Ball = [dict(b) for b in s34['bands']] + B
def union(iv):
    iv = sorted(iv); out = []
    for lo, hi in iv:
        if out and lo <= out[-1][1]: out[-1][1] = max(out[-1][1], hi)
        else: out.append([lo, hi])
    return out
def dec_in(Uv, lo, hi):
    return sum(max(0.0, log10(min(h_, hi)) - log10(max(l_, lo))) for l_, h_ in Uv if h_ > lo and l_ < hi)
mp = [b for b in Ball if b['mass_prop']]
U34 = union([(b['fmin'], b['fmax']) for b in mp if b['stage'] != '35'])
U35 = union([(b['fmin'], b['fmax']) for b in mp])
U35_le1 = union([(b['fmin'], b['fmax']) for b in mp if b['lam'] <= 1.0 and not (b['stage'] == '34')] + [(10.0, 1000.0)])  # LIGO only where lambda<=1 (<=~1 kHz)
GAP = tuple(s34['coverage']['gap']); SPAN = tuple(s34['coverage']['span'])
W = (12978178989.862732, 3.506084401023131e+18)
cov = dict(gap=GAP, gap_decades=log10(GAP[1] / GAP[0]), covered34=dec_in(U34, *GAP), covered35=dec_in(U35, *GAP), covered35_lam_le_1=dec_in(U35_le1, *GAP),
           span=SPAN, span_decades=log10(SPAN[1] / SPAN[0]), span_covered34=dec_in(U34, *SPAN), span_covered35=dec_in(U35, *SPAN),
           window=W, window_decades=log10(W[1] / W[0]), window_covered=dec_in(U35, *W), window_covered_lam_le_1=dec_in(U35_le1, *W))
cov['newly'] = cov['covered35'] - cov['covered34']; cov['open_gap'] = cov['gap_decades'] - cov['covered35']
open_w = []; prev = SPAN[0]
for lo, hi in U35:
    if lo > prev and log10(lo / prev) > 1e-3: open_w.append((prev, lo, log10(lo / prev)))
    prev = max(prev, hi)
if SPAN[1] > prev: open_w.append((prev, SPAN[1], log10(SPAN[1] / prev)))
cov['open_windows'] = open_w

res = dict(bulk_tail=dict(lam_13GHz=float(tail[0, 3]), f_cross=f_cross, f_hi_lam_le_1=float(f_tail_hi), stage33_band_bandlimited=float(lam33_bandlimited),
                          ratio_stage33=float(lam33_bandlimited / lam_bulk)),
           neutral_cutoff_keV=E_atom, lamk_at_14p5keV_nm=lamk_145,
           qubit=dict(f_pair_breaking_GHz=f_pb / 1e9, x=x_pb, log10_kernel=supp_pb_log10),
           cosmic=dict(n_b=n_b, dt_s=dt, Cdim=Cdim, I_E_per_lam=I_E_per_lam, cxb=cxb, firas=dict(Bpeak_SI=Bpk, lam=lam_firas)),
           lee=dict(alpha=alpha, P_0p2_1keV_W_per_kg=P_lee, ratio_to_bulk_residual=P_lee / av['P_resid'], phonon_quantum_13GHz_eV=phonon_quantum_eV,
                    n_quanta_for_30eV_at_tail_top=n_quanta_30eV))
dropped = [
    dict(item='Qubit quasiparticle / excess-energy rates', why=f'pair breaking in Al needs phonons >= 2Delta = {f_pb/1e9:.0f} GHz (Delta = 205 ueV, CITED); at r_C = 100 nm, x = q r_C = {x_pb:.1f}, kernel ~1e{supp_pb_log10:.0f}. No usable bound (DERIVED structural negative).'),
    dict(item='Low-threshold cryogenic detector spectra (SuperCDMS/CRESST/EDELWEISS/SENSEI/DAMIC), 30 eV - 10 keV', why=f'mass-proportional CSL at r_C = 100 nm makes no photons in neutral matter for E << {E_atom:.1f} keV (Donadi 2021) and no phonons above ~{f_cross["1"]/1e9:.0f} GHz with lambda <= 1; detector rates therefore bound nothing in this parametrization. Same reason the LEE cannot be this effect.'),
    dict(item='Optical / IR / UV spontaneous emission from lab samples', why='same neutral-atom cancellation (lambda_k >> 0.1 nm). Only free charges radiate; covered via cosmic backgrounds (very weak).'),
    dict(item='Optical/UV extragalactic background and 13.6 eV - 1 keV', why='not converted this stage (EBL numbers not fetched; 13.6 eV - 0.3 keV absorbed by Galactic H I). Would give lambda ~ O(1) at best by the same model.'),
]
flags = [
    f'Stage 33 suggested "low-threshold calorimeters staring at a cold massive sample (dGamma/dE ~ 1/E favours low energies)". That was wrong for neutral matter: below ~{E_atom:.1f} keV electrons and nucleus emit coherently and cancel (Donadi 2021). Corrected here.',
    f'IGEX band (Stage 33) starts at 14.5 keV where lambda_k = {lamk_145:.3f} nm, only just below the 0.1 nm atomic scale: the incoherent-electron assumption is marginal at its low edge.',
    f'Stage 33 assigned the white bulk-heating bound to its half-max band. Read band-limited (ccb18 Eq. 19) it is {lam33_bandlimited:.2e} ({lam33_bandlimited/lam_bulk:.2f}x weaker).',
    'LEE coincidence: NOT a candidate. (i) structural: no eV-keV deposits in this parametrization; (ii) EXCESS: rates/shapes differ between detectors, decay with time, jump after warm-ups above 10 K - a stationary universal noise cannot do that.',
    f'Energy is not the discriminator: the EDELWEISS RED20 excess carries only ~{P_lee:.1e} W/kg (0.2-1 keV), {P_lee/av["P_resid"]:.0e} of the 10 pW/kg bulk residual.',
]

with open(os.path.join(ROOT, 'results', 'stage35_sources.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for k, v in S.items(): w.writerow([k, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
with open(os.path.join(ROOT, 'results', 'stage35_bounds.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['name', 'kind', 'fmin_Hz', 'fmax_Hz', 'lambda_max_per_s_rC_1e-7m', 'label', 'mass_proportional', 'stage', 'sources', 'note'])
    for b in sorted(Ball, key=lambda x: x['fmin']):
        w.writerow([b['name'], b['kind'], f"{b['fmin']:.4g}", f"{b['fmax']:.4g}", f"{b['lam']:.3g}", b['label'], b['mass_prop'], b['stage'], ';'.join(b['sources']), b['note']])
fig_title = f"Stage 35 — 13 GHz–14.5 keV window ({cov['window_decades']:.1f} dec): formally covered {cov['window_covered']:.2f} dec, only {cov['window_covered_lam_le_1']:.2f} dec with λ ≤ 1 s$^{{-1}}$ (bulk-heating tail)"
out = dict(stage=35, title='13 GHz - 14.5 keV window', fig_title=fig_title, results=res, bands=Ball, coverage=cov, dropped=dropped, flags=flags, sources=S, runtime_s=time.time() - t0)
json.dump(out, open(os.path.join(ROOT, 'results', 'stage35.json'), 'w'), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
print(json.dumps({k: v for k, v in cov.items() if k != 'open_windows'}, indent=1, default=float))
for w_ in open_w: print(f"  open {w_[0]:.4g}-{w_[1]:.4g} Hz ({w_[2]:.2f} dec)")
for f in flags: print('FLAG', f)
