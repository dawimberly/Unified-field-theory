#!/usr/bin/env python3
"""Stage 54 - the 8.2 kHz - 4.6 MHz resonant-release gap: new narrow-line bounds from cantilevers, trampolines,
levitated nanospheres and nanostrings. Parametrization unchanged: CSL-equivalent mass-proportional lambda per nucleon
(m0 = 1 u) at r_C = 1e-7 m, all measured (or, tier B, bath-implied) heating attributed to the effect.

Formulas (DERIVED, same linearised machinery as Stage 43: eta = lam (4 pi rC^2)^{3/2}/(m0^2 (2pi)^3) Int d^3k e^{-k^2 rC^2} |F(k)|^2,
heating Gamma = hbar eta /(2 M_eff w), phonons/s):
 (P) thin plate, lateral size >> rC (Stage 36):  G_plate = 4 pi hbar rho rC^2 (1 - e^{-h^2/4rC^2}) / (m0^2 h w)
 (B) long beam of rectangular cross-section w x h, mode wavelength >> rC (NEW here; separable k-integral):
       G_beam = G_plate * g(w),  g(w) = erf(w/2rC) - (2 rC/(sqrt(pi) w)) (1 - e^{-w^2/4rC^2});  g -> 1 for w >> rC (reduces to P exactly)
     A structure made of several beams is bounded below by its smallest-g element (Gamma/lam is a mass-weighted average of local
     values; in-phase cross terms between nearby elements are >= 0), so the min is conservative.
 (S) rigid homogeneous sphere radius R, mass m (centre-of-mass mode):
       G_sph = hbar m f(R/rC) / (4 m0^2 rC^2 w),  f(x) = (12/x^6) [e^{-x^2} - 1 + (x^2/2)(e^{-x^2} + 1)]  (f -> 1: point mass = Stage 33 ion formula)
     This is Zheng et al. 2020 Eq. (3) with S_CSL = hbar^2 eta (two-sided) and Gamma = S_CSL/(2 m hbar w). It is evaluated here both in
     closed form and by direct numerical k-integration (independent check), and validated against Zheng's published 10^-6.4.
Every input value was fetched this session; url, location and quote in results/stage54_sources.csv.
"""
import json, os, csv
import numpy as np
from math import pi, exp, sqrt, erf, log10
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
S = {}
def src(key, ref, url, loc, quote, **vals): S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)

src('nist', 'NIST CODATA 2022 (fetched this session: Value?hbar, Value?ukg, Value?k)', 'https://physics.nist.gov/cgi-bin/cuu/Value?hbar',
    'Numerical value fields', '"1.054 571 817... x 10 -34 J s" ; "1.660 539 068 92 x 10 -27 kg" ; "1.380 649 x 10 -23 J K -1"',
    hbar=1.054571817e-34, u=1.66053906892e-27, kB=1.380649e-23)
# ---------- validation source ----------
src('zheng20', 'Zheng et al., "Room temperature test of the CSL model using a levitated micro-oscillator", PRR 2, 013057 (2020), arXiv:1907.06896',
    'https://arxiv.org/pdf/1907.06896', 'p.2 Eqs. (1)-(3); p.3 materials; Table I',
    'Eq.(3) "eta = 6 rC^4 lam m^2/(m0^2 R^6) [e^{-R^2/rC^2} - 1 + R^2/(2rC^2)(e^{-R^2/rC^2} + 1)]" ; "<fCSL(t)fCSL(s)> = SCSL delta(t-s), where SCSL = hbar^2 eta" ; '
    '"upper bound of the CSL force as hbar^2 eta <= dS_total" ; "has radius R = 1.0 um, corresponding to 4.7 pg" ; '
    'Table I: "sqrt(dS_total) 1.3x10^-20 N/sqrt(Hz) ; sigma 3.3x10^-20 N/sqrt(Hz) ; lambda (rC = 1e-7 m) 10^-6.4 s^-1" ("standard deviations at 95% confidence level")',
    R=1.0e-6, m=4.7e-15, sqrtdS_95=3.3e-20, lam_pub_log10=-6.4)
# ---------- tier A: measured total heating / measured mode temperature / measured force noise ----------
src('heritier18', 'Heritier et al., "Nanoladder cantilevers made from diamond and silicon", Nano Lett. 18, 1814 (2018), arXiv:1711.11358',
    'https://arxiv.org/pdf/1711.11358', 'Fig.1 caption; Table I; p.3-4 text',
    '"Devices fabricated for this study used t ~ w = 100 to 300 nm" ; Table I "(diamond) fc (kHz) 25.22 m (pg) 4.1 +- 0.6 ... Q 162,000 gamma (pg/s) 3.7 This work" ; '
    '"we calculated the mode temperature T = m x_rms^2 (2 pi fc)^2/kB by integrating the spectrum" ; "The quality factor concurrently increased from Q = 131,000 to 162,000" ; '
    '"At the lowest temperature of T ~ 140 mK, the force noise reached a value of Fth = 190 +42 -33 zN in bandwidth B = 1 Hz"',
    f=25.22e3, m=4.1e-15, dm=0.6e-15, Q_lo=131000.0, Q=162000.0, T=0.140, Fth=190e-21, dFth_up=42e-21, w_min=100e-9, t_min=100e-9, w_max=300e-9)
src('diamond_rho', 'Wikipedia "Diamond" (raw wikitext, fetched this session)', 'https://en.wikipedia.org/w/index.php?title=Diamond&action=raw', 'infobox density',
    '"|density = {{val|3.5|-|3.53|ul=g/cm3}}" (lower end used: conservative)', rho=3500.0)
src('reinhardt16', 'Reinhardt, Mueller, Bourassa & Sankey, "Ultralow-noise SiN trampoline resonators", PRX 6, 021001 (2016), arXiv:1511.01769',
    'https://arxiv.org/pdf/1511.01769', 'p.1 abstract/text; p.2-3 text and Fig.3 caption; Fig.2 caption',
    '"a thickness of 44 nm ... The fundamental mode frequency wm = 2pi x 41.4 kHz, mass meff = 2.3 ng, and ringdown time tau_m = 285 seconds of this device correspond to a thermal force noise of SF = 16.2 aN/Hz^1/2 and root-mean-squared displacement xRMS = 161 pm. Letting the system evolve without drive for 38 hours, we observe xRMS,obs = 165 +- 5 pm" ; '
    '"force noise spectrum SF,obs ... consistent (within a ~5% systematic calibration error) with the expected 16.2 +- 0.8 aN/Hz^1/2" ; "room temperature (293 K)" ; '
    '"tau_m is its amplitude ringdown time" ; "Si3N4 parameters density 2700 kg/m3" ; "a = 2.1-um-wide tether" ; "d = 100 um-wide central pad"',
    f=41.4e3, h=44e-9, tau=285.0, T=293.0, xrms_obs=165e-12, dxrms=5e-12, xrms_exp=161e-12, rho=2700.0, w_tether=2.1e-6)
src('tebbenjohanns21', 'Tebbenjohanns, Mattana, Rossi, Frimmer & Novotny, "Quantum control of a nanoparticle optically levitated in cryogenic free space", Nature 595, 378 (2021), arXiv:2103.03853',
    'https://arxiv.org/pdf/2103.03853', 'p.2 text; p.3 text; Supplementary B and symbol table',
    '"spherical silica nanoparticle (diameter 100 nm, mass m ~ 1 fg)" ; "Omega_z/(2pi) = 77.6 kHz" ; "We estimate a diameter of (106 +- 5) nm" ; '
    '"a total decoherence rate of Gamma_tot/(2pi) = (5.5 +- 0.3) kHz" ; "Gamma_th Thermal decoherence rate (phonons/s)" ; "Gamma_tot Total decoherence rate: Gamma_tot = Gamma_qba + Gamma_exc"',
    f=77.6e3, d=106e-9, dd=5e-9, m_nom=1.0e-18, Gtot_2pi=5.5e3, dG_2pi=0.3e3)
src('magrini21', 'Magrini et al., "Real-time optimal quantum control of mechanical motion at room temperature", Nature 595, 373 (2021), arXiv:2012.15188',
    'https://arxiv.org/pdf/2012.15188', 'p.2 text; Methods D',
    '"silica nanosphere of 71.5 nm radius (~ 2.8 x 10^-18 kg)" ; "Omega_z/2pi = 104 kHz" ; "decoherence rates of the particle due to backaction, Gamma_ba = 2pi * 18.8 kHz, and due to residual gas molecules, Gamma_th = 2pi * 0.6 kHz" ; '
    '"We define the decoherence rates originating from thermal force noise and measurement backaction as the average rate of phonons delivered to the particle"',
    f=104e3, R=71.5e-9, m=2.8e-18, Gba_2pi=18.8e3, Gth_2pi=0.6e3)
src('delic20', 'Delic et al., "Cooling of a levitated nanoparticle to the motional quantum ground state", Science 367, 892 (2020), arXiv:1911.04406',
    'https://arxiv.org/pdf/1911.04406', 'p.1 abstract; p.2 text; Methods',
    '"143 nm diameter silica particle" ; "Omega_x/2pi ~ 305 kHz is the axial motional frequency" ; "we extract a total heating rate as low as Gamma_x/2pi = (20.6 +- 2.3) kHz at a pressure of ~ 10^-6 mbar" ; "m = 2.83 fg is the nanoparticle mass"',
    f=305e3, d=143e-9, m=2.83e-18, G_2pi=20.6e3, dG_2pi=2.3e3)
src('jain16', 'Jain et al., "Direct measurement of photon recoil from a levitated nanoparticle", PRL 116, 243601 (2016), arXiv:1603.03420',
    'https://arxiv.org/pdf/1603.03420', 'p.2-3 text',
    '"The specific mass density of silica is rho_SiO2 = 2200 kg/m3" ; "two different particles with radii R1 = 52.7 nm and R2 = 71.6 nm ... Omega0 = 2pi x 148.8 kHz and Omega0 = 2pi x 151.2 kHz ... We obtain Gamma1 = (20.9 +- 0.2) kHz and Gamma2 = (29.4 +- 0.3) kHz" ; '
    '"a R = 68.0 nm particle ... Omega0 = 2pi x 100.5 kHz and Omega0 = 2pi x 158.8 kHz with corresponding reheating rates Gamma1 = (19.4 +- 0.1) kHz and Gamma2 = (38.0 +- 0.3) kHz" ; "the heating rate. It defines the rate at which phonons" ; "Eq. (11) predicts a reheating rate of Gamma_recoil = 13.0 kHz" (s^-1 convention, no 2pi)',
    rho=2200.0)
# ---------- tier B: bath-implied heating kB T/(hbar Q) (mode temperature NOT measured) ----------
src('norte16', 'Norte, Moura & Groeblacher, "Mechanical resonators for quantum optomechanics experiments at room temperature", PRL 116, 147202 (2016), arXiv:1511.06235',
    'https://arxiv.org/pdf/1511.06235', 'p.3 text and Fig.3 inset',
    '"ringdown for our best device at 20 nm silicon nitride thickness, 900 um Si wafer thickness, and mechanical quality factor Qm = 9.8 +- 0.2 x 10^7 at a fundamental mode frequency of 140 kHz" ; title "... at Room Temperature" (T value not stated: 300 K ASSUMED)',
    f=140e3, h=20e-9, Q=9.8e7, T_assumed=300.0)
src('gisler22', 'Gisler et al., "Soft-clamped silicon nitride string resonators at millikelvin temperatures", PRL 129, 104301 (2022), arXiv:2112.03730',
    'https://arxiv.org/pdf/2112.03730', 'abstract; Fig.2 caption; Table I; p.3 text',
    'Table I Device C "Thickness 20 nm ... floc 1.406 MHz ... Q46m (x10^6) 2300" ; "a defect at the center with a width of 500 nm" ; "wmin = 500 nm and wmax = 1200 nm" ; '
    '"thermal decoherence time of 0.38 s" ; "Tplate = 46 mK" ; "Q is found to saturate, indicating good thermalization between the sample and the stage"',
    f=1.406e6, h=20e-9, w=500e-9, Q=2.3e9, T=0.046, tdec=0.38)
src('tsaturyan17', 'Tsaturyan et al., Nature Nanotech. 12, 776 (2017), arXiv:1608.00937', 'https://arxiv.org/pdf/1608.00937', 'p.3', '"(rho = 3200 kg/m3 is the density of SiN)"', rho=3200.0)
# ---------- checked and NOT converted ----------
src('tao14', 'Tao et al., "Single-crystal diamond nanomechanical resonators with quality factors exceeding one million", Nat. Commun. 5, 3638 (2014), arXiv:1212.1347',
    'https://arxiv.org/pdf/1212.1347', 'Table I; Methods',
    'Table I "el-SCD 240 12 660(20) 32,140 ... 1,510,000 ... 6.0 (0.54c)" "c At 93 mK" ; Methods: "Cantilever temperature could also be inferred from a thermomechanical noise measurement, but this method was inaccurate in our experiments due to slight mechanical vibrations introduced by the dilution circuit."',
    f=32.14e3)
src('beccari22', 'Beccari et al., "Strained crystalline nanomechanical resonators with quality factors above 10 billion", Nat. Phys. 18, 436 (2022), arXiv:2107.02124',
    'https://arxiv.org/pdf/2107.02124', 'abstract; p.3',
    '"the heating rate is only 60 quanta-per-second" ; "Omega_m/2pi = 1.46 MHz" ; "14 nm (nominal thickness) strained silicon" ; string width not stated in the fetched text', f=1.46e6)
src('si_rho', 'Wikipedia Template:Infobox silicon (raw wikitext, fetched this session)', 'https://en.wikipedia.org/w/index.php?title=Template:Infobox_silicon&action=raw', 'density', '"|density=2.329085 g/cm3"', rho=2329.085)
src('fogliano21', 'Fogliano et al., "Ultrasensitive nano-optomechanical force sensor at dilution temperatures", Nat. Commun. 12, 4124 (2021), arXiv:2009.02912',
    'https://arxiv.org/pdf/2009.02912', 'p.2', '"oscillating at 3.94, 5.84 and 11.6 kHz" ; noise temperature reported for nanowire A (3.94 kHz) only', f=11.6e3)
src('dania24', 'Dania et al., "Ultrahigh quality factor of a levitated nanomechanical oscillator", PRL 132, 133602 (2024), arXiv:2304.02408',
    'https://arxiv.org/pdf/2304.02408', 'p.2', '"Omega_z/2pi = 1.28 kHz" ; "heating rate Gamma_tot = 3.3(2) x 10^4 phonon/s" (out of band)', f=1.28e3)

N = S['nist']['values']; HBAR, M0, KB = N['hbar'], N['u'], N['kB']; RC = 1e-7

def G_plate(rho, h, f, rC=RC): return 4 * pi * HBAR * rho * rC ** 2 * (1 - exp(-h ** 2 / (4 * rC ** 2))) / (M0 ** 2 * h * 2 * pi * f)
def g_beam(w, rC=RC): return erf(w / (2 * rC)) - (2 * rC / (sqrt(pi) * w)) * (1 - exp(-w ** 2 / (4 * rC ** 2)))
def G_beam(rho, h, w, f, rC=RC): return G_plate(rho, h, f, rC) * g_beam(w, rC)
def f_sph(x):
    if x < 0.3:   # series of the bracket (avoids cancellation): 12/x^6 * (x^6/12 - x^8/24 + x^10/80 ...)
        y = x * x; return 1 - y / 2 + 3 * y * y / 20
    return (12 / x ** 6) * (exp(-x ** 2) - 1 + (x ** 2 / 2) * (exp(-x ** 2) + 1))
def eta_sph(m, R, rC=RC):          # closed form, per unit lambda (m^-2 s^-1)
    return m ** 2 * f_sph(R / rC) / (2 * M0 ** 2 * rC ** 2)
def eta_sph_num(m, R, rC=RC):      # direct k-space integration of Eq. (2) with the sphere form factor (independent check)
    k = np.linspace(1e-3 / rC, 12 / rC, 400001); kr = k * R
    F = 3 * (np.sin(kr) - kr * np.cos(kr)) / kr ** 3
    integ = (4 * pi / 3) * np.trapezoid(k ** 4 * (m * F) ** 2 * np.exp(-(k * rC) ** 2), k)
    return (4 * pi * rC ** 2) ** 1.5 / (M0 ** 2 * (2 * pi) ** 3) * integ
def G_sph(m, R, f, rC=RC): return HBAR * eta_sph(m, R, rC) / (2 * m * 2 * pi * f)

out = {}
# ---------------- validation 1: beam formula limits ----------------
out['beam_limits'] = dict(g_1um=g_beam(1e-6), g_10um=g_beam(1e-5), g_500nm=g_beam(500e-9), g_100nm=g_beam(100e-9), g_300nm=g_beam(300e-9), g_2p1um=g_beam(2.1e-6))
# numerical check of the separable y-integral against the closed form g(w)
def g_num(w, rC=RC):
    k = np.linspace(1e-6 / rC, 40 / rC, 2000001)
    return 2 * np.trapezoid(np.exp(-(k * rC) ** 2) * 4 * np.sin(k * w / 2) ** 2 / k ** 2, k) / (2 * pi * w)
out['beam_g_numeric_check'] = {f'{w*1e9:.0f}nm': [g_beam(w), float(g_num(w))] for w in (100e-9, 500e-9, 2e-6)}
# ---------------- validation 2: sphere closed form vs numeric, and Zheng 2020 published bound ----------------
z = S['zheng20']['values']
out['sphere_closed_vs_numeric'] = {f'R={R*1e9:.0f}nm': [eta_sph(1.0, R), float(eta_sph_num(1.0, R))] for R in (20e-9, 50e-9, 71.5e-9, 1e-6)}
lam_z = z['sqrtdS_95'] ** 2 / (HBAR ** 2 * eta_sph(z['m'], z['R']))
out['zheng_validation'] = dict(lam_ours=lam_z, log10_ours=log10(lam_z), log10_pub=z['lam_pub_log10'], ratio=lam_z / 10 ** z['lam_pub_log10'],
                               note='dS_total <= (sigma_95)^2 with S two-sided = hbar^2 eta (Zheng convention); review value 6.7e-7 (Stage 33) not reproduced; paper value is')
# ---------------- tier A bounds ----------------
rows = []
def add(name, f, heat, G, tier, label, note, srcs, extra=None):
    r = dict(name=name, f=f, heat=heat, G_per_lam=G, lam=heat / G, tier=tier, label=label, note=note, sources=srcs)
    if extra: r.update(extra)
    rows.append(r)
h = S['heritier18']['values']; rho_d = S['diamond_rho']['values']['rho']; w_ = 2 * pi * h['f']
heat_T = KB * h['T'] / (HBAR * h['Q_lo'])                                        # measured mode T, lower Q (conservative)
heat_F = (h['Fth'] + h['dFth_up']) ** 2 / (4 * (h['m'] - h['dm']) * w_ * HBAR)  # one-sided S_F = 4 m w kB T/Q -> kB T/(hbar Q) = S_F/(4 m w hbar)
G_h = min(G_beam(rho_d, t, w, h['f']) for t in (h['t_min'], 200e-9, h['w_max']) for w in (h['w_min'], 200e-9, h['w_max']))
add('Diamond nanoladder cantilever 25.22 kHz (Heritier 2018)', h['f'], max(heat_T, heat_F), G_h, 'A',
    'DERIVED (beam formula, min over t,w in 100-300 nm; heating from measured mode T ~140 mK / measured F_th +1 sigma, larger used)',
    f'heat from T,Q_lo = {heat_T:.3g}; from F_th(+1sig), m(-1sig) = {heat_F:.3g} s^-1; rails treated as independent beams (conservative)', ['heritier18', 'diamond_rho', 'nist'])
r = S['reinhardt16']['values']
Q_r = pi * r['f'] * r['tau']                    # amplitude ringdown: Q = w tau/2
heat_r = KB * r['T'] / (HBAR * Q_r) * ((r['xrms_obs'] + 2 * r['dxrms']) / r['xrms_exp']) ** 2
add('SiN trampoline 41.4 kHz, 44 nm (Reinhardt 2016)', r['f'], heat_r, G_beam(r['rho'], r['h'], r['w_tether'], r['f']), 'A',
    'DERIVED (beam formula with tether width 2.1 um = smallest element; heating = thermal x (observed/expected x_rms, +2 sigma)^2)',
    f'Q = pi f tau = {Q_r:.3g}; measured x_rms 165+-5 pm vs thermal 161 pm and measured S_F spectrum within 5% of thermal', ['reinhardt16', 'nist'])
t = S['tebbenjohanns21']['values']; rho_s = S['jain16']['values']['rho']
R_t = (t['d'] - t['dd']) / 2
m_t = min(t['m_nom'], rho_s * 4 / 3 * pi * R_t ** 3)
add('Levitated silica nanosphere z 77.6 kHz, cryogenic (Tebbenjohanns 2021)', t['f'], 2 * pi * (t['Gtot_2pi'] + 2 * t['dG_2pi']), G_sph(m_t, R_t, t['f']), 'A',
    'DERIVED (sphere formula; total measured decoherence rate +2 sigma, all attributed; m = min(quoted ~1 fg, rho V(R-1 sigma)))',
    f'm used {m_t:.3g} kg, R {R_t*1e9:.1f} nm', ['tebbenjohanns21', 'jain16', 'nist'])
g = S['magrini21']['values']
add('Levitated silica nanosphere z 104 kHz, room-T (Magrini 2021)', g['f'], 2 * pi * (g['Gba_2pi'] + g['Gth_2pi']), G_sph(g['m'], g['R'], g['f']), 'A',
    'DERIVED (sphere formula; measured backaction + gas decoherence rates, all attributed)', 'no uncertainties quoted on Gamma_ba, Gamma_th', ['magrini21', 'nist'])
for (Rj, fj, Gj, dGj) in ((68.0e-9, 100.5e3, 19.4e3, 0.1e3), (52.7e-9, 148.8e3, 20.9e3, 0.2e3), (71.6e-9, 151.2e3, 29.4e3, 0.3e3), (68.0e-9, 158.8e3, 38.0e3, 0.3e3)):
    mj = rho_s * 4 / 3 * pi * Rj ** 3
    add(f'Levitated silica nanosphere {fj/1e3:.1f} kHz, R = {Rj*1e9:.1f} nm (Jain 2016)', fj, Gj + 2 * dGj, G_sph(mj, Rj, fj), 'A',
        'DERIVED (sphere formula; measured reheating rate +2 sigma, all attributed; m = rho V)', 'Gamma quoted in s^-1 (no 2pi), cf. their Eq. (11) 13.0 kHz', ['jain16', 'nist'])
dl = S['delic20']['values']
add('Levitated silica nanosphere x 305 kHz (Delic 2020)', dl['f'], 2 * pi * (dl['G_2pi'] + 2 * dl['dG_2pi']), G_sph(dl['m'], dl['d'] / 2, dl['f']), 'A',
    'DERIVED (sphere formula; measured total heating +2 sigma, all attributed)', '', ['delic20', 'nist'])
# ---------------- tier B (bath-implied, ASSUMED mode at bath temperature) ----------------
rho_sin_lo = min(S['tsaturyan17']['values']['rho'], S['reinhardt16']['values']['rho'])
n = S['norte16']['values']
add('SiN tethered membrane 140 kHz, 20 nm, room T (Norte 2016)', n['f'], KB * n['T_assumed'] / (HBAR * n['Q']), G_plate(rho_sin_lo, n['h'], n['f']), 'B',
    'ASSUMED heating kB T/(hbar Q) with T = 300 K (mode T not measured); plate formula, rho = 2700 (lower of two fetched SiN values)',
    'tether width not in fetched text: plate limit ASSUMED', ['norte16', 'tsaturyan17', 'reinhardt16', 'nist'])
gi = S['gisler22']['values']
add('Soft-clamped SiN string 1.406 MHz, 20 nm, 46 mK (Gisler 2022)', gi['f'], 1 / gi['tdec'], G_beam(rho_sin_lo, gi['h'], gi['w'], gi['f']), 'B',
    'ASSUMED heating = 1/(thermal decoherence time) quoted by authors (bath T, Q); beam formula w = 500 nm',
    f'check: kB T/(hbar Q) = {KB*gi["T"]/(HBAR*gi["Q"]):.3g} s^-1', ['gisler22', 'tsaturyan17', 'reinhardt16', 'nist'])
for rr in rows: rr['lam'] = rr['heat'] / rr['G_per_lam']

# ---------------- coverage map of 8.2 kHz - 4.6 MHz ----------------
s43 = json.load(open(os.path.join(ROOT, 'results', 'stage43.json')))
B = []
for rr in rows:
    B.append(dict(name=rr['name'], fmin=rr['f'] * 0.999, fmax=rr['f'] * 1.001, lam=rr['lam'], label=rr['label'], sources=rr['sources'],
                  kind='levitated / resonator heating', mass_prop=True, note=f"tier {rr['tier']}; " + rr['note'], curve=None, stage='54'))
Ball = [dict(b) for b in s43['bands']] + B
LO, HI = 8.2e3, 4.6e6
inband = sorted([b for b in Ball if b['mass_prop'] and b['fmax'] > LO * 0.99 and b['fmin'] < HI * 1.01], key=lambda b: b['fmin'])
iv = []
for b in inband:
    lo, hi = max(b['fmin'], LO), min(b['fmax'], HI)
    if iv and lo <= iv[-1][1]: iv[-1][1] = max(iv[-1][1], hi)
    else: iv.append([lo, hi])
covered = sum(log10(hi / lo) for lo, hi in iv)
gaps = []; prev = LO
for lo, hi in iv:
    if lo > prev: gaps.append((prev, lo, log10(lo / prev)))
    prev = max(prev, hi)
if HI > prev: gaps.append((prev, HI, log10(HI / prev)))
def strong_gaps(th):
    pts = sorted([(b['fmin'] * 0 + (b['fmin'] + b['fmax']) / 2) for b in inband if b['lam'] <= th])
    edges = [LO] + pts + [HI]; return [(a, c, log10(c / a)) for a, c in zip(edges[:-1], edges[1:])]
cov = dict(band=[LO, HI], band_decades=log10(HI / LO), covered_decades=covered, n_lines=len(inband),
           gaps_over_0p1dec=[g_ for g_ in gaps if g_[2] > 0.1], largest_gap=max(gaps, key=lambda x: x[2]),
           gaps_between_points_le_1em7=strong_gaps(1e-7), gaps_between_points_le_1em5=strong_gaps(1e-5))
out['coverage_8k_4p6M'] = cov
out['inband_lines'] = [(b['name'], (b['fmin'] + b['fmax']) / 2, b['lam'], b['stage']) for b in inband]

# ---------------- write ----------------
os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
with open(os.path.join(ROOT, 'results', 'stage54_sources.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for kk, v in S.items(): wr.writerow([kk, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
with open(os.path.join(ROOT, 'results', 'stage54_bounds.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['name', 'f_Hz', 'tier', 'heating_per_s', 'Gamma_per_lambda', 'lambda_max_per_s_rC_1e-7m', 'label', 'note', 'sources'])
    for rr in sorted(rows, key=lambda x: x['f']):
        wr.writerow([rr['name'], f"{rr['f']:.6g}", rr['tier'], f"{rr['heat']:.4g}", f"{rr['G_per_lam']:.4g}", f"{rr['lam']:.3g}", rr['label'], rr['note'], ';'.join(rr['sources'])])
with open(os.path.join(ROOT, 'results', 'stage54_coverage_8k_4p6M.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['f_center_Hz', 'lambda_max', 'stage', 'name'])
    for nm, fc, lam, st in out['inband_lines']: wr.writerow([f'{fc:.6g}', f'{lam:.3g}', st, nm])
json.dump(dict(stage=54, results=rows, validations=out, bands=Ball, coverage=cov, sources=S), open(os.path.join(ROOT, 'results', 'stage54.json'), 'w'), indent=1, default=float)
print(json.dumps({k: out[k] for k in ('beam_limits', 'beam_g_numeric_check', 'sphere_closed_vs_numeric', 'zheng_validation')}, indent=1, default=float))
for rr in sorted(rows, key=lambda x: x['f']): print(f"{rr['f']/1e3:9.2f} kHz  tier {rr['tier']}  heat {rr['heat']:.3g}  G/lam {rr['G_per_lam']:.3g}  lam <= {rr['lam']:.2g}  {rr['name']}")
print('coverage', json.dumps(cov, default=float, indent=0))
