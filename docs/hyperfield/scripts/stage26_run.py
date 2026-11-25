#!/usr/bin/env python3
"""Stage 26 - what survives of Oppenheim-type classical-quantum (CQ) Newtonian gravity, and who can reach it.

Core (DERIVED here, cross-checked against Tilloy & Diosi, PRD 96, 104045 (2017), arXiv:1706.01856):
For a CQ / measurement-and-feedback Newtonian model whose classical field couples to the mass density smeared at
length l, with decoherence kernel D0 (Lindblad -1/2 int D0 [m(x),[m(y),rho]]) and a saturated trade-off, the
white noise seen by ANY probe E (X-ray emission, force noise, interferometric contrast) is
      N_E(kappa, l) = A_E(l) * (kappa + 1/(4 kappa))          for the Diosi-Penrose (DP) kernel D0 = kappa G/(hbar r),
where A_E(l) is what the plain Diosi model (kappa = 1) predicts for that probe.  kappa is the D0 (decoherence)
strength in Diosi units; the 1/(4 kappa) term is the classical field noise (D2 at saturation, ~1/kappa).
For any other kernel (any positive D0, translation-invariant or not) Cauchy-Schwarz gives N_E >= eta_K * A_E with
eta_K >= 1 (equality only for DP).  An experiment with bound B_E therefore allows kappa in [k-, k+] with
k+- = (u +- sqrt(u^2-1))/2, u = B_E/A_E, and EXCLUDES the smearing length l for every kernel if u_E(l) < 1.
Labels: DERIVED (computed here), CITED (number from a source), ASSUMED (our choice, stated).
"""
import json, os, time
import numpy as np
from math import pi, sqrt, log10
from scipy.special import erf, erfc
from scipy.integrate import quad

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
H = 6.62607015e-34; HBAR = H / (2 * pi); G = 6.6743e-11; U = 1.66053906660e-27; KB = 1.380649e-23
GH = G * HBAR
t0 = time.time()

# ------------------------------------------------------------------ k-space machinery (DERIVED)
KLOG = np.linspace(np.log(1e-3), np.log(1e22), 60000)   # ln k, k in 1/m
K = np.exp(KLOG); DLNK = KLOG[1] - KLOG[0]

def F2_sphere(x):
    """|form factor|^2 of a uniform sphere, F = 3(sin x - x cos x)/x^3; oscillation-averaged tail for x > 60."""
    out = np.empty_like(x)
    s = x < 1e-3; out[s] = 1.0
    m = (x >= 1e-3) & (x <= 60)
    xm = x[m]; out[m] = (3 * (np.sin(xm) - xm * np.cos(xm)) / xm ** 3) ** 2
    b = x > 60; out[b] = 4.5 / x[b] ** 4
    return out

def rho2(body):
    """|rho(k)|^2 / M^2 for a body: coherent continuum sphere + incoherent nuclei above the Debye wavenumber k_c."""
    M, R, N, m2r = body['M'], body['R'], body['N'], body['m2_ratio']
    kc = body['kc']
    coh = F2_sphere(K * R)
    inc = np.where(K > kc, m2r / N, 0.0)     # N <m^2> / M^2 = <m^2>/<m>^2 / N
    return coh + inc

def A_force(body, ell):
    """DP (kappa=1) two-sided per-axis acceleration PSD of the body's c.o.m. [m^2 s^-3]:
       4 pi G hbar * (1/3) int d^3k/(2pi)^3 e^{-k^2 l^2} |rho|^2 / M^2 = (2 G hbar/(3 pi)) int k^2 e^{-k^2 l^2} |rho/M|^2 dk."""
    r2 = rho2(body)
    return (2 * GH / (3 * pi)) * np.sum(K ** 3 * np.exp(-(K * ell) ** 2) * r2) * DLNK

def A_interf(body, ell):
    """DP (kappa=1) decoherence rate of an L/R superposition dx [1/s]: (G/hbar)(2/pi) int e^{-k^2 l^2}|rho|^2 (1-sinc k dx) dk."""
    r2 = rho2(body) * body['M'] ** 2
    dx = body['dx']; x = K * dx
    one_m_sinc = np.where(x < 1e-4, x ** 2 / 6, 1 - np.sin(x) / np.where(x == 0, 1, x))
    return (G / HBAR) * (2 / pi) * np.sum(K * np.exp(-(K * ell) ** 2) * r2 * one_m_sinc) * DLNK

def A_point(ell):
    """DP (kappa=1) per-axis acceleration PSD of any single point-like nucleus smeared at l: G hbar/(6 sqrt(pi) l^3)."""
    return GH / (6 * sqrt(pi) * ell ** 3)

def make_body(name, M, rho, comp, **kw):
    """comp: list of (atomic mass in u, number fraction). Returns body dict with derived N, R, <m^2>/<m>^2, k_c."""
    mm = np.array([c[0] for c in comp]); ff = np.array([c[1] for c in comp]); ff = ff / ff.sum()
    m1 = (mm * ff).sum() * U; m2 = (mm ** 2 * ff).sum() * U ** 2
    N = M / m1; R = (3 * M / (4 * pi * rho)) ** (1 / 3); n = rho / m1
    b = dict(name=name, M=M, rho=rho, N=N, R=R, m_atom=m1, m2_ratio=m2 / m1 ** 2, n=n, kc=(6 * pi ** 2 * n) ** (1 / 3))
    b.update(kw); return b

# ------------------------------------------------------------------ inputs
SI3N4 = [(28.0855, 3), (14.007, 4)]; AUPT = [(196.967, 0.73), (195.084, 0.27)]; SIO2 = [(28.0855, 1), (15.999, 2)]
NDFEB = [(144.24, 2), (55.845, 14), (10.81, 1)]
FORCE = [  # CITED Janse et al. PRR 6, 033076 (2024) Table I (sqrt S_a, one-sided ASD), masses; rho ASSUMED where noted
    make_body('Gisler 2022 SiN nanowire', 9.30e-15, 3100.0, SI3N4, sqrtSa=1.03e-6, mode='absolute', status='existing',
              note='rho ASSUMED 3100 (Stage 24); string treated as equal-volume sphere'),
    make_body('LISA Pathfinder 2018 (Armano)', 1.928, 1.928 / 0.046 ** 3, AUPT, sqrtSa=1.74e-15, mode='differential', status='existing',
              note='differential, two independent TMs (Gamma per TM = S/4); at L1 (D0(Phi_b) caveat); cube as equal-volume sphere'),
    make_body('aLIGO 2016 (Martynov)', 40.0, 2200.0, SIO2, sqrtSa=3.92e-14, mode='absolute', status='existing', note='Janse dagger (filled by Janse)'),
    make_body('Fuchs 2024 SC-levitated magnet', 4.0e-7, 7500.0, NDFEB, sqrtSa=1.25e-9, mode='absolute', status='existing', note='rho ASSUMED 7500'),
    make_body('Westphal 2021 mm torsion balance', 2.18e-4, 19300.0, [(196.967, 1)], sqrtSa=4.16e-9, mode='absolute', status='existing', note=''),
    make_body('Monteiro 2020 levitated ng sphere', 9.42e-13, 2000.0, SIO2, sqrtSa=9.32e-7, mode='absolute', status='existing', note='rho ASSUMED 2000'),
    make_body('LISA requirement (single TM)', 1.928, 1.928 / 0.046 ** 3, AUPT, sqrtSa=3.0e-15, mode='absolute', status='planned',
              note='CITED LISA Definition Study Report ESA-SCI-DIR-RP-002 (2024): < 3 fm s^-2 Hz^-1/2 per TM'),
]
for b in FORCE:
    Sa = b['sqrtSa'] ** 2
    b['Gamma_meas'] = Sa / 4 if b['mode'] == 'differential' else Sa / 2      # DERIVED: one-sided S -> two-sided per-TM Gamma
    b['FOM_check'] = Sa * b['N']

NA = [(22.98977, 1)]
RHO_NA = 2 * 22.98977 * U / 4.2906e-10 ** 3     # bcc Na (Stage 25)
INTERF = [
    make_body('Pedalino 2025 172 kDa Na cluster', 172e3 * U, RHO_NA, NA, dx=133e-9, tau=2 * 0.983 / 160, status='existing',
              note='CITED arXiv:2507.21211; lambda <= 1/tau (Stage 15 rule)'),
    make_body('Kovachy 2015 Rb atom, 54 cm, 2T=2.08 s', 86.909 * U, 1e3, [(86.909, 1)], dx=0.54, tau=2.08, status='existing',
              note='CITED Nature 528, 530 (2015); single atom, rho irrelevant'),
    make_body('Stage-22 MUSCLE 1.415 MDa Na (our design)', 1.415e6 * U, RHO_NA, NA, dx=0.8 * 133e-9, tau=2 * 0.983 / 19.6, status='our design',
              note='Stage 22'),
    make_body('MAGIS-100 Sr atom, 1 m, 9 s', 87.906 * U, 1e3, [(87.906, 1)], dx=1.0, tau=9.0, status='planned', note='CITED arXiv:2104.02835 (Stage 17)'),
    make_body('MAQRO baseline 1e10 u silica, dx=R, 100 s', 1e10 * U, 2200.0, SIO2, dx=None, tau=100.0, status='planned',
              note='CITED MAQRO research campaign QST 8, 014006 (2023): 1e10 u, R 100-180 nm, dx ~ R, ~100 s'),
    make_body('MAQRO max 1e11 u silica, dx=R, 100 s', 1e11 * U, 2200.0, SIO2, dx=None, tau=100.0, status='planned',
              note='CITED Kaltenbaek et al. EPJ QT 3, 5 (2016): masses up to 1e11 u; dx ~ R ASSUMED as in the baseline'),
    make_body('Pino 2018 Nb 1 um, dx 0.5 um, 0.596 s', 3.3e-14, 8570.0, [(92.906, 1)], dx=0.5e-6, tau=0.596, status='proposal',
              note='CITED QST 3, 025001 (2018), arXiv:1603.01553 (Stage 16)'),
    make_body('Stage-12 BMV diamond 5.6e-15 kg, 3.85 um, 2.5 s (our design)', 5.6234e-15, 3510.0, [(12.011, 1)], dx=3.852e-6, tau=2.5,
              status='our design', note='Stage 12/21'),
    make_body('Stage-16 BMV Pb 1e-14 kg, 2.69 um, 1 s (our design)', 1e-14, 11340.0, [(207.2, 1)], dx=2.6915e-6, tau=1.0,
              status='our design', note='Stage 16/21'),
]
for b in INTERF:
    if b['dx'] is None: b['dx'] = b['R']
    b['lam_max'] = 1.0 / b['tau']

# X-ray (CITED Stage 25): R0 bounds in Diosi convention (kappa = 1)
XRAY = [dict(name='Donadi 2021 Ge (Nat Phys 17, 74)', R0=0.54e-10 / (8 * pi) ** (1 / 3), A=72.63, status='existing'),
        dict(name='XENONnT 2026 Xe (PRL 136, 120201)', R0=4.5e-10, A=131.29, status='existing')]
def u_xray(x, ell):
    """u = Gamma_X / A_point; with a nucleon-resolved correction for l ~ nuclear size (nucleons Gaussian, sigma = r_rms/sqrt 3)."""
    A = x['A']; sig = 1.2e-15 * A ** (1 / 3) * sqrt(3 / 5) / sqrt(3)
    fnuc = (A / ell ** 3 + A * (A - 1) / (ell ** 2 + sig ** 2) ** 1.5) / (A ** 2 / ell ** 3)   # Gamma_resolved / Gamma_point
    return (ell / x['R0']) ** 3 / fnuc

# cold atoms (CITED Kovachy PRL 114, 143004 (2015): spontaneous heating 20 +- 30 pK/s for 87Rb, T = 50 pK)
COLD = dict(m=86.909 * U, dTdt95=(20 + 1.645 * 30) * 1e-12, T=50e-12)
COLD['Gamma_meas'] = KB * COLD['dTdt95'] / COLD['m']                     # DERIVED per-axis d<v^2>/dt
COLD['l_diffusive_min'] = HBAR / (COLD['m'] * sqrt(KB * COLD['T'] / COLD['m']))   # kicks hbar/l << m sigma_v
L_ISL = 52e-6   # CITED Lee et al. PRL 124, 101101 (2020): 1/r^2 holds to 52 um; ASSUMED to cap the smearing length

# published (Eq. 47/151) counting for interferometers, Stage 25 normalisation: lambda_pub(kappa) = kappa 2 G N M^2/(hbar R_lam)
def A_interf_pub(b, ell): return 2 * G * b['N'] * b['M'] ** 2 / (HBAR * (ell ** 2 * b['dx']) ** (1 / 3))

# ------------------------------------------------------------------ scan
ELL = np.logspace(-16, -4, 241)
res = dict(ell=ELL.tolist(), u={}, u_pub={})
for b in FORCE:
    A = np.array([A_force(b, e) for e in ELL]); res['u'][b['name']] = (b['Gamma_meas'] / A).tolist()
    b['A_plateau'] = (4 * pi / 3) * GH * b['rho'] / b['M']
for b in INTERF:
    A = np.array([A_interf(b, e) for e in ELL]); res['u'][b['name']] = (b['lam_max'] / A).tolist()
    res['u_pub'][b['name']] = (b['lam_max'] / np.array([A_interf_pub(b, e) for e in ELL])).tolist()
for x in XRAY:
    res['u'][x['name']] = [u_xray(x, e) for e in ELL]
uc = COLD['Gamma_meas'] / np.array([A_point(e) for e in ELL])
res['u']['Kovachy 2015 cold-atom heating'] = np.where(ELL >= COLD['l_diffusive_min'], uc, np.inf).tolist()
res['u_cold_naive'] = uc.tolist()
print('scan done %.1f s' % (time.time() - t0), flush=True)

status = {b['name']: b['status'] for b in FORCE + INTERF}
status.update({x['name']: x['status'] for x in XRAY}); status['Kovachy 2015 cold-atom heating'] = 'existing'
EXIST = [k for k, s in status.items() if s == 'existing']
XN = [x['name'] for x in XRAY]

def kpm(u):
    u = np.asarray(u, float); ok = u >= 1
    s = np.sqrt(np.where(ok, u ** 2 - 1, 0)); kp = (u + s) / 2
    return np.where(ok, 1 / (4 * kp), np.nan), np.where(ok, kp, np.nan)      # k- k+ = 1/4 (avoids cancellation)

def scenario(names, xray_mode, pub=False):
    """xray_mode: 'full' (X-ray sees D0 + field noise), 'lindblad' (Stage 25: X-ray sees D0 term only), 'none'."""
    others = [n for n in names if n not in XN]
    U_o = np.min([np.array(res['u_pub'][n] if (pub and n in res['u_pub']) else res['u'][n], float) for n in others], axis=0)
    UX = np.min([np.array(res['u'][n], float) for n in XN], axis=0)
    if xray_mode == 'full':
        Uall = np.minimum(U_o, UX); km, kp = kpm(Uall)
    elif xray_mode == 'lindblad':
        km, kp = kpm(U_o); kp = np.minimum(kp, UX); bad = ~(kp >= km); km = np.where(bad, np.nan, km); kp = np.where(bad, np.nan, kp)
        Uall = U_o
    else:
        Uall = U_o; km, kp = kpm(U_o)
    W = np.log10(kp / km)
    W = np.where(ELL <= L_ISL, W, np.nan)
    allowed = np.isfinite(W) & (W >= 0)
    lmin = float(ELL[allowed][0]) if allowed.any() else None
    # refine l_min by interpolation of log u (full / none) or of W (lindblad)
    return dict(kappa_minus=km.tolist(), kappa_plus=kp.tolist(), window_orders=W.tolist(), l_min_grid=lmin, u_min=Uall.tolist())

SC = {}
SC['S1_markov_full_existing'] = scenario(EXIST, 'full')
SC['S1pub_markov_full_existing_published_counting'] = scenario(EXIST, 'full', pub=True)
SC['S2_lindblad_only_xray_existing'] = scenario(EXIST, 'lindblad')
SC['S3_no_xray_existing'] = scenario(EXIST, 'none')
ALL = list(status.keys())
LIT = [k for k, s_ in status.items() if s_ in ('existing', 'planned', 'proposal')]
SC['S4_markov_full_existing_plus_planned'] = scenario(LIT, 'full')
SC['S5_markov_full_plus_planned_plus_our_designs'] = scenario(ALL, 'full')
SC['S2b_lindblad_only_xray_onEarth_absolute'] = scenario([n for n in EXIST if 'Pathfinder' not in n], 'lindblad')

def interp_cross(name_list_u, target=1.0):
    pass

# exact l_min for S1: XENONnT u_X = 1 -> l = R0 (point nucleus); check which experiment binds at each l
def binding(names, pub=False):
    arr = np.array([np.array(res['u_pub'][n] if (pub and n in res['u_pub']) else res['u'][n], float) for n in names])
    idx = np.argmin(arr, axis=0); return [names[i] for i in idx]
bind_S1 = binding(EXIST)

def at(sc, key, ell):
    i = int(np.argmin(np.abs(np.log(ELL) - np.log(ell)))); v = SC[sc][key][i]; return v
BENCH = [1e-15, 4.5e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 5e-5]
bench = []
for e in BENCH:
    i = int(np.argmin(np.abs(np.log(ELL) - np.log(e))))
    row = dict(ell=float(ELL[i]), binding_S1=bind_S1[i])
    for s in SC: row['window_' + s] = SC[s]['window_orders'][i]; row['kappa_minus_' + s] = SC[s]['kappa_minus'][i]; row['kappa_plus_' + s] = SC[s]['kappa_plus'][i]
    row['log10u'] = {n: float(np.log10(res['u'][n][i])) for n in ALL}
    bench.append(row)

# reach table: for each experiment, min log10 u over the S1 surviving range, and the l-range it would exclude
lo = SC['S1_markov_full_existing']['l_min_grid']
surv = (ELL >= lo) & (ELL <= L_ISL)
reach = {}
for n in ALL:
    u = np.array(res['u'][n], float); lu = np.log10(u)
    j = int(np.nanargmin(np.where(surv, lu, np.inf)))
    excl = surv & (u < 1)
    reach[n] = dict(status=status[n], best_log10u_in_surviving=float(lu[j]), at_ell=float(ELL[j]),
                    excludes_from=float(ELL[excl][0]) if excl.any() else None, excludes_to=float(ELL[excl][-1]) if excl.any() else None,
                    log10u_at_R0_XENONnT=float(np.interp(np.log(4.5e-10), np.log(ELL), lu)),
                    log10u_at_1um=float(np.interp(np.log(1e-6), np.log(ELL), lu)))
    if n in res['u_pub']:
        up = np.log10(np.array(res['u_pub'][n], float))
        reach[n]['published_counting_max_log10u_in_surviving'] = float(np.max(np.where(surv, up, -np.inf)))
        reach[n]['published_counting_log10u_at_1nm'] = float(np.interp(np.log(1e-9), np.log(ELL), up))

# what remains after the planned/proposed experiments (S4)
S4W = np.array(SC['S4_markov_full_existing_plus_planned']['window_orders'], float)
remaining_S4 = [(float(ELL[i])) for i in range(len(ELL)) if np.isfinite(S4W[i])]

# ------------------------------------------------------------------ kernel inefficiency eta (DERIVED)
# eta^2 = <a,D a><b,D^-1 b>/<a,b>^2 for a single Gaussian nucleus (std l): DP = 1 exactly.
def eta_ultralocal():   # D0 ~ delta(x-y): k-space; eta^2 = (int k^4 e)(int e)/(int k^2 e)^2 with e = exp(-k^2 l^2), l = 1
    f = lambda p: quad(lambda k: k ** p * np.exp(-k * k), 0, np.inf)[0]
    return sqrt(f(4) * f(0) / f(2) ** 2)
def eta_discrete():     # D0 ~ delta(x-y)/m(x) (Oppenheim ultra-local discrete), self-sourced only (lower bound)
    g = lambda r: np.exp(-r * r / 2) / (2 * pi) ** 1.5
    gp = lambda r: -r * g(r)
    phip = lambda r: (np.exp(-r * r / 2) * sqrt(2 / pi) * r - erf(r / sqrt(2))) / r ** 2   # d/dr [erf(r/sqrt2)/r]
    w = lambda f: quad(lambda r: 4 * pi * r * r * f(r), 1e-8, 40, limit=400)[0] / 3     # (1/3) int d^3r (angular avg of x-components)
    heat = w(lambda r: r * r * g(r)); diff = w(lambda r: g(r) * phip(r) ** 2); cross = w(lambda r: gp(r) * phip(r))   # g_prime^2/g = r^2 g
    return sqrt(heat * diff) / abs(cross), dict(heat=heat, diff=diff, cross=cross)
def eta_dp_check():
    g = lambda r: np.exp(-r * r / 2) / (2 * pi) ** 1.5
    # DP: heat ~ <grad m, (1/r) grad m>, diff ~ <grad phi, (-lap) grad phi> -> both equal the cross term in k-space
    f = lambda p: quad(lambda k: k ** p * np.exp(-k * k), 0, np.inf)[0]
    return sqrt(f(0) * f(0) / f(0) ** 2)
eta_u = eta_ultralocal(); eta_d, eta_d_parts = eta_discrete()
kern = dict(DP=dict(eta=eta_dp_check(), l_min_over_R0=1.0),
            ultralocal_continuous=dict(eta=eta_u, l_min_over_R0=eta_u ** (1 / 3), note='already excluded by Oppenheim Table 1'),
            ultralocal_discrete_self=dict(eta=eta_d, l_min_over_R0=eta_d ** (1 / 3), parts=eta_d_parts,
                                          note='self-sourced field noise only; environment adds -> eta is a lower bound'))
# Oppenheim discrete model with nuclear smearing: X-ray requirement on Dhat = (l_P^3/m_P) D2  (heating side alone)
def dhat_min_xray(x, ell):
    m = x['A'] * U; GX = GH / (6 * sqrt(pi) * x['R0'] ** 3)
    return HBAR ** 2 / (ell ** 2 * m * GX)
DCEIL = 4.00e-14   # CITED Stage 15: Eq. 46 ceiling with Gisler FOM
disc = {}
for x in XRAY:
    for lab, e in (('r_N=1e-15', 1e-15), ('Gaussian sigma of the nucleus', 1.2e-15 * x['A'] ** (1 / 3) * sqrt(3 / 5) / sqrt(3))):
        d = dhat_min_xray(x, e); disc[x['name'] + ' | ' + lab] = dict(ell=e, Dhat_min_kg=d, ceiling_kg=DCEIL, excluded_orders=log10(d / DCEIL))

# ------------------------------------------------------------------ headline numbers
def u_at(n, e): return float(np.interp(np.log(e), np.log(ELL), np.log10(np.array(res['u'][n], float))))
head = dict(
    l_min_S1_XENONnT=4.5e-10 * 1.0, l_min_S1_Donadi=XRAY[0]['R0'],
    l_min_S1_grid=SC['S1_markov_full_existing']['l_min_grid'],
    l_min_S2_grid=SC['S2_lindblad_only_xray_existing']['l_min_grid'],
    l_min_S3_grid=SC['S3_no_xray_existing']['l_min_grid'],
    exclusion_orders_l1em15_XENONnT_S1=-u_at(XN[1], 1e-15), exclusion_orders_l1em15_Donadi_S1=-u_at(XN[0], 1e-15),
    gisler_alone_log10u_l1em15=u_at('Gisler 2022 SiN nanowire', 1e-15),
    lpf_alone_l_kill=None,
    window_S1_at={('%.0e' % r['ell']): r['window_S1_markov_full_existing'] for r in bench},
    window_S4_at={('%.0e' % r['ell']): r['window_S4_markov_full_existing_plus_planned'] for r in bench},
    LPF_plateau_log10u=u_at('LISA Pathfinder 2018 (Armano)', 1e-7),
)
# where LPF u = 1 in the incoherent regime and where LPF overtakes X-ray
uL = np.array(res['u']['LISA Pathfinder 2018 (Armano)']); uX = np.array(res['u'][XN[1]])
head['lpf_alone_l_kill'] = float(np.exp(np.interp(0, np.log(uL[:120]), np.log(ELL[:120]))))
head['gisler_alone_l_kill'] = float(np.exp(np.interp(0, np.log(np.array(res['u']['Gisler 2022 SiN nanowire'])[:120]), np.log(ELL[:120]))))
cross_i = np.where(uL < uX)[0]; cross_i = cross_i[ELL[cross_i] > 1e-12]
head['lpf_overtakes_xray_at'] = float(ELL[cross_i[0]]) if len(cross_i) else None
# l-range closed by each planned experiment on top of existing (S4 vs S1)
S1W = np.array(SC['S1_markov_full_existing']['window_orders'], float)
closed = np.isfinite(S1W) & ~np.isfinite(S4W)
head['S4_closes_l_from'] = float(ELL[closed][0]) if closed.any() else None
head['S4_closes_l_to'] = float(ELL[closed][-1]) if closed.any() else None
head['S4_surviving_l_ranges'] = None
if np.isfinite(S4W).any():
    fin = np.isfinite(S4W); edges = np.diff(fin.astype(int)); starts = list(np.where(edges == 1)[0] + 1); ends = list(np.where(edges == -1)[0])
    if fin[0]: starts = [0] + starts
    if fin[-1]: ends = ends + [len(ELL) - 1]
    head['S4_surviving_l_ranges'] = [(float(ELL[a]), float(ELL[b])) for a, b in zip(starts, ends)]

# consistency vs Stage 25 (Eq. 47 normalisation with dropped constants)
s25 = json.load(open(os.path.join(ROOT, 'results/stage25.json')))['B']
corr = dict(stage25_l_min_consistent_Donadi=s25['l_min_consistent_Donadi'], stage25_l_min_consistent_XENONnT=s25['l_min_consistent_XENONnT'],
            stage25_excl_Donadi=s25['exclusion_orders_Donadi'], stage25_excl_XENONnT=s25['exclusion_orders_XENONnT'],
            stage26_S2_lmin_grid=SC['S2_lindblad_only_xray_existing']['l_min_grid'])
# Stage-25 claim check: X-ray vs Pedalino interferometry gap (orders, in kappa+ terms) across l
gapP = np.log10(np.array(res['u']['Pedalino 2025 172 kDa Na cluster'], float)) - np.log10(uX)
corr['xray_minus_pedalino_orders_min_1e-16_1e-8'] = float(np.min(gapP[(ELL >= 1e-16) & (ELL <= 1e-8)]))
corr['xray_minus_pedalino_orders_at'] = {('%.0e' % e): float(np.interp(np.log(e), np.log(ELL), gapP)) for e in (1e-15, 1e-12, 1e-10, 1e-9, 1e-8)}

# validity: rare-kick check at Gisler's kill length (kick hbar/l, rate per nucleus m^2 Gamma l^2/hbar^2 at kappa = 1/2)
g = FORCE[0]; lk = head['gisler_alone_l_kill']
kick_rate_body = g['N'] * (g['m_atom'] ** 2 * A_point(lk) * lk ** 2 / HBAR ** 2)
valid = dict(gisler_kill_l=lk, gisler_kicks_per_s_at_kill=kick_rate_body, kick_momentum_kg_m_s=HBAR / lk,
             cold_atom_diffusive_l_min=COLD['l_diffusive_min'],
             cold_atom_naive_l_min=float((GH / (6 * sqrt(pi) * COLD['Gamma_meas'])) ** (1 / 3)))

out = dict(stage=26, constants=dict(G=G, hbar=HBAR), inputs=dict(
    force=[{k: v for k, v in b.items()} for b in FORCE], interf=[{k: v for k, v in b.items()} for b in INTERF], xray=XRAY, cold=COLD, L_ISL=L_ISL),
    kernel=kern, discrete_oppenheim=disc, headline=head, corrections=corr, validity=valid, benchmarks=bench, reach=reach,
    scenarios=SC, scan=res, binding_S1=bind_S1,
    D2_map=dict(note='saturated l_P^2 D2 = hbar/(16 pi G kappa) (SQL normalisation, two-sided) or hbar/(2 G kappa) (Stage 25 Eq. 47 normalisation); ratio 8 pi',
                hbar_over_16piG=HBAR / (16 * pi * G), hbar_over_2G=HBAR / (2 * G)),
    sources=dict(tilloy_diosi_2017='https://arxiv.org/abs/1706.01856', oppenheim='https://www.nature.com/articles/s41467-023-43348-2',
                 janse='https://arxiv.org/abs/2403.08912', donadi='https://arxiv.org/abs/2111.13490', xenonnt='https://arxiv.org/abs/2506.05507',
                 lisa_redbook='https://arxiv.org/abs/2402.07571', maqro_2023='https://wigner.hu/~diosi/prints/2023qst8_014006.pdf',
                 maqro_2016='https://doi.org/10.1140/epjqt/s40507-016-0043-7', pino='https://arxiv.org/abs/1603.01553',
                 kovachy_lensing='https://arxiv.org/abs/1407.6995', bilardello='https://arxiv.org/abs/1605.01891',
                 angeli='https://arxiv.org/abs/2501.13030', lee_isl='https://doi.org/10.1103/PhysRevLett.124.101101',
                 pedalino='https://arxiv.org/abs/2507.21211', magis='https://arxiv.org/abs/2104.02835'),
    runtime_s=time.time() - t0)
json.dump(out, open(os.path.join(ROOT, 'results/stage26.json'), 'w'), indent=1, default=lambda o: None if (isinstance(o, float) and not np.isfinite(o)) else float(o))
print(json.dumps(dict(head=head, kern={k: v['eta'] for k, v in kern.items()}, corr=corr, valid=valid), indent=1, default=str))
for r in bench:
    print('%.1e  S1 %.2f  S2 %s  S3 %s  S4 %s  bind %s' % (r['ell'], r['window_S1_markov_full_existing'], r['window_S2_lindblad_only_xray_existing'],
          r['window_S3_no_xray_existing'], r['window_S4_markov_full_existing_plus_planned'], r['binding_S1']))
for n, v in reach.items(): print('%-62s %-10s best %.2f @ %.1e  excl %s..%s  R0 %.2f  1um %.2f' % (n, v['status'], v['best_log10u_in_surviving'], v['at_ell'], v['excludes_from'], v['excludes_to'], v['log10u_at_R0_XENONnT'], v['log10u_at_1um']))
for k, v in disc.items(): print(k, '%.3g kg excl %.2f' % (v['Dhat_min_kg'], v['excluded_orders']))
print('done %.1f s' % (time.time() - t0))
