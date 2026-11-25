#!/usr/bin/env python3
"""Stage 29 - stress test of the three 'possibly new' Stage 27 claims.
 (1) volume law  Gamma*V = (4 pi/3) G hbar (kappa + 1/4kappa)
 (2) cross-correlated force noise  -2 G hbar/d^3 (longitudinal), +G hbar/d^3 (transverse), rho = -(3/2pi) V/d^3
 (3) source-mass modulation for the matter-sourced (Oppenheim 'discrete') kernel; LPF => l >= 2.5 um
Independent re-derivation (no Stage 27 code reused), numerical checks, literature verdicts.
Labels: DERIVED (computed here), CITED (number/statement from a source), ASSUMED (our choice, stated).
Conventions: Diosi master equation  d rho/dt = -(G/2hbar) int int [m(x),[m(y),rho]]/|x-y|  (Diosi 1989; NHH 2014 Eq. S-DP generator;
Tilloy-Diosi 2017 Eq. 10/14 with D = G/|x-y|).  CQ family (Stage 26): D0 = kappa G/(hbar r) back-action + classical-field noise
<dm dm> = D0^-1/4  -> total kernel G hbar (kappa + 1/(4 kappa))/r in acceleration units.  Gamma = two-sided acceleration PSD per axis
(<a(t)a(t')> = Gamma delta(t-t')); one-sided = 2 Gamma.
"""
import json, os, time
import numpy as np
from math import pi, sqrt, log10
from scipy.integrate import quad, dblquad
from scipy.special import erf, spherical_jn
from scipy.optimize import brentq
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
t0 = time.time()
H = 6.62607015e-34; HBAR = H / (2 * pi); G = 6.6743e-11; U = 1.66053906660e-27; KB = 1.380649e-23
GH = G * HBAR
rng = np.random.default_rng(29)
out = dict(stage=29, constants=dict(G=G, hbar=HBAR, G_hbar=GH))

# ---------------------------------------------------------------- 0. unit bookkeeping (DERIVED, exponent vectors [kg, m, s])
UN = dict(G=np.array([-1, 3, -2]), hbar=np.array([1, 2, -1]), V=np.array([0, 3, 0]), d=np.array([0, 1, 0]),
          rho=np.array([1, -3, 0]), M=np.array([1, 0, 0]), Dhat=np.array([1, 0, 1]), Gamma=np.array([0, 2, -3]))
def dim(**p): return sum(UN[k] * v for k, v in p.items())
units = dict(
    GammaV_vs_Ghbar=bool(np.all(dim(Gamma=1, V=1) == dim(G=1, hbar=1))),                 # m^5 s^-3
    Ghbar_over_d3_is_Gamma=bool(np.all(dim(G=1, hbar=1, d=-3) == UN['Gamma'])),
    discrete_self_pi_rho_G2_Dhat_over_R=bool(np.all(dim(rho=1, G=2, Dhat=1, d=-1) == UN['Gamma'])),
    discrete_source_G2_Dhat_M_over_r4=bool(np.all(dim(G=2, Dhat=1, M=1, d=-4) == UN['Gamma'])),
    Dhat_units='kg s  (from <dm dm> = (Dhat/4) rho delta^3 delta_t; same as Oppenheim et al. "(l_P^3/m_P) D2 ... kg s"). '
               'Stage 26/27 JSON keys and text label it "kg" - label error only, numbers unaffected.')
out['units'] = units

# ================================================================ 1. VOLUME LAW
# Rigid body, CM displacement X:  m~(k) -> e^{-ik.X} mu~(k).  Expanding the Diosi double commutator to 2nd order in X:
#   Lambda_ij = (G/2hbar) int d^3k/(2pi)^3 (4pi/k^2) k_i k_j |mu~|^2 ;  d<p_i p_j>/dt = 2 hbar^2 Lambda_ij
#   => Gamma_ij = 4 pi G hbar  int d^3k/(2pi)^3 (k_i k_j/k^2) |mu~/M|^2  = 4 pi G hbar N_ij  /V_eff
# with V_eff = M^2 / int rho^2 (Parseval) and N_ij = <k_i k_j/k^2>_{|mu~|^2} = the magnetostatic DEMAGNETIZATION tensor of the
# body shape (trace 1).  Sphere / cube: N = 1/3 -> Gamma V = (4 pi/3) G hbar.  Field-noise term (1/4kappa) has the same 1/k^2
# kernel (C_Phi = pi G hbar/(kappa k^2)), so the CQ total is the same structure times (kappa + 1/4kappa).
def sphere_FT(k, R): x = k * R; return np.where(x < 1e-6, 1.0, 3 * (np.sin(x) - x * np.cos(x)) / np.maximum(x, 1e-300) ** 3)
def gamma_sphere(R, ell=0.0):
    # isotropic: Gamma = 4 pi G hbar * (1/3) * (1/2pi^2) int k^2 dk |F|^2 e^{-k^2 l^2}
    f = lambda k: k * k * sphere_FT(k, R) ** 2 * np.exp(-(k * ell) ** 2)
    kmax = 3000.0 / R if ell == 0 else 9.0 / ell          # e^{-81} cut for l > 0
    edges = np.concatenate([[0], np.logspace(np.log10(0.01 / R), np.log10(kmax), 800)])
    I = sum(quad(f, a, b, limit=200)[0] for a, b in zip(edges[:-1], edges[1:]))
    if ell == 0:   # analytic tail |F|^2 ~ 9 cos^2/(kR)^4 -> int k^2 9/(2 (kR)^4) dk from K: 9/(2 R^4 K)
        I += 9 / (2 * R ** 4 * edges[-1])   # (averaged cos^2 -> 1/2)
    return 4 * pi * GH * I / (3 * 2 * pi ** 2)
R1 = 1e-6; V1 = 4 / 3 * pi * R1 ** 3
vol = {}
vol['GammaV_sphere_numeric'] = gamma_sphere(R1) * V1
vol['GammaV_claim_(4pi/3)Ghbar'] = 4 * pi / 3 * GH
vol['ratio_numeric_over_claim'] = vol['GammaV_sphere_numeric'] / vol['GammaV_claim_(4pi/3)Ghbar']
vol['Diosi_form_check'] = dict(note='Gamma = 2 hbar^2 Lambda/M^2 with Lambda = G M^2/(2 hbar R^3) (Diosi decoherence rate tau^-1 = (G M^2/2 hbar R^3) dx^2) '
                                    '-> Gamma = G hbar/R^3 = (4pi/3) G hbar/V', G_hbar_over_R3_times_V=GH / R1 ** 3 * V1)
# finite smearing l/R (sphere): fractional deficit
lr = np.array([1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0])
vol['finite_l_over_R'] = dict(l_over_R=lr.tolist(), Gamma_ratio_to_plateau=[gamma_sphere(1.0, x) / gamma_sphere(1.0) for x in lr])
# shape dependence: spheroid with semi-axes (a,a,c); mu~ = M F(q), q = sqrt(a^2 k_perp^2 + c^2 k_z^2); substitute k' -> 2D integral
def N_spheroid_numeric(c_over_a):
    a, c = 1.0, c_over_a
    # in scaled variables u = (a k_perp, c k_z) = q (sin t, cos t): k_z^2/k^2 = (cos^2/c^2)/(sin^2/a^2 + cos^2/c^2)
    # int d^3k |F(q)|^2 w = (1/(a^2 c)) int d^3u |F(|u|)|^2 w(t)  -> radial part cancels in the ratio N = <w>
    w = lambda t: (np.cos(t) ** 2 / c ** 2) / (np.sin(t) ** 2 / a ** 2 + np.cos(t) ** 2 / c ** 2)
    num = quad(lambda t: w(t) * np.sin(t), 0, pi)[0]; den = quad(lambda t: np.sin(t), 0, pi)[0]
    return num / den
def N_spheroid_analytic(c_over_a):   # standard demagnetizing factor along symmetry axis (Osborn 1945)
    m = c_over_a
    if abs(m - 1) < 1e-9: return 1 / 3
    if m > 1:
        e = sqrt(1 - 1 / m ** 2); return (1 - e ** 2) / e ** 3 * (np.arctanh(e) - e)
    e = sqrt(1 / m ** 2 - 1); return (1 + e ** 2) / e ** 3 * (e - np.arctan(e))
shape = []
for m in (0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0):
    Nz = N_spheroid_numeric(m); Nx = (1 - Nz) / 2
    shape.append(dict(c_over_a=m, N_axis_numeric=Nz, N_axis_analytic=N_spheroid_analytic(m),
                      GammaV_axis_over_claim=3 * Nz, GammaV_transverse_over_claim=3 * Nx))
vol['shape_spheroids'] = shape
vol['shape_note'] = ('Gamma_ii V = 4 pi G hbar (kappa+1/4kappa) N_ii. Sphere and cube (cubic symmetry -> N = 1/3 exactly) obey the Stage 27 law; '
                     'long rod: transverse x1.5, axial -> 0; thin plate: normal x3, in-plane -> 0. LPF (cubes), MAQRO/levitated (spheres), '
                     'Westphal (spheres) unaffected; a nanowire (Gisler) read transversely is x1.5 (+0.18 orders) above the isotropic estimate.')
# cube: cubic symmetry forces N_xx = N_yy = N_zz and trace N = 1 -> N = 1/3 exactly (no numerics needed).
vol['cube_N'] = 1 / 3
# inhomogeneous bodies: V -> V_eff = M^2/int rho^2 <= V (Cauchy-Schwarz). Example: Au shell (1 um thick) on 10 um silica sphere
def Veff_coated(Rc, t, rc, rs):
    Vc = 4 / 3 * pi * Rc ** 3; Vs = 4 / 3 * pi * ((Rc + t) ** 3 - Rc ** 3); M = rc * Vc + rs * Vs
    return M ** 2 / (rc ** 2 * Vc + rs ** 2 * Vs), Vc + Vs
ve, vt = Veff_coated(10e-6, 1e-6, 2200.0, 19300.0)
vol['inhomogeneous_example'] = dict(body='silica sphere R=10 um with 1 um Au shell (ASSUMED example)', V_eff=ve, V=vt,
                                    Gamma_over_volume_law=vt / ve,
                                    note='mass/composition independence holds only for homogeneous bodies; for composites Gamma ~ int rho^2 / M^2')
# lattice correction at l >= 4.5 A: plateau needs l >> lattice spacing. relative correction ~ sum_{G != 0} e^{-G^2 l^2} (|S_G|^2 weighted)
lat = {}
for name, a, shells in (('Au fcc a=4.078 A', 4.078e-10, [(8, sqrt(3)), (6, 2.0)]), ('Si diamond a=5.431 A', 5.431e-10, [(8, sqrt(3)), (12, sqrt(8))]),
                        ('Pt-Au LPF ~fcc a=3.95 A', 3.95e-10, [(8, sqrt(3)), (6, 2.0)])):
    for ell in (4.5e-10, 1e-9):
        corr = sum(mult * np.exp(-(2 * pi * s / a * ell) ** 2) for mult, s in shells)
        lat[f'{name} | l={ell:.1e}'] = float(corr)
vol['lattice_correction_relative'] = lat
# background scaling exponents (DERIVED from standard formulas, ASSUMED models)
def S_gas(R, P, T, rho, mgas=28 * U):   # molecular-regime gas damping: S_F = 4 kT * gamma_m, gamma_m ~ P * A * sqrt(8 m/(pi k T)) (order-unity factor dropped)
    A = pi * R ** 2; S_F = 4 * KB * T * P * A * sqrt(8 * mgas / (pi * KB * T)); M = rho * 4 / 3 * pi * R ** 3; return S_F / M ** 2
def S_bb(R, T, rho):                     # thermal photon recoil: S_F ~ (sigma T^4 4 pi R^2) * (k T)/c^2 (order-unity dropped)
    sig = 5.670374419e-8; c = 2.99792458e8; S_F = sig * T ** 4 * 4 * pi * R ** 2 * KB * T / c ** 2; M = rho * 4 / 3 * pi * R ** 3; return S_F / M ** 2
def slope(f, x1, x2): return log10(f(x2) / f(x1)) / log10(x2 / x1)
vol['scaling_exponents_numeric'] = dict(
    CQ_R=slope(lambda R: 4 * pi / 3 * GH / (4 / 3 * pi * R ** 3), 1e-6, 2e-6),
    gas_R=slope(lambda R: S_gas(R, 1e-6, 300, 2200), 1e-6, 2e-6), gas_P=slope(lambda P: S_gas(1e-6, P, 300, 2200), 1e-6, 2e-6),
    gas_T=slope(lambda T: S_gas(1e-6, 1e-6, T, 2200), 100, 200), gas_rho=slope(lambda r: S_gas(1e-6, 1e-6, 300, r), 2000, 4000),
    bb_R=slope(lambda R: S_bb(R, 300, 2200), 1e-6, 2e-6), bb_T=slope(lambda T: S_bb(1e-6, T, 2200), 100, 200),
    CQ_rho=0.0)
vol['literature'] = dict(
    diosi_2007='arXiv:quant-ph/0607110 (J. Phys. A 40, 2989), Eqs. 8-9, 12-13: "omega_G^2 = (4 pi/3M) G int f^2(r) dr ... At constant mass density '
               'omega_G^2 = (4pi/3) G fbar = G M/R^3"; decoherence time "2 hbar R^3/(G M^2) 1/|x-y|^2". => Gamma = hbar omega_G^2/M = (4pi/3) G hbar/V '
               '(the int f^2 form IS the V_eff generalisation). Original: Diosi, PRA 40, 1165 (1989).',
    diosi_2013='arXiv:1302.5365, Eq. 9-10: omega_G = sqrt(4 pi G rho/3) "does not depend on the size but on the density"; and proves the c.o.m. rate is '
               'insensitive to the spread of nuclear positions (our lattice check).',
    nhh_2014='arXiv:1405.2868 (Nimmrichter, Hornberger, Hammerer PRL 113, 020405), Eq. 11: "D_DP = (G hbar/2pi^2) int d^3k (k_x^2/k^2)|rho~(k)|^2" - '
             'the general shape-dependent form (our demagnetization-tensor reading); the G=0 term is the volume law.',
    oppenheim_2023='arXiv:2203.01982 (Nat Commun 14, 7910), Sec. V discussion: "the variation in acceleration depends on 1/sqrt(N) the number of nucleons '
                   'in the test mass, so that masses with smaller volume should yield larger uncertainty and this would be the effect to look for" - '
                   'volume scaling already proposed as the CQ signature.',
    tilloy_diosi_2017='arXiv:1706.01856 Eq. 11 D = gamma/4 + V o gamma^-1 o V^T; Eq. 13-14 PLD gamma = -2V -> D = G/|x-y| (Diosi kappa=1). '
                      'Confirms the Stage 26 (kappa + 1/4kappa) structure and min at kappa=1/2.')
out['claim1_volume_law'] = vol

# ================================================================ 2. CROSS-CORRELATION
# Cross term of the same double commutator, bodies A,B:  d<p_Ai p_Bj>/dt = 2 hbar^2 Lambda^AB_ij,
#   Lambda^AB_ij = -(1/2) int int mu_A mu_B d_i d_j D(r) ,  D = kappa G/(hbar r)  (+ field-noise term with 1/(4 kappa))
# => Gamma^AB_ij = -G hbar (kappa + 1/4kappa) d_i d_j (1/d) = G hbar (kappa+1/4kappa) (delta_ij - 3 n_i n_j)/d^3   (shell theorem, spheres)
# Fourier check: Gamma^AB_zz = 4 pi G hbar int d^3k/(2pi)^3 (k_z^2/k^2) F_A F_B cos(k_z d)
def cross_fourier(d, R, ell=0.0, axis='z'):
    # = 4 pi G hbar /(2 pi)^2 int_0^inf k^2 dk F(k)^2 e^{-k^2 l^2} int_{-1}^{1} du w(u) cos(k d u),  w = u^2 (long) or (1-u^2)/2 (trans)
    def ang(k):
        if axis == 'z':   # int u^2 cos(xu) du = 2[(x^2-2) sin x + 2x cos x]/x^3
            x = k * d; return 2 * ((x * x - 2) * np.sin(x) + 2 * x * np.cos(x)) / x ** 3
        x = k * d; s0 = 2 * np.sin(x) / x; s2 = 2 * ((x * x - 2) * np.sin(x) + 2 * x * np.cos(x)) / x ** 3; return 0.5 * (s0 - s2)
    f = lambda k: k * k * sphere_FT(k, R) ** 2 * np.exp(-(k * ell) ** 2) * ang(k)
    edges = np.linspace(1e-9, 400.0 / R, 4001)
    I = sum(quad(f, a, b, limit=100)[0] for a, b in zip(edges[:-1], edges[1:]))
    return 4 * pi * GH * I / (2 * pi) ** 2
cr = {}
for dR in (2.5, 3.0, 5.0):
    R = 1.0; d = dR * R
    gl = cross_fourier(d, R); gt = cross_fourier(d, R, axis='x')
    cr[f'd={dR}R'] = dict(long_numeric_over_Ghbar_d3=gl / (GH / d ** 3), trans_numeric_over_Ghbar_d3=gt / (GH / d ** 3),
                          rho_long_numeric=gl / gamma_sphere(R), rho_long_claim=-3 * (4 / 3 * pi * R ** 3) / (2 * pi * d ** 3))
out_cross = dict(fourier_checks=cr)
out_cross['general_kappa'] = ('Gamma_AB(long) = -2 G hbar (kappa + 1/4kappa)/d^3: |Gamma_AB| >= 2 G hbar/d^3 with equality only at kappa = 1/2. '
                              'Stage 27 / Dan-facing text "independent of the model\'s free parameters" is true only for the correlation COEFFICIENT, '
                              'and only inside the one-parameter DP family (all three TD terms share the 1/k^2 kernel).')
# point-like bodies (R << l): kernel erf(r/2l)/r ; find zero crossing and far-field independently
def phi_dd_long(r, l):   # -d^2/dz^2 [erf(r/2l)/r] along the axis, numerically by finite differences (independent of Stage 27's analytic)
    h = 1e-4 * l; f = lambda x: erf(x / (2 * l)) / x
    return -(f(r + h) - 2 * f(r) + f(r - h)) / h ** 2
self_pt = 1 / (6 * sqrt(pi))   # -d_z^2 at r->0 (l=1)
zc = brentq(lambda r: phi_dd_long(r, 1.0), 0.5, 5.0)
out_cross['point_like'] = dict(self_term_check=float(-(lambda h: (erf(h / 2) / h - 1 / sqrt(pi)) * 2 / h ** 2)(1e-3) * 1.0), self_claim=self_pt,
                               zero_crossing_d_over_l=zc, far_field_ratio_at_d_50l=phi_dd_long(50.0, 1.0) / self_pt / (-12 * sqrt(pi) / 50 ** 3),
                               rho_at_d_0p1l=phi_dd_long(0.1, 1.0) / self_pt)
# other kernels (structure only): KTM LOCC (arXiv:1401.0946 Eq. 16) has NO [x1,[x2,rho]] term -> rho = 0.
# TD with local monitoring gamma = gamma0 delta: back-action gamma/4 is local -> zero cross for non-overlapping bodies; field term
# V o gamma^-1 o V ~ int e^{ikr}/k^4 = -r/(8 pi) -> cross ~ d_i d_j r = (delta_ij - n_i n_j)/r : LONGITUDINAL = 0, transverse ~ 1/d.
out_cross['other_kernels'] = dict(
    KTM_LOCC='rho = 0 (independent noises; Eq. 16 of arXiv:1401.0946 contains only [x_k,[x_k,rho]]). The KTM self-noise magnitude hbar K = 2 G hbar m1 m2/d^3 '
             'equals |Gamma_AB| M_A M_B numerically - COINCIDENCE of magnitude (both are the tidal constant K), not the same physics.',
    TD_local_gamma='longitudinal cross = 0, transverse ~ G^2/(gamma0 d): rho depends on gamma0 and on axis differently',
    DP_family='rho_long = -(3/2pi) V/d^3 for equal spheres, kappa-independent')
out_cross['geometries'] = dict(LPF=-3 * 0.046 ** 3 / (2 * pi * 0.376 ** 3), touching=-0.25, d3R=-3 * (4 / 3 * pi) / (2 * pi * 27))
out_cross['literature'] = dict(
    gaona_reyes_2021='arXiv:2007.11980 (Gaona-Reyes, Carlesso, Bassi, PRD 103, 056011), Eq. 47: TD model for two particles "- m1 m2 (eta12 [x1,[x2,rho]] + '
                     'eta21 [x2,[x1,rho]])" with Eq. 46 kernels; text: "the corresponding term in Eq. (47) contains also position operators of different '
                     'particles ... the TD model cannot be reduced to that of KTM".',
    kryhin_sudhir_2023='arXiv:2309.09105 (PRL 134, 061501): Eq. 6 = epsilon x Diosi kernel; Eq. 8 "L_G ~ -epsilon G m1 m2/(hbar R^3)[x1,[x2,rho]]"; proposes '
                       'measuring S_x1x2 ("anti-correlated motion away from resonance, accompanied by a characteristic phase shift"). NOTE: direct expansion of '
                       'their Eq. 6 with both x along +axis gives +2 epsilon G m1 m2/(hbar R^3)[x1,[x2,rho]] (our -2 G hbar/d^3); their Eq. 8 has half the '
                       'magnitude and opposite sign (their main-text H_G = G m1 m2 x1 x2/(2R^3) is also 4x below the K = 2Gm1m2/R^3 used in their '
                       'appendix mu). Absorbed in their free epsilon; sign depends on axis convention - flagged, not resolved.',
    kafri_taylor_milburn_2014='arXiv:1401.0946 Eq. 15-16: Gamma = K/2, "-(K/2hbar) sum_k [x_k,[x_k,rho]]" - self terms only. Stage 27 statement '
                              '"this is exactly the KTM mutual noise" is WRONG: KTM noise is uncorrelated.')
out['claim2_cross'] = out_cross

# ================================================================ 3. DISCRETE (MATTER-SOURCED) KERNEL
# D0 = delta(x-y)/(Dhat rho(x)),  field noise <dm dm> = (Dhat/4) rho(x) delta  (Stage 26/27 convention; Dhat in kg s).
# (a) self term of a uniform sphere: Var F_x = (Dhat/4) int rho(z) g_x(z)^2 dz with the body's own field g = (4pi/3) G rho z  (MC check)
def self_term_MC(n=2000000):
    z = rng.random((n, 3)) * 2 - 1; z = z[np.linalg.norm(z, axis=1) <= 1]
    rho = 1.0; R = 1.0; M = 4 / 3 * pi
    gx = 4 * pi / 3 * rho * z[:, 0]                 # G = 1, Dhat = 1
    varF = 0.25 * rho * np.mean(gx ** 2) * (4 / 3 * pi)
    return varF / M ** 2 / (pi * rho / (15 * R))
disc = dict(self_term_MC_over_claim=self_term_MC())
# (b) source-mass term: point source M_s at r on x axis -> Gamma_xx = (G^2 Dhat/4) M_s / r^4 (exact for point source; transverse 0)
disc['source_term'] = 'Gamma_xx = G^2 Dhat M_s/(4 r^4) (longitudinal), 0 transverse - confirmed analytically'
# (c) DP: no modulation - Gamma_AA from the translation-invariant kernel contains no B-dependent term (cross terms only correlate). Confirmed.
# (d) X-ray heating lower bound on Dhat: Gamma_ba = hbar^2/(Dhat m^2) int (d_x mu)^2/rho_tot.  Stage 27 used rho_tot = mu (isolated nucleus).
#     In liquid Xe the smeared total density includes the neighbours once l > a (a = (m/rho)^(1/3) = 4.2 A):
mXe = 131.29 * U; rhoXe = 2940.0; aXe = (mXe / rhoXe) ** (1 / 3)
GX = GH / (6 * sqrt(pi) * 4.5e-10 ** 3)   # CITED Stage 26: XENONnT per-axis two-sided bound in the DP mapping
def I_isolated(ell): return mXe / ell ** 2                                    # int (d_x mu)^2/mu = m/l^2
def I_dense(ell): return mXe ** 2 / (16 * pi ** 1.5 * rhoXe * ell ** 5)     # rho_tot -> uniform rho
def I_lattice(ell, a=aXe, npts=48):
    """int (d_x mu)^2 / rho_tot over a simple-cubic lattice (ASSUMED stand-in for the liquid), Gaussian std l on every site."""
    L = 5 * ell; x = np.linspace(-L, L, npts); dx = x[1] - x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    g = lambda XX, YY, ZZ: np.exp(-(XX ** 2 + YY ** 2 + ZZ ** 2) / (2 * ell ** 2)) / (2 * pi * ell ** 2) ** 1.5
    mu = mXe * g(X, Y, Z); dmu = -X / ell ** 2 * mu
    if ell < 0.6 * a:
        nmax = int(np.ceil((L + 7 * ell) / a)); rt = mu.copy()
        for i in range(-nmax, nmax + 1):
            for j in range(-nmax, nmax + 1):
                for k in range(-nmax, nmax + 1):
                    if i == j == k == 0: continue
                    if sqrt(i * i + j * j + k * k) * a > L + 8 * ell: continue
                    rt += mXe * g(X - i * a, Y - j * a, Z - k * a)
    else:   # Fourier series of the periodic smeared density
        nmax = max(1, int(np.ceil(7 / (2 * pi * ell / a))))
        rt = np.zeros_like(X)
        for i in range(-nmax, nmax + 1):
            for j in range(-nmax, nmax + 1):
                for k in range(-nmax, nmax + 1):
                    Gv = 2 * pi / a * np.array([i, j, k]); w = np.exp(-0.5 * (Gv @ Gv) * ell ** 2)
                    if w < 1e-12: continue
                    rt += w * np.cos(Gv[0] * X + Gv[1] * Y + Gv[2] * Z)
        rt *= mXe / a ** 3
    rt = np.maximum(rt, 1e-300)
    return float(np.sum(dmu ** 2 / rt) * dx ** 3)
lat_ratio = []
for e_over_a in (0.15, 0.3, 0.5, 0.8, 1.0, 1.5, 2.5, 4.0):
    e = e_over_a * aXe; Il = I_lattice(e)
    lat_ratio.append(dict(l_over_a=e_over_a, I_lattice_over_isolated=Il / I_isolated(e), I_lattice_over_dense=Il / I_dense(e)))
disc['heating_integral_vs_l'] = lat_ratio
def I_best(ell):   # lattice result interpolated: isolated for l < 0.15a, dense for l > 4a, lattice between
    return min(I_isolated(ell), max(I_dense(ell), 0.0)) if ell > 4 * aXe else (I_isolated(ell) if ell < 0.15 * aXe else I_lattice(ell))
Dhat_max_LPF = 15 * (3 * 0.046 ** 3 / (4 * pi)) ** (1 / 3) * 7.57e-31 / (pi * (1.928 / 0.046 ** 3) * G ** 2)
ELL = np.logspace(-10, -4.5, 111)
Dmin_iso = np.array([HBAR ** 2 * I_isolated(e) / (mXe ** 2 * GX) for e in ELL])
Dmin_new = np.array([HBAR ** 2 * I_best(e) / (mXe ** 2 * GX) for e in ELL])
Dmax_env = 3 * ELL * GX / (sqrt(pi) * G ** 2 * rhoXe)          # field noise on a nucleus in uniform LXe, smeared at l (DERIVED)
Dmax = np.minimum(Dhat_max_LPF, Dmax_env)
# molecule interferometry lower bound (Oppenheim 2023 App. F scaling lambda ~ M_lam/Dhat for dx >~ l, CITED inputs M_lam ~ 1e-24 kg,
# lambda <~ 10 /s; DERIVED: for l > dx the rate falls as (dx/l)^2; ASSUMED dx = 1e-7 m). Isolated molecule in vacuum -> no dense correction.
M_LAM, LAM, DX = 1e-24, 10.0, 1e-7
Dmin_interf = (M_LAM / LAM) * np.minimum(1.0, (DX / ELL) ** 2)
Dmin_new = np.maximum(Dmin_new, Dmin_interf)
def first_open(dmin): ok = Dmax > dmin; return float(ELL[np.argmax(ok)]) if ok.any() else None
lmin_iso = first_open(Dmin_iso); lmin_new = brentq(lambda le: log10(HBAR ** 2 * I_best(10 ** le) / (mXe ** 2 * GX)) - log10(Dhat_max_LPF), -9.5, -6.5)
disc.update(Dhat_max_LPF_kg_s=Dhat_max_LPF, a_LXe=aXe,
            l_min_stage27_isolated_reproduced=lmin_iso, l_min_dense_corrected=10 ** lmin_new,
            Dmax_env_at_1nm=float(3e-9 * GX / (sqrt(pi) * G ** 2 * rhoXe)))
# eta of the discrete kernel in the dense limit (min over Dhat of heating + field noise, vs DP point value) - Cauchy-Schwarz check
eta_dense = (2 * sqrt(HBAR ** 2 / (16 * pi ** 1.5 * rhoXe) * rhoXe * G ** 2 * sqrt(pi) / 3)) / (GH / (6 * sqrt(pi)))
disc['eta_dense_limit'] = eta_dense
disc['eta_note'] = 'dense limit = ultralocal-continuous value sqrt(3) = %.3f; Stage 26 l >= eta^(1/3) R0 still holds (eta >= 1).' % sqrt(3)
# window and source-mass coverage (Stage 27 lock-in: 10 kg at 5 cm, LPF-class sensor; D_det CITED from Stage 27 = 1.2e-17 kg s)
s27 = json.load(open(os.path.join(ROOT, 'results/stage27.json')))
D_det = s27['tuning']['b_modulation']['discrete']['Dhat_detectable_LPFclass_kg']
D_det_tor = s27['tuning']['b_modulation']['discrete']['Dhat_detectable_torsion_1kg_1cm_kg']
win_new = np.where((Dmax > Dmin_new) & (ELL <= 52e-6), np.log10(Dmax / Dmin_new), np.nan)
cov_new = np.where(np.isfinite(win_new), np.log10(Dmax / np.maximum(Dmin_new, D_det)) / win_new, np.nan)
pick = lambda arr, e: float(np.interp(np.log(e), np.log(ELL), arr))
disc['window_orders_corrected'] = {f'{e:.0e}': pick(win_new, e) for e in (1e-8, 1e-7, 1e-6, 1e-5, 5e-5)}
disc['lockin_fraction_covered_corrected'] = {f'{e:.0e}': pick(cov_new, e) for e in (1e-8, 1e-7, 1e-6, 1e-5, 5e-5)}
disc['lockin_Dhat_detectable_kg_s'] = D_det; disc['torsion_Dhat_detectable_kg_s'] = D_det_tor
disc['interferometry_bound'] = dict(M_lambda=M_LAM, lambda_max=LAM, dx_ASSUMED=DX, Dhat_min_l_le_dx=M_LAM / LAM,
                                   crossover_l_xray_to_interf=float(ELL[np.argmax(Dmin_interf >= np.array([HBAR ** 2 * I_best(e) / (mXe ** 2 * GX) for e in ELL]))]))
disc['scan'] = dict(ell=ELL.tolist(), Dmin_isolated=Dmin_iso.tolist(), Dmin_corrected=Dmin_new.tolist(), Dmin_interf=Dmin_interf.tolist(),
                    Dmax=Dmax.tolist(), Dmax_env=Dmax_env.tolist())
disc['caveats'] = ['ASSUMED: the state-dependent kernel uses the TOTAL smeared mass density m(x) (field-theory reading, as in Oppenheim Eq. 46 where the '
                   'decoherence rate of a body ~ M_lambda m_P/(l_P^3 D2), i.e. its own bulk density). Stage 27 used the isolated nucleus.',
                   'ASSUMED: simple-cubic lattice with a = (m/rho)^(1/3) for liquid Xe in the crossover l ~ a.',
                   'OPEN: for l > a neighbouring nuclei receive strongly correlated kicks; collective X-ray emission was not modelled (same open point '
                   'applies to the Stage 26 DP X-ray mapping for l > a).',
                   'LPF self-term only (spacecraft/housing masses would add matter-sourced noise -> tighter Dhat_max); cube treated as equal-volume sphere.']
disc['literature'] = dict(
    oppenheim_2023='arXiv:2203.01982 Sec. V / App. F: "the discrete models appear less constrained due to the suppression of the noise away from the mass '
                   'density ... <D2(x,y,Phi)> = (l_P^3/m_P) D2 m(x)" and the torsion/Cavendish squeeze "10^-1 kg s >= (l_P^3/m_P) D2 >= 10^-25 kg s" '
                   '- matter-sourced noise and its test with source-mass (Cavendish) experiments already in the paper.',
    carney_taylor_2023='arXiv:2301.08378 ("Strongly incoherent gravity"): noise sourced by a source mass, Gamma_incoh = G^2 M m^2 l^2 sigma^2/(d^6 v^2), and '
                       '"the entire remaining parameter space can be tested by looking for loss of quantum coherence in small systems like atom '
                       'interferometers coupled to oscillating source masses" - modulated-source-mass test of matter-sourced classical-gravity noise.')
out['claim3_discrete'] = disc

# ================================================================ verdicts
out['verdicts'] = dict(
    claim1=dict(verdict='CONFIRMED numerically for spheres/cubes; ALREADY KNOWN (Diosi 1989/2007; NHH 2014 Eq. 11; volume signature in Oppenheim 2023); '
                        'CORRECTED in scope',
                value=vol['GammaV_sphere_numeric'], ratio=vol['ratio_numeric_over_claim'],
                corrections=['shape: Gamma_ii V = 4 pi G hbar (kappa+1/4kappa) N_ii (demagnetization tensor); claim exact only for N = 1/3 (sphere, cube)',
                             'composition/mass independence only for homogeneous bodies; general V -> V_eff = M^2/int rho^2',
                             'finite l/R deficit (see finite_l_over_R); lattice correction < 1e-30 at l >= 4.5 A']),
    claim2=dict(verdict='CONFIRMED (-2.000, +1.000 G hbar/d^3 at kappa=1/2); ALREADY KNOWN (Gaona-Reyes et al. 2021 Eq. 47; Kryhin-Sudhir 2023 Eq. 6/8); '
                        'CORRECTED wording; Stage 27 KTM attribution WRONG',
                corrections=['magnitude is -2 G hbar (kappa+1/4kappa)/d^3 >= 2 G hbar/d^3, not kappa-independent; only rho is kappa-independent and only in '
                             'the DP family', 'KTM (LOCC) predicts rho = 0; TD with local gamma predicts zero longitudinal correlation -> the sign/axis '
                             'pattern discriminates kernels (useful, and the one framing we did not find stated)']),
    claim3=dict(verdict='DP no-modulation CONFIRMED; matter-sourced modulation ALREADY KNOWN in concept (Oppenheim 2023 Eq. 46 / App. F; Carney-Taylor 2023); '
                        'Stage 27 bound l >= 2.5 um CORRECTED',
                l_min_stage27=lmin_iso, l_min_corrected=10 ** lmin_new,
                corrections=['heating integral must use total smeared density -> D-hat_min falls as l^-5 not l^-2 for l > a',
                             'Dhat units are kg s (label only)']))
out['runtime_s'] = time.time() - t0
json.dump(out, open(os.path.join(ROOT, 'results/stage29.json'), 'w'), indent=1, default=float)
show = {k: v for k, v in out.items() if k not in ('claim3_discrete',)}
print(json.dumps(show, indent=1, default=float)[:9000])
print(json.dumps({k: v for k, v in disc.items() if k != 'scan'}, indent=1, default=float))
