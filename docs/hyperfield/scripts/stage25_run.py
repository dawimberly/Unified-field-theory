#!/usr/bin/env python3
"""Stage 25 - test the two Stage 24 flags.
(A) N-counting in Oppenheim et al. Eq. 47 (non-local / Diosi-Penrose kernel floor), via the
    paper's own master equation (Nat Commun 14, 7910, Methods Eq. 135-137), evaluated as an
    explicit double sum on a bcc Na lattice cluster.
(B) Carry-over of DP spontaneous-radiation bounds (Donadi 2021, XENONnT 2026) to the CQ model with
    the same kernel and the same smearing length l.
Labels: DERIVED (computed here), CITED (number taken from a source), ASSUMED (our choice).
"""
import json, os, time
import numpy as np
from math import pi, sqrt, log10
from scipy.special import erf

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
H = 6.62607015e-34; HBAR = H / (2 * pi); C = 299792458.0; G = 6.6743e-11; U = 1.66053906660e-27
KB = 1.380649e-23; E0 = 8.8541878128e-12; QE = 1.602176634e-19
M_NA = 22.98976928 * U
A_BCC = 4.2906e-10          # CITED: Na bcc lattice constant at room temperature (~4.29 A)
RNG = np.random.default_rng(25)

# ---- inputs carried from Stage 15/22/24 (CITED there) ----
PED = dict(M=172e3 * U, tau=2 * 0.983 / 160, d=133e-9)          # Pedalino 2025 floor point
CEIL = 0.298 / 2.79e11 * 2.79e11 * (1e-15) ** 3 / G ** 2           # Gisler 2022 ceiling, Eq.47 LHS with r_N = 1e-15
FLOOR_PUB = None
CHEAP = dict(d=133e-9, xi=0.8)
SAG_MAX = 0.0317

def R_lam(dx, rN=1e-15): return (rN ** 2 * dx) ** (1 / 3)
def floor_published(M, N, tau, dx, rN=1e-15):
    """Eq. 47 RHS as published (N_lam * M_lam^2 / R_lam) * tau, M = TOTAL mass (Stage 15 reading)."""
    return N * M ** 2 * tau / R_lam(dx, rN)

# ---------------------------------------------------------------- lattice
def bcc_cluster(N):
    """N atoms of a bcc lattice nearest to an atom at the origin; integer coords in units a/2."""
    n = int(np.ceil((N / 2) ** (1 / 3) * 1.6)) + 3
    r = np.arange(-2 * n, 2 * n + 1)
    X, Y, Z = np.meshgrid(r, r, r, indexing='ij')
    X, Y, Z = X.ravel(), Y.ravel(), Z.ravel()
    ok = ((X % 2 == 0) & (Y % 2 == 0) & (Z % 2 == 0)) | ((X % 2 == 1) & (Y % 2 == 1) & (Z % 2 == 1))
    P = np.stack([X[ok], Y[ok], Z[ok]], 1)
    r2 = (P ** 2).sum(1)
    idx = np.argsort(r2, kind='stable')[:N]
    return P[idx].astype(np.int64)

def K(r, ell):
    """Pair kernel of two Gaussian-smeared (std ell) unit masses under 1/|x-y|: erf(r/2l)/r; K(0)=1/(sqrt(pi) l)."""
    r = np.asarray(r, float)
    out = np.empty_like(r)
    small = r < 1e-6 * ell
    out[small] = 1.0 / (sqrt(pi) * ell)
    rs = r[~small]
    big = rs > 14 * ell
    o = np.empty_like(rs)
    o[big] = 1.0 / rs[big]
    o[~big] = erf(rs[~big] / (2 * ell)) / rs[~big]
    out[~small] = o
    return out

def lattice_S(P, ells, dx, nrows):
    """S(l,dx) = sum_ij [2K(r_ij) - K(|r_ij+dx e|) - K(|r_ij-dx e|)], estimated by exact sums over j for
    nrows randomly chosen rows i (same rows for all three terms), scaled by N/nrows. Returns mean, s.e."""
    N = len(P); h = A_BCC / 2
    rows = RNG.choice(N, size=min(nrows, N), replace=False)
    Pf = P.astype(float) * h
    res = {e: [] for e in ells}
    for i in rows:
        dv = Pf - Pf[i]
        r0 = np.sqrt((dv ** 2).sum(1))
        rp = np.sqrt((dv[:, 0] + dx) ** 2 + dv[:, 1] ** 2 + dv[:, 2] ** 2)
        rm = np.sqrt((dv[:, 0] - dx) ** 2 + dv[:, 1] ** 2 + dv[:, 2] ** 2)
        for e in ells:
            res[e].append((2 * K(r0, e) - K(rp, e) - K(rm, e)).sum())
    out = {}
    for e in ells:
        a = np.array(res[e]) * N
        out[e] = (a.mean(), a.std(ddof=1) / sqrt(len(a)))
    return out

def pair_shells(P):
    """Exact histogram of squared pair separations (integer, units (a/2)^2) over i != j."""
    N = len(P); cnt = np.zeros(1, np.int64)
    for s in range(0, N, 400):
        d = P[s:s + 400, None, :] - P[None, :, :]
        n2 = (d * d).sum(2).ravel()
        b = np.bincount(n2)
        if len(b) > len(cnt): cnt = np.pad(cnt, (0, len(b) - len(cnt)))
        cnt[:len(b)] += b
    cnt[0] -= N
    return cnt

# ---- analytic regimes (all in "S" units: sum over ordered pairs, unit masses) ----
def sphere_mutual(dx, R):
    """sum-normalised mutual 1/r energy of two uniform unit-mass spheres (radius R) at separation dx."""
    if dx >= 2 * R: return 1.0 / dx
    u = dx / R
    return (6 / 5 - 0.5 * u ** 2 + 3 / 16 * u ** 3 - u ** 5 / 160) / R
def S_incoh(N, ell, dx): return 2 * N * (1 / (sqrt(pi) * ell) - float(K(np.array([dx]), ell)[0]))
def S_coh(N, R, dx): return 2 * N ** 2 * (6 / (5 * R) - sphere_mutual(dx, R))
def S_blob(N, ell, dx): return 2 * N ** 2 * (1 / (sqrt(pi) * ell) - float(K(np.array([dx]), ell)[0]))

def R_body(N): return (3 * N * A_BCC ** 3 / (2 * 4 * pi)) ** (1 / 3)   # bcc: 2 atoms per a^3

t0 = time.time()
ELLS = [1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 3e-12, 1e-11, 3e-11, 1e-10, 2e-10, 4e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7]
CASES = [('Pedalino 172 kDa', 7482), ('Stage-22 1.415 MDa', 61600)]
DXS = [('dx = d = 133 nm', 133e-9), ('dx = 1 A (< lattice spacing)', 1e-10)]
lat = {}
for lab, N in CASES:
    P = bcc_cluster(N)
    Rb = R_body(N)
    rmax = np.sqrt(((P * A_BCC / 2) ** 2).sum(1)).max()
    sh = pair_shells(P)
    nz = np.nonzero(sh)[0]
    rs = np.sqrt(nz) * A_BCC / 2; cs = sh[nz].astype(float)
    sum_inv_r = (cs / rs).sum()
    lat[lab] = dict(N=N, M_kg=N * M_NA, R_wigner_seitz_m=Rb, R_max_atom_m=float(rmax),
                    nearest_neighbour_m=float(rs[0]),
                    pair_sum_1_over_r_exact=float(sum_inv_r),
                    pair_sum_vs_uniform_sphere_6N2_over_5R=float(sum_inv_r / (6 * N ** 2 / (5 * Rb))), cases={})
    for dlab, dx in DXS:
        nrows = 600 if N > 20000 else 1200
        S = lattice_S(P, ELLS, dx, nrows)
        rowsout = []
        for e in ELLS:
            s, se = S[e]
            si = S_incoh(N, e, dx); sc = S_coh(N, Rb, dx); sb = S_blob(N, e, dx)
            # the shell histogram gives the l-dependent pair sum exactly; use it to cross-check large-dx case
            pair_K = (cs * erf(rs / (2 * e)) / rs).sum()
            rowsout.append(dict(ell=e, S_lattice=s, S_lattice_se=se, S_incoherent=si, S_coherent_sphere=sc,
                                S_blob_l_gg_R=sb, frac_incoh=si / s if s > 0 else float('nan'),
                                exact_pair_sum_K=float(pair_K),
                                ratio_to_NM2_form=s / (N * N ** 2 / R_lam(dx)),
                                ratio_to_N_form_incoh_count=s / (N / R_lam(dx)),
                                ratio_to_M2_form=s / (N ** 2 / R_lam(dx))))
        lat[lab]['cases'][dlab] = rowsout
    print(lab, 'done', round(time.time() - t0, 1), 's', flush=True)

# effective power of N: fit S ~ N^p across the two clusters at each (l, dx)
Npow = {}
for dlab, dx in DXS:
    a = lat[CASES[0][0]]['cases'][dlab]; b = lat[CASES[1][0]]['cases'][dlab]
    Npow[dlab] = [dict(ell=x['ell'], p_N=log10(y['S_lattice'] / x['S_lattice']) / log10(CASES[1][1] / CASES[0][1]))
                  for x, y in zip(a, b)]

# crossover lengths (DERIVED, analytic)
cross = {}
for lab, N in CASES:
    Rb = R_body(N)
    cross[lab] = dict(l_star_large_dx=Rb / (1.2 * sqrt(pi) * N),       # N/(sqrt(pi) l) = 1.2 N^2 / R
                      l_star_small_dx=Rb / N ** (1 / 3))                 # N dx^2/l^3 ~ N^2 dx^2/R^3

# ---------------------------------------------------------------- window impact (A)
def Q_derived(M, dx, ell, N=None):
    """Eq.137 rate coefficient in Eq.47 normalisation (constant factor 2 dropped as the paper does):
    Q = m_a^2 * S / 2 with the lattice result approximated by incoherent + coherent-sphere pieces."""
    N = M / M_NA if N is None else N
    ma = M / N; Rb = R_body(N)
    if ell < Rb:
        S = S_incoh(N, ell, dx) + S_coh(N, Rb, dx)
    else:
        S = S_blob(N, ell, dx)
    return ma ** 2 * S / 2

def Q_published(M, dx, ell=1e-15):
    N = M / M_NA
    return N * M ** 2 / R_lam(dx, ell)

Nped = PED['M'] / M_NA
fl_pub = floor_published(PED['M'], Nped, PED['tau'], PED['d'])
fl_der = Q_derived(PED['M'], PED['d'], 1e-15) * PED['tau']
# also from the actual lattice sum
Sl = [x for x in lat['Pedalino 172 kDa']['cases']['dx = d = 133 nm'] if x['ell'] == 1e-15][0]['S_lattice']
fl_lat = (M_NA ** 2 * Sl / 2) * PED['tau']
# nucleon-resolved variant: each Na nucleus = 23 nucleons smeared at l inside a nucleus of radius 1.2 A^1/3 fm
A_na = 23; r_nuc = 1.2e-15 * A_na ** (1 / 3); m_n = M_NA / A_na
self_atom = A_na ** 2 / (sqrt(pi) * 1e-15)
self_nucl = A_na / (sqrt(pi) * 1e-15) + A_na ** 2 * 1.2 / r_nuc
nucleon_factor = self_nucl / self_atom

def closing(Qfun, ceil):
    def g(lm):
        M = 10 ** lm * U; tau = 2 * CHEAP['xi'] * CHEAP['d'] ** 2 * M / H
        return log10(Qfun(M, CHEAP['xi'] * CHEAP['d']) * tau) - log10(ceil)
    lo, hi = 3.0, 20.0
    if g(hi) < 0: return float('inf')
    for _ in range(200):
        mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if g(mid) < 0 else (lo, mid)
    return 10 ** hi
Mc_pub = closing(lambda M, dx: Q_published(M, dx), CEIL)
Mc_der = closing(lambda M, dx: Q_derived(M, dx, 1e-15), CEIL)
def sag(Mamu): return 0.5 * 9.81 * (CHEAP['xi'] * CHEAP['d'] ** 2 * Mamu * U / H) ** 2

# Table-1 ultralocal check (paper's Eq. 152 numbers)
T1 = dict(M=1e-24, N=1e3, V=1e-37, lam=10.0)
T1['published_N_M2_over_V_lam'] = T1['N'] * T1['M'] ** 2 / (T1['V'] * T1['lam'])
T1['incoherent_M2_over_N_V_lam'] = T1['M'] ** 2 / (T1['N'] * T1['V'] * T1['lam'])
T1['incoherent_nuclear_volume_1e-45'] = T1['M'] ** 2 / (T1['N'] * 1e-45 * T1['lam'])
T1['paper_upper_bound'] = 1e-41

A_res = dict(
    floor_published_Pedalino=fl_pub, floor_derived_Pedalino_analytic=fl_der, floor_derived_Pedalino_lattice=fl_lat,
    floor_shift_orders=log10(fl_lat / fl_pub), ceiling=CEIL,
    window_published_orders=log10(CEIL / fl_pub), window_derived_orders=log10(CEIL / fl_lat),
    nucleon_resolved_factor=nucleon_factor, nucleon_resolved_shift_orders=log10(nucleon_factor),
    incoherent_fraction_at_Pedalino=S_incoh(Nped, 1e-15, PED['d']) / (S_incoh(Nped, 1e-15, PED['d']) + S_coh(Nped, R_body(Nped), PED['d'])),
    closing_mass_published_amu=Mc_pub, closing_mass_derived_amu=Mc_der,
    closing_sag_derived_m=sag(Mc_der), closing_sag_published_m=sag(Mc_pub),
    closing_v_derived_m_s=0.983 * H / (CHEAP['xi'] * CHEAP['d'] ** 2 * Mc_der * U),
    table1_ultralocal=T1, crossover_lengths=cross, N_power=Npow)
print('A: floor pub %.3g  der(lat) %.3g  shift %.2f  window %.2f -> %.2f  Mclose %.3g -> %.3g amu' % (
    fl_pub, fl_lat, A_res['floor_shift_orders'], A_res['window_published_orders'], A_res['window_derived_orders'], Mc_pub, Mc_der))

# ---------------------------------------------------------------- (B) radiation / heating carry-over
# Oppenheim Eq.135: d rho/dt = -(1/2) int D0(x,y) [m(x),[m(y),rho]].  Diosi convention: D0 = G/(hbar |x-y|).
# Write D0(x,y) = kappa * G / (hbar |x-y|) (Gaussian smearing l).  kappa = 1 <-> Diosi 1987 strength.
# Penrose/Donadi-2021 convention (4 pi G/hbar in place of G/(2 hbar)) is kappa = 8 pi.
R0_DON_PEN = 0.54e-10        # CITED Donadi et al. Nat Phys 17, 74 (2021), 95%, Penrose (4 pi G) convention
R0_DON = R0_DON_PEN / (8 * pi) ** (1 / 3)     # DERIVED: same bound in Diosi convention
R0_XE = 4.5e-10              # CITED XENONnT PRL 136, 120201 (2026), 95% CL, Diosi convention (stated in their SM)
A_DON = 1.8e-29              # CITED Donadi: z_s = a / R0^3 (Penrose convention), counts
Z_S_ALLOWED = A_DON / R0_DON_PEN ** 3    # DERIVED: signal counts allowed at the 95% bound

def kappa_max_xray(ell, R0b): return (ell / R0b) ** 3     # rate ~ kappa / l^3 (Donadi Eq. 4, l > nuclear size)
# DERIVED check of Donadi Eq. 4 normalisation from momentum diffusion: dGamma/domega (Diosi) = Z^2 e^2 G /(12 pi^2.5 eps0 c^3 R0^3 omega)
pref_diosi = 1 / (12 * pi ** 2.5); pref_donadi = 2 / (3 * pi ** 1.5)
# interferometric bound on kappa (Pedalino, lambda <= 1/tau) with the derived (lattice) counting
lam_per_kappa = (G / HBAR) * M_NA ** 2 * Sl              # Eq. 137 (no 1/2)
kappa_int = (1 / PED['tau']) / lam_per_kappa
# trade-off normalisation implied by Eq. 47 convention: floor F = m^2 S tau / 2 and kappa_int = hbar/(G m^2 S tau)
F_times_kappa = fl_lat * kappa_int                      # = hbar / (2G), experiment independent
# minimal kappa required by the torsion ceiling
def ceiling_l(ell, mode): return CEIL * ((ell / 1e-15) ** 3 if mode == 'consistent' else 1.0)
kappa_req_1e15 = F_times_kappa / CEIL
heat_W = kappa_req_1e15 * G * HBAR * M_NA / (4 * sqrt(pi) * (1e-15) ** 3)   # per Na atom, DERIVED (Diosi heating x kappa)
dTdt = 2 * heat_W / (3 * KB)
counts_pred = kappa_req_1e15 * (A_DON / (8 * pi)) / (1e-15) ** 3
B_scan = []
for e in np.logspace(-16, -8, 81):
    row = dict(ell=e)
    for nm, R0b in (('Donadi', R0_DON), ('XENONnT', R0_XE)):
        row['floor_xray_' + nm] = F_times_kappa / kappa_max_xray(e, R0b)
    row['ceiling_floor_only'] = ceiling_l(e, 'floor-only'); row['ceiling_consistent'] = ceiling_l(e, 'consistent')
    row['floor_interf_published'] = floor_published(PED['M'], Nped, PED['tau'], PED['d'], rN=e)
    row['floor_interf_derived'] = Q_derived(PED['M'], PED['d'], e) * PED['tau']
    B_scan.append(row)
def l_allowed(R0b, mode):
    # F_times_kappa (R0b/l)^3 <= CEIL * (l/1e-15)^3 [consistent]  or <= CEIL [floor-only]
    if mode == 'consistent': return (F_times_kappa * R0b ** 3 * 1e-45 / CEIL) ** (1 / 6)
    return R0b * (F_times_kappa / CEIL) ** (1 / 3)
B_res = dict(
    R0_Donadi_Penrose=R0_DON_PEN, R0_Donadi_Diosi=R0_DON, R0_XENONnT_Diosi=R0_XE,
    donadi_eq4_prefactor_over_diosi_derived=pref_donadi / pref_diosi, eight_pi=8 * pi,
    kappa_max_interferometry_Pedalino_l1em15=kappa_int,
    kappa_max_Donadi_l1em15=kappa_max_xray(1e-15, R0_DON), kappa_max_XENONnT_l1em15=kappa_max_xray(1e-15, R0_XE),
    hbar_over_2G=HBAR / (2 * G), F_times_kappa=F_times_kappa,
    kappa_required_by_ceiling_l1em15=kappa_req_1e15,
    exclusion_orders_Donadi=log10(kappa_req_1e15 / kappa_max_xray(1e-15, R0_DON)),
    exclusion_orders_XENONnT=log10(kappa_req_1e15 / kappa_max_xray(1e-15, R0_XE)),
    xray_vs_interferometry_orders_Donadi=log10(kappa_int / kappa_max_xray(1e-15, R0_DON)),
    xray_vs_interferometry_orders_XENONnT=log10(kappa_int / kappa_max_xray(1e-15, R0_XE)),
    floor_xray_Donadi_l1em15=F_times_kappa / kappa_max_xray(1e-15, R0_DON),
    floor_xray_XENONnT_l1em15=F_times_kappa / kappa_max_xray(1e-15, R0_XE),
    Ge_counts_predicted_at_min_kappa=counts_pred, Ge_counts_allowed_95=Z_S_ALLOWED,
    heating_W_per_Na_atom_at_min_kappa=heat_W, heating_dTdt_K_s=dTdt,
    nucleon_resolved_xray_suppression_Ge_A73=1 / 73, nucleon_resolved_exclusion_Donadi_orders=log10(kappa_req_1e15 / kappa_max_xray(1e-15, R0_DON)) - log10(73),
    l_min_consistent_Donadi=l_allowed(R0_DON, 'consistent'), l_min_consistent_XENONnT=l_allowed(R0_XE, 'consistent'),
    l_min_floor_only_Donadi=l_allowed(R0_DON, 'floor-only'), l_min_floor_only_XENONnT=l_allowed(R0_XE, 'floor-only'),
    scan=B_scan)
print('B: kappa_int %.3g  kappa_req %.3g  kappa_Don %.3g  kappa_Xe %.3g  excl %.2f / %.2f orders' % (
    kappa_int, kappa_req_1e15, B_res['kappa_max_Donadi_l1em15'], B_res['kappa_max_XENONnT_l1em15'],
    B_res['exclusion_orders_Donadi'], B_res['exclusion_orders_XENONnT']))
print('   l_min consistent Don %.3g Xe %.3g ; floor-only Don %.3g Xe %.3g' % (
    B_res['l_min_consistent_Donadi'], B_res['l_min_consistent_XENONnT'], B_res['l_min_floor_only_Donadi'], B_res['l_min_floor_only_XENONnT']))

out = dict(stage=25, inputs=dict(a_bcc_Na=A_BCC, PED=PED, ceiling=CEIL, cheap=CHEAP),
           lattice=lat, A=A_res, B=B_res,
           sources=dict(
               oppenheim_natcommun='https://www.nature.com/articles/s41467-023-43348-2 (Methods Eqs. 135-137, 147-152)',
               layton='https://arxiv.org/abs/2307.02557 (Sec. 7: cut-off < 1e-15 m excluded by heating, citing Tilloy-Diosi 2017)',
               tilloy_diosi='https://arxiv.org/abs/1706.01856 (PRD 96, 104045)',
               donadi='https://arxiv.org/abs/2111.13490 (Nat Phys 17, 74; Eq. 3 uses 4 pi G/hbar; Eq. 4; R0 > 0.54e-10 m)',
               xenonnt='https://arxiv.org/abs/2506.05507 (PRL 136, 120201; R0 > 4.5e-10 m 95% CL, Diosi convention)',
               piccione_smearing='https://arxiv.org/abs/2511.00644',
               piccione_gpsl='https://arxiv.org/abs/2602.19377',
               janse='https://arxiv.org/abs/2403.08912'),
           runtime_s=time.time() - t0)
json.dump(out, open(os.path.join(ROOT, 'results/stage25.json'), 'w'), indent=1, default=float)
print('written; %.1f s' % (time.time() - t0))
