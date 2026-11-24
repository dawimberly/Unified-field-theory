"""Stage 6B: 'churn of the hypersphere' -- rotating/moving compact space, Lorentz bound, crystallization (cubic invariants),
and the multi-field 'viscosity' layer picture. Mappings ASSUMED as labelled."""
import json, os, numpy as np
from math import comb, sqrt, pi
CK = 'results/stage6.json'
out = json.load(open(CK))
def save(): json.dump(out, open(CK, 'w'), indent=1, default=float)
HBAR_eVs = 6.582119569e-16; C = 2.99792458e8; HBARC = 1.973269804e-16
st5 = json.load(open('results/stage5.json'))

# ---- (1) splitting of first shells under churn
def sphere_m_multiplicities(N, l):
    """harmonics of degree l on S^{N-1}: distribution of J_12 eigenvalue m (z=x1+ix2, rest N-2 real coords)"""
    def P(L):
        dist = {}
        for a in range(L+1):
            for b in range(L+1-a):
                c = L-a-b; mult = comb(c+N-3, N-3)
                dist[a-b] = dist.get(a-b, 0) + mult
        return dist
    p = P(l); q = P(l-2) if l >= 2 else {}
    return {m: p.get(m, 0)-q.get(m, 0) for m in sorted(set(p) | set(q)) if p.get(m, 0)-q.get(m, 0)}
split = {}
for l in (1, 2, 3):
    mm = sphere_m_multiplicities(10, l); split[f'S9_l{l}'] = dict(total=sum(mm.values()), m_multiplicity=mm,
         rule='E -> E - m*hbar*Omega (DERIVED: J12 eigenvalue m couples to rotation)')
split['Z9_shell1'] = dict(total=18, pattern='velocity along e1: 2 modes shift by +/- m1*v, 16 unshifted (Doppler k.v, DERIVED)')
split['Lambda9_shell1'] = dict(total=2, pattern='both modes along the stacking normal: shift +/- m1*v_normal; no shift for v inside the E8 layer')
split['E8xS1_shell1'] = dict(total=2, pattern='circle modes: shift +/- m1*v_circle')
split['note'] = ('A flat lattice torus has NO continuous rotation isometry: rigid "churn" of a lattice is not a symmetry, '
                 'it is a time-dependent modulus (shear) of the metric. The continuous churn available on T^9 is translation '
                 '(brane sliding along a cycle, U(1)^9). On S^9 it is SO(10) rotation. A rotation of the compact space with '
                 'the brane co-moving is pure gauge (unobservable); only the brane velocity RELATIVE to the bulk geometry is physical.')
out['partB_splitting'] = split

# ---- Lorentz bound (GW170817, ApJL 848 L13: -3e-15 < (c_g-c)/c < +7e-16, verified)
# DERIVED: brane moving with internal speed v has induced metric -(1-v^2)dt^2+dx^2 -> c_gamma/c_g = sqrt(1-v^2) -> c_g/c_gamma-1 = v^2/2 > 0
vmax = sqrt(2*7e-16)
rows = []
for x in st5['lattices9']:
    R = x.get('R_sphere_m') if 'S9' in x['lattice'] else x['lam_m']
    Om = vmax*C/R
    rows.append(dict(space=x['lattice'], M_D_TeV=x['M_D_TeV'], radius_m=R, m1_GeV=x['m1_GeV'],
                     Omega_max_rad_s=Om, hbarOmega_max_eV=HBAR_eVs*Om, Omega_over_omegaKK=HBAR_eVs*Om/(x['m1_GeV']*1e9),
                     max_shell1_splitting_eV=x['m1_GeV']*1e9*vmax/(3 if 'S9' in x['lattice'] else 1)))
out['partB_lorentz'] = dict(v_max_over_c=vmax, source='GW170817/GRB170817A, arXiv:1710.05834 (upper side +7e-16 applies since c_g>c_gamma)',
    sme_note=('All SM fields share the same induced brane metric, so SM-only SME coefficients (photon, electron) vanish identically: '
              'churn is invisible to them (DERIVED). Only the gravity-vs-matter comparison sees it.'), rows=rows)

# ---- (2b) cubic Landau invariant: zero-sum triangles among shortest vectors
def tri_count(V):
    S = {tuple(np.round(v, 6)) for v in V}; n = 0
    for i in range(len(V)):
        for j in range(len(V)):
            if tuple(np.round(-V[i]-V[j], 6)) in S: n += 1
    return n//6 if n else 0          # unordered triangles
ns = {}; src = open('scripts/stage6_lattices.py').read()
exec(src.split("BEST_LAT_KISS")[0].replace("out = json.load(open(CK)) if os.path.exists(CK) else {}", "out={}"), ns)
def mv(gens):
    B = ns['lll'](ns['hnf_basis'](np.round(2*np.array(gens)).astype(int))/2); return ns['minvecs'](B)[1]
tri = {}
for d in (4, 8, 9):
    I = np.eye(d); h = np.full(d, .5); Dg = [I[i]-I[i+1] for i in range(d-1)] + [I[d-2]+I[d-1]]
    tri[f'Z{d}'] = (len(mv(list(I))), tri_count(mv(list(I))))
    tri[f'D{d}'] = (len(mv(Dg)), tri_count(mv(Dg)))
    tri[f'BCC{d}'] = (len(mv(list(I)+[h])), tri_count(mv(list(I)+[h])))
    if d == 8: tri['E8=D8+'] = (240, tri_count(mv(Dg+[h])))
E8 = np.array([[2,0,0,0,0,0,0,0],[-1,1,0,0,0,0,0,0],[0,-1,1,0,0,0,0,0],[0,0,-1,1,0,0,0,0],[0,0,0,-1,1,0,0,0],[0,0,0,0,-1,1,0,0],[0,0,0,0,0,-1,1,0],[.5]*8])
B9 = np.zeros((9, 9)); B9[:8, :8] = E8; B9[8, 0] = 1; B9[8, 8] = 1
V9 = ns['minvecs'](ns['lll'](B9))[1]; tri['Lambda9'] = (len(V9), tri_count(V9))
out['partB_cubic_invariant'] = dict(counts={k: dict(shell=v[0], zero_sum_triangles=v[1]) for k, v in tri.items()},
    landau=('ASSUMED mapping: density on S^{d-1} rho = rho0 + A sum_{k in shell} cos(k.x). F = r A^2 - w A^3 + u A^4 with w proportional to '
            '#zero-sum triangles (Landau-Alexander-McTague). w != 0 -> first-order (DERIVED): coexistence at r_c = w^2/(4u), '
            'ordered-phase spinodal r_s = 9w^2/(32u), disordered spinodal r = 0: hysteresis window 0<r<r_s, r_s/r_c = 9/8 universal.'))

# ---- (3) multi-field viscosity layers (ASSUMED mapping)
GeV_s = HBAR_eVs*1e-9
F = [  # name, mass GeV, width GeV (None = stable/confined), alpha-like coupling at mZ, note
 ('photon', 0.0, None, 1/127.95, 'm<1e-18 eV; alpha_em(mZ)'),
 ('neutrino', 0.45e-9, None, None, 'm<0.45 eV (KATRIN 2024, upper bound)'),
 ('electron', 0.51099895e-3, None, None, 'stable (tau>6.6e28 yr)'),
 ('muon', 0.1056583755, HBAR_eVs*1e-9/2.1969811e-6, None, 'tau=2.197 us'),
 ('tau', 1.77693, HBAR_eVs*1e-9/290.3e-15, None, 'tau=290.3 fs'),
 ('Higgs', 125.20, 3.7e-3, None, 'Gamma measured 3.7(+1.9-1.4) MeV (CMS 2022); SM 4.1 MeV'),
 ('W', 80.369, 2.085, 0.0338, 'alpha_2(mZ)'),
 ('Z', 91.188, 2.4955, 0.0338, ''),
 ('top', 172.57, 1.42, None, 'decays before hadronizing'),
 ('gluon/light quarks', 0.0, 0.33, 0.1180, 'confined: lock time ~1/Lambda_QCD (Lambda~0.33 GeV, ASSUMED proxy); QGP eta/s ~ 0.08-0.2 near KSS 1/4pi=0.0796'),
]
vmaxOm = out['partB_lorentz']['rows'][-2]['hbarOmega_max_eV']*1e-9   # S9 at M_D=6 TeV, in GeV
lay = []
for n, m, G, a, note in F:
    y = sqrt(2)*m/246.22 if n not in ('photon', 'gluon/light quarks', 'W', 'Z', 'Higgs') else None
    if G: tau = HBAR_eVs*1e-9/G; rate = G
    elif m > 0: tau = HBAR_eVs*1e-9/m; rate = m
    else: tau = float('inf'); rate = 0.0
    lay.append(dict(field=n, mass_GeV=m, width_GeV=G, Gamma_over_m=(G/m if G and m else None), coupling=a, yukawa=y,
                    tau_s=tau, response_rate_GeV=rate, Omega_tau_at_max_churn=(vmaxOm/rate if rate else float('inf')),
                    state=('locked (follows churn)' if rate and vmaxOm/rate < 1 else 'free (lags churn)'), note=note))
out['partB_multifield'] = dict(hbarOmega_max_GeV=vmaxOm, rule='Omega*tau<1 -> fluid co-rotates (locked); >1 -> lags (free). tau=1/Gamma (unstable) or 1/m (stable, ASSUMED)', rows=lay,
    fermi_check='for weak decays Gamma ~ G_F^2 m^5 -> Gamma/m ~ (m/M_W)^4: muon/tau ratio (m_tau/m_mu)^4*BR-correction; layer order = mass hierarchy (DERIVED, not new)')
mu, ta = lay[3], lay[4]
out['partB_multifield']['fermi_ratio_check'] = dict(measured=(ta['Gamma_over_m']/mu['Gamma_over_m']), predicted_m4_over_BR=(1.77693/0.1056584)**4/0.1782)
save()
print(json.dumps({k: out[k] for k in ['partB_splitting', 'partB_lorentz', 'partB_cubic_invariant']}, indent=1, default=float)[:6000])
for r in lay: print(r['field'], r['Gamma_over_m'], r['tau_s'], r['Omega_tau_at_max_churn'], r['state'])
print(out['partB_multifield']['fermi_ratio_check'])
