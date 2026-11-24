"""Stage 5b: bulk gravity G_N(r), topology vs alpha, 'crystalline' lattice compactifications, FFZ term size."""
import json, os, numpy as np
from math import gamma, pi, lgamma, exp, log
from scipy.special import jn_zeros
CK = 'results/stage5.json'
out = json.load(open(CK)) if os.path.exists(CK) else {}
def save(): json.dump(out, open(CK, 'w'), indent=1, default=float)
HBARC = 1.973269804e-16; MPLR = 2.435e18
def R_of(n, MD_TeV): return (MPLR**2/(MD_TeV*1e3)**(n+2))**(1/n)*HBARC
LIM = {4: (7.1, 'ATLAS'), 5: (6.4, 'ATLAS'), 6: (5.9, 'ATLAS'), 7: (5.2, 'CMS PRD 104, 092008? arXiv:2107.13021 Table 3')}

# ---------------- shell counts of Z^n (theta series by convolution)
def zn_shells(n, Nmax):
    one = np.zeros(Nmax+1); k = np.arange(0, int(np.sqrt(Nmax))+1)
    one[k**2] += 1; one[k[1:]**2] += 1
    c = one.copy()
    for _ in range(n-1):
        c = np.convolve(c, one)[:Nmax+1]
    return c   # c[N] = #k in Z^n with |k|^2 = N
def geff_add(n, x):
    """G_eff/G at r = x R for Z^n torus: 1 + (4/3) sum_{k!=0} exp(-|k| x); continuum asymptote for small x"""
    if x >= 0.5:
        Nmax = int((60/x)**2)+1
        c = zn_shells(n, Nmax); N = np.arange(Nmax+1)
        return 1 + (4/3)*np.sum(c[1:]*np.exp(-np.sqrt(N[1:])*x))
    S = 2*pi**(n/2)/gamma(n/2)
    return 1 + (4/3)*(S*gamma(n)/x**n - 1)
if 'geff' not in out:
    res = {}
    r = np.logspace(-20, 0, 201)
    for n in [4, 5, 6, 7, 9]:
        MD = LIM.get(n, (6.0, 'n=9 benchmark (no verified limit)'))[0]
        R = R_of(n, MD)
        g = np.array([geff_add(n, ri/R) for ri in r])
        # continuity check at x=0.5
        S = 2*pi**(n/2)/gamma(n/2)
        cont = 1 + (4/3)*(S*gamma(n)/0.5**n - 1)
        r_dev3 = float(r[np.where(g-1 > 1e-3)[0].max()]) if np.any(g-1 > 1e-3) else None
        r_dev1 = float(r[np.where(g-1 > 1)[0].max()]) if np.any(g-1 > 1) else None
        res[f'ADD_n{n}'] = dict(M_D_TeV=MD, R_m=R, r_dev_1em3_m=r_dev3, r_order_one_m=r_dev1,
                                check_x05_direct_vs_continuum=[geff_add(n, 0.5), cont], curve=g.tolist())
    # RS1 remnant: Lambda_pi = 11.74 TeV (k/M=0.1, m1=4.5 TeV); couplings 1/Lambda_pi per vertex (ASSUMED standard)
    Lam = 11.744e3; kM = 0.1
    xs = jn_zeros(1, 40); m = xs*kM*Lam
    a_i = (4/3)*(MPLR/Lam)**2
    gr = np.array([1 + a_i*np.sum(np.exp(-m*ri/HBARC)) for ri in r])
    res['RS1'] = dict(Lambda_pi_TeV=Lam/1e3, m1_TeV=m[0]/1e3, alpha_per_mode=a_i, range_m=HBARC/m[0],
                      r_dev_1em3_m=float(r[np.where(gr-1 > 1e-3)[0].max()]),
                      r_dev_analytic_m=HBARC/m[0]*np.log(a_i/1e-3), curve=gr.tolist())
    res['r_grid_m'] = r.tolist()
    out['geff'] = res; save()

# ---------------- topology vs alpha (n = 4..7), equal volume (2 pi R)^n
if 'topology_alpha' not in out:
    rows = []
    for n in range(4, 8):
        MD = LIM[n][0]; R = R_of(n, MD); Vol = (2*pi*R)**n
        Om = 2*pi**((n+1)/2)/gamma((n+1)/2)            # area of unit S^n
        Rs = (Vol/Om)**(1/n)
        rows.append(dict(n=n, M_D_TeV=MD, R_torus_m=R,
            torus_Zn=dict(N1=2*n, alpha=8*n/3, lam_m=R),
            orbifold_TnZ2_fixed_point=dict(N1=n, alpha=8*n/3, lam_m=R, note='half the modes, each with |psi(0)|^2=2'),
            orbifold_TnZ2_generic_point_avg=dict(alpha=4*n/3, note='alpha = (8/3) sum_i cos^2(y_i/R): ranges 0..8n/3 with brane position'),
            one_large_dim_anisotropic=dict(alpha=8/3, note='only the largest radius sets the first shell'),
            sphere_Sn_flux=dict(N1=n+1, alpha=4*(n+1)/3, R_sphere_m=Rs, lam_m=Rs/np.sqrt(n), note='m1^2=n/Rs^2 (l=1), Freund-Rubin type, ASSUMED stabilized'),
            twisted_antiperiodic_bulk_fermion=dict(zero_mode=False, m1_over_inv_R=0.5, note='not allowed for the graviton itself (metric must be periodic); relevant for bulk nu_R')))
    out['topology_alpha'] = rows; save()

# ---------------- crystalline: 9D lattices (13D = 3+1 brane + 9 compact)
def lll_free_enum(B, C):
    """all nonzero integer combos x with |x B|^2 <= C (Fincke-Pohst on Gram)"""
    G = B @ B.T; n = len(G)
    Q = np.zeros((n, n)); A = G.copy().astype(float)
    # LDL-type decomposition: x^T G x = sum_i q_ii (x_i + sum_{j>i} q_ij x_j)^2
    for i in range(n):
        Q[i, i] = A[i, i] - sum(Q[k, k]*Q[k, i]**2 for k in range(i))
        for j in range(i+1, n):
            Q[i, j] = (A[i, j] - sum(Q[k, k]*Q[k, i]*Q[k, j] for k in range(i))) / Q[i, i]
    res = []
    x = np.zeros(n, dtype=int)
    def rec(i, rem):
        if i < 0:
            if np.any(x): res.append(x.copy())
            return
        c = -sum(Q[i, j]*x[j] for j in range(i+1, n))
        w = np.sqrt(max(rem, 0)/Q[i, i])
        for xi in range(int(np.ceil(c - w - 1e-9)), int(np.floor(c + w + 1e-9))+1):
            x[i] = xi
            t = Q[i, i]*(xi - c)**2
            if t <= rem + 1e-9: rec(i-1, rem - t)
        x[i] = 0
    rec(n-1, C)
    return np.array(res)
def shells(B, k=3):
    G = B @ B.T
    m = min(np.diag(G))
    C = 1.6*m
    while True:
        vs = lll_free_enum(B, C)
        nr = np.round(np.einsum('ij,jk,ik->i', vs, G, vs), 6)
        u, c = np.unique(nr, return_counts=True)
        if len(u) >= k: break
        C *= 1.6
    return [(float(a), int(b)) for a, b in zip(u[:k], c[:k])]
E8 = np.array([[2,0,0,0,0,0,0,0],[-1,1,0,0,0,0,0,0],[0,-1,1,0,0,0,0,0],[0,0,-1,1,0,0,0,0],
               [0,0,0,-1,1,0,0,0],[0,0,0,0,-1,1,0,0],[0,0,0,0,0,-1,1,0],[.5]*8])
L9 = np.zeros((9, 9)); L9[:8, :8] = E8; L9[8, 0] = 1.0; L9[8, 8] = 1.0          # E8 layers shifted by deep hole h=(1,0..0), spacing 1
E8C = np.zeros((9, 9)); E8C[:8, :8] = E8; E8C[8, 8] = np.sqrt(2)               # E8 + circle of length sqrt2
Z9 = np.eye(9)
if 'lattices9' not in out:
    rows = []
    for name, B in [('Z9 hypercubic', Z9), ('Lambda9 laminated', L9), ('E8 x circle', E8C)]:
        n = 9; detV = abs(np.linalg.det(B))
        direct = shells(B)
        Bd = np.linalg.inv(B).T
        dual = shells(Bd)
        for MD in [6.0, 10.0]:
            R = R_of(9, MD); L = 2*pi*R/detV**(1/9)        # equal volume (2 pi R)^9
            m1 = 2*pi*np.sqrt(dual[0][0])/L * HBARC       # GeV
            rows.append(dict(lattice=name, det_volume=detV, direct_shells_norm_count=direct, dual_shells_norm_count=dual,
                             M_D_TeV=MD, R_ref_m=R, m1_GeV=m1, lam_m=HBARC/m1, N1=dual[0][1], alpha=4/3*dual[0][1],
                             next_shells_mass_ratio=[float(np.sqrt(s[0]/dual[0][0])) for s in dual[1:]],
                             next_shells_count=[s[1] for s in dual[1:]]))
    # sphere S^9 at equal volume
    for MD in [6.0, 10.0]:
        R = R_of(9, MD); Vol = (2*pi*R)**9; Om = 2*pi**5/gamma(5); Rs = (Vol/Om)**(1/9)
        rows.append(dict(lattice='S9 sphere (flux-stabilized, ASSUMED)', M_D_TeV=MD, R_sphere_m=Rs,
                         m1_GeV=HBARC*np.sqrt(9)/Rs, lam_m=Rs/3, N1=10, alpha=40/3, next_shells_count=[54, 210],
                         next_shells_mass_ratio=[np.sqrt(2*10/9), np.sqrt(3*11/9)]))
    out['lattices9'] = rows
    # Weyl-law check: number of KK modes below E is shape independent at leading order
    E = 1e3
    out['lattice_weyl'] = dict(note='N(E) ~ Vol * Omega_9 E^9/(9 (2pi)^9) for every lattice; shape enters only at relative order (m1/E)^2',
                               rel_shape_effect_at_1TeV=[(r['m1_GeV']/E)**2 for r in rows if r['M_D_TeV'] == 6.0])
    save()

# ---------------- FFZ brane-localized curvature term in the bulk track
if 'ffz_brane_term' not in out:
    S = 265.0; gapf = 0.009651/(1/24)
    rows = []
    for n in [4, 5, 6, 7, 9]:
        MD = LIM.get(n, (6.0,))[0]*1e3
        Mind2 = S*MD**2/(16*pi**2)
        r_ind = (Mind2/MD**(n+2))**(1/n)*HBARC
        rows.append(dict(n=n, M_D_TeV=MD/1e3, M_ind2_over_Mbar2=Mind2/MPLR**2, r_ind_m=r_ind, r_ind_times_MD=r_ind/(HBARC/MD),
                         gap_shift_factor_on_r_ind=(1+gapf)**(1/n)))
    out['ffz_brane_term'] = dict(rows=rows, gap_fraction=gapf,
        atlas_theory_unc='ATLAS 2021: PDF 11-43%, scale 23-36% on ADD signal (verified, PRD 103 112006)',
        kappa_sign='repo: 1/G_ind proportional to C (ffz_induced_collapse.py, from Gamma > -R/(16 pi G)); kappa enters additively dC=kappa (ffz_scheme_invariant.py); positive prefactor ASSUMED (ffz_running_gind.py). So kappa=-0.00965 lowers 1/G by 23% of the 3-loop piece; a bulk volume term only raises 1/G.')
    save()

# ---------------- figures
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
g = out['geff']; r = np.array(g['r_grid_m'])
fig, ax = plt.subplots(figsize=(9, 6))
for k, v in g.items():
    if k == 'r_grid_m': continue
    lab = k + (f" (M_D={v['M_D_TeV']} TeV, R={v['R_m']:.1e} m)" if 'R_m' in v else f" (m1={v['m1_TeV']:.1f} TeV)")
    ax.loglog(r, np.array(v['curve'])-1+1e-30, label=lab)
ax.axhline(1e-3, color='gray', ls=':'); ax.axvline(52e-6, color='k', ls='--'); ax.text(6e-5, 1e3, 'torsion 52 um', fontsize=8)
ax.set_ylim(1e-6, 1e35); ax.set_xlabel('r [m]'); ax.set_ylabel('G_eff(r)/G_N - 1')
ax.set_title('Stage 5: bulk-gravity track, deviation from Newton (at the largest allowed size)'); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig('figures/stage5_geff.png', dpi=110)
fig, ax = plt.subplots(figsize=(9, 5))
rows = [x for x in out['lattices9'] if x['M_D_TeV'] == 6.0]
for i, x in enumerate(rows):
    ms = [1.0] + x['next_shells_mass_ratio']; cs = [x['N1']] + x['next_shells_count']
    ax.bar(np.array(ms)*x['m1_GeV'] + 0.03*i, cs, width=0.05, label=f"{x['lattice']}: alpha={x['alpha']:.1f}, lambda={x['lam_m']:.1e} m")
ax.set_yscale('log'); ax.set_xlabel('KK graviton mass [GeV] (equal volume, M_D=6 TeV, n=9)'); ax.set_ylabel('modes in shell')
ax.set_title("'Crystalline' 9D compactifications: first KK shells"); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig('figures/stage5_lattice_shells.png', dpi=110)
print(json.dumps({k: v for k, v in out.items() if k not in ('geff',)}, indent=1, default=float)[:9000])
print({k: {kk: vv for kk, vv in v.items() if kk != 'curve'} for k, v in out['geff'].items() if k != 'r_grid_m'})
