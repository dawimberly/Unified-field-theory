"""Stage 2: D-dimensional hyperfields. D = number of SPACE dimensions; hidden sphere is S^D in R^{D+1}.
(a) 1/r potential in D dims (Fock-Alliluev SO(D+1) case)   (b) Gauss-law potential 1/r^{D-2}
(c) 3D slice / projection tests against observed hydrogen   (d) compact extra dimensions: how small?"""
import json, sys, os, time
T0=time.time()
def log(*a): print(f'[{time.time()-T0:7.1f}s]', *a, flush=True)
import numpy as np
from scipy.integrate import quad
from scipy.special import eval_gegenbauer, gamma, gammaln
sys.path.insert(0, os.path.dirname(__file__))
from hflib import *
from stage1_fock_s3 import kernel_eig_S  # re-uses the quadrature (stage1 re-runs on import; cheap)

DS = list(range(3, 14)) + [25, 26]
out = {'D_list': DS}

# ---------------- (a) 1/r in D dims: three independent routes to E_n --------------------------
def p0_fock(D, k, Z=1.0):
    cD = 2**(D-1) * np.pi**((D-1)/2) * gamma((D-1)/2)       # FT of 1/r in D dims = cD/|q|^{D-1}
    return Z * cD * kernel_eig_S(D, k) / (2*np.pi)**D       # derived in HYPERFIELD_RESULTS.md
a = []
for D in DS:
    rmax = 60*(4 + (D-3)/2)**2 + 200
    rec = dict(D=D, n=[], E_fock=[], E_formula=[], E_radial_l0=[], deg_numeric=[], deg_formula=[], l_content=[])
    ev0 = radial_levels(D, 0, lambda r: -1/r, nlev=4, rmax=rmax, N=12000, sigma=-1.0)
    for n in range(1, 5):
        rec['n'].append(n)
        rec['E_fock'].append(float(-p0_fock(D, n-1)**2/2))
        rec['E_formula'].append(E_coulomb(D, n))
        rec['E_radial_l0'].append(float(ev0[n-1]))
        rec['deg_numeric'].append(int(harmonic_dim_numeric(D+1, n-1)))
        rec['deg_formula'].append(int(harmonic_dim_formula(D+1, n-1)))
        rec['l_content'].append({f'l={l}': int(harmonic_dim_formula(D, l)) for l in range(n)})
    # check l-independence (hidden SO(D+1)) : 4f-like state l=3, n=4
    if D <= 13:
        rec['E_radial_n4_l3'] = float(radial_levels(D, 3, lambda r: -1/r, nlev=1, rmax=rmax, N=12000, sigma=-1.0)[0])
    a.append(rec); log('coulomb D', D)
out['coulomb_1_over_r'] = a

# observed scale-free ratio f(1S-2S)/f(1S-3S)
f12 = 2466061413187.035e3   # Hz, Parthey et al. PRL 107, 203001 (2011), +-10 Hz
f13 = 2922743278665.79e3    # Hz, Grinin et al. Science 370, 1061 (2020), +-0.72 kHz
R_obs = f12 / f13
cmp = []
for D in DS:
    nu = lambda n: n + (D-3)/2
    R_D = (1/nu(1)**2 - 1/nu(2)**2) / (1/nu(1)**2 - 1/nu(3)**2)
    E1 = -RYD_INF_EV * MU_RATIO / nu(1)**2      # no fit: m_e, alpha fixed
    cmp.append(dict(D=D, ratio_12_over_13=R_D, rel_dev_from_obs=R_D/R_obs - 1,
                    E1_eV=E1, f12_pred_Hz=(1/nu(1)**2 - 1/nu(2)**2)*RYD_INF_EV*MU_RATIO/4.135667696e-15,
                    shell2_states=int(harmonic_dim_formula(D+1, 1)), shell3_states=int(harmonic_dim_formula(D+1, 2))))
out['observed_ratio'] = R_obs
out['vs_observed'] = cmp

# ---------------- (c) 3D slice and 3D projection of D-dim clouds ------------------------------
def norm3(f):
    return np.sqrt(quad(lambda r: f(r)**2 * r**2, 0, np.inf, limit=400)[0])
sl = []
for D in DS:
    for (n, l) in [(1,0),(2,0),(2,1),(3,2),(4,3)]:
        f3 = lambda r: R_coulomb(3, n, l, r)
        fD = lambda r: R_coulomb(D, n, l, r)     # psi_D = fD(rho) * Y_lm(x1,x2,x3/r3)  (valid D-dim harmonic)
        # slice x4..xD = 0
        ov = quad(lambda r: f3(r)*fD(r)*r**2, 0, np.inf, limit=400)[0] / (norm3(f3)*norm3(fD))
        rmean_slice = quad(lambda r: fD(r)**2*r**3, 0, np.inf, limit=400)[0] / norm3(fD)**2
        rmean_3d = quad(lambda r: f3(r)**2*r**3, 0, np.inf, limit=400)[0] / norm3(f3)**2
        # projection (marginal density) : rho3(r) = |r^l Y|^2 * Int |g(sqrt(r^2+y^2))|^2 y^{D-4} dy  where fD = r^l g
        if D > 3:
            k = D - 3
            Sk = 2*np.pi**(k/2)/gamma(k/2)
            g = lambda rho: R_coulomb(D, n, l, rho) / rho**l
            tg, wt = np.polynomial.legendre.leggauss(400)
            t = (tg + 1)/2; wt = wt/2
            sc = 2*(n + (D-3)/2)                      # natural length of the D-dim state
            y = sc*t/(1-t); wy = sc*wt/(1-t)**2        # y in (0, inf)
            rr = np.linspace(1e-6, 60*n + 20*(D-3), 800)
            RHO = np.sqrt(rr[:, None]**2 + y[None, :]**2)
            p3 = rr**(2*l) * Sk * np.sum(g(RHO)**2 * y**(k-1) * wy, axis=1)
            q3 = f3(rr)**2
            wr = rr**2
            bc = np.trapezoid(np.sqrt(p3*q3)*wr, rr) / np.sqrt(np.trapezoid(p3*wr, rr)*np.trapezoid(q3*wr, rr))
            proj = float(bc**2)
        else:
            proj = 1.0
        sl.append(dict(D=D, orbital=f"{n}{'spdf'[l]}", slice_fidelity=float(ov**2), projection_fidelity=proj,
                       rmean_slice=float(rmean_slice), rmean_3d=float(rmean_3d)))
        log('slice', D, n, l)
out['slice_projection'] = sl

# ---------------- (b) Gauss-law potential -g/r^{D-2}, regularized at r_c -------------------------
GLC = 'results/stage2_gauss_cache.json'
gl = json.load(open(GLC)) if os.path.exists(GLC) else []
for D in ([] if gl else [4, 5, 6, 7, 9, 10, 11, 13]):
    for g in ([0.4, 1.0] if D == 4 else [1.0]):
        row = dict(D=D, g=g, rc=[], E0=[])
        for rc in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3]:
            V = lambda r, rc=rc: -g / np.maximum(r, rc)**(D-2)
            sig = 1.3 * (-g / rc**(D-2)) - 1.0
            try:
                sig = -1.0 if (D == 4 and g < 0.5) else sig
                e0 = radial_levels(D, 0, V, nlev=1, rmin=rc*1e-3, rmax=400, N=12000, sigma=sig)[0]
            except Exception as ex:
                e0 = float('nan')
            row['rc'].append(rc); row['E0'].append(float(e0))
        gl.append(row); log('gauss', D, g)
json.dump(gl, open(GLC, 'w'))
out['gauss_law'] = gl

# ---------------- (d) compact extra dimensions ----------------------------------------------------
# ASSUMED model: photon (and electron) see N extra flat circles of radius R; lowest KK shell has 2N modes:
#   V(r) = -1/r * (1 + 2N exp(-r/R))   [a.u.]
psi0sq = lambda n: 1/(np.pi*n**3)          # |psi_nS(0)|^2 in a.u.
def dE_first_order(n, R, N):
    f = lambda r: R_coulomb(3, n, 0, r)
    nn = quad(lambda r: f(r)**2*r**2, 0, np.inf, limit=400)[0]
    return -2*N*quad(lambda r: f(r)**2*r*np.exp(-r/R), 0, 50*n, points=[R, 10*R], limit=800)[0]/nn
comp = dict(checks=[], bounds=[])
for R in [1e-2, 1e-3]:
    for N in [1]:
        V = lambda r: -1/r*(1 + 2*N*np.exp(-r/R))
        ev = radial_levels(3, 0, V, nlev=2, rmin=1e-8, rmax=300, N=16000, sigma=-1.5)
        ev_ref = radial_levels(3, 0, lambda r: -1/r, nlev=2, rmin=1e-8, rmax=300, N=16000, sigma=-1.5)
        comp['checks'].append(dict(R_a0=R, N=N, dE1S_nonpert=float(ev[0]-ev_ref[0]), dE1S_first=dE_first_order(1, R, N),
                                   dE1S_contact=-2*N*4*np.pi*R**2*psi0sq(1),
                                   dE2S_nonpert=float(ev[1]-ev_ref[1]), dE2S_first=dE_first_order(2, R, N)))
f12_au = 0.5*(1 - 1/4)   # 1S-2S in hartree (Bohr)
for N in [1, 2, 6, 7, 8, 22]:
    # contact limit: d(f12)/f12 = 2N*4pi R^2 (|psi1|^2-|psi2|^2)/f12 = 2N*4pi R^2*(7/8)/pi / (3/8) = (56/3) N R^2
    coef = 2*N*4*np.pi*(psi0sq(1)-psi0sq(2))/f12_au
    Rmax_naive = np.sqrt(4.2e-15/coef) * A0_M
    # realistic: shift equals delta(r_p^2) = -12 N R^2 ; tolerance on r_p^2 from muonic-H vs e-p scattering spread
    drp2_fm2 = 2*0.8409*0.014     # PRad 0.831(14) fm vs muonic 0.84087 fm -> ~0.024 fm^2 (ASSUMED tolerance)
    Rmax_real = np.sqrt(drp2_fm2/(12*N))*1e-15
    comp['bounds'].append(dict(N_extra=N, coef=coef, Rmax_m_naive_1S2S=Rmax_naive, KK_mass_GeV_naive=HBARC_EVM/Rmax_naive/1e9,
                               Rmax_m_proton_radius_limited=Rmax_real, KK_mass_GeV_real=HBARC_EVM/Rmax_real/1e9))
# brane world: only gravity in the bulk (ADD). Gravity/Coulomb in H, enhanced by (R/r)^N for r < R
G, me, mp, ke2 = 6.67430e-11, 9.1093837015e-31, 1.67262192369e-27, 2.307077e-28  # SI, ke2 = e^2/(4 pi eps0)
grav_ratio = G*me*mp/ke2
add = []
for N, R in [(2, 30e-6), (3, 1e-9), (6, 1e-13)]:
    # 1S: <(R/r)^N / r> / <1/r> ; for N>=2 diverges at small r -> cut at proton radius 0.84 fm (ASSUMED)
    rp = 0.84e-15 / A0_M
    Ra = R / A0_M
    num = quad(lambda r: 4*r**2*np.exp(-2*r)*(Ra/r)**N/r, rp, 1e3, limit=400, points=[1e-4, 1e-2, 1])[0]
    add.append(dict(N=N, R_m=R, fractional_shift_1S=grav_ratio*num/1.0))
comp['add_gravity_only'] = dict(grav_to_coulomb=grav_ratio, cases=add)
out['compactification'] = comp

json.dump(out, open('results/stage2.json', 'w'), indent=1)

# ---------------- figure ---------------------------------------------------------------------------
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(17, 5))
for j, (n, l) in enumerate([(1,0),(2,1),(3,2)]):
    rr = np.linspace(0, 40*n, 2000)
    for D, c in zip([3,4,5,9,10,11,13], plt.cm.viridis(np.linspace(0,1,7))):
        f = R_coulomb(D, n, l, rr); f = f/np.sqrt(np.trapezoid(f**2*rr**2, rr))
        ax[j].plot(rr, f**2*rr**2, color=c, lw=2.5 if D == 3 else 1.3, label=f'D={D}' + (' (actual H)' if D == 3 else ' slice'))
    ax[j].set_title(f'{n}{"spdf"[l]}: radial density of 3D slice of D-dim 1/r hyperfield'); ax[j].set_xlabel('r [Bohr radii]'); ax[j].set_ylabel('r$^2$|R|$^2$ (3D-normalized)')
    ax[j].legend(fontsize=8)
fig.tight_layout(); fig.savefig('figures/stage2_Dslices_vs_actual.png', dpi=110)

fig, ax = plt.subplots(figsize=(9, 5))
for i, D in enumerate(DS):
    for n in range(1, 5):
        E = E_coulomb(D, n) * HARTREE_EV * MU_RATIO
        ax.hlines(E, i-0.35, i+0.35, color='k' if D == 3 else 'C0')
        ax.text(i+0.37, E, str(int(harmonic_dim_formula(D+1, n-1))), fontsize=6, va='center')
ax.set_xticks(range(len(DS))); ax.set_xticklabels([str(D) for D in DS]); ax.set_xlabel('space dimensions D (hidden sphere S$^D$)')
ax.set_ylabel('E$_n$ [eV]'); ax.set_title('1/r hyperfield spectra n=1..4 (numbers = states per shell). Observed H = D=3 (black)')
fig.tight_layout(); fig.savefig('figures/stage2_spectra_vs_D.png', dpi=110)
print('done')
