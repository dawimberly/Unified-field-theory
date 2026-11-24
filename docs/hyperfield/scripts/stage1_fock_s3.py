"""Stage 1 (D=3 space, hidden S^3): Fock's stereographic construction, verified numerically.
Outputs results/stage1.json and figures/stage1_slices_vs_actual.png"""
import json, sys, os
import numpy as np
from scipy.integrate import quad, simpson
from scipy.special import eval_gegenbauer, spherical_jn, sph_harm_y, gamma, factorial
sys.path.insert(0, os.path.dirname(__file__))
from hflib import *

out = {}
# ---------------------------------------------------------------- 1. Fock kernel eigenvalues on S^3
# Fock: (p^2/2 - E) phi(p) = (Z/2pi^2) Int phi(p')/|p-p'|^2 d^3p'.  With p0^2=-2E and stereographic map
# the equation becomes  Y(xi) = Z/(2 pi^2 p0) Int_{S^3} Y(xi')/|xi-xi'|^2 dOmega'.
# Funk-Hecke gives the kernel eigenvalue on degree-k harmonics; we compute it by quadrature (no formula).
def kernel_eig_S(Dsph, k):
    """eigenvalue of 1/|xi-xi'|^{Dsph-1} on degree-k harmonics on S^{Dsph} (Dsph = space dim D)."""
    a = (Dsph - 1) / 2                 # Gegenbauer index for S^Dsph
    volS = 2 * np.pi ** (Dsph / 2) / gamma(Dsph / 2)  # vol(S^{Dsph-1})
    f = lambda t: (eval_gegenbauer(k, a, t) / eval_gegenbauer(k, a, 1.0)) * (1 + t) ** ((Dsph - 2) / 2) / 2 ** ((Dsph - 1) / 2)
    # integrand (1-t^2)^{(D-2)/2} / (2(1-t))^{(D-1)/2} = (1+t)^{(D-2)/2} (1-t)^{-1/2} / 2^{(D-1)/2}
    val, err = quad(f, -1, 1, weight='alg', wvar=(0.0, -0.5), limit=200, epsabs=1e-14, epsrel=1e-13)
    return volS * val
def main():
    global out
    Z = 1.0
    rows = []
    for k in range(0, 4):
        lam = kernel_eig_S(3, k)
        p0 = Z * lam / (2 * np.pi**2)       # from 1 = Z lam/(2 pi^2 p0)
        E = -p0**2 / 2
        n = k + 1
        rows.append(dict(n=n, kernel_eig=lam, kernel_eig_exact=2*np.pi**2/n, p0=p0,
                         E_hartree=E, E_eV_Rinf=E*HARTREE_EV, E_eV_reduced_mass=E*HARTREE_EV*MU_RATIO,
                         target_eV=-13.6057/n**2))
    out['fock_kernel'] = rows

    # ---------------------------------------------------------------- 2. degeneracy n^2 by counting harmonics in 4 variables
    out['degeneracy'] = [dict(n=n, harmonics_deg_n_minus_1_in_R4=int(harmonic_dim_numeric(4, n-1)),
                              sum_2l_plus_1=sum(2*l+1 for l in range(n)), n2=n*n) for n in range(1, 6)]

    # ---------------------------------------------------------------- 3. independent position-space solve (no SO(4) used)
    fd = []
    for l in range(4):
        ev = radial_levels(3, l, lambda r: -1/r, nlev=4-l if l < 4 else 1, rmax=600, N=12000)
        for i, e in enumerate(ev):
            n = l + 1 + i
            if n <= 4:
                fd.append(dict(n=n, l=l, E_numeric=float(e), E_bohr=-0.5/n**2, rel_err=float(e/(-0.5/n**2)-1)))
    out['position_space_fd'] = fd

    # ---------------------------------------------------------------- 4. momentum-space fidelities
    def R3(n, l, r):
        """normalized 3D hydrogen radial function"""
        R = R_coulomb(3, n, l, r)
        nrm = np.sqrt(quad(lambda s: (R_coulomb(3, n, l, s)*s)**2, 0, np.inf, limit=200)[0])
        return R / nrm
    def phi_FT(n, l, p):
        """radial momentum wavefunction by numerical Hankel transform of the position orbital"""
        rgrid = np.linspace(0, 30*n + 30, 100001)
        Rr = R3(n, l, rgrid)
        out = np.empty(len(p))
        for i0 in range(0, len(p), 20):
            J = spherical_jn(l, np.outer(p[i0:i0+20], rgrid))
            out[i0:i0+20] = np.sqrt(2/np.pi) * simpson(J * Rr * rgrid**2, x=rgrid, axis=1)
        return out

    def realY(l, m, th, ph):
        if m == 0:
            return np.real(sph_harm_y(l, 0, th, ph))
        Y = sph_harm_y(l, abs(m), th, ph)
        return np.sqrt(2) * (-1)**m * (np.real(Y) if m > 0 else np.imag(Y))

    # quadrature in p-space: p = p0 tan(chi/2), truncated at p = PMAX*p0 (norm fraction beyond ~1e-12; same domain for both sides)
    PMAX = 200.0
    CHIMAX = 2*np.arctan(PMAX)
    def pgrid(p0, Nchi=160):
        xg, wg = np.polynomial.legendre.leggauss(Nchi)
        chi = (xg + 1) * CHIMAX / 2; wchi = wg * CHIMAX / 2
        p = p0 * np.tan(chi / 2); dp = p0 / (2*np.cos(chi/2)**2) * wchi
        return p, dp

    radial_fid = []
    for n in range(1, 5):
        p0 = 1.0 / n            # DERIVED above from the kernel eigenvalue
        p, dp = pgrid(p0)
        w = p**2 * dp
        for l in range(n):
            cchi = (p0**2 - p**2) / (p0**2 + p**2); schi = 2*p0*p / (p0**2 + p**2)
            fock = (p0**2 + p**2)**-2 * schi**l * eval_gegenbauer(n-l-1, l+1, cchi)
            ft = phi_FT(n, l, p)
            F = np.sum(w*fock*ft)**2 / (np.sum(w*fock**2) * np.sum(w*ft**2))
            radial_fid.append(dict(n=n, l=l, orbital=f"{n}{'spdf'[l]}", fidelity=float(F), one_minus_F=float(1-F)))
    out['radial_fidelity_fock_vs_FT'] = radial_fid

    # unlabeled-subspace test: span{ (p0^2+p^2)^-2 H(xi) : H harmonic poly deg n-1 in R^4 } vs span{FT hydrogen orbitals}
    def harmonic_basis(nvar, deg):
        src = monomials(nvar, deg)
        if deg < 2:
            return src, np.eye(len(src))
        tgt = {m: i for i, m in enumerate(monomials(nvar, deg - 2))}
        L = np.zeros((len(tgt), len(src)))
        for j, m in enumerate(src):
            for k in range(nvar):
                if m[k] >= 2:
                    mm = list(m); mm[k] -= 2
                    L[tgt[tuple(mm)], j] += m[k]*(m[k]-1)
        u, s, vt = np.linalg.svd(L)
        rank = np.sum(s > 1e-10 * s.max())
        return src, vt[rank:].T          # columns = coefficient vectors of harmonic polys

    xc, wc = np.polynomial.legendre.leggauss(20); th = np.arccos(xc)
    ph = np.linspace(0, 2*np.pi, 40, endpoint=False); wph = 2*np.pi/40
    sub = []
    for n in range(1, 5):
        p0 = 1.0/n
        p, dp = pgrid(p0, 100)
        P, TH, PH = np.meshgrid(p, th, ph, indexing='ij')
        W = (p**2*dp)[:, None, None] * wc[None, :, None] * wph
        px, py, pz = P*np.sin(TH)*np.cos(PH), P*np.sin(TH)*np.sin(PH), P*np.cos(TH)
        den = p0**2 + P**2
        xi = [(p0**2 - P**2)/den, 2*p0*px/den, 2*p0*py/den, 2*p0*pz/den]
        mons, C = harmonic_basis(4, n-1)
        monvals = np.stack([np.prod([xi[k]**m[k] for k in range(4)], axis=0) for m in mons], -1)
        A = (monvals @ C) * (den**-2)[..., None]          # Fock-mapped hyperfield functions
        Bl = []
        for l in range(n):
            rad = phi_FT(n, l, p)[:, None, None]
            for m in range(-l, l+1):
                Bl.append(rad * realY(l, m, TH, PH))
        B = np.stack(Bl, -1)
        A = A.reshape(-1, A.shape[-1]); B = B.reshape(-1, B.shape[-1]); W = np.broadcast_to(W, P.shape).reshape(-1)
        def orth(M):
            G = (M * W[:, None]).T @ M
            ev, U = np.linalg.eigh(G)
            return M @ (U / np.sqrt(ev))
        QA, QB = orth(A), orth(B)
        S = (QA * W[:, None]).T @ QB
        sv = np.linalg.svd(S, compute_uv=False)
        # per-orbital captured weight of each actual orbital inside the hyperfield span
        Bn = B / np.sqrt(np.sum(B*B*W[:, None], axis=0))
        cap = np.sum(((QA * W[:, None]).T @ Bn)**2, axis=0)
        sub.append(dict(n=n, dim_hyperfield_span=int(A.shape[-1]), dim_hydrogen_shell=int(B.shape[-1]),
                        min_principal_cos2=float(sv.min()**2), max_principal_cos2=float(sv.max()**2),
                        min_orbital_captured=float(cap.min())))
    out['subspace_fidelity'] = sub

    os.makedirs('results', exist_ok=True)
    json.dump(out, open('results/stage1.json', 'w'), indent=1)

    # ---------------------------------------------------------------- 5. figure: hyperfield slice vs actual (p_x-p_z plane)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    orbs = [(1,0,0),(2,1,0),(3,2,0),(4,3,0),(3,1,0),(4,2,1)]
    fig, ax = plt.subplots(3, len(orbs), figsize=(3*len(orbs), 9))
    for j, (n, l, m) in enumerate(orbs):
        p0 = 1/n; L = 3.0*p0
        g = np.linspace(-L, L, 201); PX, PZ = np.meshgrid(g, g)
        Pm = np.hypot(PX, PZ); TH = np.arccos(np.clip(PZ/np.maximum(Pm,1e-12), -1, 1)); PH = np.where(PX >= 0, 0.0, np.pi)
        den = p0**2 + Pm**2; cchi = (p0**2-Pm**2)/den; schi = 2*p0*Pm/den
        fock = den**-2 * schi**l * eval_gegenbauer(n-l-1, l+1, cchi) * realY(l, m, TH, PH)
        gf = np.linspace(0, np.sqrt(2)*L, 3001); ft = np.interp(Pm, gf, phi_FT(n, l, gf)) * realY(l, m, TH, PH)
        a, b = fock**2/np.max(fock**2), ft**2/np.max(ft**2)
        ax[0,j].imshow(a, extent=[-L,L,-L,L], origin='lower', cmap='magma'); ax[0,j].set_title(f'{n}{"spdf"[l]} (m={m})\nS$^3$ harmonic, stereographic')
        ax[1,j].imshow(b, extent=[-L,L,-L,L], origin='lower', cmap='magma'); ax[1,j].set_title('actual H orbital (Fourier)')
        ax[2,j].imshow(a-b, extent=[-L,L,-L,L], origin='lower', cmap='coolwarm', vmin=-1e-3, vmax=1e-3); ax[2,j].set_title(f'difference (max {np.abs(a-b).max():.1e})')
        for i in range(3): ax[i,j].set_xlabel('$p_x$ [a.u.]'); ax[i,j].set_ylabel('$p_z$')
    fig.suptitle('Stage 1: 4D hyperspherical harmonics on S$^3$ (Fock map, $p_0=1/n$, nothing fitted) vs actual hydrogen orbitals, |phi|$^2$ in the $p_x$-$p_z$ plane')
    fig.tight_layout(); os.makedirs('figures', exist_ok=True); fig.savefig('figures/stage1_slices_vs_actual.png', dpi=110)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
