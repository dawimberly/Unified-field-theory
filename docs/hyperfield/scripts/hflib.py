"""Shared helpers for the hyperfield tests. Atomic units (hbar = m_e = e = 4 pi eps0 = 1)."""
import numpy as np
from math import comb
from itertools import combinations_with_replacement
from scipy.linalg import eigh_tridiagonal
from scipy.special import eval_genlaguerre, gammaln

HARTREE_EV = 27.211386245988          # CODATA 2018
RYD_INF_EV = HARTREE_EV / 2           # 13.605693122994 eV  (infinite nuclear mass)
MU_RATIO = 1 / (1 + 1 / 1836.15267343)  # reduced-mass factor m_r/m_e for hydrogen
A0_M = 5.29177210903e-11              # Bohr radius [m]
HBARC_EVM = 197.3269804e-9            # hbar c [eV m]

# ---------- harmonic-polynomial counting (pure linear algebra, no formula used) ----------
def monomials(nvar, deg):
    return [tuple(np.bincount(c, minlength=nvar)) for c in combinations_with_replacement(range(nvar), deg)] if deg >= 0 else []

def harmonic_dim_numeric(nvar, deg):
    """dim of kernel of the Laplacian P_deg -> P_{deg-2} in nvar variables, computed by matrix rank."""
    src = monomials(nvar, deg)
    if deg < 2:
        return len(src)
    tgt = {m: i for i, m in enumerate(monomials(nvar, deg - 2))}
    L = np.zeros((len(tgt), len(src)))
    for j, m in enumerate(src):
        for k in range(nvar):
            if m[k] >= 2:
                mm = list(m); mm[k] -= 2
                L[tgt[tuple(mm)], j] += m[k] * (m[k] - 1)
    return len(src) - np.linalg.matrix_rank(L)

def harmonic_dim_formula(nvar, deg):
    """dim of degree-deg harmonics on S^{nvar-1}."""
    if deg < 0:
        return 0
    return comb(deg + nvar - 1, nvar - 1) - (comb(deg + nvar - 3, nvar - 1) if deg >= 2 else 0)

# ---------- radial solver on a logarithmic grid ----------
import scipy.sparse as _sp
import scipy.sparse.linalg as _sl
def radial_levels(D, l, V, nlev=4, rmin=1e-7, rmax=400.0, N=8000, sigma=None, dense=False):
    """Lowest eigenvalues of -1/2 Lap_D + V(r) at hyper-angular momentum l in D space dims.
    psi = r^{-(D-1)/2} u ;  r = e^x, u = e^{x/2} v  =>  -v''/2 + (l+(D-2)/2)^2/2 v + r^2 V v = E r^2 v.
    Generalized sparse problem, shift-invert about sigma (must lie below the ground state)."""
    x = np.linspace(np.log(rmin), np.log(rmax), N)
    h = x[1] - x[0]
    r = np.exp(x)
    k2 = (l + (D - 2) / 2) ** 2
    Vr = V(r)
    if sigma is None:
        sigma = -1.0
    o = -0.5 / h**2 * np.ones(N - 1)
    A = _sp.diags([o, 1 / h**2 + k2 / 2 + r**2 * Vr, o], [-1, 0, 1], format='csc')
    B = _sp.diags(r**2, format='csc')
    if dense:
        from scipy.linalg import eigh
        ev = eigh(A.toarray(), B.toarray(), eigvals_only=True, subset_by_index=[0, nlev - 1])
        return np.sort(ev)
    ev = _sl.eigsh(A, k=nlev, M=B, sigma=sigma, which='LM', return_eigenvectors=False, maxiter=3000)
    return np.sort(ev)

# ---------- exact D-dim Coulomb (1/r) radial functions ----------
def R_coulomb(D, n, l, r, Z=1.0):
    """Unnormalized radial factor f(r) multiplying the solid harmonic? No: full radial R(r) with psi = R(r) Y_l(angles).
    R = r^l exp(-Z r/nu) L_{n-l-1}^{(2l+D-2)}(2 Z r/nu), nu = n + (D-3)/2."""
    nu = n + (D - 3) / 2
    return r**l * np.exp(-Z * r / nu) * eval_genlaguerre(n - l - 1, 2 * l + D - 2, 2 * Z * r / nu)

def E_coulomb(D, n, Z=1.0):
    return -Z**2 / (2 * (n + (D - 3) / 2) ** 2)
