"""Talbot-Lau visibility model of Pedalino et al. (arXiv:2507.21211, Methods Eqs. 1-9), quantum and classical.
Shared by stage22_run.py. CITED formulas; nothing tuned."""
import numpy as np
from scipy.special import jv
H = 6.62607015e-34; HBAR = H / (2 * np.pi); C = 299792458.0; EPS0 = 8.8541878128e-12; U = 1.66053906660e-27
LAM_L = 266e-9          # CITED grating laser wavelength, d = 133 nm
ALPHA_ATOM = -4 * np.pi * EPS0 * 4.5e-30   # CITED alpha_266 per Na atom = -4 pi eps0 (4.5 +- 0.5) A^3
M_NA = 22.98976928 * U

def sigma_ion(m_kDa):
    """CITED sigma_ion,266 = (0.537 m/kDa - 1.5) 1e-20 m^2 (fit to 100-200 kDa; extrapolated above)."""
    return max(0.537 * m_kDa - 1.5, 0.0) * 1e-20

def grating_params(m_kg, v, P, w=600e-6, alpha_scale=1.0, sigma_scale=1.0):
    """Eqs. (3),(4): n0 and phi0 for incident power P, waist w along the beam."""
    m_kDa = m_kg / U / 1e3
    n0 = 8 * sigma_scale * sigma_ion(m_kDa) * P * LAM_L / (np.sqrt(2 * np.pi) * H * C * w * v)
    alpha = alpha_scale * ALPHA_ATOM * m_kg / M_NA
    phi0 = np.sqrt(8 / np.pi) * alpha * P / (HBAR * C * EPS0 * w * v)
    return n0, phi0

def Bn(n, xi, n0, phi0, quantum=True):
    if quantum:
        zc = phi0 * np.sin(np.pi * xi); zi = 0.5 * n0 * np.cos(np.pi * xi)
    else:
        zc = phi0 * np.pi * xi; zi = 0.5 * n0 + 0 * xi
    a = zc + zi + 1e-12; b = zc - zi
    ratio = (b / a).astype(complex) if isinstance(b, np.ndarray) else complex(b / a)
    arg = np.sign(a) * np.sqrt(np.asarray(zc ** 2 - zi ** 2, dtype=complex))
    return np.exp(-n0 / 2) * np.power(ratio, n / 2) * jv(n, arg)

def signal(xi, g1, g2, g3, quantum=True):
    """Eq. (2): S_l = B_-l^(1)(0) B_2l^(2)(l xi) B_l^(3)(0); returns S0, S1 (real)."""
    S0 = Bn(0, 0.0, *g1, quantum) * Bn(0, 0.0, *g2, quantum) * Bn(0, 0.0, *g3, quantum)
    S1 = Bn(-1, 0.0, *g1, quantum) * Bn(2, xi, *g2, quantum) * Bn(1, 0.0, *g3, quantum)
    return np.real(S0), np.real(S1)

def visibility(m_kg, v_mean, dv_sigma, L, d, P1=None, P2=None, P3=None, n0_13=None, phi2=None, w=600e-6,
               dm_frac=0.0, alpha_scale=1.0, sigma_scale=1.0, quantum=True, extra_coh=1.0, nv=41, nm=9):
    """Velocity- (Gaussian sigma dv_sigma) and mass- (uniform +-dm_frac) averaged V = 2|<S1>|/<S0>.
    Either powers (P1,P2,P3) or design parameters at the mean (n0 of G1=G3, phi0 of G2) are given;
    in the latter case the powers are fixed and n0, phi0 scale with 1/v and with m as the CITED cross sections do."""
    vs = v_mean + dv_sigma * np.linspace(-3, 3, nv) if dv_sigma > 0 else np.array([v_mean])
    wv = np.exp(-0.5 * ((vs - v_mean) / max(dv_sigma, 1e-30)) ** 2) if dv_sigma > 0 else np.ones(1)
    ms = m_kg * (1 + dm_frac * np.linspace(-1, 1, nm)) if dm_frac > 0 else np.array([m_kg])
    if P1 is None:   # convert design (n0_13, phi2) at mean m, v into powers
        n_unit, p_unit = grating_params(m_kg, v_mean, 1.0, w, 1.0, 1.0)
        P1 = P3 = n0_13 / n_unit; P2 = abs(phi2 / p_unit)
    S0t = 0.0; S1t = 0.0
    for m in ms:
        for v, ww in zip(vs, wv):
            xi = L * H / (m * v * d ** 2)
            g1 = grating_params(m, v, P1, w, alpha_scale, sigma_scale); g2 = grating_params(m, v, P2, w, alpha_scale, sigma_scale)
            g3 = grating_params(m, v, P3, w, alpha_scale, sigma_scale)
            s0, s1 = signal(xi, g1, g2, g3, quantum)
            S0t += ww * s0; S1t += ww * s1 * extra_coh
    return 2 * abs(S1t) / S0t, S0t / (wv.sum() * len(ms)), (P1, P2, P3)
