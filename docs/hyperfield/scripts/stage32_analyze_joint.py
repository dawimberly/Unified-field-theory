#!/usr/bin/env python3
"""Stage 32 -- ready-to-run analysis of the JOINT (emission angle x bond length) Grundmann et al. 2020 H2 data.

Usage (when the Frankfurt group sends data):
    python stage32_analyze_joint.py JOINT.csv [--out results/stage32_joint_result.json] [--nbins 4] [--edges 0.60,0.70,0.76,0.82,0.96]
    python stage32_analyze_joint.py --selftest          # runs the pipeline on a simulated file

Accepted input (auto-detected from the header, comma or whitespace separated, Zenodo-style 'Variables =' preamble allowed):
  (a) binned histogram: columns containing 'cos(alpha)'/'cos_alpha', 'cos(beta)'/'cos_beta', one of 'KER'/'Kinetic Energy' (eV)
      or 'R' (Angstrom), and 'count';
  (b) event list: same columns without 'count' (each row one event); cos(alpha) is histogrammed with --dcos (default 0.0444,
      the Fig. 3 bin width) and cos(beta) into 15 rows (Fig. 3).
Method (identical to Stage 31, which reproduced the paper's Fig. 3C):
  per (R bin, cos beta row): Poisson-likelihood Gaussian+constant fit of the central fringe within |cos a| <= w_R,
  w_R = 0.4 * sigma_R/sigma_S (window scaled with the fringe width; 0.4 at subset S, Stage 31 headline);
  error on cos(alpha0) from the Fisher matrix (cross-checked against a Poisson bootstrap in stage32_run.py);
  tau = cos(alpha0) * R / v_ph with v_ph = E_e/p_e (paper Eq. 2), R per bin from KER = e^2/(4 pi eps0 R) (paper, reflection approx.);
  weighted line through tau vs cos beta per R bin -> slope s(R).
Hypotheses on s(R) (all in zs per unit cos beta):
  H0     : s = R/c + nd(R)                         (light travel time, paper Eq. 3, + nondipole term of the paper's red line)
  H0free : s = a*R + nd(R)                         (anything that scales with R: travel time, momentum-zero offset, R-scale error)
  H1     : s = a*R + nd(R) + d1                    (do-si-do: extra release term NOT scaling with R)   <- test (B), systematic-free
  H1abs  : s = R/c + nd(R) + d1 + delta*R          (delta ~ N(0, sig_delta): momentum-zero + R-scale priors) <- test (A)
"""
import sys, os, json, argparse
import numpy as np
from math import sqrt, pi
from scipy.optimize import minimize

# ---- constants: CODATA 2022 as fetched this session (physics.nist.gov/cuu/Constants/Table/allascii.txt); see stage32_sources.csv
C = 299792458.0                    # m/s (exact)
BOHR = 5.29177210544e-11           # m
HARTREE_EV = 27.211386245981       # eV
AU_T = 2.4188843265864e-17         # s
E_CH = 1.602176634e-19             # C (exact)
EPS0 = 8.8541878188e-12            # F/m
KE_EVA = E_CH / (4 * pi * EPS0) / 1e-10   # e^2/(4 pi eps0) in eV*Angstrom (DERIVED) = 14.39964...
# ---- paper inputs (CITED, Grundmann et al. 2020, arXiv:2010.08298v1, fetched this session)
EE0 = 735.0          # eV fast electron  "[Ee = (735+-15) eV]"
EPH = 800.0          # eV photon         "photons with an energy of 800 eV"
IP2 = 31.03          # eV                "adiabatic double ionization energy of H2 (31.03 eV)"
R_S = 0.74           # Angstrom          "R = (0.74 +- 0.02) A [selected subset S]"
TAU_ND_S = 271.7 - 246.9   # zs at R_S: red (nondipole+travel) minus blue (travel) line, DIGITIZED this session (fig3C 600 dpi)
W_S = 0.4            # ASSUMED headline window in cos(alpha) at subset S (Stage 31)

DP0_AU = 0.04        # a.u.              "zero point of the momentum distribution ... approximately +-0.04 atomic units"
DR_REL = 0.02 / 0.74 # R scale           "R = (0.74 +- 0.02) A" (used as R-scale uncertainty, ASSUMED interpretation as in Stage 31)

def Ee_of_KER(KER, Eph=EPH, Ee0=EE0, KER0=19.19):
    """Fast-electron energy vs KER and photon energy: E_e = f*(Eph - IP2 - KER), with f = Ee0/(800 - IP2 - KER0) fixed at the measured
    800 eV value (DERIVED from CITED Ee0, IP2; ASSUMED the same sharing fraction at other KER / photon energies); KER0 = median KER of Fig. 2."""
    f = Ee0 / (EPH - IP2 - KER0)
    return f * (Eph - IP2 - np.asarray(KER, float))

def zs_per_cos(R_A, Ee_eV):
    """tau_b [zs] per unit cos(alpha0): R/v_ph, v_ph = E_e/p_e (paper Eq. 2), nonrelativistic p_e (0.14% correction neglected)."""
    pe = np.sqrt(2 * np.asarray(Ee_eV) / HARTREE_EV); vph = (np.asarray(Ee_eV) / HARTREE_EV) / pe
    return (np.asarray(R_A) * 1e-10 / BOHR) / vph * AU_T * 1e21

def tau_Rc(R_A):
    return np.asarray(R_A) * 1e-10 / C * 1e21

def nd_term(R_A, scaling=-1.0):
    """nondipole contribution to the slope (zs). Magnitude DIGITIZED at R_S; R-scaling ASSUMED: -1 (envelope-tilt toy, DERIVED in
    stage32_run.py: cos-shift ~ 1/(p_e R)^2 -> tau ~ 1/R), 0 or +1 as variants."""
    return TAU_ND_S * (np.asarray(R_A) / R_S) ** scaling

# DERIVED default prior for test (A): momentum-zero offset dp0 gives tau = dp0*R/E_e (a.u.), R-scale gives DR_REL * slope (~R)
_tau_dp = DP0_AU * (R_S * 1e-10 / BOHR) / (EE0 / HARTREE_EV) * AU_T * 1e21
SIG_DELTA_DEFAULT = float(sqrt(_tau_dp ** 2 + (DR_REL * 504.5) ** 2))

# ---------------------------------------------------------------- fitting of one fringe
def gmodel(q, x):
    return np.clip(np.exp(q[0]) * np.exp(-0.5 * ((x - q[1]) / np.exp(q[2])) ** 2) + np.exp(q[3]), 1e-9, None)

def fit_fringe(x, y, sig0=0.107):
    def nll(q):
        lam = gmodel(q, x); return np.sum(lam - y * np.log(lam))
    q0 = [np.log(max(y.max() - y.min(), 1.0)), 0.0, np.log(sig0), np.log(max(y.min(), 1.0))]
    r = minimize(nll, q0, method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-9, maxiter=20000))
    return r.x

def fisher_mu_err(A, mu, sg, B, x):
    """sigma(mu) from the Poisson Fisher matrix of A*exp(-(x-mu)^2/2sg^2)+B (params A, mu, sg, B)."""
    g = np.exp(-0.5 * ((x - mu) / sg) ** 2); lam = np.clip(A * g + B, 1e-12, None)
    J = np.vstack([g, A * g * (x - mu) / sg ** 2, A * g * (x - mu) ** 2 / sg ** 3, np.ones_like(x)])
    F = (J / lam) @ J.T
    return float(sqrt(np.linalg.inv(F)[1, 1]))

def fit_row(x, y, win, sig0=0.107):
    m = np.abs(x) <= win; xx, yy = x[m], y[m]
    if yy.sum() < 30: return np.nan, np.nan, None
    q = fit_fringe(xx, yy, sig0)
    A, mu, sg, B = np.exp(q[0]), q[1], np.exp(q[2]), np.exp(q[3])
    return mu, fisher_mu_err(A, mu, sg, B, xx), (A, mu, sg, B)

def wline(x, y, e):
    """weighted fit y = a + b x; returns b, sigma_b."""
    w = 1 / e ** 2; X = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(X.T @ (X * w[:, None])); beta = cov @ X.T @ (w * y)
    return beta[1], sqrt(cov[1, 1]), beta, cov

# ---------------------------------------------------------------- hypothesis fits on slopes s(R)
def gls(X, y, e, prior=None):
    W = np.diag(1 / e ** 2); P = X.T @ W @ X; v = X.T @ W @ y
    if prior is not None:                 # list of (param index, sigma) Gaussian priors centred on 0
        for i, s in prior: P[i, i] += 1 / s ** 2
    cov = np.linalg.inv(P); b = cov @ v
    chi2 = float(((y - X @ b) / e) @ ((y - X @ b) / e))
    return b, cov, chi2

def hypothesis_tests(R, s, e, nd_scaling=-1.0, sig_delta_zs_at_RS=None):
    """R [A], slopes s [zs], errors e [zs]. Returns dict with test (B) (systematic-free shape) and (A) (absolute with priors)."""
    R = np.asarray(R, float); s = np.asarray(s, float) - nd_term(R, nd_scaling); e = np.asarray(e, float)
    out = {}
    # H0 fixed: s = R/c
    r0 = s - tau_Rc(R); out['H0_fixed_chi2'] = float((r0 / e) @ (r0 / e)); out['dof_H0_fixed'] = len(R)
    # (B) H0free: a*R ; H1: a*R + d1
    b0, c0, chi0 = gls(R[:, None], s, e); out['B_H0free_a_zs_per_A'] = float(b0[0]); out['B_H0free_chi2'] = chi0
    b1, c1, chi1 = gls(np.vstack([R, np.ones_like(R)]).T, s, e)
    out.update(B_a_zs_per_A=float(b1[0]), B_a_err=float(sqrt(c1[0, 0])), B_d1_zs=float(b1[1]), B_d1_err=float(sqrt(c1[1, 1])),
               B_Z=float(b1[1] / sqrt(c1[1, 1])), B_H1_chi2=chi1, B_dchi2=chi0 - chi1, B_a_over_inv_c=float(b1[0] / tau_Rc(1.0)))
    # (A) a fixed at 1/c, delta*R nuisance with prior (zs at R_S -> zs per A)
    if sig_delta_zs_at_RS is None: sig_delta_zs_at_RS = SIG_DELTA_DEFAULT
    y = s - tau_Rc(R)
    bA, cA, chiA = gls(np.vstack([np.ones_like(R), R]).T, y, e, prior=[(1, sig_delta_zs_at_RS / R_S)])
    out.update(A_d1_zs=float(bA[0]), A_d1_err=float(sqrt(cA[0, 0])), A_Z=float(bA[0] / sqrt(cA[0, 0])), A_delta_zs_at_RS=float(bA[1] * R_S),
               A_sig_delta_prior_zs_at_RS=sig_delta_zs_at_RS)
    return out

# ---------------------------------------------------------------- full pipeline on a 3-D histogram
def analyse_hist(Hc, ca, cb, Rc, KERc=None, Eph=EPH, sig_S=0.107, sigR_scale=None, return_rows=False):
    """Hc[iR, ib, ia] counts; ca cos(alpha) centres; cb cos(beta) centres; Rc R-bin centres (count-weighted means)."""
    if KERc is None: KERc = KE_EVA / np.asarray(Rc)
    Ee = Ee_of_KER(KERc, Eph=Eph)
    slopes, errs, rows = [], [], []
    for iR, R in enumerate(Rc):
        scale = (R_S / R) if sigR_scale is None else sigR_scale[iR]
        win = W_S * scale
        mus, ers = [], []
        for ib in range(len(cb)):
            mu, er, _ = fit_row(ca, Hc[iR, ib].astype(float), win, sig_S * scale)
            mus.append(mu); ers.append(er)
        mus, ers = np.array(mus), np.array(ers); ok = np.isfinite(mus) & np.isfinite(ers)
        k = zs_per_cos(R, Ee[iR])
        b, db, _, _ = wline(cb[ok], mus[ok] * k, ers[ok] * k)
        slopes.append(b); errs.append(db); rows.append(dict(R=float(R), mu=mus.tolist(), dmu=ers.tolist(), zs_per_cos=float(k)))
    res = dict(R_A=[float(r) for r in Rc], slope_zs=[float(v) for v in slopes], slope_err_zs=[float(v) for v in errs])
    if return_rows: res['rows'] = rows
    return res

# ---------------------------------------------------------------- input handling
def _read_table(path):
    with open(path) as f: lines = [l.strip() for l in f if l.strip()]
    lines = [l for l in lines if not l.lower().startswith('variables')]
    hdr = [h.strip().lower() for h in (lines[0].split(',') if ',' in lines[0] else lines[0].split())]
    data = np.array([[float(v) for v in (l.split(',') if ',' in l else l.split())] for l in lines[1:]])
    def col(*keys):
        for i, h in enumerate(hdr):
            if any(k in h for k in keys): return i
        return None
    return hdr, data, col

def load_joint(path, edges=None, nbins=4, dcos=0.0444444, nbeta=15):
    hdr, d, col = _read_table(path)
    ia, ib = col('cos(alpha)', 'cos_alpha', 'cosalpha'), col('cos(beta)', 'cos_beta', 'cosbeta')
    ik, ir, ic = col('ker', 'kinetic energy'), col('r (', 'r_a', 'r[a', 'internuclear', 'bond'), col('count')
    if ia is None or ib is None or (ik is None and ir is None): raise SystemExit(f'cannot identify columns in header {hdr}')
    R = KE_EVA / d[:, ik] if ik is not None else d[:, ir]
    w = d[:, ic] if ic is not None else np.ones(len(d))
    if ic is not None:   # binned: keep the provided cos grids
        ca = np.unique(d[:, ia]); cb = np.unique(d[:, ib])
    else:
        ea = np.arange(-1, 1 + 1e-9, dcos); ca = 0.5 * (ea[1:] + ea[:-1])
        eb = np.linspace(-1, 1, nbeta + 1); cb = 0.5 * (eb[1:] + eb[:-1])
    if edges is None:     # equal-count R bins
        order = np.argsort(R); cw = np.cumsum(w[order]) / w.sum()
        edges = [R.min() - 1e-9] + [R[order][np.searchsorted(cw, q)] for q in np.arange(1, nbins) / nbins] + [R.max() + 1e-9]
    edges = np.asarray(edges, float); nR = len(edges) - 1
    Hc = np.zeros((nR, len(cb), len(ca))); Rsum = np.zeros(nR); Ksum = np.zeros(nR); W = np.zeros(nR)
    jR = np.clip(np.searchsorted(edges, R, side='right') - 1, 0, nR - 1); inR = (R >= edges[0]) & (R < edges[-1])
    if ic is not None:
        ja = np.searchsorted(ca, d[:, ia]); jb = np.searchsorted(cb, d[:, ib])
    else:
        ja = np.clip(np.searchsorted(ea, d[:, ia], side='right') - 1, 0, len(ca) - 1); jb = np.clip(np.searchsorted(eb, d[:, ib], side='right') - 1, 0, nbeta - 1)
    for k in np.where(inR)[0]:
        Hc[jR[k], jb[k], ja[k]] += w[k]; Rsum[jR[k]] += w[k] * R[k]; Ksum[jR[k]] += w[k] * KE_EVA / R[k]; W[jR[k]] += w[k]
    return Hc, ca, cb, Rsum / W, Ksum / W, edges, W

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('path', nargs='?'); ap.add_argument('--out', default=None); ap.add_argument('--nbins', type=int, default=4)
    ap.add_argument('--edges', default=None); ap.add_argument('--nd-scaling', type=float, default=-1.0)
    ap.add_argument('--eph', type=float, default=EPH, help='photon energy of this data set (eV)')
    ap.add_argument('--selftest', action='store_true')
    a = ap.parse_args(argv)
    if a.selftest:
        here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
        import stage32_run as S32
        path = S32.write_selftest_file(); print('self-test file:', path); a.path = path
    if not a.path: ap.error('give the joint data file or --selftest')
    edges = [float(v) for v in a.edges.split(',')] if a.edges else None
    Hc, ca, cb, Rc, Kc, edges, W = load_joint(a.path, edges, a.nbins)
    res = analyse_hist(Hc, ca, cb, Rc, Kc, Eph=a.eph)
    res.update(edges_A=list(map(float, edges)), counts=list(map(float, W)))
    res['tests'] = {f'nd_scaling_{s:+.0f}': hypothesis_tests(Rc, res['slope_zs'], res['slope_err_zs'], s) for s in (-1.0, 0.0, 1.0)}
    t = res['tests'][f'nd_scaling_{a.nd_scaling:+.0f}']
    print('R bins (A):', np.round(Rc, 3)); print('slopes (zs):', np.round(res['slope_zs'], 1), '+-', np.round(res['slope_err_zs'], 1))
    print('travel time R/c:', np.round(tau_Rc(Rc), 1))
    print(f"(B) R-independent extra d1 = {t['B_d1_zs']:.1f} +- {t['B_d1_err']:.1f} zs  (Z = {t['B_Z']:.2f}); R-scaling coefficient = {t['B_a_over_inv_c']:.2f} x (1/c)")
    print(f"(A) with travel time fixed + priors: d1 = {t['A_d1_zs']:.1f} +- {t['A_d1_err']:.1f} zs (Z = {t['A_Z']:.2f})")
    print(f"chi2 H0 fixed (R/c+nd) = {t['H0_fixed_chi2']:.1f} / {t['dof_H0_fixed']}")
    if a.out:
        with open(a.out, 'w') as f: json.dump(res, f, indent=1)
        print('wrote', a.out)
    return res

if __name__ == '__main__':
    main()
