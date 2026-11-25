#!/usr/bin/env python3
"""Stage 28 - LISA Pathfinder pressure-free white-noise floor -> bound on the CQ volume-law noise.

DATA ACTUALLY USED (stated plainly):
  * Per-run Brownian (frequency-independent) Delta-g PSD floors S_Brown,k with 1-sigma bars, 13 noise runs + the cold run,
    DIGITIZED from Fig. 3 of Armano et al. 2024, PRD 110, 042004, arXiv:2405.05207 (file BrownianEvolution_paper.png
    from https://arxiv.org/html/2405.05207v3). Pixel calibration to the printed log grid is in this script.
  * Run temperatures and mean times from that paper's Table 1 (CITED); t_v = launch + 62.4 d (CITED).
  * NOT used: raw LPF Legacy Archive time series. The ESA archive (lpf.esac.esa.int/lpfsa) is a GWT web app with
    CAS login and no documented scripted interface; the HEASARC mirror (heasarc.gsfc.nasa.gov/FTP/lpf/data) holds only
    NASA DRS telemetry (capacitive sensing, no LTP interferometer Delta-g). See results/stage28.json 'data_access'.
  * LPF has NO pressure gauge in the GRS: "pressure" is itself inferred from S_Brown (P = S/(2 kappa_gas)). So the
    zero-pressure extrapolation must use an outgassing MODEL for P(t, T) (ASSUMED, the paper's own form):
        S_k = a (t_v/(t_k - t_v))^gamma exp[-T_a (1/T_k - 1/T0)] + S_0 ,   T0 = 293 K
    S_0 (pressure-free white floor) contains any CQ noise. Temperature lever arm (runs 10, 11 at ~284 K) breaks the
    a <-> S_0 degeneracy.
"""
import json, os, urllib.request
import numpy as np
from math import pi, sqrt, log10
from scipy.optimize import least_squares
from scipy import ndimage
from PIL import Image
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
IMG = os.path.join(ROOT, 'data/armano2024_fig3_BrownianEvolution.png')
URL = 'https://arxiv.org/html/2405.05207v3/BrownianEvolution_paper.png'
if not os.path.exists(IMG):
    os.makedirs(os.path.dirname(IMG), exist_ok=True); urllib.request.urlretrieve(URL, IMG)
HBAR = 6.62607015e-34 / (2 * pi); G = 6.6743e-11; GH = G * HBAR

# ---------------------------------------------------------------- digitization (DERIVED from the published figure)
im = np.array(Image.open(IMG).convert('RGB')).astype(int)
# axis calibration from the major grid lines found in the image (pixel centres): x: 50 d -> 664.5, 100 d -> 1482.5,
# 200 d -> 2301.5, 500 d -> 3384.5 ; y: 10 -> 747.0, 2 -> 1435.5 (fm^2 s^-4 Hz^-1)
xp = np.array([664.5, 1482.5, 2301.5, 2780.5, 3120.5, 3384.5]); xv = np.log10([50, 100, 200, 300, 400, 500])
yp = np.array([58.5, 153.5, 276.5, 450.5, 747.0, 1043.5, 1262.5, 1435.5]); yv = np.log10([50, 40, 30, 20, 10, 5, 3, 2])
cx = np.polyfit(xp, xv, 1); cy = np.polyfit(yp, yv, 1)
calib = dict(x_fit=cx.tolist(), y_fit=cy.tolist(), x_resid_dec=float(np.max(np.abs(np.polyval(cx, xp) - xv))),
             y_resid_dec=float(np.max(np.abs(np.polyval(cy, yp) - yv))))
X = lambda p: 10 ** np.polyval(cx, p); Y = lambda p: 10 ** np.polyval(cy, p)
def blobs(mask):
    lab, n = ndimage.label(mask); out = []
    for s in ndimage.find_objects(lab):
        ys, xs = s
        if (xs.stop - xs.start) * (ys.stop - ys.start) > 300 and xs.start > 500: out.append((xs.start, xs.stop, ys.start, ys.stop))
    return sorted(out)
MARK = 25.0   # marker diameter in px (width of every blob)
def read(b):
    x0, x1, y0, y1 = b; xc = (x0 + x1 - 1) / 2
    top, bot = Y(y0), Y(y1 - 1)
    if (y1 - y0) <= MARK + 4:        # error bar hidden inside marker: centre = blob centre, sigma <= half-marker
        yc = (y0 + y1 - 1) / 2; s = Y(yc); sig = (Y(yc - MARK / 2) - Y(yc + MARK / 2)) / 2; hidden = True
    else:
        s = (top + bot) / 2; sig = (top - bot) / 2; hidden = False
    return dict(t_minus_tv=float(X(xc)), S=float(s), sigma=float(sig), errbar_hidden=hidden, px=(int(x0), int(x1), int(y0), int(y1)))
blue = [read(b) for b in blobs(np.abs(im - np.array([1, 26, 128])).sum(2) < 60)]
cold = [read(b) for b in blobs((im[:, :, 2] > 200) & (im[:, :, 0] < 100) & (im[:, :, 1] > 120))]
# Table 1 (CITED): mean time from launch (d), temperature (K)
T1 = [(112, 295.37, .04), (127, 295.30, .03), (166, 294.97, .06), (202, 294.93, .01), (229, 296.62, .05), (237, 296.50, .02),
      (302, 296.50, .02), (354, 296.86, .04), (398, 295.38, .04), (446, 284.72, .03), (535, 284.2, .1), (547, 295.78, .01), (558, 295.91, .03)]
TV = 62.4
assert len(blue) == 13, len(blue)
runs = []
for k, (b, (tm, T, dT)) in enumerate(zip(blue, T1)):
    runs.append(dict(run=k + 1, t_minus_tv_table=tm - TV, t_minus_tv_digitized=b['t_minus_tv'], T=T, S=b['S'], sigma=b['sigma'], errbar_hidden=b['errbar_hidden']))
colddat = dict(T=274.8, T_err=0.5, t_minus_tv=cold[0]['t_minus_tv'], S=cold[0]['S'], sigma=cold[0]['sigma'])
# digitization accuracy check: x positions vs Table-1 mean times
dx = np.array([r['t_minus_tv_digitized'] / r['t_minus_tv_table'] - 1 for r in runs])
DIG_REL = 0.02   # ASSUMED 2% digitization uncertainty on S (grid residuals + marker centring), added in quadrature

t = np.array([r['t_minus_tv_table'] for r in runs]); T = np.array([r['T'] for r in runs]); S = np.array([r['S'] for r in runs])
sg = np.sqrt(np.array([r['sigma'] for r in runs]) ** 2 + (DIG_REL * S) ** 2)
T0 = 293.0
def model(p, t, T):
    a, gam, Ta, S0 = p
    return a * (TV / t) ** gam * np.exp(-Ta * (1 / T - 1 / T0)) + S0
def fit(mask, S0_fixed=None, extra=None):
    tt, TT, SS, ss = t[mask], T[mask], S[mask], sg[mask]
    if extra is not None:
        tt = np.append(tt, extra[0]); TT = np.append(TT, extra[1]); SS = np.append(SS, extra[2]); ss = np.append(ss, extra[3])
    def res(q):
        p = list(q) + ([S0_fixed] if S0_fixed is not None else [])
        return (model(p, tt, TT) - SS) / ss
    x0 = [27.0, 0.8, 7000.0] + ([0.5] if S0_fixed is None else [])
    lb = [0, 0, 0] + ([0] if S0_fixed is None else []); ub = [1e4, 5, 5e4] + ([50] if S0_fixed is None else [])
    r = least_squares(res, x0, bounds=(lb, ub))
    return r, float((r.fun ** 2).sum()), len(SS) - len(x0)
def profile(mask, extra=None, grid=None):
    grid = np.linspace(0, 3.0, 301) if grid is None else grid
    chi = np.array([fit(mask, S0_fixed=s, extra=extra)[1] for s in grid])
    i = int(np.argmin(chi)); cmin = chi[i]
    up95 = float(np.interp(cmin + 2.71, chi[i:], grid[i:])) if chi[-1] > cmin + 2.71 else float('nan')   # one-sided 95% (Delta chi2 = 2.71)
    lo68 = float(np.interp(cmin + 1, chi[:i + 1][::-1], grid[:i + 1][::-1])) if i > 0 and chi[0] > cmin + 1 else 0.0
    hi68 = float(np.interp(cmin + 1, chi[i:], grid[i:]))
    return dict(grid=grid.tolist(), chi2=chi.tolist(), S0_best=float(grid[i]), chi2_min=float(cmin), S0_68=(lo68, hi68), S0_95_upper=up95,
                delta_chi2_at_S0_zero=float(chi[0] - cmin), significance_sigma_vs_zero=float(sqrt(max(chi[0] - cmin, 0))))
allm = np.ones(13, bool)
variants = {}
ex_cold = (colddat['t_minus_tv'], colddat['T'], colddat['S'], sqrt(colddat['sigma'] ** 2 + (DIG_REL * colddat['S']) ** 2))
for name, mask, extra in (('A: all 13 runs (baseline)', allm, None),
                          ('B: drop runs 1-2 (early fast-decay phase)', np.arange(13) >= 2, None),
                          ('C: 13 runs + cold run (T = 274.8 K)', allm, ex_cold),
                          ('D: only runs 8-13 (late, includes both cold runs)', np.arange(13) >= 7, None)):
    r, chi, dof = fit(mask, extra=extra)
    pr = profile(mask, extra)
    J = r.jac; cov = np.linalg.pinv(J.T @ J)
    variants[name] = dict(a=float(r.x[0]), gamma=float(r.x[1]), T_a=float(r.x[2]), S0_fit=float(r.x[3]), S0_sigma_lin=float(sqrt(cov[3, 3])),
                          chi2=chi, dof=dof, **{k: v for k, v in pr.items() if k not in ('grid', 'chi2')}, profile=dict(grid=pr['grid'], chi2=pr['chi2']))
# validation: S0 fixed to 0 should reproduce the paper's own fit (a = 27.0 +- 0.7, gamma = 0.80 +- 0.02, T_a = 7.0 +- 0.2 kK)
r0, chi0, dof0 = fit(allm, S0_fixed=0.0)
paper_check = dict(a=float(r0.x[0]), gamma=float(r0.x[1]), T_a=float(r0.x[2]), chi2=chi0, dof=dof0)
# model-independent floor: S_0 <= lowest measured white level (cold run) at 95%
cold_up = colddat['S'] + 1.645 * sqrt(colddat['sigma'] ** 2 + (DIG_REL * colddat['S']) ** 2)
run10_up = runs[9]['S'] + 1.645 * sg[9]

# ---------------------------------------------------------------- convert to CQ bound (Stage 26 conventions)
s26 = json.load(open(os.path.join(ROOT, 'results/stage26.json')))
M_LPF = 1.928; V = M_LPF / 19800.0
A_pl = (4 * pi / 3) * GH / V                      # plateau DP (kappa = 1) Gamma, two-sided per axis
GAM26 = (1.74e-15) ** 2 / 4                        # Stage 26: Gamma = S_dg/4 (differential, one-sided -> two-sided)
FEMTO2 = 1e-30                                     # fm^2 s^-4 Hz^-1 -> m^2 s^-4 Hz^-1
def bound(S0_up):
    Gam = S0_up * FEMTO2 / 4; u = Gam / A_pl
    kp = (u + sqrt(max(u * u - 1, 0))) / 2; km = 1 / (4 * kp)
    return dict(S0_upper_fm2=S0_up, Gamma_bound=Gam, GammaV_bound=Gam * V, log10u=log10(u), gain_vs_stage26=log10(GAM26 / Gam),
                kappa_max=kp, kappa_min=km)
B = {k: bound(v['S0_95_upper']) for k, v in variants.items() if np.isfinite(v['S0_95_upper'])}
B['model-independent: cold run S <= ' + '%.2f' % cold_up] = bound(cold_up)
B['model-independent: run 10 (Stage 26 input) S <= ' + '%.2f' % run10_up] = bound(run10_up)
# kappa vs l curve: LPF u(l) from Stage 26 scaled by Gamma_new/Gamma_26
E26 = np.array(s26['scan']['ell']); u26 = np.array(s26['scan']['u']['LISA Pathfinder 2018 (Armano)'], float)
base = B['A: all 13 runs (baseline)']
u_new = u26 * base['Gamma_bound'] / GAM26
kp = (u_new + np.sqrt(np.clip(u_new ** 2 - 1, 0, None))) / 2
kcurve = dict(ell=E26.tolist(), u26=u26.tolist(), u_new=u_new.tolist(), kappa_max=kp.tolist(), kappa_min=(1 / (4 * kp)).tolist(),
              excluded_new=[float(e) for e in E26[u_new < 1]][:1] + [float(e) for e in E26[u_new < 1]][-1:])
out = dict(stage=28,
           data_access=dict(
               esa_archive='https://lpf.esac.esa.int/lpfsa/ : GWT single-page app, auth_method=cas (js/config.js); anonymous role exists but no documented scripted/TAP interface (tap_adql_interface_active=false; /lpf-tap/ 404). Raw Delta-g time series not obtained.',
               heasarc='https://heasarc.gsfc.nasa.gov/FTP/lpf/data/ : public, downloads work, but it is the NASA DRS (ST7) telemetry (1 Hz/10 Hz science + housekeeping, mostly mode Off during LTP runs); no LTP interferometric Delta-g -> not usable for a fm s^-2 floor.',
               used='published per-run S_Brown (13 runs + cold run) digitized from Armano et al. 2024 Fig. 3 (arXiv:2405.05207v3) + Table 1 temperatures'),
           calibration=calib, digitization_x_check_rel=dict(max=float(np.max(np.abs(dx))), rms=float(np.sqrt(np.mean(dx ** 2)))),
           digitization_rel_uncertainty_assumed=DIG_REL, runs=runs, cold_run=colddat,
           paper_fit_check=dict(published=dict(a=27.0, gamma=0.80, T_a=7000.0), ours_with_S0_zero=paper_check, note='paper fits with S_0 = 0 (no floor)'),
           variants=variants, bounds=B, A_plateau_DP=A_pl, Gamma_stage26=GAM26, V=V, kappa_curve=kcurve,
           sources=dict(armano2024='https://arxiv.org/abs/2405.05207', armano2018='PRL 120, 061101 (2018)', figure=URL))
json.dump(out, open(os.path.join(ROOT, 'results/stage28.json'), 'w'), indent=1, default=float)
print(json.dumps(dict(paper_check=paper_check, calib=calib, dx=out['digitization_x_check_rel'], runs=[(r['run'], round(r['t_minus_tv_table'], 1), round(r['t_minus_tv_digitized'], 1), r['T'], round(r['S'], 3), round(r['sigma'], 3), r['errbar_hidden']) for r in runs],
                      cold=colddat, variants={k: {kk: vv for kk, vv in v.items() if kk != 'profile'} for k, v in variants.items()}, bounds=B,
                      excl=kcurve['excluded_new']), indent=1, default=str))
