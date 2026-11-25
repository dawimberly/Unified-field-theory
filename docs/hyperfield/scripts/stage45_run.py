#!/usr/bin/env python3
"""Stage 45 - short-time gate (dead time t_g) from a single-photon decay trace that contains its own excitation pulse.

MODEL (ASSUMED). After excitation at t = 0 (centre of the 15 ps pi pulse):
 Gate-a 'clock gate': the hook only forms at t_g; the whole emission pattern is delayed:  S(t) = Th(t-t_g) e^{-(t-t_g)/tau} sin^2(a(t-t_g))
 Gate-b 'release gate': the state evolves from t = 0 but release is blocked until t_g:      S(t) = Th(t-t_g) e^{-(t-t_g)/tau} sin^2(a t)
 a = Delta_FSS/(2 hbar). Standard QM: t_g = 0, the sin^2 onset (fine-structure precession into the H-polarised collection mode,
 Ollivier et al. Eq. 4) is the only ps-scale non-exponential start in this trace. The generic short-time (Zeno) quadratic regime lasts about
 the inverse bandwidth of the coupled reservoir; for optical emission this is ASSUMED << 1 ps, so it cannot mimic a 5-50 ps gate.
Measured = A_L Gauss(t - t_L; s_L) [laser, same trace] + A_X (S conv Gauss jitter s_j) + B.  t_g is measured relative to t_L in the same
trace, so the axis zero and cable delays cancel.
Data: Ollivier et al., "Reproducibility of high-performance quantum dot single-photon sources", ACS Photonics 7, 1050 (2020), arXiv:1910.08863,
Fig. 2b linear panel, DIGITIZED from the embedded raster (917x632 px) -> data/stage45_ollivier2020_fig2b_digitized.csv.
"""
import json, os, csv
import numpy as np
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
S = {}
def src(key, ref, url, loc, quote, **vals): S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)
src('ollivier20', 'Ollivier, Maillette de Buy Wenniger, ... Senellart et al., ACS Photonics 7, 1050 (2020), arXiv:1910.08863', 'https://arxiv.org/pdf/1910.08863', 'Sec. III.A, Eq. (4), Fig. 2b/2d',
    '"For visualization purposes, the light of the laser was not perfectly extinguished for these data so that we detect both the laser pulse and the exciton emission dynamics to compare the relative timescales" ; '
    '"The maximum of the single-photon emission is delayed by approximately 200 ps from the laser excitation pulse" ; "accurately reproduced by eq. (4) using the following parameters: tau = 252+-3 ps and Delta_FSS = 8.58+-0.03 ueV" ; '
    '"consistent with a finite gaussian detector jitter time with a full width at half maximum (FWHM) of 53 ps, which also dominated the observed width of the 15 ps Gaussian excitation pulse" ; '
    'Eq. (4): "|<H|psi(t)>|^2 = e^{-t/tau} sin^2(t Delta_FSS/2hbar) sin^2(2 theta)" ; trion S11: "mono-exponential decay with a lifetime of 164.9 +- 0.9 ps ... the finite detector response of roughly 50 ps is consistent with the observed rise time"',
    tau=0.252, dtau=0.003, D=8.58, dD=0.03, jitter_fwhm=0.053, pulse=0.015)
src('digitization', 'DIGITIZED (this stage) from arXiv:1910.08863 page-4 raster image', 'data/stage45_ollivier2020_fig2b_digitized.csv', 'Fig. 2b top (linear) panel',
    'x: tick-label centres 0.0/0.5/1.0/1.5 ns at px 120.5/206/292.5/378.5 (172 px/ns, 5.8 ps/px, linear to 1 px); y: labels 0.0..1.0 at px 412..304.5 (107.5 px/unit); '
    'value = mean row of black (data) pixels per column; legend box columns restricted to y < 0.13', px_per_ns=172.0, x0=120.5)
src('nist', 'NIST CODATA 2022 hbar (fetched this session) and exact SI e', 'https://physics.nist.gov/cgi-bin/cuu/Value?hbar', 'value field', '"1.054 571 817... x 10 -34 J s"; e = 1.602 176 634e-19 C exact (SI 2019) -> hbar = 6.5821e-16 eV s (DERIVED)', hbar_eVs=1.054571817e-34 / 1.602176634e-19)
src('stage44', 'HYPERFIELD_RESULTS.md Stage 44', 'docs/hyperfield/HYPERFIELD_RESULTS.md', 'D4', 'Norman 60Co: t_g < 1.3e4 s = 5.3e-5 tau')
d = np.loadtxt(os.path.join(ROOT, 'data', 'stage45_ollivier2020_fig2b_digitized.csv'), delimiter=','); t, y = d[:, 0], d[:, 1]
HB = S['nist']['values']['hbar_eVs']; tt = np.arange(-0.6, 2.6, 0.0005); kx = np.arange(-0.15, 0.15001, 0.0005)
gauss = lambda x, s: np.exp(-x ** 2 / (2 * s * s))
def model(t, AL, tL, sL, AX, tau, D, tg, B, sj, v):
    a = D * 1e-6 / (2 * HB) * 1e-9
    ph = tt - tL - tg if v == 'a' else tt - tL
    Sx = np.where(tt >= tL + tg, np.exp(-(tt - tL - tg) / tau) * np.sin(a * ph) ** 2, 0)
    k = gauss(kx, sj); k /= k.sum()
    return AL * gauss(t - tL, sL) + AX * np.interp(t, tt, np.convolve(Sx, k, mode='same')) + B
P0 = [0.23, 0.003, 0.02, 2.4, 0.255, 8.5, 0.003, 0.0225]; LB = [0, -.05, .005, .1, .1, 7, -.05, .005]; UB = [1, .05, .06, 10, .5, 10, .05, .06]
def fit_at(tg, v):
    g = lambda t, AL, tL, sL, AX, tau, D, B, sj: model(t, AL, tL, sL, AX, tau, D, tg, B, sj, v)
    q, _ = curve_fit(g, t, y, p0=P0, bounds=(LB, UB), maxfev=40000)
    return float(np.sum((y - g(t, *q)) ** 2)), q
out = {}
grids = {'a': np.round(np.arange(-0.02, 0.02001, 0.0005), 5), 'b': np.round(np.arange(0, 0.07001, 0.0025), 5)}
for v in ('a', 'b'):
    prof = [(float(tg),) + (lambda r: (r[0], r[1].tolist()))(fit_at(tg, v)) for tg in grids[v]]
    out[v] = dict(profile=[(p[0], p[1]) for p in prof], best=min(prof, key=lambda p: p[1]))
# residual correlation -> inflation of sigma^2 (digitized points are not independent)
sb, qb = out['a']['best'][1], out['a']['best'][2]
r = y - model(t, *qb[:6], out['a']['best'][0], qb[6], qb[7], 'a'); r -= r.mean()
ac = [np.sum(r[:-k] * r[k:]) / np.sum(r * r) for k in range(1, 30)]; tau_int = 1 + 2 * sum(x for x in ac if x > 0)
s2 = np.var(r, ddof=9) * tau_int; out['sigma2_eff'] = s2; out['tau_int'] = tau_int
for v in ('a', 'b'):
    pr = np.array(out[v]['profile']); dchi = (pr[:, 1] - pr[:, 1].min()) / s2; out[v]['dchi2'] = dchi.tolist()
    ok = pr[dchi <= 3.84, 0]; out[v]['tg_95'] = [float(ok.min()), float(ok.max())]
    out[v]['dchi2_at_0'] = float(dchi[np.argmin(abs(pr[:, 0]))])
    out[v]['best_tg_ps'] = out[v]['best'][0] * 1e3; out[v]['best_params'] = dict(zip(['AL', 'tL', 'sL', 'AX', 'tau', 'D', 'B', 'sj'], out[v]['best'][2]))
SYS = 0.005   # ASSUMED: excitation instant within +-5 ps of the 15 ps pulse centre (pi-pulse build-up) + 1-px rendering offset between curves
out['a']['tg_95_with_sys_ps'] = [(out['a']['tg_95'][0] - SYS) * 1e3, (out['a']['tg_95'][1] + SYS) * 1e3]
out['b']['tg_95_with_sys_ps'] = [0.0, (out['b']['tg_95'][1] + SYS) * 1e3]
tau = out['a']['best_params']['tau']
out['ratios'] = dict(a=max(abs(x) for x in out['a']['tg_95_with_sys_ps']) * 1e-3 / tau, b=out['b']['tg_95_with_sys_ps'][1] * 1e-3 / tau)
# ---- figure ----
fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
ax[0].plot(t, y, 'k.', ms=3, label='Ollivier 2020 Fig. 2b (DIGITIZED)')
qa = out['a']['best'][2]; ax[0].plot(t, model(t, *qa[:6], out['a']['best'][0], qa[6], qa[7], 'a'), 'r-', lw=1, label=f"gate-a best t_g = {out['a']['best_tg_ps']:.1f} ps")
ax[0].plot(t, qa[0] * gauss(t - qa[1], qa[2]), color='orange', lw=1, label='laser (same trace) → t = 0')
fb, qbb = fit_at(0.05, 'b'); ax[0].plot(t, model(t, *qbb[:6], 0.05, qbb[6], qbb[7], 'b'), 'b--', lw=1, label='gate-b forced t_g = 50 ps')
ax[0].set_xlim(-0.15, 0.8); ax[0].set_xlabel('time after excitation (ns)'); ax[0].set_ylabel('norm. counts'); ax[0].legend(fontsize=7); ax[0].set_title('QD single-photon trace incl. its own laser pulse', fontsize=9)
for v, c in (('a', 'r'), ('b', 'b')):
    pr = np.array(out[v]['profile']); ax[1].plot(pr[:, 0] * 1e3, out[v]['dchi2'], c + '-o', ms=3, label=f'gate-{v}')
ax[1].axhline(3.84, color='k', ls=':', label=f'95% (Δχ² = 3.84, σ² inflated ×{tau_int:.1f})'); ax[1].set_ylim(0, 30)
ax[1].set_xlabel('t_g (ps)'); ax[1].set_ylabel('Δχ² (profile, other params free)'); ax[1].legend(fontsize=7); ax[1].set_title('profile likelihood of the dead time', fontsize=9)
names = ['Stage 44: ⁶⁰Co\n(created ≤3.5 h)', 'Stage 45 gate-b\n(QD, release blocked)', 'Stage 45 gate-a\n(QD, clock delayed)']
vals = [1.26e4, out['b']['tg_95_with_sys_ps'][1] * 1e-12, max(abs(x) for x in out['a']['tg_95_with_sys_ps']) * 1e-12]
rel = [5.3e-5, out['ratios']['b'], out['ratios']['a']]
ax[2].barh(names, vals, color=['tab:orange', 'tab:blue', 'tab:red']); ax[2].set_xscale('log'); ax[2].set_xlim(1e-12, 1e6)
for i, (v_, r_) in enumerate(zip(vals, rel)): ax[2].text(v_, i, f' {v_:.1e} s  (t_g/τ < {r_:.1e})', va='center', fontsize=8)
ax[2].set_xlabel('95% upper limit on dead time t_g (s)'); ax[2].set_title('gate (dead-time) limits', fontsize=9)
fig.suptitle('Stage 45 — short-time gate: emission starts with the excitation pulse; dead time < ~10 ps (clock) / < ~53 ps (release)', fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(ROOT, 'figures', 'stage45_short_time_gate.png'), dpi=110)
with open(os.path.join(ROOT, 'results', 'stage45_sources.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for kk, v in S.items(): wr.writerow([kk, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
for v in ('a', 'b'): out[v]['best'] = [out[v]['best'][0], out[v]['best'][1]]
json.dump(dict(stage=45, results=out, sources=S), open(os.path.join(ROOT, 'results', 'stage45.json'), 'w'), indent=1, default=float)
print(json.dumps({k: (v if k in ('sigma2_eff', 'tau_int', 'ratios') else {kk: vv for kk, vv in v.items() if kk not in ('profile', 'dchi2')}) for k, v in out.items()}, indent=1, default=float))
