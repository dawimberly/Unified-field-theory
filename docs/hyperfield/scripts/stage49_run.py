"""Stage 49: eta at x ~ 1-5 from transmon primary thermometry of propagating microwaves.
Scigliuzzo et al., PRX 10, 041054 (2020), Fig. 4(d): steady-state radiation temperature T_r (from thermal occupation n_r, converted
with Bose-Einstein, paper Eq. E3) vs independently measured mixing-chamber temperature T_MX (calibrated RuO2).
Stage 31 definition: n = 1/((1+eta) e^x - 1), x = hbar w / k T.  If the paper converts n_r -> T_r with plain BE:
   T_r = hbar w / (k ln(1 + 1/n_r)).
Model fitted (ASSUMED): n_r = 1/((1+eta) e^{x(T_MX)} - 1) + n0  (n0 = residual occupation from hotter stages, T-independent).
"""
import json, csv, numpy as np
import pymupdf as fitz
from scipy.optimize import least_squares
OUT = '/workspace/uft-stage33/docs/hyperfield/'; SRC = '/workspace/stage49_src/scig.pdf'
hbar = 1.054571817e-34; kB = 1.380649e-23
f_ge = 5.332e9; TQ = 2*np.pi*hbar*f_ge/kB          # K
# ---- DIGITIZED from the vector graphics of Fig. 4(d) (page 6): purple filled markers inside the panel frame ----
pg = fitz.open(SRC)[5]
P = (0.4941176474094391, 0.24705882370471954, 0.5960784554481506)
pts = []
for d in pg.get_drawings():
    r = d['rect']
    if r.x0 > 444 and r.x1 < 530 and r.y0 > 290 and r.y1 < 380 and d.get('fill') and abs(d['fill'][0]-P[0]) < 1e-3 and 2 < r.width < 3:
        pts.append(((r.x0+r.x1)/2, (r.y0+r.y1)/2))
pts = sorted(pts)
# axis calibration from tick marks / labels: x ticks 0,50,100,150,200 mK at 445.94, 462.51, 479.07, 495.63, 512.19 pt ; y labels 50..200 mK at 356.96..307.27, frame bottom 373.52 = 0 mK
sx = (512.19-445.94)/200; sy = (356.96-307.27)/150
TMX = np.array([(x-445.94)/sx for x, y in pts])*1e-3; Tr = np.array([(373.52-y)/sy for x, y in pts])*1e-3
np.savetxt(OUT+'data/stage49_scigliuzzo2020_fig4d_digitized.csv', np.c_[TMX*1e3, Tr*1e3], delimiter=',', header='T_MX_mK,T_r_mK (DIGITIZED from PDF vector coordinates of Fig. 4d)', fmt='%.3f', comments='')
nBE = lambda T, eta=0.0: 1/((1+eta)*np.exp(TQ/T) - 1)
n_r = 1/np.expm1(TQ/Tr)
def model(p, T): eta, n0 = p[0], p[1]; return nBE(T, eta) + n0
def Tr_model(p, T):                       # p = eta, n0, dT  (dT: offset of the RuO2 reading, nuisance)
    return TQ/np.log1p(1/model(p, T + (p[2] if len(p) > 2 else 0.0)))
R = {}
fits = {}
for name, p0, lo, hi in [('eta_n0', [0.0, 1e-3], [-0.5, 0], [5, 0.05]), ('eta_n0_dT', [0.0, 1e-3, 0.0], [-0.5, 0, -0.02], [5, 0.05, 0.02]),
                         ('n0_dT_eta0', [1e-3, 0.0], [0, -0.02], [0.05, 0.02])]:
    if name == 'n0_dT_eta0':
        f = lambda q: Tr_model([0.0, q[0], q[1]], TMX) - Tr
    else:
        f = lambda q: Tr_model(q, TMX) - Tr
    fit = least_squares(f, p0, bounds=(lo, hi), x_scale='jac')
    dof = len(TMX)-len(p0); s2 = np.sum(fit.fun**2)/dof; cov = np.linalg.inv(fit.jac.T@fit.jac)*s2
    fits[name] = fit
    R[name] = dict(params=list(fit.x), se=list(np.sqrt(np.diag(cov))), rms_mK=np.sqrt(s2)*1e3, chi2_unit=np.sum(fit.fun**2), dof=dof)
R['all'] = dict(eta=R['eta_n0_dT']['params'][0], eta_se=R['eta_n0_dT']['se'][0], n0=R['eta_n0_dT']['params'][1], dT_mK=R['eta_n0_dT']['params'][2]*1e3, dT_se_mK=R['eta_n0_dT']['se'][2]*1e3)
fr = R['all']
# systematics (DERIVED mapping; inputs CITED/ASSUMED):
#  occupation-scale error s: ln n shifts by s; d ln n / d eta = -(1+n) -> eta_equiv = s/(1+n) at the hottest points (n ~ 0.4)
#  T_MX calibration error eps: d ln n = x (1+n) eps -> eta_equiv = x eps
hot = TMX > 0.15; nh = np.mean(nBE(TMX[hot])); xh = np.mean(TQ/TMX[hot])
s_scale = 0.01      # CITED: Fig. 3 slope agrees with prediction "within 1%"
eps_T = 0.01        # ASSUMED: RuO2 calibration 1% (paper gives no number)
sys_scale = s_scale/(1+nh); sys_T = xh*eps_T
tot = np.sqrt(fr['eta_se']**2 + sys_scale**2 + sys_T**2)
R['systematics'] = dict(n_hot=nh, x_hot=xh, scale_1pct=sys_scale, T_cal_1pct_ASSUMED=sys_T)
R['eta95'] = abs(fr['eta']) + 1.96*tot
R['eta95_stat_only'] = abs(fr['eta']) + 1.96*fr['eta_se']
R['x_range_used'] = [float(TQ/TMX.max()), float(TQ/0.060)]
R['firas'] = 9e-5; R['vs_firas'] = R['eta95']/9e-5
R['TQ_K'] = TQ
R['points'] = dict(TMX_mK=list(np.round(TMX*1e3, 2)), Tr_mK=list(np.round(Tr*1e3, 2)))
R['note'] = 'Bound is set by the overall agreement T_r ~ T_MX at the 1-2% level around 150-215 mK (x ~ 1.2-1.7). eta and an occupation-scale error are partly degenerate.'
json.dump(R, open(OUT+'results/stage49.json', 'w'), indent=1, default=float)
S = [dict(key='scigliuzzo2020', ref='Scigliuzzo, Bengtsson, Besse, Wallraff, Delsing, Gasparinetti, PRX 10, 041054 (2020), arXiv:2003.13522; open-access full text at Chalmers',
          url='https://research.chalmers.se/publication/521793/file/521793_Fulltext.pdf', location='Sec. II A, II D, Fig. 4(d), Appendix E Eq. (E3)',
          quote='"We operate the transmon at its highest fundamental frequency wge = 2pi x 5.332 GHz" ; "we extract the thermal occupation, convert it into radiation temperature Tr, and compare it to the temperature of the base (mixing chamber) plate, TMX, measured independently by a calibrated Ruthenium oxide thermometer" ; '
                '"The two values essentially agree over more than 2 orders of magnitude in thermal occupation number, from nr = 0.001 to 0.4 (corresponding to a temperature range of 35 to 210 mK)" ; "(d) Steady-state values of Tr versus TMX (circles); the solid line indicates Tr = TMX." ; "n(w, Tr) = 1/(exp(hbar w/kB Tr) - 1)"',
          values=json.dumps(dict(f_ge=5.332e9))),
     dict(key='scigliuzzo2020_fig3', ref='same', url='https://research.chalmers.se/publication/521793/file/521793_Fulltext.pdf', location='Sec. II C',
          quote='"a best fit of the linear slope of the low-power trace for dnr << 1 agrees with the predicted one within 1%. This agreement confirms the validity of the calibration scheme"', values=json.dumps(dict(scale=0.01))),
     dict(key='scigliuzzo2020_fig4d_digitized', ref='same, Fig. 4(d)', url='https://research.chalmers.se/publication/521793/file/521793_Fulltext.pdf', location='page 6, Fig. 4(d) vector paths',
          quote='DIGITIZED: 14 filled purple markers (2.7 pt) read from PDF drawing coordinates with PyMuPDF; axes calibrated from tick marks (0-200 mK, 0.33125 pt/mK both axes) and tick labels', values='{}'),
     dict(key='stage31_def', ref='This repo, Stage 31 (7)', url='docs/hyperfield/HYPERFIELD_RESULTS.md', location='Stage 31 (7)',
          quote='"If E/A = e^(-x)/(1+eta), then n = 1/((1+eta)e^x - 1)" ; "|eta| < 9.0x10^-5"', values='{"eta_firas":9e-5}')]
with open(OUT+'results/stage49_sources.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(S[0].keys())); w.writeheader(); [w.writerow(s) for s in S]
print(json.dumps({k: v for k, v in R.items() if k != 'points'}, indent=1, default=float))
print(R['points'])
# ---- figure ----
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))
a = ax[0]
TT = np.linspace(0.005, 0.23, 300)
a.plot(TMX*1e3, Tr*1e3, 'o', color='purple', label='data, Fig. 4(d) (DIGITIZED, vector)')
a.plot(TT*1e3, TT*1e3, '-', color='orange', label='T_r = T_MX (paper line)')
pf = fits['eta_n0_dT'].x; Tfit = Tr_model(pf, TT)
a.plot(TT*1e3, Tfit*1e3, 'k--', label='fit: eta = %.3f, n0 = %.1e, dT = %.1f mK' % (pf[0], pf[1], pf[2]*1e3))
for e, c in [(0.05, 'tab:red'), (-0.05, 'tab:blue')]:
    a.plot(TT*1e3, Tr_model([e, pf[1], pf[2]], TT)*1e3, ':', color=c, label='eta = %+.2f' % e)
a.set_xlabel('T_MX (mK), RuO2 thermometer'); a.set_ylabel('T_r (mK), transmon primary thermometer'); a.legend(fontsize=7); a.grid(alpha=0.3)
a.set_title('Scigliuzzo 2020: radiation vs plate temperature')
b = ax[1]
xx = TQ/TMX
b.semilogx(xx, Tr/TMX - 1, 'o', color='purple', label='data: T_r/T_MX - 1')
xs = TQ/TT
b.semilogx(xs, Tfit/TT - 1, 'k--', label='fit (eta + residual n0 + RuO2 offset)')
for e, c in [(9e-5, 'k'), (0.02, 'tab:red'), (-0.02, 'tab:blue')]:
    b.semilogx(xs, Tr_model([e, pf[1], pf[2]], TT)/TT - 1, ':' if e != 9e-5 else '-', color=c, lw=1.2, label='eta = %+.0e%s' % (e, ' (FIRAS limit)' if e == 9e-5 else ''))
b.set_ylim(-0.15, 0.3); b.set_xlim(1, 30); b.axvspan(1.05, 11.1, color='grey', alpha=0.12, label='FIRAS x range')
b.set_xlabel('x = hbar w / k T_MX'); b.set_ylabel('fractional deviation'); b.legend(fontsize=7, loc='upper left'); b.grid(alpha=0.3, which='both')
b.set_title('deviation vs x (low-x end carries the eta information)')
c = ax[2]
labs = ['transmon\n(this stage)', 'CO ladder\n(Stage 48)', 'FIRAS\n(Stage 31)']
vals = [R['eta95'], 2.3e-5, 9e-5]; xr = ['x 1.2-4.3', 'x 0.02-0.5', 'x 1.05-11']
c.bar(range(3), vals, color=['purple', 'tab:red', 'k']); c.set_yscale('log'); c.set_xticks(range(3)); c.set_xticklabels(labs, fontsize=8)
for i, v in enumerate(vals): c.text(i, v*1.2, '%.1e\n%s' % (v, xr[i]), ha='center', fontsize=8)
c.set_ylim(1e-5, 1); c.set_ylabel('|eta| 95% upper limit'); c.set_title('eta bounds'); c.grid(alpha=0.3, axis='y')
fig.suptitle('Stage 49 - eta at x ~ 1-4 from transmon microwave thermometry: |eta| < %.1e (95%%), %.0f x weaker than FIRAS; FIRAS remains best at x >= 1' % (R['eta95'], R['vs_firas']), fontsize=11)
fig.tight_layout(); fig.savefig(OUT+'figures/stage49_transmon_eta.png', dpi=110)
