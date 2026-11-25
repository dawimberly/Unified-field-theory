"""Stage 48: path dependence of the hysteresis eta, using Li et al. 2025 (Sci. Adv. 11, eadz6560) Table S1.

Measured: line-intensity ratio I1/I2 of CO (3-0) lines -> Einstein-A ratio via Boltzmann factor at the sensor T (paper's Eq. S5).
If populations deviate from Boltzmann by N_J = N_J^B exp(-g(J)), the paper's inferred A-ratio carries exp(-(g(J1)-g(J2))).
So r = ln(A_exp / A_UCL2025) = -(g(J1)-g(J2)) + theory error.
  * R(J)/P(J) share the lower level J -> g cancels exactly -> control: r = theory (+exp) error only.
  * R(J1)/R(J2), different lower levels -> sensitive to g.
Readings of the Stage 31 eta (per quantum: emission/absorption = e^{-x}/(1+eta)):
  L  ladder: eta acts on each Delta J = 1 rung (rotational dipole selection rule) -> g(J) = J ln(1+eta)
  G  ground-only: all J >= 1 share one factor -> g(J) = ln(1+eta) for J >= 1, 0 for J = 0
  S  pair-contextual (Stage 47 reading): every different-level ratio carries one (1+eta) factor; not a population vector
     -> predicts a closure defect: (R1/R3)(R3/R5)(R5/R10) / (R1/R10) = (1+eta)^{-2}
"""
import json, csv, numpy as np
OUT = '/workspace/uft-stage33/docs/hyperfield/'
# Table S1 (CITED): pair, P(Pa), T(K), I_exp, sig_I, A_exp, A_UCL2025
T1 = [
 ('R0/R1',15.88,298.429,0.503012,20e-6,0.822015,0.822050), ('R1/R3',8.14,298.430,0.534955,15e-6,0.875766,0.875781),
 ('R3/R5',6.95,298.426,0.767994,18e-6,0.936956,0.937037), ('R5/R10',7.12,298.431,1.073886,17e-6,0.901332,0.901365),
 ('R1/R10',7.18,298.431,0.441206,10e-6,0.739608,0.739695), ('R10/R13',7.18,298.430,1.472357,19e-6,0.950862,0.950942),
 ('R13/R15',7.18,298.429,1.458483,58e-6,0.968856,0.968894), ('R15/R17',14.33,298.426,1.594931,48e-6,0.969778,0.969828),
 ('R17/R19',26.12,298.430,1.738204,35e-6,0.970451,0.970496), ('R19/R21',50.62,298.430,1.889597,38e-6,0.970908,0.970993),
 ('R21/R22',50.62,298.432,1.417658,30e-6,0.985592,0.985542), ('R22/R23',98.36,298.434,1.446427,52e-6,0.985549,0.985628),
 ('R23/R24',98.36,298.437,1.475842,22e-6,0.985708,0.985706), ('R24/R25',206.0,298.434,1.505442,45e-6,0.985731,0.985776),
 ('R25/R26',422.2,298.429,1.535529,49e-6,0.985812,0.985840), ('R26/R27',423.9,298.428,1.566004,63e-6,0.985883,0.985899),
 ('R1/P1',12.88,298.428,2.076074,125e-6,0.416684,0.416691),
 ('R5/P5',3.29,298.430,1.376124,28e-6,0.965122,0.965077), ('R5/P5',7.58,298.429,1.376132,41e-6,0.965128,0.965077), ('R5/P5',11.21,298.928,1.376095,43e-6,0.965102,0.965077),
 ('R10/P10',6.42,298.428,1.428264,49e-6,1.209416,1.209413), ('R15/P15',12.88,298.428,1.567518,47e-6,1.428752,1.428804),
 ('R17/P17',12.88,298.427,1.634831,65e-6,1.519477,1.519466),
 ('R20/P20',49.12,298.428,1.745511,70e-6,1.661507,1.661486),
 ('R21/P21',49.12,298.427,1.784819,93e-6,1.710702,1.710718), ('R21/P21',77.07,298.430,1.784829,36e-6,1.710712,1.710718), ('R21/P21',97.11,298.428,1.784802,45e-6,1.710686,1.710718),
 ('R23/P23',121.8,298.431,1.867087,56e-6,1.812346,1.812367), ('R23/P23',125.7,298.427,1.867099,75e-6,1.812358,1.812367),
 ('R23/P23',139.6,298.427,1.867157,50e-6,1.812414,1.812367), ('R23/P23',153.5,298.429,1.867182,56e-6,1.812439,1.812367),
 ('R26/P26',225.8,298.427,1.999426,130e-6,1.973568,1.973515), ('R26/P26',334.6,298.428,1.999467,80e-6,1.973609,1.973515),
 ('R26/P26',417.6,298.431,1.999478,99e-6,1.973620,1.973515), ('R26/P26',440.7,298.428,1.999611,64e-6,1.973751,1.973515),
 ('R26/P26',529.8,298.429,1.999612,52e-6,1.973752,1.973515), ('R26/P26',543.5,298.428,1.999598,56e-6,1.973738,1.973515),
 ('R26/P26',647.6,298.427,1.999587,56e-6,1.973727,1.973515)]
# lower-state energies: E = B J(J+1) - D J^2(J+1)^2, B, D DERIVED from the paper's cited E(R1 lower)=3.8450, E(R17 lower)=587.7209 cm-1; check vs E(J=20)=806.3828
M = np.array([[2, -4], [306, -306**2]]); B, D = np.linalg.solve(M, [3.8450, 587.7209])
E = lambda J: B*J*(J+1) - D*(J*(J+1))**2
c2 = 1.438776877
def parse(p):
    a, b = p.split('/'); return a[0], int(a[1:]), b[0], int(b[1:])
rows = []
for p, P, T, I, sI, Ae, Au in T1:
    t1, J1, t2, J2 = parse(p)
    same = (J1 == J2)
    sT = 0.0 if same else c2*abs(E(J1)-E(J2))/T**2 * 5e-3   # 5 mK probe calibration (CITED), relative effect on ratio
    rows.append(dict(pair=p, P=P, J1=J1, J2=J2, same=same, r=np.log(Ae/Au), s_stat=sI/I, s_T=sT))
ctrl = [x for x in rows if x['same']]; diff = [x for x in rows if not x['same']]
def wstats(v, s):
    w = 1/np.asarray(s)**2; m = np.sum(w*v)/np.sum(w); return m, 1/np.sqrt(np.sum(w)), np.sum(w*(v-m)**2)/(len(v)-1)
rc = np.array([x['r'] for x in ctrl]); sc = np.array([x['s_stat'] for x in ctrl])
mc, smc, chi2c = wstats(rc, sc)
# theory+systematic floor from the control: extra scatter needed for chi2/dof = 1
def floor(r, s):
    lo, hi = 0, 1e-3
    for _ in range(60):
        mid = (lo+hi)/2; _, _, c = wstats(r, np.sqrt(s**2+mid**2)); lo, hi = (mid, hi) if c > 1 else (lo, mid)
    return hi
s_th = floor(rc, sc)
rd = np.array([x['r'] for x in diff]); sd = np.sqrt(np.array([x['s_stat'] for x in diff])**2 + np.array([x['s_T'] for x in diff])**2)
_, _, chi2d_raw = wstats(rd, sd)
s_th_d = floor(rd, sd)
sig_d = np.sqrt(sd**2 + s_th**2)                      # theory floor taken from the eta-free control (ASSUMED transferable)
J1 = np.array([x['J1'] for x in diff]); J2 = np.array([x['J2'] for x in diff])
# --- Model L: r = -(J1-J2) eta ; Model G: r = -([J1>=1]-[J2>=1]) eta ---
def fit(X):
    w = 1/sig_d**2; eta = np.sum(w*X*rd)/np.sum(w*X*X); se = 1/np.sqrt(np.sum(w*X*X))
    chi2 = np.sum(w*(rd-eta*X)**2); return eta, se, chi2
XL = -(J1-J2).astype(float); XG = -((J1 >= 1).astype(float) - (J2 >= 1).astype(float))
eL, sL, cL = fit(XL); eG, sG, cG = fit(XG)
chi0 = np.sum(rd**2/sig_d**2)
# --- closure (model S) ---
d = {x['pair']: x for x in rows}
Aexp = {p: Ae for p, P, T, I, sI, Ae, Au in T1}; sIr = {p: sI/I for p, P, T, I, sI, Ae, Au in T1}
prod = Aexp['R1/R3']*Aexp['R3/R5']*Aexp['R5/R10']; direct = Aexp['R1/R10']
s_prod = np.sqrt(sum(sIr[k]**2 for k in ['R1/R3', 'R3/R5', 'R5/R10', 'R1/R10']))
lnclos = np.log(prod/direct)       # model S predicts -2 ln(1+eta)
eta_S = -lnclos/2; seS = s_prod/2
# --- free g(J) along the measured R ladder (path reconstruction), g(0)=0 ---
ladder = ['R0/R1', 'R1/R3', 'R3/R5', 'R5/R10', 'R10/R13', 'R13/R15', 'R15/R17', 'R17/R19', 'R19/R21', 'R21/R22', 'R22/R23', 'R23/R24', 'R24/R25', 'R25/R26', 'R26/R27']
gJ = [0]; gJs = [0]; gs2 = 0; Js = [0]
for p in ladder:
    x = d[p]; gJ.append(gJ[-1] + x['r']); gs2 += x['s_stat']**2 + x['s_T']**2 + s_th**2; gJs.append(np.sqrt(gs2)); Js.append(x['J2'])
# r = -(g(J1)-g(J2)) = g(J2)-g(J1) -> cumulative sum gives g(J) - g(0)
# direct R1/R10 vs ladder path:
g_path_1_10 = gJ[Js.index(10)] - gJ[Js.index(1)]; g_dir = d['R1/R10']['r']
eta_firas = 9e-5
R = dict(
  B_cm=B, D_cm=D, E20_check=E(20), E20_cited=806.3828,
  control=dict(n=len(ctrl), mean_r=mc, sem=smc, chi2_dof_stat_only=chi2c, theory_floor=s_th,
               note='R(J)/P(J) share the lower level, so any eta/g(J) cancels; their residuals vs UCL2025 measure theory+systematics only'),
  different_level=dict(n=len(diff), chi2_dof_stat_T_only=chi2d_raw, floor_needed=s_th_d, chi2_with_control_floor=chi0, dof=len(diff)),
  ladder_L=dict(eta=eL, se=sL, eta95=abs(eL)+1.96*sL, chi2=cL, dchi2=chi0-cL, vs_firas=(abs(eL)+1.96*sL)/eta_firas),
  ground_G=dict(eta=eG, se=sG, eta95=abs(eG)+1.96*sG, chi2=cG, dchi2=chi0-cG, vs_firas=(abs(eG)+1.96*sG)/eta_firas),
  closure_S=dict(product=prod, direct=direct, ln_ratio=lnclos, se=s_prod, eta=eta_S, se_eta=seS, eta95=abs(eta_S)+1.96*seS, vs_firas=(abs(eta_S)+1.96*seS)/eta_firas,
                 paper_values='0.739584(31) vs 0.739599(21)',
                 note='any model with well-defined level populations passes closure identically; this bounds only a pair-contextual factor (model S)'),
  path_vs_direct=dict(g_ladder_1_10=g_path_1_10, g_direct_1_10=g_dir),
  gJ=dict(J=Js, g=gJ, sigma=gJs),
  firas_eta95=eta_firas)
json.dump(R, open(OUT+'results/stage48.json', 'w'), indent=1, default=float)
S = [dict(key='li2025_tableS1', ref='Li et al., Sci. Adv. 11, eadz6560 (2025), Supplementary Materials (PMC12442875, pmc-oa-opendata S3)',
          url='https://pmc-oa-opendata.s3.amazonaws.com/PMC12442875.1/sciadv.adz6560_sm.pdf', location='Table S1 (pp. 4-5), all 38 rows transcribed into stage48_run.py',
          quote='"Table S1: Line-intensity ratios of CO lines in the (3-0) band. Uncertainties in parentheses are 1 sigma statistical errors." ; e.g. "R(1) / R(10) 7.18 298.431 0.441206(10) 0.441258 0.441384 0.739608 0.739695 0.739907"', values='{}'),
     dict(key='li2025_closure', ref='same, Supplementary S2', url='https://pmc-oa-opendata.s3.amazonaws.com/PMC12442875.1/sciadv.adz6560_sm.pdf', location='S2, after Eq. S5',
          quote='"A(1)/A(3) x A(3)/A(5) x A(5)/A(10) = 0.875765(25) x 0.936955(22) x 0.901324(18) = 0.739584(31). Furthermore, a separate measurement was conducted for the intensity ratio between R(1) and R(10), yielding 0.739599(21)."', values='{"prod":0.739584,"direct":0.739599}'),
     dict(key='li2025_T', ref='same, Supplementary S1/S2', url='https://pmc-oa-opendata.s3.amazonaws.com/PMC12442875.1/sciadv.adz6560_sm.pdf', location='S1 (line 75), S2',
          quote='"temperature measurement uncertainties were predominantly due to probe calibration error of 5 mK" ; "Note that for a line pair sharing the same lower state, the R(J)/P(J) ratio is temperature-independent." ; "The standard deviation of the obs.-calc. residues is 0.003%" ; "0.007% for R(J)/R(J + 1) and about 0.006% for R(J)/R(J + 2)"', values='{"dT_K":0.005}'),
     dict(key='li2025_E', ref='Li et al. 2025 main text (science.org, fetched via WebFetch)', url='https://www.science.org/doi/full/10.1126/sciadv.adz6560', location='LRT section',
          quote='"For the R(1)/R(17) line pair, the lower-state energies are 3.8450 cm-1 and 587.7209 cm-1" ; "For the R(0)/R(20) pair, the lower-state energies are 0 and 806.3828 cm-1"', values='{"E1":3.8450,"E17":587.7209,"E20":806.3828}'),
     dict(key='stage31_def', ref='This repo, Stage 31 (7) / Stage 47', url='docs/hyperfield/HYPERFIELD_RESULTS.md', location='Stage 31 (7)',
          quote='"If E/A = e^(-x)/(1+eta), then n = 1/((1+eta)e^x - 1)" ; "|eta| < 9.0x10^-5"', values='{"eta_firas":9e-5}')]
with open(OUT+'results/stage48_sources.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(S[0].keys())); w.writeheader(); [w.writerow(s) for s in S]
print(json.dumps({k: v for k, v in R.items() if k != 'gJ'}, indent=1, default=float))
# ---------- robustness of the ladder trend ----------
# (1) does the eta-free control show a J trend? (theory J-dependence)
Jc = np.array([x['J1'] for x in ctrl], float); scf = np.sqrt(sc**2 + s_th**2)
Wc = 1/scf**2; Ac = np.vstack([np.ones_like(Jc), Jc]).T
cov = np.linalg.inv(Ac.T @ (Wc[:, None]*Ac)); pc = cov @ (Ac.T @ (Wc*rc))
# (2) joint fit: ladder eta + temperature offset dT (r_T = -c2 (E1-E2) dT / T^2)
Tm = 298.43; XT = -c2*(E(J1)-E(J2))/Tm**2
Ad = np.vstack([XL, XT]).T; Wd = 1/sig_d**2
covd = np.linalg.inv(Ad.T @ (Wd[:, None]*Ad)); pd = covd @ (Ad.T @ (Wd*rd)); chi_j = np.sum(Wd*(rd - Ad@pd)**2)
# (3) leave-one-out on the ladder fit
loo = []
for i in range(len(rd)):
    m = np.ones(len(rd), bool); m[i] = False; w = Wd[m]
    loo.append(np.sum(w*XL[m]*rd[m])/np.sum(w*XL[m]**2))
R['robust'] = dict(control_J_slope=pc[1], control_J_slope_se=np.sqrt(cov[1, 1]),
                   joint_eta=pd[0], joint_eta_se=np.sqrt(covd[0, 0]), joint_dT_K=pd[1], joint_dT_se=np.sqrt(covd[1, 1]), joint_chi2=chi_j,
                   joint_eta95=abs(pd[0]) + 1.96*np.sqrt(covd[0, 0]), loo_min=min(loo), loo_max=max(loo))
json.dump(R, open(OUT+'results/stage48.json', 'w'), indent=1, default=float)
print(json.dumps(R['robust'], indent=1, default=float))

# ---------------- figure ----------------
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))
a = ax[0]
for i, x in enumerate(ctrl): a.errorbar(x['P'], x['r']*1e6, x['s_stat']*1e6, fmt='o', color='tab:green', ms=4, label='R(J)/P(J): same lower level (eta-free control)' if i == 0 else None)
for i, (x, s) in enumerate(zip(diff, sd)): a.errorbar(x['P'], x['r']*1e6, s*1e6, fmt='s', color='tab:red', ms=4, label='R(J1)/R(J2): different levels (eta-sensitive)' if i == 0 else None)
a.axhline(0, color='k', lw=0.8); a.set_xscale('log'); a.set_xlabel('pressure (Pa)'); a.set_ylabel('ln(A_exp / A_UCL2025)  (ppm)')
a.set_title('residuals vs ab initio (Table S1)'); a.legend(fontsize=7, loc='lower left'); a.grid(alpha=0.3)
b = ax[1]
Jg = np.array(Js); gg = np.array(gJ)*1e6; ss = np.array(gJs)*1e6
b.fill_between(Jg, gg-1.96*ss, gg+1.96*ss, color='tab:blue', alpha=0.2, label='path-summed g(J), 95%')
b.plot(Jg, gg, 'o-', color='tab:blue')
jj = np.linspace(0, 27, 50)
b.plot(jj, jj*eta_firas*1e6, 'k--', label='ladder eta = FIRAS 9e-5')
b.plot(jj, -jj*eta_firas*1e6, 'k--')
b.plot(jj, jj*R['robust']['joint_eta95']*1e6, 'm:', label='ladder 95%% (T offset free): |eta| < %.1e' % R['robust']['joint_eta95'])
b.plot(jj, -jj*R['robust']['joint_eta95']*1e6, 'm:')
b.errorbar([10], [(gJ[Js.index(1)] + g_dir)*1e6], [np.sqrt(d['R1/R10']['s_stat']**2 + d['R1/R10']['s_T']**2 + s_th**2)*1e6], fmt='D', color='tab:orange', label='J=10 via direct R(1)/R(10)')
b.set_xlabel('rotational level J (lower state)'); b.set_ylabel('g(J) = -ln(N_J / N_J,Boltzmann)  (ppm)')
b.set_title('population deviation summed along the measured path'); b.legend(fontsize=7, loc='upper left'); b.grid(alpha=0.3)
c = ax[2]
names = ['ladder (L), T free\nper Delta J=1', 'ground-only (G)', 'pair-contextual (S)\nclosure', 'FIRAS\n(Stage 31)']
vals = [R['robust']['joint_eta95'], R['ground_G']['eta95'], R['closure_S']['eta95'], eta_firas]
cols = ['tab:purple', 'tab:cyan', 'tab:orange', 'k']
c.bar(range(4), vals, color=cols); c.set_yscale('log'); c.set_xticks(range(4)); c.set_xticklabels(names, fontsize=8)
for i, v in enumerate(vals): c.text(i, v*1.15, '%.1e' % v, ha='center', fontsize=8)
c.set_ylabel('|eta| 95% upper limit'); c.set_title('readings of eta: lab CO vs FIRAS'); c.set_ylim(1e-6, 1e-3); c.grid(alpha=0.3, axis='y')
fig.suptitle('Stage 48 - path dependence of eta (CO, Li 2025): ladder |eta| < %.1e, ground-only < %.1e, contextual < %.1e; FIRAS 9e-5' % (vals[0], vals[1], vals[2]), fontsize=11)
fig.tight_layout(); fig.savefig(OUT+'figures/stage48_eta_path.png', dpi=110)

