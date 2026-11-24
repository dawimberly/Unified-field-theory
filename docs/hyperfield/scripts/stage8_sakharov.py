"""Stage 8: Sakharov induced gravity from Standard-Model field content (one loop, heat kernel a1)."""
import json, os, numpy as np
from math import pi, sqrt, log
from scipy.optimize import least_squares
CK = 'results/stage8.json'
out = json.load(open(CK)) if os.path.exists(CK) else {}
def save(): json.dump(out, open(CK, 'w'), indent=1, default=float)
MBAR = 2.435e18; MPL = 1.22089e19; MGUT = 2e16

# ---- 1. per-species induced coefficient. Proper-time cutoff:
#  1/(16 pi G) = (1/(192 pi^2)) * sum_i u_i * J(m_i),  J = int_{1/L^2} ds e^{-m^2 s}/s^2 = L^2 - m^2 ln(L^2/m^2) + O(m^2)
#  "unit" u: real minimally-coupled scalar = 1  (DERIVED from Gamma_E = -1/2 int ds/s Tr e^{-sD}, a1 = (1/6 - xi) R)
units = {
 'real scalar, coupling xi': '6*(1/6 - xi)  (xi=0: +1, conformal xi=1/6: 0)',
 'Weyl fermion': '+1  (Dirac = +2; D^2 = -box + R/4, tr a1 = 4(1/6-1/4)R, fermion loop sign)',
 'massless gauge vector + 2 FP ghosts (Feynman gauge)': '-4  (tr a1 = 4/6 - 1 = -1/3 for A_mu, ghosts -2*(1/6))',
 'massive vector (R_xi: vector + ghosts + eaten Goldstone with the doublet xi)': '-4 + 6(1/6 - xi)  (= -3 at xi=0)',
}
def u_scalar(xi): return 6*(1/6-xi)
def sm_units(xi=0.0, nuR=False, phase='above'):
    nW = 48 if nuR else 45
    if phase == 'above': return dict(higgs=4*u_scalar(xi), fermions=nW, vectors=12*(-4))
    return dict(higgs=1*u_scalar(xi), fermions=nW, vectors=3*(-4+u_scalar(xi))+9*(-4))   # eaten Goldstones carry the doublet xi
tab = []
for xi in (0.0, 1/6, -1.0, 1.0):
    for nuR in (False, True):
        for ph in ('above', 'below'):
            u = sm_units(xi, nuR, ph); tab.append(dict(xi=xi, nuR=nuR, phase=ph, **u, N_eff=sum(u.values())))
# vector ambiguity: (a) heat-kernel -4 per massless vector (Larsen-Wilczek 1996 / Visser 2002 'N0+2N1/2-4N1');
# (b) 'physical polarizations only' (2 minimal scalars each, Kabat 1995 contact-term debate): +2
alt = dict(vectors_as_2_scalars_N_eff_xi0=4+45+12*2, vectors_as_2_scalars_N_eff_xi0_nuR=4+48+24)
out['species'] = dict(units=units, table=tab, vector_ambiguity=alt,
    refs='Larsen & Wilczek, NPB 458, 249 (1996); Kabat, NPB 453, 281 (1995); Visser, MPLA 17, 977 (2002) (from memory, UNVERIFIED citations)',
    note='EWSB conserves N_eff (Goldstones are eaten 1:1): above = below, checked in table')
save()

# ---- 2. quadratic coefficient -> cutoff needed
cut = []
for row in tab:
    if row['phase'] != 'above': continue
    N = row['N_eff']
    if N > 0:
        L = MBAR*sqrt(96*pi**2/N)
        cut.append(dict(xi=row['xi'], nuR=row['nuR'], N_eff=N, sign='correct (attractive)', Lambda_GeV=L, Lambda_over_Mbar=L/MBAR, Lambda_over_MPl=L/MPL,
                        Mbar2_from_GUT_cutoff_fraction=N*MGUT**2/(96*pi**2)/MBAR**2))
    else:
        cut.append(dict(xi=row['xi'], nuR=row['nuR'], N_eff=N, sign='wrong (repulsive) or zero', Lambda_GeV=None))
for k, N in alt.items():
    cut.append(dict(case=k, N_eff=N, Lambda_GeV=MBAR*sqrt(96*pi**2/N), Lambda_over_MPl=MBAR*sqrt(96*pi**2/N)/MPL))
out['quadratic'] = dict(formula='Mbar^2 = 1/(8 pi G) = N_eff Lambda^2/(96 pi^2) (proper-time cutoff; other regulators change the O(1) prefactor)', rows=cut)
save()

# ---- 3. FFZ finiteness sums with SM masses (GeV); d = signed dof (fermion -), u = induced-G units
vt = 246.22
sm = [  # name, multiplicity(Weyl or dof count unit), type, mass
 ('top', 6, 'D', 172.57), ('bottom', 6, 'D', 4.18), ('charm', 6, 'D', 1.27), ('strange', 6, 'D', 0.0935), ('up', 6, 'D', 2.16e-3), ('down', 6, 'D', 4.70e-3),
 ('tau', 2, 'D', 1.77693), ('muon', 2, 'D', 0.10566), ('electron', 2, 'D', 0.000511), ('neutrinos', 3, 'W', 0.0),
 ('W', 2, 'V', 80.369), ('Z', 1, 'V', 91.188), ('photon', 1, 'g', 0.0), ('gluon', 8, 'g', 0.0), ('Higgs', 1, 'S', 125.20)]
def du(t, xi=0.0):
    return {'D': (-2, 1), 'W': (-2, 1), 'V': (3, -3), 'g': (2, -4), 'S': (1, u_scalar(xi))}[t]   # (signed dof per unit, units per unit); D counted per Weyl
S = dict(sum_d=0.0, sum_d_m2=0.0, sum_d_m4=0.0, sum_d_m4lnm2=0.0, sum_u=0.0, sum_u_m2=0.0)
for n, mult, t, m in sm:
    d, u = du(t)
    S['sum_d'] += mult*d; S['sum_u'] += mult*u
    S['sum_d_m2'] += mult*d*m**2; S['sum_u_m2'] += mult*u*m**2
    S['sum_d_m4'] += mult*d*m**4; S['sum_d_m4lnm2'] += mult*d*m**4*(log(m**2) if m > 0 else 0)
ccobs = (2.24e-12)**4   # GeV^4, observed vacuum energy (2.24 meV)^4 (ASSUMED Planck-2018-like)
out['ffz_sm'] = dict(sums=S, conditions='finite & Lambda_cc-free induced G requires sum_d = sum_d m^2 = sum_d m^4 = 0 (vacuum energy), sum_u = sum_u m^2 = 0 (G finite);'
    ' then 1/(16 pi G) = -(1/(192 pi^2)) sum_u m^2 ln m^2  (DERIVED from J(m) expansion)',
    induced_if_finite_GeV2=-S['sum_u_m2']*0/1, sm_finite_part_Mbar2_fraction=abs(sum(mult*du(t)[1]*m**2*log(m**2) for n, mult, t, m in sm if m > 0))/(96*pi**2)/MBAR**2,
    cc_note='SM sum_d m^4 = %.3g GeV^4 vs observed vacuum energy %.2g GeV^4: ratio %.1e (the cosmological-constant problem is inherited)' % (S['sum_d_m4'], ccobs, abs(S['sum_d_m4'])/ccobs))
# repo relation
out['ffz_repo'] = dict(C_ledger=0.05131719, target=1/24, gap=0.05131719-1/24, kappa=-0.00965,
    G_ind_over_GN='(1/24)/C_ledger = %.4f with G_N<->1/24 identification ASSUMED in repo (ffz_distortion_phase.py)' % ((1/24)/0.05131719),
    relation=('The repo coefficient is a multi-loop residue in its own h^2 z/K model whose normalisation to G_N is an ASSUMED identification. '
              'It does not fix the one-loop N_eff found here; at one loop the SM gives N_eff=+1 (xi=0) and needs Lambda=7.5e19 GeV. '
              'kappa=-0.00965 would move C to 1/24, i.e. lower 1/G_ind by 18.8% -- a percent-level statement about a quantity whose overall size is not derived.'))
save()

# ---- 4. minimal finite extensions: heavy scalars (xi free) + heavy Weyl fermions, physical statistics only
# Reduction (DERIVED): d-equations are homogeneous in x=m^2 -> fix 3 mass RATIOS; given masses, the two scalar unit-weights
# U1=a(1-6xi1), U2=b(1-6xi2) follow linearly from sum_u=0 and sum_u x=0; then F=-sum u x ln x is scale covariant (x->lam x: F->lam F),
# so the target 1/(8piG)=Mbar^2 only FIXES the overall heavy scale lam = 96 pi^2/F (needs F>0). G is input, never output.
from scipy.optimize import fsolve
def solve_ext(a, b, cW, eW, seed):
    rng = np.random.default_rng(seed)
    d = np.array([a, b, -2*cW, -2*eW], float)
    def eqs(z):
        x = np.exp(np.concatenate([[0.0], z]))
        return [np.dot(d, x)/np.dot(np.abs(d), x), np.dot(d, x**2)/np.dot(np.abs(d), x**2), np.dot(d, x**2*np.log(x))/np.dot(np.abs(d), x**2*(1+np.abs(np.log(x))))]
    sols = []
    for k in range(60):
        z, info, ier, msg = fsolve(eqs, rng.normal(0, 2, 3), full_output=True, xtol=1e-14)
        if ier == 1 and np.max(np.abs(eqs(z))) < 1e-11:
            x = np.exp(np.concatenate([[0.0], z]))
            if min(abs(x[0]-x[1]), abs(x[2]-x[3])) < 1e-6*x.max() or x.max()/x.min() > 1e8: continue   # degenerate / runaway
            # linear solve for U1, U2: U1+U2 = -(SMu + cW + eW); U1 x1 + U2 x2 = -(cW x3 + eW x4)
            A = np.array([[1, 1], [x[0], x[1]]]); rhs = np.array([-(S['sum_u'] + cW + eW), -(cW*x[2] + eW*x[3])])
            U = np.linalg.solve(A, rhs); u = np.array([U[0], U[1], cW, eW])
            F = -np.dot(u, x*np.log(x))
            key = tuple(np.round(np.sort(x), 6))
            if any(key == s_['key'] for s_ in sols): continue
            sol = dict(key=key, x=x, U=U, F=F)
            if F > 0:
                lam = 96*pi**2/F; m = np.sqrt(lam*x)
                sol.update(masses_over_Mbar=m.tolist(), xi=[(1-U[0]/a)/6, (1-U[1]/b)/6])
            sols.append(sol)
    return [dict(a=a, b=b, cW=cW, eW=eW, mass_ratios_sq=s_['x'].tolist(), F=float(s_['F']), sign_ok=bool(s_['F'] > 0),
                 masses_over_Mbar=s_.get('masses_over_Mbar'), xi=s_.get('xi')) for s_ in sols]
if 'extensions' not in out:
    res = []
    for cW, eW in [(1, 1), (1, 2), (2, 2), (2, 3), (3, 3), (1, 5), (4, 4)]:
        for a in (1, 5, 10, 20, 32, 50):
            b = 2*(cW+eW) - int(S['sum_d']) - a
            if b < 1: continue
            res += solve_ext(a, b, cW, eW, seed=a+7*cW+13*eW)
    good = [r for r in res if r['sign_ok']]
    out['extensions'] = dict(found=len(res), positive_G=len(good), examples=good[:6], wrong_sign_examples=[r for r in res if not r['sign_ok']][:3],
        parameter_count=('Per choice of integer multiplicities: 4 masses + 2 xi = 6 continuous parameters vs 6 equations '
                         '(3 vacuum-energy sums, sum_u=0, sum_u m^2=0, target G) -> isolated solutions. The target equation only sets the overall '
                         'heavy scale, so 1/G is an input, not an output; what is "predicted" are heavy mass ratios and xi values at ~M_Pl, all unobservable. '
                         'Multiplicities are free discrete choices (each choice gives different solutions).'),
        tuning=('sum_d m^4 ln m^2 = 0 is the cosmological-constant condition, needed to (rho_vac,obs)/Mbar^4 ~ 1e-120 relative accuracy; '
                'imposed by hand, not by a symmetry -> standard tuning. Without it, sum_d m^4 ln m^2 ~ Mbar^4 gives rho_vac ~ 1e120 too large.'),
        dof_constraint='heavy sector must supply net bosonic dof a+b-2(cW+eW) = %d to cancel the SM fermion excess' % (-int(S['sum_d'])))
    save()
print(json.dumps({k: out[k] for k in ['quadratic', 'ffz_sm', 'ffz_repo']}, indent=1, default=float)[:5000])
print('ext', out['extensions']['found'], out['extensions']['positive_G'])
for e in out['extensions']['examples'][:5]: print(e)
for r in tab: print(r)
