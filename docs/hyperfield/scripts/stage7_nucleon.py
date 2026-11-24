"""Stage 7: nucleon mass budget, gravitational form factors/pressure, brane gravity at fm, equivalence principle."""
import json, os, numpy as np
from math import gamma, pi, sqrt
CK = 'results/stage7.json'
out = json.load(open(CK)) if os.path.exists(CK) else {}
def save(): json.dump(out, open(CK, 'w'), indent=1, default=float)
HBARC = 0.1973269804  # GeV fm
GeVfm3_Pa = 1.602176634e-10/1e-45   # 1 GeV/fm^3 in Pa
MN = 0.93827

# 1. mass budget (Yang et al. PRL 121, 212001 (2018), arXiv:1808.08677, verified; MSbar 2 GeV)
mb = dict(source='chiQCD, Yang et al. PRL 121 212001 (2018), verified abstract',
          quark_condensate_Hm=(0.09, 0.022), quark_energy_HE=(0.33, 0.057), glue_field_HG=(0.37, 0.064), quarter_anomaly=(0.23, 0.014))
mb['sum'] = sum(v[0] for k, v in mb.items() if isinstance(v, tuple))
mb['field_not_Higgs_fraction'] = 1 - mb['quark_condensate_Hm'][0]
mb['naive_valence_current_mass_fraction'] = (2*2.16e-3+4.70e-3)/MN    # PDG m_u, m_d (MSbar 2 GeV)
mb['mapping'] = "ASSUMED: Dan's 'locked field energy' = everything except H_m: 91(2)% of the nucleon mass is field energy (E=mc^2 of confined gluon+quark kinetic fields)"
# other systems: fraction of mass NOT from Higgs/rest masses
mb['systems'] = {
  'electron': dict(field_fraction=0.0, note='elementary, mass = Yukawa (plus ~alpha/pi radiative)'),
  'pion': dict(quark_mass_fraction=0.5, field_fraction=0.5, note='DERIVED (Feynman-Hellmann + GMOR m_pi^2 ~ m_q): m_q dm_pi/dm_q = m_pi/2'),
  'nucleon': dict(quark_mass_fraction=0.09, field_fraction=0.91),
  'glueball 0++ (pure YM)': dict(field_fraction=1.0, note='~1.7 GeV, lattice; not cleanly observed (UNVERIFIED mass)'),
  'nucleus (Fe56)': dict(binding_fraction=8.79/931.5, note='nuclear binding = -0.94% of mass on top of nucleon field energy'),
  'neutron star 1.4 Msun': dict(note='gravitational binding fraction computed in stage7_ns.py'),
}
out['mass_budget'] = mb

# 2. gravitational form factors -> pressure/shear (dipoles from Hackett, Pefkou, Shanahan PRL 132 251904 (2024), arXiv:2310.08484 Table 3, verified; m_pi=170 MeV)
fits = dict(A_g=(0.501, 1.262), A_q=(0.510, 1.477), J_g=(0.255, 1.399), J_q=(0.251, 1.62), D_g=(-2.57, 0.538), D_q=(-1.30, 0.81))
r = np.linspace(0.005, 2.5, 5000)       # fm
def Dtil(D0, L, r):                      # FT of dipole D0/(1+Q^2/L^2)^2 ; L in fm^-1
    return D0*L**3/(8*pi)*np.exp(-L*r)
def p_s(D0, LGeV, r):
    L = LGeV/HBARC; M = MN/HBARC                     # fm^-1
    # p = 1/(6M) 1/r^2 d/dr r^2 d/dr Dt ;  s = -1/(4M) r d/dr (1/r d/dr Dt)
    pre = D0*L**3/(8*pi)
    p = pre/(6*M)*(L**2 - 2*L/r)*np.exp(-L*r)          # fm^-4
    s = -pre/(4*M)*(L**2 + L/r)*np.exp(-L*r)           # fm^-4  (positive for D0<0)
    return p*HBARC, s*HBARC                             # GeV/fm^3
pg, sg = p_s(*fits['D_g'], r); pq, sq = p_s(*fits['D_q'], r); pt, st = pg+pq, sg+sq
def zero(r, f):
    i = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]; return [float(r[k]) for k in i]
vonlaue = float(np.trapezoid(r**2*pt, r)*4*pi)
pk = int(np.argmax(r**2*pt)); nk = int(np.argmin(r**2*pt))
# radii
M = MN
A0 = fits['A_g'][0]+fits['A_q'][0]; Ap = 2*fits['A_g'][0]/fits['A_g'][1]**2 + 2*fits['A_q'][0]/fits['A_q'][1]**2   # GeV^-2
D0 = fits['D_g'][0]+fits['D_q'][0]
r2_mass = (6*Ap - 3*D0/(2*M**2))/A0*HBARC**2
r2_scalar = 6*Ap/A0*HBARC**2
intD = fits['D_g'][0]*fits['D_g'][1]**2 + fits['D_q'][0]*fits['D_q'][1]**2
r2_mech = 6*D0/intD*HBARC**2
gfrac_D = fits['D_g'][0]/D0
out['gff'] = dict(source='dipole fits Hackett-Pefkou-Shanahan 2024 (m_pi~170 MeV, a=0.091 fm, single ensemble)',
    exp_ref='Burkert-Elouadrhiri-Girod Nature 557 396 (2018): peak ~1e35 Pa near centre, zero crossing ~0.6 fm, binding peak ~0.8 fm, d1(0)=-2.04(14)(33) (verified abstract/JLab); Matters Arising 2019 questions uncertainties',
    r_zero_total_fm=zero(r, pt), r_zero_gluon_fm=zero(r, pg), r_zero_quark_fm=zero(r, pq),
    analytic_zero='DERIVED: dipole D gives p ~ (L^2-2L/r) e^{-Lr} -> node at r0 = 2/Lambda: gluon %.3f fm, quark %.3f fm' % (2*HBARC/0.538, 2*HBARC/0.81),
    r2p_max_at_fm=float(r[pk]), r2p_min_at_fm=float(r[nk]),
    p_at_r_fm={str(x): float(np.interp(x, r, pt))*GeVfm3_Pa for x in (0.1, 0.25, 0.5, 0.8, 1.0)},
    p_units='Pa', mean_core_pressure_Pa=float(np.trapezoid(r[r < zero(r, pt)[0]]**2*pt[r < zero(r, pt)[0]], r[r < zero(r, pt)[0]])/np.trapezoid(r[r < zero(r, pt)[0]]**2, r[r < zero(r, pt)[0]]))*GeVfm3_Pa,
    von_laue_integral=vonlaue, shear_positive_everywhere=bool(np.all(st > 0)),
    D0_total=D0, gluon_share_of_D=gfrac_D, r_mass_fm=sqrt(r2_mass), r_scalar_fm=sqrt(r2_scalar), r_mech_fm=sqrt(r2_mech),
    r_charge_fm=0.84075, r_charge_src='CODATA 2022 (from memory, UNVERIFIED digits)',
    note_dipole='dipole p(r) ~ 1/r at r->0 (fit artefact): "peak pressure" only meaningful as r^2 p or shell average')
# 'pressure makes gravity repel?' -> Tolman active mass
eps_dens = None
out['gff']['tolman'] = ('DERIVED: active gravitational mass = int (eps + 3p) d^3r (Tolman). von Laue: int p d^3r = 0 for any stable system -> '
                        'pressure adds exactly zero net gravitating mass; locally the repulsive core (p>0) gravitates MORE (eps+3p > eps), '
                        'the confining shell (p<0) gravitates less. So positive pressure does NOT make gravity repel; it is the negative (binding) '
                        'pressure that reduces active mass, and only locally.')
# local 3p/eps at r=0.25 fm using energy density of A-dipole (approx eps ~ M * FT[A]) -- order of magnitude
def Atil(A0, LGeV, r):
    L = LGeV/HBARC; return A0*L**3/(8*pi)*np.exp(-L*r)
eps = MN*(Atil(*fits['A_g'], r)+Atil(*fits['A_q'], r))     # GeV/fm^3 (leading term)
out['gff']['three_p_over_eps'] = {str(x): float(3*np.interp(x, r, pt)/np.interp(x, r, eps)) for x in (0.1, 0.25, 0.5, 0.8)}
save()
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for f, lab in [(pg, 'gluon'), (pq, 'quark'), (pt, 'total')]:
    ax[0].plot(r, r**2*f, label=lab)
ax[0].axhline(0, c='k', lw=.5); ax[0].set_xlabel('r [fm]'); ax[0].set_ylabel('r^2 p(r) [GeV/fm]'); ax[0].legend(); ax[0].set_xlim(0, 2)
ax[0].set_title('proton pressure from lattice GFF dipoles (HPS 2024)')
for f, lab in [(sg, 'gluon'), (sq, 'quark'), (st, 'total')]:
    ax[1].plot(r, r**2*f, label=lab)
ax[1].set_xlabel('r [fm]'); ax[1].set_ylabel('r^2 s(r) [GeV/fm]'); ax[1].set_xlim(0, 2); ax[1].legend(); ax[1].set_title('shear force')
fig.tight_layout(); fig.savefig('figures/stage7_proton_pressure.png', dpi=110)

# 3. brane gravity at nucleon scale
MPLR = 2.435e18; HB = 1.973269804e-16
def R_of(n, MD): return (MPLR**2/(MD*1e3)**(n+2))**(1/n)*HB
rows = []
for n, MD in [(4, 7.1), (5, 6.4), (6, 5.9), (7, 5.2), (8, 6.0), (8, 10.0)]:
    R = R_of(n, MD); S = 2*pi**(n/2)/gamma(n/2)
    rr = 1e-15
    G_ratio_torus = 1 + (4/3)*(S*gamma(n)*(R/rr)**n - 1) if rr < 0.5*R else None
    row = dict(n=n, M_D_TeV=MD, R_m=R, Geff_over_G_at_1fm=G_ratio_torus,
               grav_over_strong_at_1fm=G_ratio_torus*5.906e-39/1.0 if G_ratio_torus else None)
    if n == 8:
        lam = R/np.sqrt(2)   # E8 self-dual: dual min norm 2 at equal volume (2 pi R)^8 -> m1 = sqrt2/R
        row.update(E8_N1=240, E8_alpha=320.0, E8_lambda_m=lam, Z8_alpha=64/3, Z8_lambda_m=R,
                   E8_note='alpha=(4/3)*240=320 at lambda=R/sqrt2; short-distance (r<<R) G_eff identical to torus (Weyl law, shape-independent)')
    rows.append(row)
# RS1
Lam, kM = 11.744e3, 0.1; m1 = 3.8317*kM*Lam; aRS = (4/3)*(MPLR/Lam)**2
rows.append(dict(model='RS1', m1_GeV=m1, lambda_m=HB/m1, alpha_per_mode=aRS, Geff_minus1_at_1fm=float(aRS*np.exp(-1e-15/(HB/m1)))))
out['brane_fm'] = dict(rows=rows,
    neutron_bound=dict(src='Haddock et al. PRD 97 062002 (2018), arXiv:1712.02984 (verified): alpha(0.1 nm) = (1.6 +/- 16.2)e22 at 95% CL -> |alpha| < ~1.8e22 at lambda=0.1 nm; Kamiya PRL 114 161101 (2015) similar above 0.1 nm',
                       comparison='ADD/E8 alpha = 10-320 at lambda <= 5e-13 m; at lambda = 0.1 nm their Yukawa is e^{-1e-10/5e-13}= e^{-200}: >20 orders below the neutron limit. E8 alpha=320 vs Z8 alpha=21.3 changes nothing observable.'),
    nucleon_EMT=('KK gravitons couple to the same T_mu_nu whose matrix elements are the GFFs; the correction to graviton-proton amplitudes is ~ (G_eff/G - 1) * G m^2 ~ '
                 '1e-25 of the strong force at 1 fm even for n=4 -- no effect on GFF extractions (DVCS measures them via quark GPDs, not via gravity).'))
save()

# 5. equivalence principle (MICROSCOPE final, Touboul et al. PRL 129 121102 (2022), arXiv:2209.15487, verified)
eta, eta_err = -1.5e-15, np.hypot(2.3e-15, 1.5e-15)
eta_lim1 = abs(eta)+eta_err; eta_lim2 = abs(eta)+2*eta_err
mats = {'Ti48': dict(Z=22, A=48, BA=8.723), 'Pt195': dict(Z=78, A=195, BA=7.915)}
def comp(m):
    Z, A = m['Z'], m['A']; N = A-Z
    coul = 0.711*Z**2/A**(1/3)/A/0.9315      # Coulomb energy fraction (SEMF a_C=0.711 MeV), per nucleon mass 931.5 MeV
    qmass = (N-Z)/A*2.52/938.9/2              # strong-isospin quark-mass part: (m_n-m_p)_QCD=2.52 MeV (BMW 2015, UNVERIFIED digits); N-Z excess /2 relative to symmetric
    return dict(binding_fraction=m['BA']/931.5, coulomb_fraction=coul/1000, isospin_quarkmass_fraction=qmass, sigma_fraction=0.09)
c = {k: comp(v) for k, v in mats.items()}
d = {k: abs(c['Pt195'][k]-c['Ti48'][k]) for k in c['Ti48']}
out['equivalence'] = dict(eta_TiPt=eta, eta_err_1sigma=eta_err, eta_abs_1sigma=eta_lim1, eta_abs_2sigma=eta_lim2,
    compositions=c, differences_PtminusTi=d,
    coupling_difference_bounds_1sigma={k: (eta_lim1/v if v > 0 else None) for k, v in d.items()},
    coupling_difference_bounds_2sigma={k: (eta_lim2/v if v > 0 else None) for k, v in d.items()},
    interpretation=('ASSUMED: gravity (or a brane mode) couples to energy component X with strength (1+eps_X). eta = sum eps_X * Delta f_X. '
                    'The sigma-term fraction (9%) is the same for both nuclei to O(1e-4), so the bulk "lock" fraction itself is untested; '
                    'MICROSCOPE tests only composition differences ~1e-4-1e-3.'))
save()
print(json.dumps({k: out[k] for k in ['gff', 'brane_fm', 'equivalence']}, indent=1, default=float)[:7000])
