"""Stage 4: brane-world placements. Which fields live where, and what survives the data.
Writes results/stage4.json (checkpointed per section) and figures/stage4_*.png"""
import json, os, sys
import numpy as np
import sympy as sp
from scipy.special import jn_zeros
from math import gamma, pi
CK = 'results/stage4.json'
out = json.load(open(CK)) if os.path.exists(CK) else {}
def save(): json.dump(out, open(CK, 'w'), indent=1)

HBARC = 1.973269804e-16      # GeV m
MPLR = 2.435e18              # reduced Planck mass [GeV] (PDG)  -- convention used by ATLAS/PDG: Mbar^2 = M_D^{n+2} R^n
MPL = 1.22089e19             # non-reduced
MPC = 3.0857e22              # m
H0 = 67.4 * 1e3 / MPC        # s^-1  (Planck 2018, ASSUMED)
C = 2.99792458e8
LH = C / H0                  # Hubble length [m]

# ---------------- 1. ADD: derivation and R(n, M_D) -----------------------------------------------
if 'add_derivation' not in out:
    Ms, n, R, Mb = sp.symbols('M_star n R Mbar', positive=True)
    # S = (M*^{n+2}/2) Int d^4x d^ny sqrt(G) R_{4+n};  factorized flat torus: R_{4+n}=R_4, Int d^ny = (2 pi R)^n
    Vn = (2*sp.pi*R)**n
    rel = sp.Eq(Mb**2, Ms**(n+2)*Vn)
    out['add_derivation'] = dict(volume_convention=str(rel), pdg_convention='Mbar^2 = M_D^(n+2) R^n  =>  M_D = (2 pi)^(n/(n+2)) M_star',
                                 note='DERIVED: dimensional reduction of the bulk Einstein-Hilbert term over a flat n-torus of radius R')
    save()

def R_of(nn, MD_TeV):
    MD = MD_TeV * 1e3
    return (MPLR**2 / MD**(nn+2))**(1/nn) * HBARC
def MD_of(nn, R_m):
    return (MPLR**2 / (R_m/HBARC)**nn)**(1/(nn+2)) / 1e3   # TeV

atlas = {2: 11.2, 3: 8.5, 4: 7.1, 5: 6.4, 6: 5.9}       # ATLAS monojet 139 fb-1, PRD 103, 112006 (2021), Table 10 (verified)
cms_n2 = 10.7                                            # CMS monojet 137 fb-1 (via PDG 2025 review)
# Hannestad & Raffelt PRD 67, 125008 (2003) Table V: upper limits on R [m] (verified from arXiv hep-ph/0304029)
hr_sn = {1: 4.9e2, 2: 1.0e-6, 3: 1.1e-9, 4: 3.8e-11, 5: 4.9e-12, 6: 1.2e-12, 7: 4.4e-13}
hr_ns = {1: 8.3, 2: 5.9e-8, 3: 1.3e-10, 4: 5.9e-12, 5: 9.4e-13, 6: 2.8e-13, 7: 1.2e-13}
torsion = {1: 30e-6, 2: 37e-6}  # Lee+2020 (n=1 torus radius <30 um); Kapner+2007 via PDG listing (n=2: <37 um)
solar = 1e11 * 0.01              # ASSUMED: 1/r^2 verified at ~AU; n=1 with R ~ 1e12 m excluded outright

if 'add_table' not in out:
    rows = []
    for nn in range(1, 8):
        r = dict(n=nn, R_for_MD_1TeV=R_of(nn, 1), R_for_MD_3TeV=R_of(nn, 3), R_for_MD_10TeV=R_of(nn, 10),
                 alpha_naive_2n=2*nn, alpha_spin2_8n_over_3=8*nn/3)
        lims = {}
        if nn in atlas: lims['LHC monojet (ATLAS)'] = R_of(nn, atlas[nn])
        lims['SN1987A cooling (HR03)'] = hr_sn[nn]
        lims['NS heating (HR03)'] = hr_ns[nn]
        if nn in torsion: lims['torsion balance'] = torsion[nn]
        best = min(lims, key=lims.get)
        r['limits_R_m'] = lims
        r['R_max_m'] = lims[best]; r['strongest'] = best
        r['MD_min_TeV_implied'] = MD_of(nn, lims[best])
        r['KK_mass_min_eV'] = HBARC*1e9 / lims[best]
        r['hierarchy_solved_MD_le_10TeV'] = bool(r['MD_min_TeV_implied'] <= 10)
        rows.append(r)
    out['add_table'] = rows
    out['add_notes'] = dict(
        n1='R(1 TeV)=%.2e m: solar-system scale, excluded by planetary 1/r^2 tests (DERIVED)' % R_of(1, 1),
        n7='ATLAS table stops at n=6; for n=7 only astrophysical limits are used (collider limit not verified)',
        pdg_ns_n2='PDG listing quotes R<1.6e-10 m for n=2 from NS heating (Hannestad-Raffelt 2002 M>1700 TeV); the 2003 paper Table V gives 5.9e-8 m. Unresolved convention/revision difference; the weaker 2003 value is used.',
        astro_caveat='SN/NS bounds assume KK gravitons decay to SM photons on our brane; they weaken if KK gravitons decay to other branes (ASSUMED).')
    save()

# Yukawa coefficient check: sum the full KK lattice of a torus vs the lowest shell
if 'yukawa_check' not in out:
    chk = []
    for nn in [1, 2, 3, 6]:
        K = {1: 40, 2: 30, 3: 15, 6: 5}[nn]
        grids = np.meshgrid(*[np.arange(-K, K+1)]*nn, indexing='ij')
        kk = np.sqrt(sum(g.astype(float)**2 for g in grids)).ravel(); kk = kk[kk > 0]
        for x in [3.0, 10.0]:
            full = (4/3)*np.sum(np.exp(-kk*x))
            lowest = (8*nn/3)*np.exp(-x)
            chk.append(dict(n=nn, r_over_R=x, full_sum=float(full), lowest_shell=float(lowest), ratio=float(full/lowest)))
    out['yukawa_check'] = chk
    save()

# ---------------- 2. RS1 / RS2 ----------------------------------------------------------------------
if 'rs' not in out:
    x1 = float(jn_zeros(1, 1)[0])
    rs1 = []
    for Lam in [0.246, 1, 3, 10, 11.75, 30]:   # TeV, Lambda_pi = Mbar e^{-k pi r_c}
        krc = np.log(MPLR/(Lam*1e3))/np.pi
        rs1.append(dict(Lambda_pi_TeV=Lam, k_rc=krc, m1_TeV_k0p1=x1*0.1*Lam, m1_TeV_k0p01=x1*0.01*Lam,
                        tuning_v2_over_Lam2=(0.246/Lam)**2))
    lhc = dict(atlas_diphoton={'0.01': 2.2, '0.05': 3.9, '0.1': 4.5},  # TeV, ATLAS PLB 822 (2021) 136651, verified
               cms_diphoton='2.2-5.6 TeV for k/Mbar 0.01-0.2 (PDG 2026 review)', cms_dilepton='2.47-4.78 TeV for k/Mbar 0.01-0.1 (PDG 2026 review)')
    Lmin = {c: m/(x1*float(c)) for c, m in lhc['atlas_diphoton'].items()}
    ell = 44e-6   # Kapner+2007 RS2 bound (secondary citation)
    k_GeV = HBARC/ell
    M5_rs2 = (MPLR**2 * k_GeV)**(1/3)
    rs2_dev = {f'{r*1e6:.0f}um': 2*ell**2/(3*r**2) for r in [52e-6, 1e-4, 1e-3]}
    out['rs'] = dict(x1_first_J1_zero=x1, rs1=rs1, lhc=lhc, Lambda_pi_min_TeV=Lmin,
                     krc_max_allowed={c: float(np.log(MPLR/(L*1e3))/np.pi) for c, L in Lmin.items()},
                     rs2=dict(ell_max_m=ell, M5_min_GeV=M5_rs2, deviation_at=rs2_dev,
                              note='V=-GMm/r(1+2 l^2/(3 r^2)); Mbar^2 = M5^3/k (ASSUMED convention)'))
    save()

# ---------------- 3. split placements: what hydrogen can see ----------------------------------------
if 'split' not in out:
    alpha = 1/137.035999
    a0 = 5.29177210903e-11
    def rp_equiv_bound(Nextra=1, drp2_fm2=0.024):
        return np.sqrt(drp2_fm2/(12*Nextra))*1e-15
    v = 246.22
    mnu_for_Mstar = lambda Mstar_GeV, y=1.0: y*v/np.sqrt(2)*Mstar_GeV/MPLR*1e9   # eV  (n-independent volume suppression)
    Mstar_for_mnu = lambda m_eV, y=1.0: m_eV*1e-9*MPLR/(y*v/np.sqrt(2))
    out['split'] = dict(
        photon_in_bulk=dict(hydrogen_R_max_m=rp_equiv_bound(), source='stage 2 (DERIVED, r_p-degeneracy limited)'),
        UED_all_in_bulk=dict(hydrogen_R_max_m=rp_equiv_bound(), LHC_mUED_inv_R_min_TeV=1.86, LHC_R_max_m=HBARC/1.86e3,
                             LHC_source='arXiv:2110.00500 recast of ATLAS Run 2 (phenomenology paper)',
                             hydrogen_weaker_by=rp_equiv_bound()/(HBARC/1.86e3)),
        gluons_only=dict(hydrogen='no tree-level effect: Coulomb potential unchanged (DERIVED)'),
        higgs_only=dict(hydrogen_ratio_higgs_to_coulomb_at_a0=float((2.94e-6*1.0)/(4*np.pi*alpha) * np.exp(-125.25/HBARC*1e0*a0)) ,
                        note='Higgs-exchange e-p force is Yukawa with range 1.6e-18 m: exp(-a0 m_h) underflows to 0'),
        nuR_in_bulk=dict(mnu_eV_for_Mstar_1TeV=mnu_for_Mstar(1e3), mnu_eV_for_Mstar_10TeV=mnu_for_Mstar(1e4),
                         Mstar_TeV_for_mnu_0p05eV=Mstar_for_mnu(0.05)/1e3,
                         note='m_nu = y v/sqrt2 * M*/Mbar (volume suppression, any n); ASSUMED y=1'),
        gravity_only=dict(hydrogen_fractional_shift_n2_R30um=5.5e-27, source='stage 2'))
    save()

# ---------------- 4. DGP and induced gravity -----------------------------------------------------------
if 'dgp' not in out:
    def M5_for_rc(rc_m):
        rc = rc_m/HBARC
        return (MPLR**2/(2*rc))**(1/3)
    bounds = {'LLR r_c>120 Mpc': 120*MPC, 'normal-branch cosmology H0 r_c>3.5 (Lombriser+2009)': 3.5*LH,
              'SPT+Planck PR4 H0 r_c>0.71 (2026)': 0.71*LH, 'self-accelerating needs r_c~1/H0 (excluded)': LH}
    out['dgp'] = dict(Hubble_length_m=LH, rows=[dict(bound=k, r_c_m=v, M5_max_GeV=M5_for_rc(v)) for k, v in bounds.items()])
    # induced fraction of Mbar^2 in ADD / RS1 if brane matter induces M_ind^2 ~ |S| Lambda^2/(16 pi^2) with Lambda=M_D
    S_dan = 265.0   # |S| from earlier rounds (UV gauge-Kabat dominated), ASSUMED to carry over
    ind = []
    for MD in [1, 3, 10]:
        Mind2 = S_dan*(MD*1e3)**2/(16*np.pi**2)
        ind.append(dict(M_D_TeV=MD, M_ind2_GeV2=Mind2, induced_fraction_of_Mbar2=Mind2/MPLR**2,
                        N_KK_below_MD_n2=np.pi*(MPLR/(MD*1e3))**2, N_KK_below_MD_n6=(np.pi**3/6)*(MPLR**2/(MD*1e3)**2)))
    out['induced'] = dict(ADD_or_RS1=ind,
                          ffz_gap=dict(gap=0.009651, kappa_bare=-0.00965, target=1/24,
                                       gap_fraction_of_target=0.009651/(1/24)))
    # coincidence: TeV^2/Mbar vs dark-energy scale
    rho_de_quarter = 2.25e-12  # GeV  (~2.25 meV, ASSUMED from Omega_L=0.69, H0=67.4)
    out['coincidences'] = dict(TeV2_over_Mbar_GeV=1e6/MPLR, dark_energy_quarter_GeV=rho_de_quarter,
                               ratio=rho_de_quarter/(1e6/MPLR), R_n2_MD1TeV_m=R_of(2, 1),
                               dark_energy_length_m=HBARC/rho_de_quarter)
    save()

# ---------------- figures ------------------------------------------------------------------------------
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(9, 6))
ns = np.arange(1, 8)
for MD, c in [(1, 'C0'), (3, 'C1'), (10, 'C2')]:
    ax.semilogy(ns, [R_of(n, MD) for n in ns], 'o-', color=c, label=f'R needed, M_D = {MD} TeV')
ax.semilogy(list(atlas), [R_of(n, m) for n, m in atlas.items()], 'kv', ms=9, label='LHC monojet max R (ATLAS 2021)')
ax.semilogy(ns, [hr_sn[n] for n in ns], 'm^', label='SN1987A max R (HR 2003)')
ax.semilogy(ns, [hr_ns[n] for n in ns], 'rs', label='NS heating max R (HR 2003)')
ax.semilogy([1, 2], [torsion[1], torsion[2]], 'b*', ms=12, label='torsion balance max R')
ax.axhline(1.5e11, color='gray', ls=':'); ax.text(1.1, 3e11, '1 AU', color='gray')
ax.set_xlabel('number of extra dimensions n (gravity only in bulk; n=6: 10D string, n=7: 11D M-theory)')
ax.set_ylabel('compact radius R [m]'); ax.set_title('Stage 4 ADD: radius required vs largest radius allowed')
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig('figures/stage4_add_R_vs_n.png', dpi=110)

fig, ax = plt.subplots(figsize=(9, 6))
r = np.logspace(-9, -2, 400)
for n, MD in [(2, 11.2), (3, 8.5), (6, 5.9)]:
    Rn = R_of(n, MD)
    ax.loglog(r, (8*n/3)*np.exp(-r/Rn), label=f'ADD n={n}, M_D={MD} TeV (R={Rn:.1e} m), alpha=8n/3')
ax.loglog(r, 2*44e-6**2/(3*r**2), label='RS2, l = 44 um (Kapner 2007 bound)')
ax.axvline(52e-6, color='k', ls='--'); ax.text(55e-6, 1e-6, 'Lee 2020 closest\nseparation 52 um', fontsize=8)
ax.set_ylim(1e-8, 1e3); ax.set_xlabel('separation r [m]'); ax.set_ylabel('fractional deviation from Newton  dV/V')
ax.set_title('Stage 4: largest still-allowed short-range gravity deviations'); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig('figures/stage4_gravity_deviation.png', dpi=110)
print(json.dumps(out, indent=1, default=float))
