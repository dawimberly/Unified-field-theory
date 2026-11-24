"""Stage 5a: brane localization as Dan's 'lock' (goals 1-2). Checkpoints into results/stage5.json"""
import json, os, numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import quad
CK = 'results/stage5.json'
out = json.load(open(CK)) if os.path.exists(CK) else {}
def save(): json.dump(out, open(CK, 'w'), indent=1)
HBARC = 1.973269804e-16  # GeV m

def fd_levels(Vfun, L=40.0, N=8001, k=6):
    y = np.linspace(-L, L, N); h = y[1]-y[0]
    d = 2/h**2 + Vfun(y); e = -np.ones(N-1)/h**2
    return eigh_tridiagonal(d, e, select='i', select_range=(0, k-1), eigvals_only=True)

# ---- (1) Jackiw-Rebbi / Rubakov-Shaposhnikov fermion on a phi^4 kink, units mu = 1/thickness = 1
if 'fermion_zero_mode' not in out:
    rows = []
    for s in [0.5, 1.0, 2.0, 3.5, 6.0]:      # s = g v / mu  (Yukawa x vev / kink inverse width)
        V = lambda y: s**2 - s*(s+1)/np.cosh(y)**2     # H_- = D^dagger D, exact Poschl-Teller
        ev = fd_levels(V, k=4)
        exact = [s**2 - (s-j)**2 for j in range(int(np.ceil(s))) if s - j > 0]
        f2 = lambda y: np.cosh(y)**(-2*s)
        nrm = quad(f2, -np.inf, np.inf)[0]
        width = np.sqrt(quad(lambda y: y*y*f2(y), -np.inf, np.inf)[0]/nrm)
        tail = 2*quad(f2, 3, np.inf)[0]/nrm
        rows.append(dict(s=s, numeric_E2=[float(x) for x in ev], exact_bound_E2=exact, continuum_E2=s**2,
                         gap_to_unlocked_over_mu=s, first_excited_locked_over_mu=(np.sqrt(2*s-1) if s > 1 else None),
                         zero_mode_rms_width_over_thickness=width, zero_mode_prob_beyond_3_thickness=tail))
    out['fermion_zero_mode'] = rows
    # topological charge of the kink and index: one normalizable chiral zero mode per unit winding
    y = np.linspace(-30, 30, 20001); phi = np.tanh(y)
    out['kink_topology'] = dict(Q=float((phi[-1]-phi[0])/2), index_theorem='n_L - n_R = Q * sign(g) = 1 (Jackiw-Rebbi 1976; Callias)',
                                note='The zero mode exists for any s>0 and any kink shape with Q=1: protected by winding, not by parameters (DERIVED numerically for 5 s values)')
    # scalar fluctuations: branon (translation zero mode) and shape mode
    ev = fd_levels(lambda y: 4 - 6/np.cosh(y)**2, k=3)
    out['kink_scalar_modes'] = dict(numeric_omega2_over_mu2=[float(x) for x in ev], exact=[0, 3, 4],
                                    labels=['translation zero mode (branon)', 'shape mode sqrt(3) mu', 'continuum 2 mu = bulk scalar mass'])
    save()

# ---- bounds on thickness / unlocked gap, and leakage
if 'bounds' not in out:
    CI = 35.8e3; mUED = 1.86e3
    YR = 3.156e7; HBAR_eVs = 6.582119569e-16
    tau_e = 6.6e28*YR
    G_oPs = HBAR_eVs/142.05e-9
    out['bounds'] = dict(
        thickness_from_contact_interactions_m=HBARC/CI, CI_TeV=35.8, CI_source='ATLAS dilepton 139 fb-1, JHEP 11 (2020) 005 (LL constructive, verified)',
        thickness_from_KK_SM_m=HBARC/mUED, KK_source='mUED recast 1/R>1.86 TeV (arXiv:2110.00500)',
        mapping='ASSUMED: Lambda_CI ~ 1/delta and m_KK ~ 1/delta up to O(1)',
        unlocked_gap_min_TeV_if_s1=1.86, unlocked_gap_min_TeV_CI=35.8,
        leakage=dict(
            electron_Gamma_max_eV=HBAR_eVs/tau_e, electron_Gamma_over_me=HBAR_eVs/tau_e/0.51099895e6,
            electron_source='Borexino e->nu gamma, tau>6.6e28 yr (PRL 115, 231802), used as charge-nonconservation proxy (ASSUMED)',
            oPs_invisible_BR_max=4.2e-7, oPs_Gamma_inv_max_eV=4.2e-7*G_oPs,
            ee_to_invisible_BR_max=2.1e-8, oPs_source='Badertscher et al. PRD 75, 032004 (2007), verified'),
        model_prediction='Exact kink model: every state with E < g v is a true bound state; leakage rate = 0 identically (DERIVED). Leakage needs a gapless bulk channel (RS2-like) or E > g v.')
    save()

# ---- (2) hysteresis: phi^4 cannot, phi^6 can
if 'hysteresis' not in out:
    # phi^4: s = g v / mu with mu = v sqrt(lambda/2)  ->  s = g sqrt(2/lambda), independent of v (DERIVED)
    res = dict(phi4='s = g*sqrt(2/lambda) is independent of v: number of locked states cannot change as v(rho) drops; the wall simply vanishes at v=0 (continuous, no hysteresis)')
    # phi^6: V = (phi^2-1)^2 (phi^2 + a)/4 ; scan a
    def minima(a):
        p = np.linspace(-1.6, 1.6, 32001); V = (p**2-1)**2*(p**2+a)/4
        i = np.where((V[1:-1] < V[:-2]) & (V[1:-1] < V[2:]))[0]+1
        return p[i], V[i]
    scan = []
    for a in [1.0, 0.5, 0.2, 0.1, 0.03, 0.01, 1e-3, 0.0, -0.1, -0.3, -0.5, -0.9, -1.1]:
        pm, Vm = minima(a)
        sym_min = bool(np.any(np.abs(pm) < 1e-3)); brk = bool(np.any(np.abs(np.abs(pm)-1) < 5e-3))
        row = dict(a=a, symmetric_is_local_min=sym_min, broken_is_local_min=brk)
        if brk and a > -1:
            # BPS kink -1 -> +1: dy = dphi / sqrt(2 V);  thickness = y(0.9)-y(-0.9)
            f = lambda p: 1/np.sqrt(2*max((p**2-1)**2*(p**2+a)/4, 1e-300))
            if a > 0:
                row['wall_thickness_10_90'] = 2*quad(f, 0, 0.8)[0]
            else:
                row['wall_thickness_10_90'] = float('inf')
        scan.append(row)
    res['phi6_scan'] = scan
    res['coexistence_window'] = 'symmetric and broken both local minima for -1 < a < 0 ... and a > 0 symmetric is local max? see scan'
    out['hysteresis'] = res
    save()
print(json.dumps(out, indent=1, default=str)[:6000])
