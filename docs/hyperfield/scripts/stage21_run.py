"""Stage 21 - design changes that LOOSEN the acceleration-noise requirement (not meet it).

Builds on Stage 18 (exact correlated-phase PPT model: per-mass requirement 0.44 fm s^-2/rtHz at the trapped
diamond cheap point; ground residual +6.0 orders seismic, +8.0 with plate surface E-field noise; space-grade +4.0).
Every lever is scored with the same exact model:  rho(Phi, gamma t) x Gaussian which-path phases, N5 = 25/lambda^2.
Residual for a noise spectrum S(f) under protocol P:  orders = 0.5 log10( Var_P[S] / sigma^2_allowed(P) ),
Var_P[S] = kappa^2 int S(f) |W_P(f)|^2 df  (1 mHz..300 Hz), W_P = Fourier transform of the branch profile w(t).
Labels: DERIVED / CITED (URL) / ASSUMED.  Nothing is tuned.
"""
from __future__ import annotations
import json, os, sys
from math import pi, sqrt, log10
import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import stage18_run as s18                                                    # noqa: E402
from stage18_run import rho_noisy, lmin, nlnm_asd                            # noqa: E402
from stage12_lib import G, HBAR, K_B                                         # noqa: E402

M, DX, D, R, PHI = s18.M, s18.DX, s18.D, s18.R, s18.PHI
TAU = s18.TAU_ACC
T_INT = s18.T_INT
GT_RATE = s18.GT / s18.T_TOT          # gas+BB exponent per second (gas-saturated, ~Dx independent; DERIVED from Stage 16)
# Shot criterion for ALL designs (ASSUMED, Stage 16 N_SHOTS_MAX): N5 = 25/lambda^2 <= 1e4  ->  lambda <= -0.05.
# (Stage 18's 0.44 fm used N5 = 5049, which leaves noise only 4% of the slack; with a fixed 5049 every protocol
#  that adds any duration is infeasible, so that criterion cannot rank levers.  Both are reported.)
N_CAP = 1e4
L_TARGET = -5.0 / sqrt(N_CAP)
L_TARGET_S18 = s18.L_TARGET
N5_BASE = N_CAP
KAPPA = s18.KAPPA
F_LO, F_HI = 1e-3, 300.0              # same band as Stage 18 (drifts < 1 mHz tracked; ASSUMED)
FGRID = np.geomspace(F_LO, F_HI, 2500)
C_LIGHT = 2.99792458e8
G0 = 9.81


# ---------------------------------------------------------------------------
# protocols: branch profile w(t) (linear ramps of TAU, swaps of 2 TAU; ASSUMED shape as in Stage 18)
# ---------------------------------------------------------------------------
def profile(n_swap=0, t_hold=T_INT, dt=5e-4):
    segs = [(TAU, 0.0, 1.0)]
    hseg = t_hold / (n_swap + 1)
    sgn = 1.0
    for k in range(n_swap + 1):
        segs.append((hseg, sgn, sgn))
        if k < n_swap:
            segs.append((2 * TAU, sgn, -sgn))
            sgn = -sgn
    segs.append((TAU, sgn, 0.0))
    ts, ws, t0 = [], [], 0.0
    for L, a, b in segs:
        n = max(2, int(round(L / dt)))
        tt = np.linspace(0, L, n, endpoint=False)
        ts.append(t0 + tt); ws.append(a + (b - a) * tt / L); t0 += L
    t = np.concatenate(ts); w = np.concatenate(ws)
    return t, w, t0


def W2_profile(t, w, f=FGRID):
    dt = t[1] - t[0]
    out = np.empty(len(f))
    for i0 in range(0, len(f), 200):
        ff = f[i0:i0 + 200, None]
        out[i0:i0 + 200] = np.abs((w[None, :] * np.exp(-2j * pi * ff * t[None, :])).sum(1) * dt) ** 2
    return out


def var_allowed(Phi, gt, lam_target=L_TARGET):
    if lmin(rho_noisy(Phi, gt, 0, 0, 0)) >= lam_target:
        return 0.0
    return brentq(lambda v: lmin(rho_noisy(Phi, gt, sqrt(v), sqrt(v), 0.0)) - lam_target, 1e-14, 50.0)


def proto(n_swap=0, t_hold=T_INT, kappa=KAPPA, Phi_hold_rate=PHI / T_INT, gt_rate=GT_RATE, lam_target=L_TARGET):
    t, w, ttot = profile(n_swap, t_hold)
    Phi = Phi_hold_rate * t_hold          # ASSUMED conservative: ramps and swaps add no entangling phase
    gt = gt_rate * ttot
    v = var_allowed(Phi, gt, lam_target)
    W2 = W2_profile(t, w)
    intw2 = float(np.sum(w ** 2) * (t[1] - t[0]))
    req_white = sqrt(v / (kappa ** 2 * 0.5 * intw2)) if v > 0 else 0.0
    return dict(n_swap=n_swap, t_hold=t_hold, t_tot=ttot, Phi=Phi, gt=gt, var_allowed=v, W2=W2, intw2=intw2,
                req_white=req_white, kappa=kappa)


def residual_orders(P, asd):
    S = np.array([asd(f) for f in FGRID]) ** 2
    var = P['kappa'] ** 2 * np.trapezoid(S * P['W2'], FGRID)
    if P['var_allowed'] <= 0:
        return np.inf
    return 0.5 * log10(var / P['var_allowed'])


# ---------------------------------------------------------------------------
# spectra (per mass, along the split axis)
# ---------------------------------------------------------------------------
A_ISO = (2 * pi) ** 2 * s18.ISO_DISP_1HZ                       # 3.9e-10, Stage 18
SPECTRA = {
    'seismic_nlnm': (nlnm_asd, 'NLNM quiet site, no isolation (CITED Peterson 1993)'),
    'seismic_iso_flat': (lambda f: A_ISO, 'aLIGO-class isolation, flat in acceleration (Stage 18 shape, ASSUMED)'),
    'seismic_iso_passive': (lambda f: A_ISO * min(1.0, (1.0 / f) ** 2),
                            'aLIGO-class isolation + passive stage: accel falls f^-2 above 1 Hz (ASSUMED)'),
    'tilt': (lambda f: G0 * s18.OMEGA_RL / (2 * pi * f), 'tilt x g, ring-laser bound (CITED, Stage 18)'),
    'efield_white': (lambda f: s18.P_DIP * sqrt(s18.S_E_REF * (s18.Z_REF / s18.Z_C) ** 4) / s18.Z_C / M,
                     'plate E-field noise on dipole, white (Stage 18)'),
    'efield_1overf': (lambda f: s18.P_DIP * sqrt(s18.S_E_REF * (s18.Z_REF / s18.Z_C) ** 4 * 1e6 / f) / s18.Z_C / M,
                      'plate E-field noise on dipole, 1/f (Stage 18)'),
    'johnson_plate': (lambda f: 5.876e-9, 'Au-plate magnetic Johnson noise, white (Stage 18)'),
    # ISS spec SSP 41000E: 1/3-octave rms a <= 16 ug x f (0.1-100 Hz), 1.6 ug (0.01-0.1 Hz)  -> ASD = a/sqrt(0.2316 f)
    'iss_spec': (lambda f: (16e-6 * G0 * f if f > 0.1 else 1.6e-6 * G0) / sqrt(0.2316 * f),
                 'ISS vibratory requirement (CITED SSP 41000E)'),
    # ZARM drop tower / Einstein-Elevator: residual 1e-6 g; ASSUMED spread flat over 0.1-10 Hz -> 3.1e-6 /rtHz
    'drop_tower': (lambda f: 1e-6 * G0 / sqrt(10.0), 'drop tower / Einstein-Elevator residual 1e-6 g (CITED), flat 0.1-10 Hz (ASSUMED)'),
}


def main():
    res = dict(stage=21, baseline=dict(m=M, dx=DX, D=D, Phi=PHI, kappa=KAPPA, N5=N5_BASE, gt_rate=GT_RATE))
    base = proto(0)
    res['baseline']['req_white'] = base['req_white']
    res['baseline']['req_white_stage18_criterion'] = proto(0, lam_target=L_TARGET_S18)['req_white']
    base_res = {k: residual_orders(base, a) for k, (a, _) in SPECTRA.items()}
    res['baseline']['residual_orders'] = base_res

    # ---- Lever 1: echo / branch-swap filter functions ------------------------------------------------
    L1 = []
    for n in (0, 1, 2, 3, 4, 6, 8):
        P = proto(n)
        L1.append(dict(n_swap=n, t_tot=P['t_tot'], gt=P['gt'], req_white=P['req_white'],
                       residual_orders={k: residual_orders(P, a) for k, (a, _) in SPECTRA.items()},
                       W2=P['W2'].tolist()))
    # same total time 3 s (hold shrinks)
    L1b = []
    for n in (1, 2):
        th = s18.T_TOT - 2 * TAU - 2 * TAU * n
        P = proto(n, th)
        L1b.append(dict(n_swap=n, t_hold=th, Phi=P['Phi'], req_white=P['req_white'],
                        residual_orders={k: residual_orders(P, a) for k, (a, _) in SPECTRA.items()}))
    res['lever1_echo'] = dict(same_hold=L1, same_total=L1b, f=FGRID.tolist())

    # ---- Lever 1b: hold time (N5 fixed) ---------------------------------------------------------------
    holds = [0.25, 0.5, 1.0, 2.5, 5.0, 10.0, 20.0]
    L1c = []
    for th in holds:
        for n in (0, 2):
            P = proto(n, th)
            L1c.append(dict(t_hold=th, n_swap=n, req_white=P['req_white'],
                            res_iso=residual_orders(P, SPECTRA['seismic_iso_flat'][0]),
                            res_nlnm=residual_orders(P, SPECTRA['seismic_nlnm'][0]),
                            res_tilt=residual_orders(P, SPECTRA['tilt'][0])))
    res['lever1_hold'] = L1c
    # analytic scaling (DERIVED): sigma^2 ~ Phi ~ G m^2 dx^2 t/(hbar D^3) and sigma^2 = kappa^2 S t/2
    #   => S_req ~ G hbar / D^3, independent of m, dx, t (Grossardt 2020 form 2 G hbar/d^3)
    res['scaling'] = dict(grossardt_sqrt_2Ghbar_D3=sqrt(2 * G * HBAR / D ** 3),
                          ratio_req_to_grossardt=base['req_white'] / sqrt(2 * G * HBAR / D ** 3))

    # ---- Lever 2: shorter shots, more repetitions at fixed total wall-clock time -----------------------
    T_OH = 2.0   # ASSUMED per-shot overhead (load, cool, readout), s
    budget = N5_BASE * (s18.T_TOT + T_OH)
    L2 = []
    for th in (0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0, 20.0):
        ttot = th + 2 * TAU
        N = budget / (ttot + T_OH)
        lam_need = -5.0 / sqrt(N)
        P = proto(0, th, lam_target=lam_need)
        L2.append(dict(t_hold=th, N_shots=N, lam_need=lam_need, req_white=P['req_white'],
                       res_iso=residual_orders(P, SPECTRA['seismic_iso_flat'][0]) if P['var_allowed'] > 0 else None))
    res['lever2_shots'] = dict(T_overhead=T_OH, wallclock_s=budget, rows=L2)

    # ---- Lever 3: remove the plate ---------------------------------------------------------------------
    # unshielded sphere-sphere Casimir-Polder (retarded, D >> R): U = -(23 hbar c / 4 pi) a1 a2 / D^7,
    # a = R^3 (eps-1)/(eps+2)  (CITED Casimir & Polder 1948; same form used by Bose et al. PRL 119, 240401 (2017)).
    # Curvature ratio |U_CP''| / |U_G''| = 23*56 hbar c a^2 / (8 pi G m^2 D^6); a ~ R^3 ~ m => D independent of m.
    eps = 5.7
    a_pol = R ** 3 * (eps - 1) / (eps + 2)
    D_free = (23 * 56 * HBAR * C_LIGHT * a_pol ** 2 / (8 * pi * G * M ** 2 * 0.1)) ** (1 / 6)
    scale = D_free / D
    m_free = M * sqrt(scale)                        # keeps Phi at fixed t and u = (dx/D)^2 (DERIVED)
    dx_free = DX * scale
    kappa_free = m_free * dx_free / HBAR
    # gas: saturated rate ~ R^2 ~ m^(2/3): keep gamma t by lowering P by that factor (cost in vacuum orders; Stage 20: not binding)
    P_factor = (m_free / M) ** (2 / 3)
    Pfree = proto(0, kappa=kappa_free)
    tighten = log10(base['req_white'] / Pfree['req_white'])
    # neutral-dipole engineering instead (keep plate, p -> 0, q = 0, no normal metal):
    # residual = induced dipole in static patch field: F ~ alpha_pol E0 dE / z  (E0 = V_patch/z, V_patch ASSUMED 10 mV)
    alpha_SI = 4 * pi * 8.8541878128e-12 * a_pol
    E0 = 0.010 / s18.Z_C
    dE_white = sqrt(s18.S_E_REF * (s18.Z_REF / s18.Z_C) ** 4)
    a_induced = alpha_SI * E0 * dE_white / s18.Z_C / M
    # trap surfaces when the plate is gone (ground): nearest surface z_trap ASSUMED 100 um
    z_trap = 100e-6
    aE_trap_white = s18.P_DIP * sqrt(s18.S_E_REF * (s18.Z_REF / z_trap) ** 4) / z_trap / m_free
    aJ_trap = 5.876e-9 * (s18.Z_C / z_trap) ** 2 * (m_free / M) ** 0   # thin-film z^-2 (ASSUMED), V/m ratio cancels at fixed B0
    res['lever3_plate'] = dict(D_free=D_free, m_free=m_free, dx_free=dx_free, P_factor_lower=P_factor,
                               req_white_free=Pfree['req_white'], tighten_orders=tighten,
                               neutral_dipole_induced_a=a_induced,
                               neutral_dipole_orders=log10(a_induced / base['req_white']),
                               z_trap=z_trap, aE_trap_white=aE_trap_white,
                               aE_trap_orders=log10(aE_trap_white / Pfree['req_white']),
                               aJ_trap=aJ_trap, aJ_trap_orders=log10(aJ_trap / Pfree['req_white']))

    # ---- Lever 4: trap stiffness (DERIVED: no gain); separate SG-source mount penalty ------------------
    w_trap = 2 * pi * 100.0
    res['lever4_trap'] = dict(gain_orders=0.0, f_trap=100.0,
                              X_rel_allowed=base['req_white'] / w_trap ** 2,
                              note='trap-frame inertial force -m*a is independent of omega; if the Stern-Gerlach '
                                   'source and the trap are on different mounts, their relative motion X couples as '
                                   'm*omega^2*X (enhanced, not suppressed)')

    # ---- Lever 5: free fall -----------------------------------------------------------------------------
    ff = {}
    # (a) trapped particles in a falling capsule: capsule residual couples like ground platform noise
    for key in ('drop_tower', 'iss_spec'):
        ff[key + '_trapped'] = dict(orders=base_res[key])
    # (b) released (free-floating) particles: platform noise decouples; residuals are direct forces.
    # spacecraft self-gravity jitter: 2 G M_sc dr / r^3 (M_sc 500 kg at 0.5 m, dr 1e-9 m/rtHz drag-free jitter: ASSUMED)
    M_sc, r_sc, dr_sc = 500.0, 0.5, 1e-9
    a_sg = 2 * G * M_sc * dr_sc / r_sc ** 3
    # distant electrode E-field on dipole at z_el = 1 mm (ASSUMED), white and 1/f, for the no-plate mass
    z_el = 1e-3
    aE_el_w = s18.P_DIP * sqrt(s18.S_E_REF * (s18.Z_REF / z_el) ** 4) / z_el / m_free
    Pspace = Pfree
    sp_rows = dict(self_gravity=a_sg, efield_1mm_white=aE_el_w)
    sp = {k: log10(v / Pspace['req_white']) for k, v in sp_rows.items()}
    # 1/f case through the filter
    Pspace_dict = dict(Pspace)
    sp['efield_1mm_1overf'] = residual_orders(Pspace_dict, lambda f: s18.P_DIP * sqrt(
        s18.S_E_REF * (s18.Z_REF / z_el) ** 4 * 1e6 / f) / z_el / m_free)
    # shots: N5 at facility rates
    ff['shots'] = dict(N5=N5_BASE, zarm_days=N5_BASE / 3.0, einstein_elevator_days=N5_BASE / 300.0,
                       zarm_src='https://www.zarm.uni-bremen.de/fileadmin/user_upload/Research/Testing_and_Facilities/ZARM_BDT_PUG.pdf',
                       ee_src='https://www.hcmr.uni-hannover.de/en/details/einstein-elevator')
    ff['space_free_floating'] = dict(a=sp_rows, orders=sp, req=Pspace['req_white'], M_sc=M_sc, r_sc=r_sc, dr_sc=dr_sc,
                                     z_el=z_el, dr_allowed=Pspace['req_white'] * r_sc ** 3 / (2 * G * M_sc),
                                     p_needed_ecm_1mm=s18.P_DIP / s18.E_CM * Pspace['req_white'] / aE_el_w)
    # GGN tidal only for a released pair on the ground (drop tower): person at 2 m, Stage 18 a_d -> no-plate D
    a_d2 = 6 * G * 70.0 * 0.3 * D_free / 2.0 ** 4
    ff['drop_tower_released_ggn_2m_orders'] = log10(s18.line_white_equivalent(a_d2, 1.0)[0] / base['req_white']) + tighten
    res['lever5_freefall'] = ff

    # ---- Lever 6: relative-coordinate / one-trap protocols ---------------------------------------------
    # (i) collective dephasing: the only decoherence-free subspace for two qubits is span{LR,RL}; its only product
    # states |LR>,|RL> acquire just a global phase from gravity -> cannot witness gravity-mediated entanglement.
    dfs = rho_noisy(PHI, 0.0, 3.0, 3.0, 1.0)
    psi = np.zeros(4, complex); psi[1] = psi[2] = 1 / sqrt(2)
    rho_dfs = np.outer(psi, psi.conj())
    # common noise on the DFS Bell state (entangled at preparation) leaves it unchanged:
    ph = np.array([0, 1, 1, 2.0]); th = 3.0
    lam_dfs = lmin(rho_dfs)
    # (ii) Gaussian two-mode version: CM mode thermal (n from common noise) x relative mode squeezed (r):
    # particles 1,2 entangled iff smallest PT symplectic eigenvalue < 1/2.  DERIVED numerically.
    def nu_min(n, r):
        Vc = (n + 0.5) * np.eye(2)
        Vr = 0.5 * np.diag([np.exp(-2 * r), np.exp(2 * r)])
        S = np.zeros((4, 4))           # (x1,p1,x2,p2) from (X,P,x,p): x1=(X+x)/sqrt2, x2=(X-x)/sqrt2
        S[0, 0] = S[0, 2] = S[1, 1] = S[1, 3] = 1 / sqrt(2)
        S[2, 0] = S[3, 1] = 1 / sqrt(2); S[2, 2] = S[3, 3] = -1 / sqrt(2)
        V = np.zeros((4, 4)); V[:2, :2] = Vc; V[2:, 2:] = Vr
        V12 = S @ V @ S.T
        Tp = np.diag([1, 1, 1, -1]); Vpt = Tp @ V12 @ Tp
        Om = np.array([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])
        return float(np.min(np.abs(np.linalg.eigvals(1j * Om @ Vpt))))
    gauss = []
    for r in (0.1, 0.5, 1.0):
        nth = (np.exp(2 * r) - 1) / 2       # predicted threshold e^{2r} = 2n+1
        gauss.append(dict(r=r, n_threshold_pred=nth, nu_below=nu_min(0.9 * nth, r), nu_above=nu_min(1.1 * nth, r)))
    res['lever6_relative'] = dict(dfs_bell_lambda=lam_dfs, gaussian=gauss,
                                  stage18_anti_parallel_cmrr=1.006)

    # ---- combined designs ---------------------------------------------------------------------------------
    best_n = min(L1, key=lambda r: max(r['residual_orders']['seismic_iso_flat'], r['residual_orders']['tilt']))
    comb = {}
    def tot(orders):
        return 0.5 * log10(sum(10 ** (2 * o) for o in orders))
    ro = best_n['residual_orders']
    comb['ground_plate_echo'] = dict(n_swap=best_n['n_swap'], inertial=tot([ro['seismic_iso_flat'], ro['tilt']]),
                                     inertial_passive=tot([ro['seismic_iso_passive'], ro['tilt']]),
                                     with_surface=tot([ro['seismic_iso_flat'], ro['tilt'], ro['efield_white'], ro['johnson_plate']]))
    # ground, plate kept but neutral-dipole engineered (p=0, q=0) and non-metal (no Johnson): induced-dipole floor
    ind_echo = res['lever3_plate']['neutral_dipole_orders'] + (ro['efield_white'] - base_res['efield_white'])
    comb['ground_plate_echo_neutral'] = dict(induced=ind_echo,
                                             total=tot([ro['seismic_iso_flat'], ro['tilt'], ind_echo]))
    comb['ground_baseline_same_criterion'] = dict(inertial=tot([base_res['seismic_iso_flat'], base_res['tilt']]),
                                                  with_surface=tot([base_res['seismic_iso_flat'], base_res['tilt'],
                                                                    base_res['efield_white'], base_res['johnson_plate']]))
    # ground, no plate, echo: requirement tightened, plate noise -> trap surface noise at 100 um
    comb['ground_noplate_echo'] = dict(inertial=comb['ground_plate_echo']['inertial'] + tighten,
                                       inertial_passive=comb['ground_plate_echo']['inertial_passive'] + tighten,
                                       with_trap_surface=tot([comb['ground_plate_echo']['inertial'] + tighten,
                                                              res['lever3_plate']['aE_trap_orders'],
                                                              res['lever3_plate']['aJ_trap_orders']]))
    comb['space_trapped_lpf'] = dict(orders=None)
    Plpf = proto(0)
    comb['space_trapped_lpf']['orders'] = residual_orders(Plpf, lambda f: (2 * pi * f) ** 2 * s18.LPF_READ)
    comb['space_free_floating_noplate'] = dict(orders=tot([sp['self_gravity'], sp['efield_1mm_white']]),
                                               orders_1overf=tot([sp['self_gravity'], sp['efield_1mm_1overf']]),
                                               binding=max(sp, key=lambda k: sp[k] if k != 'efield_1mm_1overf' else -99))
    res['combined'] = comb
    res['spectra_labels'] = {k: v[1] for k, v in SPECTRA.items()}
    with open(os.path.join(ROOT, 'results', 'stage21.json'), 'w') as fh:
        json.dump(res, fh, indent=1, default=float)

    # ---- print ----
    print('baseline req %.3f fm (Stage 18 0.442); N5 %.0f; Grossardt sqrt(2Ghbar/D^3)=%.2f fm ratio %.3f' %
          (base['req_white'] * 1e15, N5_BASE, res['scaling']['grossardt_sqrt_2Ghbar_D3'] * 1e15, res['scaling']['ratio_req_to_grossardt']))
    print('baseline residual orders:', {k: round(v, 2) for k, v in base_res.items()})
    for r in L1:
        print('echo N=%d ttot=%.1f req=%.3f fm ' % (r['n_swap'], r['t_tot'], r['req_white'] * 1e15),
              {k: round(v, 2) for k, v in r['residual_orders'].items()})
    for r in L1b:
        print('echo same-total N=%d hold=%.2f Phi=%.3f req=%.3f' % (r['n_swap'], r['t_hold'], r['Phi'], r['req_white'] * 1e15),
              {k: round(v, 2) for k, v in r['residual_orders'].items()})
    for r in L1c:
        print('hold %5.2f N=%d req=%.3f fm iso %+.2f nlnm %+.2f tilt %+.2f' % (r['t_hold'], r['n_swap'], r['req_white'] * 1e15,
              r['res_iso'], r['res_nlnm'], r['res_tilt']))
    for r in L2:
        print('lever2 hold %5.2f N=%.3g lam_need=%.4f req=%.3f fm iso=%s' % (r['t_hold'], r['N_shots'], r['lam_need'],
              r['req_white'] * 1e15, r['res_iso']))
    print('lever3:', {k: (('%.3g' % v) if isinstance(v, float) else v) for k, v in res['lever3_plate'].items()})
    print('lever4:', res['lever4_trap']['X_rel_allowed'])
    print('lever5:', json.dumps(ff, default=float)[:1500])
    print('lever6:', res['lever6_relative'])
    print('combined:', json.dumps(comb, default=float))


if __name__ == '__main__':
    main()
