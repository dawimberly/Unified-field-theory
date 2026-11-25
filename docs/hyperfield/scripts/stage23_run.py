"""Stage 23 - firm up the Stage 21 space design (released, drag-free, no plate) with MEASURED inputs.

Stage 21 space design: diamond, D = 273 um (Casimir-Polder 10% cut), m = 2.8e-14 kg, dx = 95.8 um, 3 s shot,
N5 <= 1e4, per-mass requirement 1.03e-17 m s^-2/rtHz.  Residual +2.0 orders came from ASSUMED inputs:
p = 1e-2 e cm, ion-trap-extrapolated E-field noise at 1 mm, spacecraft jitter 1 nm/rtHz at 500 kg/0.5 m.
Here every input is replaced by a measured value (CITED, URL) or flagged ASSUMED.  Same exact model and filter
as Stage 21 (stage21_run.proto / residual_orders).
"""
from __future__ import annotations
import json, os, sys
from math import pi, sqrt, log10
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import stage21_run as s21                                                     # noqa: E402
import stage18_run as s18                                                     # noqa: E402
from stage12_lib import G, HBAR, K_B                                          # noqa: E402

E = 1.602176634e-19
EUM = E * 1e-6               # 1 e um in C m  (= 1e-4 e cm)
EPS0 = 8.8541878128e-12
HBARC = HBAR * 2.99792458e8
M0, DX0, D0 = s18.M, s18.DX, s18.D

# ---------------------------------------------------------------------------
# MEASURED inputs
# ---------------------------------------------------------------------------
SRC = dict(
    rider2016='Rider et al., PRL 117, 101101 (2016): silica r=2.5 um, permanent dipoles 151+-6, 89+-10, 192+-30 e um; '
              'UV discharge with charge quantisation. https://web.stanford.edu/group/grattalab/Publications/PhysRevLett.117.101101.pdf',
    afek2021='Afek et al., PRA 104, 053512 (2021), arXiv:2108.04406: 15 um silica p0=(-200+-360) e um, '
             'alpha=401+-8 e um/(V/mm); precessive control cancels dipole backgrounds by x120. https://arxiv.org/abs/2108.04406',
    priel2022='Priel et al., Sci. Adv. (2022), PMC9565793: O(100 e um) dipole signals seen in all levitated-MS techniques '
              '(7.5 um silica). https://pmc.ncbi.nlm.nih.gov/articles/PMC9565793/',
    diamond='No measured permanent dipole for levitated (micro)diamond found; rotating-diamond experiments report no '
            'observable dipole coupling (arXiv:2309.01545) -> silica surface density used (ASSUMED).',
    lpf_patch='Armano et al., PRL 118, 171101 (2017): stray-potential fluctuations 34+-2 uV/rtHz at 1 mHz, 190+-30 at 0.1 mHz; '
              'static stray potentials up to 100 mV in lab; charging lambda_eff1+lambda_eff2 = 2220+-260 /s; '
              'UV discharge keeps V_tm within +-80 mV. https://eprints.gla.ac.uk/141213/1/141213.pdf',
    lpf_patch_shape='Armano et al. arXiv:2404.08340: Delta_x PSD ~ 1/f (ASD ~ f^-1/2); 2-sigma limit on an extra 1/f ASD '
                    'term 65 uV/rtHz x (0.1 mHz/f). https://arxiv.org/abs/2404.08340',
    lpf_jitter='Armano et al., PRL 116, 231101 (2016): spacecraft motion bulge 20-200 mHz, max jitter 8-30 nm/rtHz on the '
               'different axes; stiffness omega_2^2 = (-525+-30)e-9 s^-2, gravitation-dominated, agrees with self-gravity '
               'model to 10%. https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.116.231101',
    lpf_release='LPF GPRM: release requirement <5 um/s; in flight typically ~10 um/s linear (Dalla Ricca et al. ESMATS 2021). '
                'https://www.esmats.eu/esmatspapers/pastpapers/pdfs/2021/dalla-ricca.pdf ; '
                'https://sci.esa.int/web/lisa-pathfinder/-/56868-lisa-pathfinder-test-mass',
    gravitower='Prakash et al., npj Microgravity (2026): first release-recapture of optically levitated silica nanoparticles '
               'in microgravity (GraviTower Bremen); usable free flights 5-6 us, 9 us not separable. '
               'https://doi.org/10.1038/s41526-026-00596-y',
    goce='GOCE accelerometers 2e-12 m s^-2/rtHz spec (5-100 mHz), in-flight 3.1-6.7e-12. '
         'https://sf2a.eu/proceedings/2010/2010sf2a.conf..0113C.pdf',
)
P_SILICA_R25 = 151.0 * EUM            # median of Rider 2016 (CITED)
R_SILICA = 2.5e-6
ALPHA_AFEK = 401.0 * EUM / 1e3        # C m per (V/m) for 15 um silica (CITED) ; R = 7.5 um
R_AFEK = 7.5e-6
SPIN_CANCEL = 120.0                   # CITED Afek 2021
DELTA_1MHZ = 34e-6                    # V/rtHz at 1 mHz (CITED)
DELTA_EXTRA_01MHZ = 65e-6             # V/rtHz at 0.1 mHz, extra 1/f ASD 2-sigma limit (CITED)
V_STATIC = 0.100                      # V, static stray potential (CITED upper level)
LAMBDA_LPF = 2220.0 / 2               # /s per TM (CITED)
A_LPF = 6 * 0.046 ** 2                # 46 mm cube area, m^2 (CITED LPF TM size)
OMEGA2_SC = 525e-9                    # s^-2 (CITED)
JITTER = (8e-9, 30e-9)                # m/rtHz (CITED, 20-200 mHz max)
JITTER_1MHZ = 0.2e-9                  # m/rtHz simulated at 1 mHz (CITED Vitale & Weber, arXiv:1012.5968)
V_RELEASE_LPF = 10e-6                 # m/s typical in flight (CITED)
T_FREE_DEMO = 6e-6                    # s usable free flight demonstrated (CITED)
GOCE = 3.1e-12

MATS = dict(
    diamond=dict(rho=3510.0, cm=(5.7 - 1) / (5.7 + 2), coef=23 / (4 * pi), conductor=False),
    silica=dict(rho=2200.0, cm=(3.8 - 1) / (3.8 + 2), coef=23 / (4 * pi), conductor=False),
    osmium=dict(rho=22590.0, cm=1.0, coef=143 / (16 * pi), conductor=True),
)
# Casimir-Polder, retarded, D >> R (CITED Casimir & Polder 1948; dielectric 23/(4pi); perfect conductor e+m 143/(16 pi)).


def design(mat):
    """D from the 10% curvature cut (Stage 16/21 rule): U = -coef hbar c (cm R^3)^2 / D^7, U_G = -G m^2/D.
    D is independent of m (DERIVED); m then chosen to keep Phi at the Stage 21 value (m ~ sqrt(D))."""
    p = MATS[mat]
    # ratio = coef*56 hbar c cm^2 R^6 / (2 G m^2 D^6) with m = 4/3 pi rho R^3  -> R cancels
    k = p['coef'] * 56 * HBARC * p['cm'] ** 2 / (2 * G * (4 / 3 * pi * p['rho']) ** 2)
    D = (k / 0.1) ** (1 / 6)
    m = M0 * sqrt(D / D0)
    R = (3 * m / (4 * pi * p['rho'])) ** (1 / 3)
    dx = DX0 * D / D0
    kappa = m * dx / HBAR
    P = s21.proto(0, kappa=kappa)
    return dict(mat=mat, D=D, m=m, R=R, dx=dx, kappa=kappa, P=P, req=P['req_white'])


def jitter_asd(f, J):
    """Spacecraft-TM jitter: flat J in 20-200 mHz (CITED max), f^-2 above 0.2 Hz (ASSUMED, force-noise limited),
    log-linear down to 0.2 nm at 1 mHz below 20 mHz (simulated value CITED)."""
    if f >= 0.2:
        return J * (0.2 / f) ** 2
    if f >= 0.02:
        return J
    s = np.log(J / JITTER_1MHZ) / np.log(0.02 / 1e-3)
    return JITTER_1MHZ * (f / 1e-3) ** s


def delta_asd(f, extra=False):
    a = DELTA_1MHZ * sqrt(1e-3 / f)
    if extra:
        a = sqrt(a ** 2 + (DELTA_EXTRA_01MHZ * 1e-4 / f) ** 2)
    return a


def dipole(mat, R, law):
    if MATS[mat]['conductor']:
        # metal sphere: no bulk dipole; patch-potential dipole p ~ 4 pi eps0 R^2 V_patch, V_patch = 10 mV (ASSUMED)
        return 4 * pi * EPS0 * R ** 2 * 0.010
    return P_SILICA_R25 * (R / R_SILICA) ** law


def evaluate(d, z_el, law=2, spin=False, J=JITTER[0], extra=False, q_e=0):
    m, R, P = d['m'], d['R'], d['P']
    p = dipole(d['mat'], R, law) / (SPIN_CANCEL if spin else 1.0)
    alpha = ALPHA_AFEK * (R / R_AFEK) ** 3
    if MATS[d['mat']]['conductor']:
        alpha = 4 * pi * EPS0 * R ** 3
    E0 = V_STATIC / z_el
    terms = {}
    terms['permanent_dipole'] = s21.residual_orders(P, lambda f: p * delta_asd(f, extra) / z_el ** 2 / m)
    terms['induced_dipole'] = s21.residual_orders(P, lambda f: alpha * E0 * delta_asd(f, extra) / z_el ** 2 / m)
    if q_e:
        terms['charge'] = s21.residual_orders(P, lambda f: q_e * E * delta_asd(f, extra) / z_el / m)
    terms['sc_jitter'] = s21.residual_orders(P, lambda f: OMEGA2_SC * jitter_asd(f, J))
    tot = 0.5 * log10(sum(10 ** (2 * v) for v in terms.values()))
    return dict(p_eum=p / EUM, terms=terms, total=tot, binding=max(terms, key=terms.get))


def main():
    res = dict(stage=23, sources=SRC)
    designs = {k: design(k) for k in MATS}
    res['designs'] = {k: dict(D=v['D'], m=v['m'], R=v['R'], dx=v['dx'], req=v['req'],
                              req_vs_diamond_orders=log10(v['req'] / designs['diamond']['req'])) for k, v in designs.items()}
    dd = designs['diamond']
    # ---- (1) dipole: measured vs required -----------------------------------------------------------
    p2 = dipole('diamond', dd['R'], 2); p3 = dipole('diamond', dd['R'], 3)
    # required p at z_el for permanent-dipole term alone to reach 0 orders (DERIVED: residual ~ linear in p)
    req_p = {}
    for z in (1e-3, 4e-3, 1e-2):
        o = s21.residual_orders(dd['P'], lambda f: EUM * delta_asd(f) / z ** 2 / dd['m'])
        req_p[str(z)] = EUM * 10 ** (-o) / EUM
    res['dipole'] = dict(silica_meas_R25_eum=[151, 89, 192], diamond_scaled_R2_eum=p2 / EUM, diamond_scaled_R3_eum=p3 / EUM,
                         stage21_assumed_eum=1e-2 * 1e4 * 1.0, required_eum=req_p,
                         stage21_required_eum_1mm=1.0, spin_cancel=SPIN_CANCEL,
                         afek_alpha_scaled=ALPHA_AFEK * (dd['R'] / R_AFEK) ** 3,
                         cm_alpha=4 * pi * EPS0 * dd['R'] ** 3 * MATS['diamond']['cm'])
    # ---- (2) E-field noise: measured spectrum at the filter, vs Stage 18/21 ion-trap extrapolation ------
    dE_1mm_1Hz_lpf = delta_asd(1.0) / 1e-3
    dE_1mm_iontrap = sqrt(s18.S_E_REF * (s18.Z_REF / 1e-3) ** 4)
    res['efield'] = dict(delta_1Hz_uV=delta_asd(1.0) * 1e6, delta_extra_1Hz_uV=delta_asd(1.0, True) * 1e6,
                         dE_1mm_1Hz_lpf=dE_1mm_1Hz_lpf, dE_1mm_iontrap_white=dE_1mm_iontrap,
                         lpf_over_iontrap_orders=log10(dE_1mm_1Hz_lpf / dE_1mm_iontrap),
                         note='LPF Delta_x includes actuation-electronics noise, so it is an upper bound for a passive '
                              'housing; it is the only measured space value; extrapolated from 0.1-20 mHz to ~0.3 Hz (ASSUMED shape)')
    # ---- (3) jitter ---------------------------------------------------------------------------------------
    jit = {}
    for J in JITTER:
        o = s21.residual_orders(dd['P'], lambda f: OMEGA2_SC * jitter_asd(f, J))
        jit['%.0fnm' % (J * 1e9)] = o
    # required flat jitter (in-band) for 0 orders, DERIVED (residual linear in J)
    J_req = JITTER[0] * 10 ** (-jit['8nm'])
    stage21_model_stiff = 2 * G * 500.0 / 0.5 ** 3
    res['jitter'] = dict(orders=jit, J_required=J_req, J_measured=JITTER, omega2_measured=OMEGA2_SC,
                         omega2_stage21_model=stage21_model_stiff, stage21_assumed_J=1e-9, stage21_required_J=1.9e-11,
                         omega2_needed_at_8nm=OMEGA2_SC * 10 ** (-jit['8nm']))
    # ---- (4) release / D stability ------------------------------------------------------------------------
    t = s18.T_TOT
    dD_lpf = V_RELEASE_LPF * t
    cas_factor = (dd['D'] / (dd['D'] - dD_lpf)) ** 6
    v_req = 0.01 * dd['D'] / t
    T_eff_req = dd['m'] * v_req ** 2 / K_B
    dPhi = 3 * dD_lpf / dd['D']
    res['release'] = dict(v_lpf=V_RELEASE_LPF, dD_lpf=dD_lpf, dD_over_D=dD_lpf / dd['D'], casimir_ratio_at_closest=0.1 * cas_factor,
                          dPhi_over_Phi=dPhi, v_required_1pct=v_req, T_eff_required=T_eff_req,
                          t_free_demo=T_FREE_DEMO, t_free_needed=t, free_flight_gap_orders=log10(t / T_FREE_DEMO),
                          dc_selfgravity_drift=0.5 * OMEGA2_SC * dd['D'] * t ** 2)
    # ---- (5) charge ---------------------------------------------------------------------------------------
    rate = LAMBDA_LPF * (4 * pi * dd['R'] ** 2) / A_LPF
    o_q1 = s21.residual_orders(dd['P'], lambda f: E * delta_asd(f) / 1e-3 / dd['m'])
    res['charge'] = dict(rate_per_s=rate, mean_days_between=1 / rate / 86400, prob_per_shot=rate * t,
                         one_electron_1mm_orders=o_q1,
                         note='area scaling of the LPF cosmic-ray charging rate is ASSUMED; discharge to exactly 0 e by UV '
                              'with charge quantisation is demonstrated in levitated microspheres (Rider 2016)')
    # ---- scenarios ------------------------------------------------------------------------------------------
    sc = {}
    sc['diamond_measured_1mm'] = evaluate(dd, 1e-3, law=2)
    sc['diamond_measured_1mm_R3'] = evaluate(dd, 1e-3, law=3)
    sc['diamond_measured_1mm_extra1f'] = evaluate(dd, 1e-3, law=2, extra=True)
    sc['diamond_1cm'] = evaluate(dd, 1e-2, law=2)
    sc['diamond_1cm_spin'] = evaluate(dd, 1e-2, law=2, spin=True)
    sc['diamond_1cm_spin_30nm'] = evaluate(dd, 1e-2, law=2, spin=True, J=JITTER[1])
    sc['diamond_1mm_one_electron'] = evaluate(dd, 1e-3, law=2, q_e=1)
    for mat in ('silica', 'osmium'):
        d = designs[mat]
        sc[mat + '_1mm'] = evaluate(d, 1e-3, law=2)
        sc[mat + '_1cm_spin' if mat != 'osmium' else mat + '_1cm'] = evaluate(d, 1e-2, law=2, spin=(mat != 'osmium'))
    sc['osmium_1cm_spin'] = evaluate(designs['osmium'], 1e-2, law=2, spin=True)   # spin-averaging of a patch dipole: ASSUMED applicable
    # electrode distance needed for the permanent-dipole term to reach 0 orders (residual ~ z^-2, DERIVED)
    res['z_needed_m'] = {k: 1e-2 * 10 ** (sc[k]['terms']['permanent_dipole'] / 2) for k in ('diamond_1cm_spin', 'diamond_1cm', 'osmium_1cm_spin')}
    # fraction of the dipole-term variance from the measured LPF band (1-20 mHz), DERIVED
    f = s21.FGRID; W2 = dd['P']['W2']; S = np.array([delta_asd(x) ** 2 for x in f])
    tot = np.trapezoid(S * W2, f); m = f <= 0.02
    res['measured_band_fraction'] = float(np.trapezoid((S * W2)[m], f[m]) / tot)
    res['scenarios'] = sc
    best = min(sc, key=lambda k: sc[k]['total'] if 'electron' not in k and 'extra' not in k else 99)
    res['best'] = best
    res['comparators'] = dict(goce_inflight=GOCE, goce_orders=log10(GOCE / dd['req']),
                              lpf_1p74=1.74e-15, lpf_orders=log10(1.74e-15 / dd['req']))
    with open(os.path.join(ROOT, 'results', 'stage23.json'), 'w') as fh:
        json.dump(res, fh, indent=1, default=float)
    # ---- print ----
    for k, v in res['designs'].items():
        print('%-8s D=%.0f um m=%.3g R=%.3g req=%.3g (%+.2f vs diamond)' % (k, v['D'] * 1e6, v['m'], v['R'], v['req'], v['req_vs_diamond_orders']))
    print('dipole:', {k: v for k, v in res['dipole'].items()})
    print('efield:', {k: v for k, v in res['efield'].items() if k != 'note'})
    print('jitter:', res['jitter'])
    print('release:', res['release'])
    print('charge:', {k: v for k, v in res['charge'].items() if k != 'note'})
    for k, v in sc.items():
        print('%-30s p=%.3g e um total %+.2f binding %s  %s' % (k, v['p_eum'], v['total'], v['binding'],
              {a: round(b, 2) for a, b in v['terms'].items()}))
    print('z_needed', res['z_needed_m'], 'measured-band fraction', res['measured_band_fraction'])
    print('best:', best, 'comparators', res['comparators'])


if __name__ == '__main__':
    main()
