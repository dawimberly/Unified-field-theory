#!/usr/bin/env python3
"""
Stage 19 — Creating and reading out a spatial superposition of a massive object WITHOUT an embedded spin.

Question: Stage 16's best heavier test mass (Pb Meissner sphere, 1e-14 kg, 15 mK) needs dx = 2.69 um,
i.e. 4.65e5 zero-point widths (2 z_zpf at 100 Hz), with no demonstrated which-path tag. Which spin-free
routes exist, what delocalisation can each give (demonstrated vs proposed), and what gap do they leave
inside the unchanged Stage 16 budget (stage16_run.evaluate / min_dx_point)?

Metric (same as Stage 16): F = dx / (2 z_zpf), z_zpf = sqrt(hbar / (2 m w_ref)), w_ref = 2 pi 100 Hz.
ASSUMED: a coherent Gaussian of position spread sigma entangles like a two-branch cat of dx = 2 sigma
(O(1) factor, ~0.3 orders uncertainty). For transient expansions (free expansion, loop protocols) the
time-averaged dx^2 is ASSUMED to be dx_max^2/3 (DERIVED for linear growth), i.e. a sqrt(3) penalty.
Writes results/stage19.json.
"""
from __future__ import annotations
import json, os, sys
from math import pi, sqrt, log10
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import stage16_run as s16  # noqa: E402

HBAR, G, K_B, U = s16.HBAR, s16.G, s16.K_B, s16.U_DA
W_REF = 2 * pi * s16.F_TRAP
CAND = {c['key']: c for c in s16.CANDIDATES}
TRANSIENT = sqrt(3.0)  # DERIVED for linear growth of dx(t): <dx^2> = dx_max^2 / 3


def zref(m):
    return sqrt(HBAR / (2 * m * W_REF))


def req(m, t, key):
    """Smallest healthy dx at mass m, hold t (Stage 16 rules, unchanged)."""
    r = s16.min_dx_point(m, t, CAND[key])
    if r is None:
        return None
    return dict(m=m, t=t, dx=r['dx'], F=r['dx'] / (2 * zref(m)), N5=r['N5'], sqrt_Sa_fm=r['sqrt_Sa_req'] * 1e15,
                P_req_Pa=r['P_req'], Q_min_100Hz=r['Q_min_100Hz'])


def hold_for(t):
    """Stage 16 hold grid; take the largest grid hold not above t (conservative)."""
    grid = [0.5, 1.0, 2.5, 5.0]
    ok = [g for g in grid if g <= t + 1e-12]
    return ok[-1] if ok else None


def main():
    out = dict(stage=19, w_ref_Hz=s16.F_TRAP, transient_penalty=TRANSIENT, routes=[], checks={})

    # ---------------- requirement curves (DERIVED, Stage 16 code) ----------------
    lead_pt = req(1e-14, 1.0, 'lead_meissner')
    curve = []
    for m in np.geomspace(1e-14, 1e-6, 17):
        for t in (0.5, 1.0, 2.5, 5.0):
            r = req(m, t, 'lead_meissner')
            if r:
                curve.append(r)
    Fmin_1s = min((c for c in curve if c['t'] == 1.0), key=lambda c: c['F'])
    Fmin_all = min(curve, key=lambda c: c['F'])
    base = s16.evaluate(5.6234e-15, 3.852e-6, 2.5, CAND['diamond_spin'])
    cheap = dict(m=5.6234e-15, dx=3.852e-6, t=2.5, F=3.852e-6 / (2 * zref(5.6234e-15)), N5=base['N5'],
                 sqrt_Sa_fm=base['sqrt_Sa_req'] * 1e15)
    out['requirements'] = dict(stage16_lead_1e14kg_1s=lead_pt, lead_min_F_at_1s=Fmin_1s, lead_min_F_any_hold=Fmin_all,
                               stage12_cheap_point=cheap, lead_curve=curve)

    R = out['routes']

    def add(**k):
        R.append(k)
        return k

    # ---------------- (a) superconducting / SQUID / flux-qubit routes ----------------
    # Romero-Isart et al. PRL 109, 147205 (2012), arXiv:1112.5609: Pb R = 2 um, w_t = 2pi 28 kHz,
    # g0 = 2pi 1.3 kHz (flux qubit), T1 = T2 = 10 us; l_s = 8 x_zp g0 / w_t (their Eq. 4 text).
    m_ri = 11360 * 4 / 3 * pi * (2e-6) ** 3
    wt = 2 * pi * 28e3
    xzp_ri = sqrt(HBAR / (2 * m_ri * wt))
    ls = 8 * xzp_ri * (2 * pi * 1.3e3) / wt
    tstar = pi / wt
    r_ri = req(m_ri, 0.5, 'lead_meissner')
    add(key='romero_isart_2012', cls='a', status='proposed', name='Flux-qubit magnetomechanics (Romero-Isart 2012)',
        tag='yes (flux qubit, external to the sphere)', m=m_ri, dx=ls, F=ls / (2 * zref(m_ri)), t_coh=None,
        readout='flux qubit state (which-path tag) / SQUID',
        gap_dx_orders=log10(r_ri['dx'] / ls), req=r_ri,
        notes=f'DERIVED l_s/x_zp = 8 g0/w_t = {8*1.3/28:.3f}; l_s = {ls:.2e} m. 2t*/T2 = {2*tstar/10e-6:.1f} with the paper\'s own '
              f'T2 = 10 us (condition 2t* << T2 violated; needs ~1 ms qubits). Superposition << zero-point width; '
              f'any BMV use needs amplification (route b).',
        src='https://arxiv.org/abs/1112.5609')

    # Pino et al. QST 3, 025001 (2018), arXiv:1603.01553: R = 1 um, ~2e13 amu, slit d = 0.5 um, sigma_d = 11.6 nm,
    # free expansion t3 = 483 ms to coherence length xi = 607 nm; split -> measure takes t4..t7 = 96 ms (cat lifetime
    # after split ~ t5+t6+t7 = 83.5 ms); 50 mK, 1e-17 mbar. Readout: SQUID position measurement of fringes.
    m_p = 2e13 * U
    r_p05 = req(m_p, 0.5, 'niobium_meissner')
    xi = 607e-9
    gap_exp = log10(r_p05['dx'] * TRANSIENT / xi)
    # cat d = 0.5 um held 83.5 ms: Phi ~ dx^2 t at fixed D -> dx_req(t) = dx_req(0.5) sqrt(0.5/t)
    dx_req_cat = r_p05['dx'] * sqrt(0.5 / 0.0835)
    add(key='pino_2018', cls='a', status='proposed', name='On-chip SC double slit, magnetic skatepark (Pino 2018)',
        tag='no (position fringes via SQUID)', m=m_p, dx=xi, F=xi / (2 * zref(m_p)), t_coh=0.483,
        readout='SQUID position measurement of interference fringes',
        gap_dx_orders=gap_exp, gap_cat_orders=log10(dx_req_cat / 0.5e-6), req=r_p05,
        notes='Expansion state (xi = 607 nm over 483 ms) vs Stage 16 requirement at the same mass and a 0.5 s hold, '
              'with the sqrt(3) transient penalty. The d = 0.5 um cat exists only ~84 ms. Pino propose a DP-collapse '
              'test of ONE sphere; a two-sphere gravity witness from fringe correlations is NOT in the paper (ASSUMED possible).',
        src='https://arxiv.org/abs/1603.01553')

    # Hexapole / double-well Meissner traps, arXiv:2607.03622, QST (2026) doi:10.1088/2058-9565/aeae6f:
    # 1 um particle, 100 Hz, 10 mK, D_DWP/z_zpf < 100, Q >= 1e11 (current 2e7), ppb field stability.
    m_h = 11340 * 4 / 3 * pi * (1e-6) ** 3
    F_h = 100 / 2.0
    r_h = req(m_h, 0.5, 'lead_meissner')
    add(key='hexapole_2026', cls='a', status='proposed', name='Magnetic higher-order (hexapole) double-well traps (2026)',
        tag='no (position statistics)', m=m_h, dx=F_h * 2 * zref(m_h), F=F_h, t_coh=None,
        readout='pulsed cavity-magnetomechanical position measurement',
        gap_dx_orders=log10(r_h['F'] / F_h), req=r_h,
        gap_Q_orders=log10(1e11 / 2.6e7),
        notes='CITED: nonlinearity ~100 z_zpf, Q >= 1e11 needed vs 2e7 demonstrated, ppb field stability. A 1-um particle is ASSUMED Pb here.',
        src='https://arxiv.org/abs/2607.03622')

    add(key='hofer_2023', cls='a', status='demonstrated (classical only)', name='Meissner-levitated PbSn sphere (Hofer 2023)',
        tag='no', m=5.6e-9, dx=0.0, F=0.0, t_coh=None, readout='SQUID (thermal motion)', gap_dx_orders=None,
        notes='CITED: 100 um, 5.6 ug, Q up to 2.6e7, 15 mK. Thermal state, no ground state, no superposition.',
        src='https://doi.org/10.1103/PhysRevLett.131.043603')
    add(key='chalmers_2026', cls='a', status='demonstrated (components)', name='Flux-tunable resonator + chip levitation (Chalmers thesis 2026)',
        tag='no', m=None, dx=0.0, F=0.0, t_coh=None, readout='dispersive SQUID resonator, 20 GHz/Phi0',
        gap_dx_orders=None,
        notes='CITED: chip levitation at 40 mK, 20 GHz/Phi0 responsivity, 1.6% flux transfer. Readout/cooling hardware, no quantum state yet.',
        src='https://doi.org/10.63959/chalmers.dt/5848')

    # ---------------- (b) potential-shaping / loop protocols ----------------
    for key, name, F, status, note, src in [
        ('rossi_2025', 'Stiffness-modulation delocalisation (Rossi 2025)', 3.5, 'demonstrated',
         'CITED: xi = 73 pm = 3.5 x 21 pm initial, 1.2 fg silica, optical (not SC). The only demonstrated coherent CoM expansion.',
         'https://arxiv.org/abs/2408.01264'),
        ('weiss_2021', 'Bang-bang inverted-potential loop (Weiss 2021)', 1e3, 'proposed',
         'CITED: eta = 1e3 example; gravity coupling g/w0 ~ 1e-23 for 100 nm silica, "a formidable task".',
         'https://arxiv.org/abs/2012.12260'),
        ('roda_llordes_2024_XL', 'Wide double-well loop, size XL (Roda-Llordes 2024)', 1e4 / sqrt(2), 'proposed',
         'CITED Table 1: sqrt(2) eta = 1e4, w/W = 1e-3, D/X_W = 1e6, 77 dB squeezing.', 'https://arxiv.org/abs/2303.07959'),
        ('roda_llordes_2024_XXL', 'Wide double-well loop, size XXL (Roda-Llordes 2024)', 1e5 / sqrt(2), 'proposed',
         'CITED Table 1: sqrt(2) eta = 1e5, w/W = 1e-4, D/X_W = 1e8, 97 dB squeezing. At W = 2pi 100 Hz: w = 2pi 0.01 Hz '
         '(DERIVED), loop time ~ 1/w ~ 16 s >> 1 s hold, and D = 1e8 X_W = 0.29 mm at 1e-14 kg.', 'https://arxiv.org/abs/2303.07959'),
    ]:
        gap_lead = log10(lead_pt['F'] * TRANSIENT / F)
        gap_min = log10(Fmin_1s['F'] * TRANSIENT / F)
        gap_cheap = log10(cheap['F'] * TRANSIENT / F)
        add(key=key, cls='b', status=status, name=name, tag='no', F=F, readout='position / quadratures after recompression',
            gap_dx_orders=gap_lead, gap_dx_best_mass_1s=gap_min, gap_dx_stage12_cheap=gap_cheap, notes=note, src=src,
            m=None, dx=None, t_coh=None)

    # ---------------- (c) gravity-only / position-only witnesses ----------------
    # Lami, Pedernales, Plenio PRX 14, 021022 (2024): coherent states only; gamma = (8 pi/3) G rho / w_I;
    # Table 1: m = 1e-4 kg, R = 2 mm, w_I = 7e-3, gamma = 1.6e-4 Hz, t = 1e2 s, T/Q < 1e-13 K.
    rho_lpp = 1e-4 / (4 / 3 * pi * (2e-3) ** 3)
    gam_lpp = 8 * pi / 3 * G * rho_lpp / 7e-3
    gam_pb100 = 8 * pi / 3 * G * 11340 / W_REF
    t_pb100 = 0.005 / gam_pb100
    TQ_hofer = 0.015 / 2.6e7
    TQ_shin = 300.0 / 1.4e7
    N_runs = 25 * 0.99 * 0.01 / 0.01 ** 2  # 5 sigma on a 0.01 fidelity drop, binomial (DERIVED)
    add(key='lami_pedernales_plenio_2024', cls='c', status='proposed', name='LOCC-inequality test with coherent states (LPP 2024)',
        tag='none needed', F=None, m=1e-4, dx=None, t_coh=100.0, readout='displacement + ground-state (zero-phonon) detection',
        gap_dx_orders=0.0,
        gap_TQ_vs_hofer_orders=log10(TQ_hofer / 1e-13), gap_TQ_vs_torsion_orders=log10(TQ_shin / 1e-13),
        gap_ground_state_torsion_orders=log10(6000 / 1.0),
        checks=dict(rho_implied=rho_lpp, gamma_recomputed=gam_lpp, gamma_paper=1.6e-4,
                    gamma_Pb_100Hz=gam_pb100, t_needed_Pb_100Hz_s=t_pb100, runs_5sigma=N_runs,
                    total_time_days=N_runs * 100 / 86400),
        notes='No delocalisation, but needs w ~ mHz, T/Q < 1e-13 K and ground state of a mHz mode. In the Stage 16 '
              f'100 Hz Pb trap the coupling is {gam_pb100:.2e} s^-1, so a 1% fidelity drop needs {t_pb100:.1e} s (DERIVED): does not rescue it.',
        src='https://arxiv.org/abs/2302.03075')
    add(key='krisnanda_2020', cls='c', status='proposed', name='CV entanglement of two released Os spheres (Krisnanda 2020)',
        tag='none (Gaussian, quadratures)', m=7.9e-8, dx=33e-12, F=None, t_coh=1.0, readout='position/momentum quadratures',
        gap_dx_orders=None, notes='Stage 16 (c\'\'): works only in space at 2.7 K; joint gap 15.6, needs sqrt(S_a) = 0.0087 fm/s^2/rtHz.',
        src='https://arxiv.org/abs/1906.08808')
    add(key='hanif_2024', cls='c', status='proposed', name='No-disturbance test (Hanif, Das, Ulbricht, Bose 2024)',
        tag='REQUIRES embedded spins (Stern-Gerlach)', m=None, dx=None, F=None, t_coh=None, readout='spin',
        gap_dx_orders=None, notes='Honest negative: needs no entanglement but still needs spin-tagged superpositions. Not a spin-free route.',
        src='https://arxiv.org/abs/2307.08133')
    add(key='grossardt_2020', cls='c', status='analysis', name='Acceleration-noise constraint (Grossardt 2020)',
        tag='n/a', m=None, dx=None, F=None, t_coh=None, readout='n/a', gap_dx_orders=None,
        notes='S_a << 2 G hbar / d^3, mass-independent; reproduced in Stage 16. Applies to every route above.',
        src='https://arxiv.org/abs/2007.05762')

    # ---------------- (d) other current ingredients ----------------
    add(key='kotler_2021', cls='d', status='demonstrated', name='Deterministic CV entanglement of two 70 pg drums (Kotler 2021)',
        tag='none (quadrature tomography)', m=7e-14, dx=None, F=None, t_coh=None, readout='near-quantum-limited quadratures',
        gap_dx_orders=None, notes='Shows a spin-free entanglement READOUT exists for massive oscillators. Coupling is electromechanical; dx ~ zero-point.',
        src='https://www.science.org/doi/10.1126/science.abf2998')
    add(key='youssefi_2023', cls='d', status='demonstrated', name='Squeezed drum with 7.7 ms quantum lifetime (Youssefi 2023)',
        tag='none', m=2.3e-12, dx=1.42e-15, F=None, t_coh=7.7e-3, readout='microwave quadratures', gap_dx_orders=None,
        notes='Coherence ingredient only; x_zpf 1.4 fm (Stage 17).', src='https://arxiv.org/abs/2208.13082')
    add(key='torsion_cooling', cls='d', status='demonstrated', name='Torsion-pendulum cooling (Shin 2025; Agafonova 2024)',
        tag='none', m=1e-6, dx=None, F=None, t_coh=None, readout='optical lever',
        gap_dx_orders=None, notes='CITED: 6000 phonons (10 mK, Q 1.4e7, cm-scale); 1 mg at 18 Hz to 240 uK (2.8e5 phonons). No torsion ground state yet.',
        src='https://arxiv.org/abs/2409.02275 ; https://arxiv.org/abs/2408.09445')

    # ---------------- combinations (DERIVED) ----------------
    ri = next(r for r in R if r['key'] == 'romero_isart_2012')
    xxl = next(r for r in R if r['key'] == 'roda_llordes_2024_XXL')
    amp_needed = r_ri['F'] * TRANSIENT / ri['F']
    out['combinations'] = dict(
        tag_then_inflate=dict(desc='flux-qubit tag (RI 2012) amplified by an XXL inverted/double-well loop; both branches '
                                   'and widths scale by the same cosh factor (DERIVED for an inverted harmonic potential)',
                              amplification_needed=amp_needed, amplification_XXL=xxl['F'],
                              gap_orders=log10(amp_needed / xxl['F'])),
        inflate_only_CV=dict(desc='XXL loop + CV witness (no tag) at the Stage 16 Pb point',
                             gap_orders=xxl['gap_dx_orders'], gap_best_mass_1s=xxl['gap_dx_best_mass_1s'],
                             best_mass=Fmin_1s['m'], sqrt_Sa_needed_fm=Fmin_1s['sqrt_Sa_fm']),
    )

    # demonstrated-only gap (what exists today)
    out['demonstrated_gap'] = dict(
        lead_1e14=log10(lead_pt['F'] / 3.5), cheap_point=log10(cheap['F'] / 3.5),
        note='vs Rossi 2025 (3.5, optical silica). For superconductors no coherent CoM state has been shown at all.')

    # ranking: spin-free routes that give a number at the lead point / own mass
    LPF_FM, GROUND_FM = 1.74, 1.5e4   # Stage 16 constants (LPF 2018; Fuchs 2024 thermal force noise)
    for r in R:
        q = r.get('req')
        if q:
            r['gap_noise_vs_LPF_orders'] = log10(LPF_FM / q['sqrt_Sa_fm'])
            r['gap_noise_vs_ground_orders'] = log10(GROUND_FM / q['sqrt_Sa_fm'])
    for r in R:
        if r['cls'] == 'b':
            r['gap_noise_vs_LPF_lead'] = log10(LPF_FM / lead_pt['sqrt_Sa_fm'])
            r['gap_noise_vs_LPF_best_mass_1s'] = log10(LPF_FM / Fmin_1s['sqrt_Sa_fm'])
    pino = next(r for r in R if r['key'] == 'pino_2018')
    pino['gap_vacuum_orders'] = log10(1e-15 / pino['req']['P_req_Pa'])  # Pino assume 1e-17 mbar = 1e-15 Pa
    cand = [r for r in R if r.get('gap_dx_orders') is not None and r['status'] == 'proposed' and r['cls'] in 'ab']
    cand.sort(key=lambda r: r['gap_dx_orders'])
    out['ranking_proposed'] = [(r['key'], r['gap_dx_orders']) for r in cand]

    out['checks'] = dict(ri_ls_over_xzp=8 * 1.3 / 28, ri_mass_kg=m_ri, ri_2tstar_over_T2=2 * tstar / 10e-6,
                         lpp_gamma_ratio_recomputed_over_paper=gam_lpp / 1.6e-4,
                         pino_mass_kg=m_p, pino_F=xi / (2 * zref(m_p)), lead_zpf_m=zref(1e-14),
                         free_expansion_F_formula='F = sqrt(w_ref * W) * t (ground state of trap W released, long t; DERIVED)',
                         free_expansion_W_needed_for_lead_1s_Hz=(lead_pt['F'] * TRANSIENT) ** 2 / W_REF / (2 * pi))

    with open(os.path.join(ROOT, 'results/stage19.json'), 'w') as f:
        json.dump(out, f, indent=1, default=lambda o: float(o) if isinstance(o, (np.floating,)) else str(o))
    # print summary
    print('lead point', {k: lead_pt[k] for k in ('dx', 'F', 'sqrt_Sa_fm')})
    print('min F at 1 s', {k: Fmin_1s[k] for k in ('m', 'dx', 'F', 'sqrt_Sa_fm')}, 'any', {k: Fmin_all[k] for k in ('m', 't', 'F', 'sqrt_Sa_fm')})
    print('cheap', cheap)
    for r in R:
        extra = {k: r[k] for k in r if k.startswith('gap_')}
        print(f"{r['key']:30s} {r['status']:28s} F={r.get('F')} ", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in extra.items()})
    print('combinations', json.dumps(out['combinations'], indent=1, default=str))
    print('demonstrated gap', out['demonstrated_gap'])
    print('checks', json.dumps(out['checks'], indent=1, default=str))
    print('LPP checks', R[[r['key'] for r in R].index('lami_pedernales_plenio_2024')]['checks'])
    print('ranking', out['ranking_proposed'])


if __name__ == '__main__':
    main()
