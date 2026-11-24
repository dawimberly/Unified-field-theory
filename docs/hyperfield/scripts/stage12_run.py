#!/usr/bin/env python3
"""Stage 12: decoherence budget + design optimisation for the BMV/QGEM test.

Reproduces published working points, builds a constraint-aware grid, and
writes results/stage12.json. Short, vectorised-enough nested loops (~few 1e5
point evaluations, a few seconds).
"""
from __future__ import annotations
import sys, json, os
from math import pi, log10, sqrt, inf
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from stage12_lib import (
    G, HBAR, K_B, RHO_DIAMOND, E_CM, T2_ELECTRONIC_S, T2_NUCLEAR_S,
    sphere_radius, phases_linear, phases_linear_centers, phases_parallel,
    r_min_cp_unshielded, concurrence_dephased, dephased_ppt_lmin,
    n_shots_witness_5sigma, n_shots_phase_5sigma,
    Point, evaluate, required_pressure, asdict,
    gamma_air_vdk, gamma_air_geometric, gamma_bb,
    dBdx_for_split, drop_height, wavepacket_sigma, humpty_rel_gradient,
    F_cp_sphere_plate, F_dipole_plate, p_schut_default,
)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
os.chdir(ROOT)
os.makedirs('results', exist_ok=True)

# ---------- literature points (CITED inputs) ----------
R15 = sphere_radius(1e-15)
R14 = sphere_radius(1e-14)
POINTS = [
    Point('Bose_2017', 'linear', False, False,
          m=1e-14, d=200e-6, dx=250e-6, t_int=2.5, tau_acc=0.5,
          T_env=0.15, T_int_temp=0.15, T_gas=300.0, P=1e-15,
          notes='Original protocol, centre-to-centre 450 um, closest 200 um'),
    Point('VDK_2020', 'linear', True, False,
          m=1e-15, d=57 * R15, dx=23e-6, t_int=1.0, tau_acc=0.25,
          T_env=1.0, T_int_temp=1.0, T_gas=4.0, P=5.5e-16,
          W_plate=1e-6, notes='van de Kamp et al. PRA 102, 062807; N=57, 10^4 T/m, Phi~0.01 in 1 s'),
    Point('Schut_2023_ff_14', 'parallel', True, False,
          m=1e-14, d=41e-6, dx=30e-6, t_int=0.5, tau_acc=0.25,
          T_env=1.0, T_int_temp=1.0, T_gas=4.0, P=1e-15,
          z_plate=None, notes='Schut et al. PRR 5, 043170 fig.8-like: m=1e-14, d~41 um, dx~30 um, free-fall + plate'),
    Point('Schut_2024_trap_d61', 'parallel', True, True,
          m=1e-14, d=61e-6, dx=8.5e-6, t_int=1.0, tau_acc=0.1,
          T_env=1.0, T_int_temp=1.0, T_gas=4.0, P=1e-15,
          z_plate=30e-6, notes='Schut et al. PRR 6, 013199 Table I, z=30 um, |phi1|=0.01'),
    Point('Schut_2024_trap_d21', 'parallel', True, True,
          m=1e-14, d=21e-6, dx=1.7e-6, t_int=1.0, tau_acc=0.1,
          T_env=1.0, T_int_temp=1.0, T_gas=4.0, P=1e-15,
          z_plate=10e-6, notes='Table I, z=10 um, |phi1|=0.01'),
    Point('Schut_2024_trap_d11', 'parallel', True, True,
          m=1e-14, d=11e-6, dx=0.65e-6, t_int=1.0, tau_acc=0.1,
          T_env=1.0, T_int_temp=1.0, T_gas=4.0, P=1e-15,
          z_plate=5e-6, notes='Table I, z=5 um, |phi1|=0.01'),
]

# ---------- state of the art (CITED) ----------
SOA = dict(
    m_superposition_kg=172000 * 1.66053906660e-27,  # Pedalino 2025, 172 kDa
    dx_superposition_m=133e-9,
    t_superposition_s=0.01,
    m_cooled_levitated_kg=1e-18,  # Tebbenjohanns 2021; not a cat
    P_typical_UHV_Pa=1e-10,
    P_XHV_claimed_Pa=1e-15,  # cryogenic pumping, claimed possible
    T_cryostat_K=4.0,
    dBdx_routine_T_per_m=1e4,
    dBdx_ambitious_T_per_m=1e6,
    hd_rel_ambitious=1e-6,
    N_shots_two_weeks=10_000,
    drop_tower_m=110.0,  # ZARM Bremen class
    NV_T2e_s=T2_ELECTRONIC_S,
    P_QGEM_design_Pa=1e-15,  # Bose/VDK design point; claimed possible with cryo pumping
)


def orders_over(req, soa, higher_is_harder=True):
    if req is None or soa is None or req <= 0 or soa <= 0 or not np.isfinite(req):
        return None
    if higher_is_harder:
        return log10(req / soa)
    return log10(soa / req)


def difficulty(ev, P_used, soa=SOA, mass_ref=None):
    """Sum of log10 gaps vs SOA. mass_ref=None uses interference mass;
    set mass_ref to cooled-NP mass to ask 'once we can levitate this mass'."""
    mref = soa['m_superposition_kg'] if mass_ref is None else mass_ref
    terms = {
        'mass': orders_over(ev['m'], mref, True),
        'dx': orders_over(ev['dx'], soa['dx_superposition_m'], True),
        'time': orders_over(ev['t_tot'], soa['t_superposition_s'], True),
        'vacuum_vs_UHV': orders_over(P_used, soa['P_typical_UHV_Pa'], False) if P_used and P_used > 0 else None,
        'vacuum_vs_XHV': orders_over(P_used, soa['P_QGEM_design_Pa'], False) if P_used and P_used > 0 else None,
        'gradient': orders_over(ev['dBdx'], soa['dBdx_routine_T_per_m'], True),
        'humpty_gs': orders_over(soa['hd_rel_ambitious'], ev['hd_rel_gs'], True) if ev['hd_rel_gs'] < soa['hd_rel_ambitious'] else 0.0,
        'shots': orders_over(ev['N_shots_5sig_witness'], soa['N_shots_two_weeks'], True)
                 if np.isfinite(ev['N_shots_5sig_witness']) else 3.0,
        'drop': orders_over(ev['drop_height_m'], soa['drop_tower_m'], True)
                if ev['drop_height_m'] > 0 else 0.0,
    }
    def _sum(keys):
        s = 0.0
        n = 0
        for k in keys:
            v = terms.get(k)
            if v is None or not np.isfinite(v):
                continue
            s += max(0.0, v)
            if v > 0.05:
                n += 1
        return s, n
    routine_keys = ('mass', 'dx', 'time', 'vacuum_vs_UHV', 'gradient', 'humpty_gs', 'shots', 'drop')
    qgem_keys = ('mass', 'dx', 'time', 'vacuum_vs_XHV', 'gradient', 'humpty_gs', 'shots', 'drop')
    tot, n = _sum(routine_keys)
    return dict(terms={k: (None if v is None else float(v)) for k, v in terms.items()},
                total=float(tot), n_hard=int(n),
                total_qgem=float(_sum(qgem_keys)[0]))


def literature_block():
    out = {'r_min_unshielded_diamond_0p1_m': float(r_min_cp_unshielded()),
           'r_min_note': 'VDK eq. 5 quotes 157 um; we recompute 156.9 um (DERIVED, match). Bose 200 um is conservative.',
           'points': []}
    for pt in POINTS:
        ev = evaluate(pt)
        ev['difficulty_vs_interference'] = difficulty(ev, pt.P)
        ev['difficulty_vs_cooled_mass'] = difficulty(ev, pt.P, mass_ref=SOA['m_cooled_levitated_kg'])
        ev['notes'] = pt.notes
        # JSON-safe
        clean = {}
        for k, v in ev.items():
            if isinstance(v, (float, int, str, bool)) or v is None:
                clean[k] = v if not isinstance(v, float) else (None if not np.isfinite(v) else float(v))
            elif isinstance(v, dict):
                clean[k] = v
        out['points'].append(clean)
    # Schut 3.4 um claim check
    claim = []
    for d in (61e-6, 41e-6, 33.2e-6, 21e-6, 11e-6):
        ph = phases_parallel(1e-14, 1.0, d, 3.4e-6)
        claim.append(dict(d_um=d * 1e6, phi1=float(ph['phi1']), Phi=float(ph['Phi'])))
    out['schut_3p4um_claim'] = dict(
        rows=claim,
        verdict=('PRR 6, 013199 says "approximately 3.4 um" at m=1e-14, 0.01 rad, 1 s. '
                 'Their Table I is internally consistent with |phi1|=0.01: 8.5, 1.7, 0.65 um at '
                 'd=61, 21, 11 um. |phi1|=0.01 at dx=3.4 um wants d≈33 um, which is not a table row. '
                 'If they meant Phi=phi1+phi2=0.01, then dx=3.4 um wants d≈41 um (and |phi1|=0.005). '
                 'The 3.4 um sentence does not match Table I under a single consistent Phi convention. '
                 'We treat Table I as the source of record and flag the 3.4 um sentence as inconsistent.')
    )
    # VDK Phi check
    ph = phases_linear(1e-15, 1.0, 57 * sphere_radius(1e-15), 23e-6)
    out['vdk_phi_check'] = dict(Phi=float(ph['Phi']), phi1=float(ph['phi1']),
                                expect_step2='~0.01', match=bool(abs(abs(ph['Phi']) - 0.01) < 0.003))
    # air-rate factor
    out['air_rate_vdk_over_geometric'] = 8 * pi / 3
    return out


# ---------- optimisation grid ----------
def needed_P(pt):
    """Choose P so air uses at most half of the remaining PPT budget after BB+GGN.

    Budget for total γ t is ~ |Phi|/2 (the small-phase PPT edge). We keep a
    factor 0.4 of |Phi|/2 for air,  after subtracting γ_bb t and γ_d t.
    Floor P at 1e-18 Pa (beyond any lab).
    """
    # evaluate at tiny P to get non-air contributions, then set P
    pt0 = Point(**{**asdict(pt), 'P': 1e-18})
    ev0 = evaluate(pt0)
    Phi = abs(ev0['Phi'])
    budget = 0.4 * (0.5 * Phi)  # 0.4 of the PPT edge
    already = ev0['gamma_bb'] * ev0['t_tot'] + ev0['gamma_d_t']
    remain = budget - already
    if remain <= 0:
        return 1e-18, ev0, False
    P = required_pressure(pt, remain, air_model='vdk')
    P = float(min(max(P, 1e-18), 1e-6))
    ev = evaluate(Point(**{**asdict(pt), 'P': P}))
    return P, ev, True


CONSTRAINTS_NOTES = []


def feasible(pt, ev, P):
    """Hard cuts for a 'workable' point. Returns (ok, reasons_failed)."""
    fail = []
    if ev['clash']:
        fail.append('geometry_clash')
    if not ev['entangled']:
        fail.append('not_entangled')
    if ev['N_shots_5sig_witness'] > 1e5:
        fail.append('too_many_shots')
    if not ev['electronic_ok']:
        fail.append('electronic_T2')
    if not ev['nuclear_ok']:
        fail.append('nuclear_T2')
    if not pt.trapped and ev['drop_height_m'] > 150.0:
        fail.append('drop_too_tall')
    if ev['dBdx'] > 1e7:
        fail.append('gradient_unreal')
    # EM constraints
    if not pt.shielded:
        if ev['cp_ss_ratio'] > 0.1:
            fail.append('unshielded_casimir')
    else:
        if ev['z_plate'] is None or ev['z_plate'] <= 0:
            fail.append('plate_gap_invalid')
        elif ev['z_plate'] < ev['R'] + 0.5e-6:
            fail.append('too_close_to_plate')
        if not pt.trapped:
            # VDK: a_cp / a_mag ≤ 0.1; we also cap dipole the same way
            if ev['ratio_cp_mag'] > 0.1:
                fail.append('plate_casimir_vs_magnet')
            if ev['ratio_dip_mag'] > 0.1:
                fail.append('plate_dipole_vs_magnet')
        # trapped: trap is assumed to hold if we passed the gap cut
    return (len(fail) == 0), fail


def run_grid():
    """Scan geometries × (m, d, dx, t, tau_acc, dipole scenario). Keep workable + some metadata."""
    geoms = [
        dict(geometry='linear', shielded=False, trapped=False, tau_accs=(0.25, 0.5), p=p_schut_default()),
        dict(geometry='linear', shielded=True, trapped=False, tau_accs=(0.25, 0.5), p=p_schut_default()),
        dict(geometry='parallel', shielded=False, trapped=False, tau_accs=(0.25, 0.5), p=p_schut_default()),
        dict(geometry='parallel', shielded=True, trapped=False, tau_accs=(0.25, 0.5), p=p_schut_default()),
        dict(geometry='parallel', shielded=True, trapped=True, tau_accs=(0.1, 0.25), p=p_schut_default()),
        # same trapped geometry with 100× better dipole control
        dict(geometry='parallel', shielded=True, trapped=True, tau_accs=(0.1, 0.25), p=1e-4 * E_CM),
    ]
    masses = np.logspace(-16, -13, 13)
    dxs = np.logspace(-7, -3.3, 15)      # 0.1 um .. 0.5 mm
    ds = np.logspace(np.log10(11e-6), np.log10(5e-4), 14)  # >= Schut closest d=11 um
    times = (0.5, 1.0, 2.5, 5.0)
    T_env = T_int = 1.0
    T_gas = 4.0

    workable = []
    n_tried = 0
    fail_counts = {}
    best = None
    best_lev = None
    best_qgem = None
    best_product = None

    for g in geoms:
        for m in masses:
            R = sphere_radius(m)
            for d in ds:
                for dx in dxs:
                    # crude geometry pre-cuts
                    if g['geometry'] == 'linear' and d <= 2 * R:
                        continue
                    if g['geometry'] == 'parallel' and d <= 2 * R + (1e-6 if g['shielded'] else 0):
                        continue
                    z_guess = 0.5 * d - 0.5e-6
                    if g['shielded'] and z_guess < 5e-6:
                        # Schut 2024 treats z=5 um as the closest studied point (needs more trap analysis)
                        continue
                    if g['shielded'] and z_guess <= R + 0.5e-6:
                        continue
                    for t in times:
                        for ta in g['tau_accs']:
                            n_tried += 1
                            z_plate = 0.5 * d - 0.5e-6 if g['shielded'] else None
                            pt = Point(
                                name='grid', geometry=g['geometry'],
                                shielded=g['shielded'], trapped=g['trapped'],
                                m=float(m), d=float(d), dx=float(dx),
                                t_int=float(t), tau_acc=float(ta),
                                T_env=T_env, T_int_temp=T_int, T_gas=T_gas,
                                P=1e-15, p_dipole=g['p'],
                                z_plate=z_plate, W_plate=1e-6,
                            )
                            P, ev, okP = needed_P(pt)
                            ok, fail = feasible(pt, ev, P)
                            if not ok:
                                for f in fail:
                                    fail_counts[f] = fail_counts.get(f, 0) + 1
                                continue
                            ev = evaluate(Point(**{**asdict(pt), 'P': P}))
                            dif = difficulty(ev, P)
                            dif_lev = difficulty(ev, P, mass_ref=SOA['m_cooled_levitated_kg'])
                            rec = dict(
                                geometry=g['geometry'], shielded=g['shielded'], trapped=g['trapped'],
                                p_dipole_ecm=g['p'] / E_CM,
                                m=float(m), R=ev['R'], d=float(d), dx=float(dx),
                                t_int=float(t), tau_acc=float(ta), t_tot=ev['t_tot'],
                                P=float(P), Phi=ev['Phi'], phi1=ev['phi1'],
                                concurrence=ev['concurrence'], lambda_min=ev['lambda_min'],
                                N_shots=ev['N_shots_5sig_witness'], N_shots_phase=ev['N_shots_5sig_phase'],
                                dBdx=ev['dBdx'], drop_m=ev['drop_height_m'],
                                gamma_t=ev['gamma_t'], gamma_air=ev['gamma_air'], gamma_bb=ev['gamma_bb'],
                                hd_rel_gs=ev['hd_rel_gs'], hd_rel_1mK=ev['hd_rel_1mK'],
                                ratio_cp_mag=ev['ratio_cp_mag'], ratio_dip_mag=ev['ratio_dip_mag'],
                                cp_ss_ratio=ev['cp_ss_ratio'], z_plate=ev['z_plate'],
                                difficulty=dif['total'], difficulty_terms=dif['terms'],
                                difficulty_lev=dif_lev['total'], difficulty_lev_terms=dif_lev['terms'],
                                # QGEM-design difficulty: vacuum measured against 1e-15 Pa, not routine UHV
                                difficulty_qgem=float(dif.get('total_qgem', dif['total'])),
                                product_mdt=float(m * dx * ev['t_tot']),
                            )
                            workable.append(rec)
                            if best is None or rec['difficulty'] < best['difficulty']:
                                best = rec
                            if best_lev is None or rec['difficulty_lev'] < best_lev['difficulty_lev']:
                                best_lev = rec
                            if best_qgem is None or rec['difficulty_qgem'] < best_qgem['difficulty_qgem']:
                                best_qgem = rec
                            if best_product is None or rec['product_mdt'] < best_product['product_mdt']:
                                best_product = rec
    return dict(n_tried=n_tried, n_workable=len(workable), fail_counts=fail_counts,
                best=best, best_levitated_conditional=best_lev,
                best_qgem=best_qgem, best_product=best_product, points=workable)


def pareto(points, keys=('m', 'dx', 't_tot', 'difficulty')):
    """Minimal Pareto front: keep points not dominated in all of keys (lower is better)."""
    arr = [{k: p[k] for k in list(keys) + ['geometry', 'shielded', 'trapped', 'P', 'Phi',
                                           'concurrence', 'N_shots', 'dBdx', 'd', 't_int',
                                           'p_dipole_ecm', 'hd_rel_gs', 'drop_m']}
           for p in points]
    front = []
    for i, a in enumerate(arr):
        dominated = False
        for j, b in enumerate(arr):
            if i == j:
                continue
            if all(b[k] <= a[k] * (1 + 1e-12) for k in keys) and any(b[k] < a[k] for k in keys):
                dominated = True
                break
        if not dominated:
            front.append(a)
    front.sort(key=lambda r: r['difficulty'])
    return front


def mass_time_trade():
    """At fixed |Phi|=0.05 (Schut 10k-shot 5σ phase) and fixed geometries, t vs m."""
    target = 0.05
    rows = []
    for label, geom, d, dx, shielded, trapped in [
        ('linear_noshield_Bose_like', 'linear', 200e-6, 250e-6, False, False),
        ('linear_shield_VDK_like', 'linear', 57 * sphere_radius(1e-15), 23e-6, True, False),
        ('parallel_trap_d21', 'parallel', 21e-6, 1.7e-6, True, True),
        ('parallel_trap_d61', 'parallel', 61e-6, 8.5e-6, True, True),
    ]:
        for m in np.logspace(-16, -13, 25):
            # scale t so |Phi|=target (DERIVED, gravity only, no decoherence)
            if geom == 'linear':
                ph1 = phases_linear(m, 1.0, d if d > 2 * sphere_radius(m) else d, dx)
            else:
                ph1 = phases_parallel(m, 1.0, d, dx)
            if not ph1['ok'] or abs(ph1['Phi']) < 1e-30:
                continue
            t = target / abs(ph1['Phi'])
            rows.append(dict(label=label, geometry=geom, shielded=shielded, trapped=trapped,
                             m=float(m), d=float(d), dx=float(dx), t_for_Phi_0p05=float(t),
                             Phi_per_second=float(ph1['Phi'])))
    return rows


def pedernales_note(points):
    """Does a smaller superposition help? Only if d can shrink with it (shield+trap)."""
    # compare, at m=1e-14, t=1, the |Phi| vs dx at fixed d vs at fixed d/dx
    m, t = 1e-14, 1.0
    rows = []
    for d in (200e-6, 61e-6, 21e-6, 11e-6):
        for dx in np.logspace(-7, -4, 12):
            ph = phases_parallel(m, t, d, dx)
            rows.append(dict(d=float(d), dx=float(dx), Phi=float(ph['Phi']),
                             phi1=float(ph['phi1'])))
    return dict(
        statement=('Pedernales & Plenio, arXiv:2311.04745: the entangling rate is set by spatial '
                   'delocalisation, not by the state being a cat rather than a squeezed Gaussian '
                   'of the same size. Shrinking Δx at fixed d therefore always hurts (Phi ~ Δx^2/d^3 '
                   'in parallel). Shrinking Δx *and* d, which the plate+trap geometry allows, can '
                   'keep Phi fixed: a 10× smaller cat at a 10^{2/3}≈4.6× smaller d has the same phase. '
                   'That is the real win of Schut 2024, not a loophole around delocalisation.'),
        grid=rows,
    )


def binding_constraint(best):
    """Which difficulty term is largest at the cheapest point?"""
    if best is None:
        return None
    terms = {k: (v if v is not None else -inf) for k, v in best['difficulty_terms'].items()}
    k = max(terms, key=terms.get)
    terms_lev = {k2: (v if v is not None else -inf) for k2, v in best['difficulty_lev_terms'].items()}
    klev = max(terms_lev, key=terms_lev.get)
    return dict(vs_interference=k, orders=terms[k],
                vs_cooled_mass=klev, orders_lev=terms_lev[klev],
                all=best['difficulty_terms'], all_lev=best['difficulty_lev_terms'])


def main():
    lit = literature_block()
    print('literature points:')
    for p in lit['points']:
        print(f"  {p['name']}: Phi={p['Phi']:.4f} C={p['concurrence']:.3f} "
              f"gt={p['gamma_t']:.3f} N={p['N_shots_5sig_witness']} "
              f"dBdx={p['dBdx']:.3e} drop={p['drop_height_m']:.1f} "
              f"cp={p['cp_ss_ratio']:.3f} rCPmag={p['ratio_cp_mag']:.3f} "
              f"rDIPmag={p['ratio_dip_mag']:.3f} hd={p['hd_rel_gs']:.2e} "
              f"dif={p['difficulty_vs_interference']['total']:.1f}")
    print('r_min', lit['r_min_unshielded_diamond_0p1_m'])
    print('scanning grid...')
    grid = run_grid()
    print('tried', grid['n_tried'], 'workable', grid['n_workable'], 'fails', grid['fail_counts'])
    if grid['best']:
        b = grid['best']
        print('best:', {k: b[k] for k in
              ['geometry', 'shielded', 'trapped', 'm', 'd', 'dx', 't_int', 'tau_acc',
               'P', 'Phi', 'concurrence', 'N_shots', 'dBdx', 'drop_m', 'difficulty',
               'difficulty_lev', 'p_dipole_ecm', 'hd_rel_gs']})
    # don't dump every workable point into json (can be thousands); store a
    # thinned copy + the Pareto front
    pts = grid['points']
    front = pareto(pts) if pts else []
    # thin: keep best 200 by difficulty and a random-ish stride of the rest
    pts_sorted = sorted(pts, key=lambda r: r['difficulty'])
    thin = pts_sorted[:200]
    if len(pts_sorted) > 200:
        thin += pts_sorted[200::max(1, len(pts_sorted) // 200)]
    trade = mass_time_trade()
    ped = pedernales_note(pts)
    bind = binding_constraint(grid['best']) if grid['best'] else None
    bind_lev = binding_constraint(grid['best_levitated_conditional']) if grid['best_levitated_conditional'] else None

    out = dict(
        soa=SOA,
        literature=lit,
        grid=dict(n_tried=grid['n_tried'], n_workable=grid['n_workable'],
                  fail_counts=grid['fail_counts'],
                  T_env_K=1.0, T_int_K=1.0, T_gas_K=4.0,
                  note='P is set per point so air uses <= 0.4 of the PPT budget |Phi|/2 after BB+GGN.'),
        cheapest=grid['best'],
        cheapest_levitated_conditional=grid['best_levitated_conditional'],
        cheapest_qgem_design=grid['best_qgem'],
        cheapest_product_mdt=grid['best_product'],
        binding=bind,
        binding_levitated_conditional=bind_lev,
        pareto=front,
        workable_thin=thin,
        mass_time_trade=trade,
        pedernales=ped,
        errors_found=[
            'Bose et al. 2017: written witness W=|⟨XZ⟩-⟨YZ⟩| and the W>1 test do not hold (Stage 11). Later QGEM papers switched to the PPT witness, which is the correct tool.',
            'van de Kamp et al. 2020: unshielded CP bound of 157 um is correct (we recompute 156.9 um). Bose 200 um is a conservative loose bound, not an error.',
            'van de Kamp P~5e-16 is 5.5e-16 Pa at T=4 K and n=1e7 m^{-3} (unit omitted in the text but implied Pa).',
            'Schut et al. PRR 6, 013199 (2024): the conclusion sentence "approximately 3.4 um" at m=1e-14, 0.01 rad, 1 s does not match Table I under one Phi convention. Table I is consistent with |phi1|=0.01 (8.5/1.7/0.65 um at d=61/21/11 um). 3.4 um gives |phi1|=0.01 only at d≈33 um, or Phi=phi1+phi2=0.01 at d≈41 um. Worth a clarification; we treat Table I as authoritative.',
            'Phi convention splits the literature: Bose/VDK quote Phi=phi1+phi2; the 2024 trap paper\'s |ω_ent|=0.01 Hz is |phi1|/t = |Phi|/(2t) in parallel. A factor-of-two trap for comparisons.',
            'The QGEM air-molecule rate (VDK B6) is 8π/3 ≈ 8.38 times the geometric πR^2 ⟨v⟩ collision rate. Using the geometric rate understates the vacuum requirement by that factor.',
            'Stage 11 prose said Casimir-safe r_min scales as m^{1/3}. That is wrong: at fixed density r_min is mass-independent (the Stage 11 *calculation* already gave a constant 157 um; only the sentence was wrong).',
        ],
        worth_sending_to_authors=(
            'The 3.4 um vs Table I mismatch in Schut et al. PRR 6, 013199 is the one concrete clarification worth sending, if Dan wants. '
            'The Bose 2017 witness wording is already superseded by later papers from the same groups (Chevalier/Tilly/Schut) and may not need a letter. '
            'Do not contact anyone unless Dan says so.'
        ),
    )
    # make json safe
    def sanitise(x):
        if isinstance(x, dict):
            return {k: sanitise(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [sanitise(v) for v in x]
        if isinstance(x, (np.floating, float)):
            return None if not np.isfinite(x) else float(x)
        if isinstance(x, (np.integer,)):
            return int(x)
        if isinstance(x, (np.bool_, bool)):
            return bool(x)
        return x
    json.dump(sanitise(out), open('results/stage12.json', 'w'), indent=1)
    print('wrote results/stage12.json, pareto', len(front), 'thin', len(thin))


if __name__ == '__main__':
    main()
