#!/usr/bin/env python3
"""Exact Stage 18/21 noise model applied to the Pb Meissner (Stage 16) and Pino Nb (Stage 19/20) points.

Uses stage21_run.proto / residual_orders unchanged, with each point's own kappa = m dx / hbar, Phi, gas+BB
exponent rate and hold from stage16_run.evaluate, shot cap N5 <= 1e4 (Stage 21 rule), ramps tau = 0.25 s.
Only platform (inertial) spectra are scored: plate E-field / Johnson entries in Stage 18/21 use the diamond's
permanent dipole and Au plate and do not transfer to superconducting spheres (flagged, not computed).
Writes results/summary21_pb_pino.json.
"""
import json, os, sys
from math import sqrt, log10, pi
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import stage16_run as s16
import stage21_run as s21
from stage12_lib import G, HBAR

C = {c['key']: c for c in s16.CANDIDATES}
pts = dict(cheap_check=('diamond_spin', 5.6234e-15, 3.852e-6, 2.5),
           pb=('lead_meissner', 1e-14, 2.6915348e-6, 1.0),
           pino=('niobium_meissner', 2e13 * s16.U_DA, 1.2302687708e-06, 0.5))
tot = lambda o: 0.5 * log10(sum(10 ** (2 * x) for x in o))
out = {}
for k, (ck, m, dx, t) in pts.items():
    r = s16.evaluate(m, dx, t, C[ck])
    kappa = m * dx / HBAR
    gt_rate = r['gamma_t'] / r['t_tot']
    row = dict(m=m, dx=dx, t_hold=t, D=r['D'], Phi=r['Phi'], stage16_rule_req=r['sqrt_Sa_req'],
               ceiling_sqrt_2Ghbar_D3=sqrt(2 * G * HBAR / r['D'] ** 3))
    for n in (0, 1):
        P = s21.proto(n, t, kappa=kappa, Phi_hold_rate=r['Phi'] / t, gt_rate=gt_rate)
        res = {s: s21.residual_orders(P, a) for s, (a, _) in s21.SPECTRA.items()
               if s in ('seismic_iso_flat', 'seismic_nlnm', 'tilt', 'seismic_iso_passive')}
        row[f'swap{n}'] = dict(req_white=P['req_white'], residual=res,
                               inertial_total=tot([res['seismic_iso_flat'], res['tilt']]),
                               ratio_to_ceiling=P['req_white'] / row['ceiling_sqrt_2Ghbar_D3'])
    row['best_ground_inertial'] = min(row['swap0']['inertial_total'], row['swap1']['inertial_total'])
    out[k] = row
json.dump(out, open(os.path.join(ROOT, 'results/summary21_pb_pino.json'), 'w'), indent=1, default=float)
for k, r in out.items():
    print(k, 'D=%.3g um ceil=%.3g fm s16=%.3g fm' % (r['D'] * 1e6, r['ceiling_sqrt_2Ghbar_D3'] * 1e15, r['stage16_rule_req'] * 1e15))
    for n in (0, 1):
        s = r[f'swap{n}']
        print('   swap%d req=%.3g fm ratio=%.3f inertial=%.2f' % (n, s['req_white'] * 1e15, s['ratio_to_ceiling'], s['inertial_total']),
              {a: round(b, 2) for a, b in s['residual'].items()})
