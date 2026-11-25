#!/usr/bin/env python3
"""Stage 24 - robustness of the Stage 22 non-local window / closing mass to the conventions Eq. 47 leaves free.

Published Eq. 47 (Oppenheim et al. Nat Commun 14, 7910 (2023)):
    sigma_a^2 N r^3 dT / G^2  >=  l_P^2 D2  >=  N_lam M_lam^2 / (R_lam lam),   R_lam = (r^2 dx)^(1/3)
Conventions swept (one at a time from the Stage 15/22 baseline, and as physically paired sets):
  smearing length l (= r in both sides when 'consistent'; floor only when 'floor-only'),
  superposition size dx in {d, d/2, xi d (TL path separation at G2)},
  constituent counting N in {atoms, nucleons, rigid body (N=1, r = body radius)}, plus the independent-constituent
  Diosi-Penrose sum M^2/N (DERIVED alternative, flagged).
Regularisation when l > dx (ASSUMED): smeared-DP small-separation expansion, rate ~ (dx/l)^2 / l, i.e.
  floor = coeff * tau * dx^2 / l^3 (continuous with the R=(l^2 dx)^(1/3) form only up to the 1/3-power; flagged).
Ceiling: FOM = S_a N (Janse), N counted the same way as the floor; Gisler 2022 sensor (9.3e-15 kg SiN, N_atoms 2.79e11).
Closing mass: Stage 22 cheapest MUSCLE geometry, L = 0.983 m, d = 133 nm, xi = 0.8 (tau = 2 xi d^2 M / h).
"""
import json, os
from math import log10, pi, sqrt
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
H = 6.62607015e-34; HBAR = H / (2 * pi); C = 299792458.0; G = 6.6743e-11; U = 1.66053906660e-27
M_NA = 22.98976928 * U; RS_NA = 2.08e-10
S22 = json.load(open(os.path.join(ROOT, 'results/stage22.json')))
# floor experiment (Pedalino 2025) - CITED via Stage 15
PED = dict(M=172e3 * U, tau=2 * 0.983 / 160, d=133e-9, v=160.0, L=0.983)
PED['xi'] = PED['L'] * H / (PED['M'] * PED['v'] * PED['d'] ** 2)
# ceiling experiment (Gisler 2022) - CITED via Stage 15
GIS = dict(M=9.3e-15, N_atoms=2.79e11, Sa=0.298 / 2.79e11, rho=3100.0)   # Sa from CITED FOM 0.298 / N   # rho ASSUMED (Si3N4 ~3.1 g/cc)
GIS['m_avg_u'] = GIS['M'] / GIS['N_atoms'] / U
GIS['r_body'] = (3 * GIS['M'] / (4 * pi * GIS['rho'])) ** (1 / 3)        # equivalent sphere (string is not a sphere: flagged)
CHEAP = dict(L=0.983, d=133e-9, xi=0.8)
SAG_MAX = 0.0317   # largest sag among Stage 22 designs that passed all constraints (xi=1.5); reach criterion ASSUMED

def radius_cluster(M): return RS_NA * (M / M_NA) ** (1 / 3)

def floor_val(M, tau, dx, ell, Nmode):
    a = radius_cluster(M)
    if Nmode == 'atoms': coeff = (M / M_NA) * M ** 2
    elif Nmode == 'nucleons': coeff = (M / U) * M ** 2
    elif Nmode == 'indepDP': coeff = M ** 2 / (M / M_NA)
    elif Nmode == 'v1_noN': coeff = M ** 2
    elif Nmode == 'rigid': coeff = M ** 2; ell = a if ell == 'body' else ell
    if ell == 'body': ell = a
    if Nmode == 'v1_noN': return coeff * tau / 1e-9            # arXiv v1 stated R ~ 1e-9 m
    if ell >= dx: return coeff * tau * dx ** 2 / ell ** 3
    return coeff * tau / (ell ** 2 * dx) ** (1 / 3)

def ceiling_val(ell, Nmode):
    if Nmode in ('atoms', 'indepDP', 'v1_noN'): N = GIS['N_atoms']
    elif Nmode == 'nucleons': N = GIS['M'] / U
    elif Nmode == 'rigid': N = 1.0
    r = GIS['r_body'] if ell == 'body' else ell
    return GIS['Sa'] * N * r ** 3 / G ** 2

CEIL_BASE = ceiling_val(1e-15, 'atoms')

def dx_of(mode, d, xi):
    return dict(d=d, d_half=d / 2, xi_d=xi * d)[mode]

def evaluate(ell_floor, ell_ceil, dxmode, Nmode):
    fl = floor_val(PED['M'], PED['tau'], dx_of(dxmode, PED['d'], PED['xi']), ell_floor, Nmode)
    ce = ceiling_val(ell_ceil, Nmode)
    width = log10(ce / fl)
    # closing mass at the cheapest MUSCLE point (xi fixed => dx_of uses xi=0.8)
    def g(lm):
        M = 10 ** lm * U; tau = 2 * CHEAP['xi'] * CHEAP['d'] ** 2 * M / H
        return log10(floor_val(M, tau, dx_of(dxmode, CHEAP['d'], CHEAP['xi']), ell_floor, Nmode)) - log10(ce)
    lo, hi = 0.0, 16.0
    if g(lo) > 0: Mc = 10 ** lo
    elif g(hi) < 0: Mc = float('inf')
    else:
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if g(mid) < 0 else (lo, mid)
        Mc = 10 ** hi
    v = CHEAP['L'] * H / (CHEAP['xi'] * CHEAP['d'] ** 2 * Mc * U) if np.isfinite(Mc) else float('nan')
    sag = 0.5 * 9.81 * (CHEAP['xi'] * CHEAP['d'] ** 2 * Mc * U / H) ** 2 if np.isfinite(Mc) else float('inf')
    reach = bool(np.isfinite(Mc) and sag <= SAG_MAX)
    return dict(floor=fl, ceiling=ce, sag_m=sag, reachable_in_stage22_constraints=reach, window_orders=width, status='open' if width > 0 else 'closed',
                closing_mass_amu=Mc, closing_mass_MDa=Mc / 1e6, closing_v_m_s=v,
                ceiling_shift_orders=log10(ce / CEIL_BASE))

ELLS = [  # (label, value, basis)
    ('cluster Compton hbar/(Mc), 1.4 MDa', HBAR / (1.415e6 * U * C), 'DERIVED'),
    ('Na atom Compton hbar/(m_Na c)', HBAR / (M_NA * C), 'DERIVED'),
    ('nucleon Compton hbar/(m_p c)', HBAR / (1.67262192e-27 * C), 'DERIVED'),
    ('nuclear r_N = 1e-15 m (Oppenheim 2023; Stage 15)', 1e-15, 'CITED'),
    ('Oppenheim R_lam ~ V^(1/3) ~ 1e-12 m used as l', 1e-12, 'CITED (published value of R, used here as l)'),
    ('Donadi 2020 DP lower limit R0 > 0.54 A', 0.54e-10, 'CITED'),
    ('atomic radius 1e-10 m', 1e-10, 'ASSUMED class value'),
    ('cluster radius (own body: 4.1 nm Pedalino, 8.2 nm at 1.4 MDa)', 'body', 'DERIVED (Wigner-Seitz)'),
    ('CSL/DP-style r_C = 1e-7 m', 1e-7, 'CITED (Bassi et al. RMP 85, 471)'),
]
rows = []
base = evaluate(1e-15, 1e-15, 'd', 'atoms')
rows.append(dict(group='baseline', label='Stage 15/22: l=1e-15, dx=d, N=atoms', **base))
for lab, ell, basis in ELLS:
    if ell == 1e-15: continue
    for mode in ('floor-only', 'consistent'):
        r = evaluate(ell, 1e-15 if mode == 'floor-only' else ell, 'd', 'atoms')
        rows.append(dict(group='smearing l', label=f'{lab} [{mode}]', l=ell if ell != 'body' else 'body', basis=basis, mode=mode, **r))
for dm in ('d_half', 'xi_d'):
    for ell in (1e-15, 1e-7):
        r = evaluate(ell, 1e-15, dm, 'atoms')
        rows.append(dict(group='dx', label=f'dx = {dm}, l = {ell:g} [floor-only]', **r))
for nm in ('nucleons', 'rigid', 'indepDP', 'v1_noN'):
    ell = 'body' if nm == 'rigid' else 1e-15
    r = evaluate(ell, ell, 'd', nm)
    rows.append(dict(group='N counting', label=f'N = {nm} (l = {"body radius" if ell=="body" else "1e-15"}, consistent)', **r))
# physically paired sets
pairs = [('Oppenheim-literal: nucleons + nuclear r_N (their N~1e3 for a 602 u molecule)', 1e-15, 'nucleons'),
         ('atoms + atomic radius 1e-10 m', 1e-10, 'atoms'),
         ('atoms + Donadi DP limit 0.54 A', 0.54e-10, 'atoms'),
         ('rigid body + body radius', 'body', 'rigid'),
         ('independent-constituent DP sum + atomic radius', 1e-10, 'indepDP')]
for lab, ell, nm in pairs:
    r = evaluate(ell, ell, 'd', nm)
    rows.append(dict(group='paired', label=lab, **r))

# full grid for the figure: l from 1e-23 to 1e-6, floor-only and consistent, three dx modes, atoms
grid = {}
ls = np.logspace(-23, -6, 120)
for mode in ('floor-only', 'consistent'):
    for dm in ('d', 'd_half', 'xi_d'):
        w = []; mc = []
        for l in ls:
            r = evaluate(l, 1e-15 if mode == 'floor-only' else l, dm, 'atoms')
            w.append(r['window_orders']); mc.append(r['closing_mass_MDa'])
        grid[f'{mode}|{dm}'] = dict(l=ls.tolist(), window=w, closing_MDa=mc)

def crossing(key, field, target):
    gg = grid[key]; x = np.log10(gg['l']); y = np.array(gg[field], dtype=float)
    for i in range(len(x) - 1):
        if np.isfinite(y[i]) and np.isfinite(y[i + 1]) and (y[i] - target) * (y[i + 1] - target) <= 0:
            return 10 ** (x[i] + (target - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]))
    return None
M_reach = 1.415 * sqrt(SAG_MAX / 0.012359)
thresholds = {f'{m}': dict(l_window_closes=crossing(f'{m}|d', 'window', 0.0), l_closing_mass_leaves_reach=crossing(f'{m}|d', 'closing_MDa', M_reach))
              for m in ('floor-only', 'consistent')}
thresholds['M_reach_MDa'] = M_reach
# robustness summary (DERIVED)
def spread(group):
    ws = [r['window_orders'] for r in rows if r['group'] == group]
    return (min(ws), max(ws))
summary = dict(
    exponent_window_vs_l_consistent='window ~ l^(3+2/3) for l < dx: 3.67 orders per decade of l',
    exponent_window_vs_l_floor_only='window ~ l^(2/3) for l < dx: 0.67 orders per decade',
    closing_mass_vs_l='M_close ~ l^(1/6) (floor-only) or l^(11/12) (consistent) at fixed xi, l < dx',
    dx_spread_orders=spread('dx'), N_spread_orders=spread('N counting'), l_spread_orders=spread('smearing l'),
    paired_spread_orders=spread('paired'))
out = dict(stage=24, thresholds=thresholds, inputs=dict(pedalino=PED, gisler=GIS, cheapest=CHEAP, ceiling_base=CEIL_BASE),
           rows=rows, grid=grid, summary=summary,
           notes=dict(
               regularisation='ASSUMED: for l >= dx the smeared-DP rate is suppressed as (dx/l)^2 with R -> l.',
               ceiling_dependence='The Eq. 47 upper bound carries r^3 and N; if the same smearing length sets both sides the window '
                                  'moves as l^(11/3). Stage 15 kept the ceiling at r_N = 1e-15 m.',
               indepDP='DERIVED alternative: for N constituents each of mass M/N separated by more than l, the DP self-energy sum is '
                       'N (M/N)^2 = M^2/N, and the ceiling sigma_a^2 ~ 1/N uses the same independent-constituent picture. The published '
                       'floor uses N M^2 (total mass). We cannot tell from the text which is intended; flagged, not adopted.',
               gisler_body='Gisler string treated as an equivalent sphere, r = %.2f um (rho ASSUMED 3100 kg/m^3).' % (GIS['r_body'] * 1e6)))
json.dump(out, open(os.path.join(ROOT, 'results/stage24.json'), 'w'), indent=1, default=float)
for r in rows:
    print('%-11s %-78s floor=%.2e ceil=%.2e window=%+7.2f %-6s Mclose=%.3g MDa v=%.3g sag=%.3g reach=%s ceil_shift=%+.1f' % (
        r['group'], r['label'][:78], r['floor'], r['ceiling'], r['window_orders'], r['status'], r['closing_mass_MDa'], r['closing_v_m_s'], r['sag_m'], r['reachable_in_stage22_constraints'], r['ceiling_shift_orders']))
print(summary); print(thresholds); print(PED['xi'], GIS['m_avg_u'], GIS['r_body'])
