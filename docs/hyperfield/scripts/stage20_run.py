#!/usr/bin/env python3
"""
Stage 20 — The real vacuum gap: species, temperature, and demonstrated number densities.

Key point (DERIVED): in the short-wavelength (saturated) regime every residual-gas collision fully decoheres,
and the Stage 12/16 rate Gamma = (16 pi/3) n R^2 sqrt(2 pi k T / m_a) equals (32 pi^2/3) R^2 J with
J = n vbar/4 the impingement flux.  The requirement is therefore a FLUX J_req, independent of species and
gas temperature.  A pressure quoted at one temperature cannot be compared with a pressure needed at another.

Steps
 1. Required flux / density per design from stage16_run.evaluate (unchanged budget, P solved by the 40% rule).
 2. Regime check: thermal de Broglie wavelength of He / H2 at each T vs R and dx (Stage 16 gas_rate interpolation).
 3. Demonstrated limits, converted to n and J (lifetime limits via Langevin capture, velocity-independent).
 4. He-3 / He-4 saturated vapour pressure (ITS-90 Table 3 where valid, zero-entropy extrapolation below),
    H2 SVP (Souers 1977), and 2D-adsorbed He on graphite equilibrium (Cole-Frankl-Goodstein E_b = 12.2 meV).
 5. Gap per design in orders vs each demonstrated limit.
Writes results/stage20.json.
"""
from __future__ import annotations
import json, os, sys
from math import pi, sqrt, log, log10, exp
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import stage16_run as s16  # noqa: E402

HBAR, K_B, U = s16.HBAR, s16.K_B, s16.U_DA
E = 1.602176634e-19
EPS0 = 8.8541878128e-12
M = dict(He4=4.002602 * U, He3=3.016029 * U, H2=2.01588 * U, pbar=1.67262192e-27)
ALPHA_VOL = dict(He4=0.2050522e-30, H2=0.787e-30)  # polarizability volumes (m^3), CRC handbook values
CAND = {c['key']: c for c in s16.CANDIDATES}
YEAR = 365.25 * 86400


def vbar(T, m):
    return sqrt(8 * K_B * T / (pi * m))


def lam_th(T, m):
    return 2 * pi * HBAR / sqrt(2 * pi * m * K_B * T)


def langevin(species):
    """Langevin capture rate coefficient for pbar + neutral (velocity independent), m^3/s (DERIVED):
       k = 2 pi sqrt(alpha_vol e^2 / (4 pi eps0 mu)). Every capture ASSUMED to annihilate."""
    mu = M['pbar'] * M[species] / (M['pbar'] + M[species])
    return 2 * pi * sqrt(ALPHA_VOL[species] * E * E / (4 * pi * EPS0 * mu))


# ---------------- helium vapour pressure ----------------
ITS90 = {  # Table 3 (Preston-Thomas, Metrologia 27, 3 (1990))
    'He3': dict(rng=(0.65, 3.2), A=[1.053447, 0.980106, 0.676380, 0.372692, 0.151656, -0.002263, 0.006596, 0.088966, -0.004770, -0.054943], B=7.3, C=4.3),
    'He4_low': dict(rng=(1.25, 2.1768), A=[1.392408, 0.527153, 0.166756, 0.050988, 0.026514, 0.001975, -0.017976, 0.005409, 0.013259, 0.0], B=5.6, C=2.9),
    'He4_high': dict(rng=(2.1768, 5.0), A=[3.146631, 1.357655, 0.413923, 0.091159, 0.016349, 0.001826, -0.004325, -0.004973, 0.0, 0.0], B=10.3, C=1.9),
}


def its90_T(p, key):
    c = ITS90[key]
    x = (log(p) - c['B']) / c['C']
    return sum(a * x ** i for i, a in enumerate(c['A']))


BRACKET = {'He3': (10.0, 2e5), 'He4_low': (20.0, 6e3), 'He4_high': (4e3, 3e5)}  # p ranges inside each ITS-90 sub-range


def its90_p(T, key):
    lo, hi = BRACKET[key]
    for _ in range(200):
        mid = sqrt(lo * hi)
        if its90_T(mid, key) < T:
            lo = mid
        else:
            hi = mid
    return sqrt(lo * hi)


def p_zero_entropy(T, m, L0_K):
    """DERIVED low-T limit (condensed-phase entropy -> 0): p = (m kT / 2 pi hbar^2)^{3/2} kT exp(-L0/T)."""
    nq = (m * K_B * T / (2 * pi * HBAR ** 2)) ** 1.5
    return nq * K_B * T * exp(-L0_K / T)


def fit_L0(T, p, m):
    nq = (m * K_B * T / (2 * pi * HBAR ** 2)) ** 1.5
    return -T * log(p / (nq * K_B * T))


def svp(T, iso):
    if iso == 'He3':
        if T >= 0.65:
            return its90_p(T, 'He3'), 'ITS-90'
        return p_zero_entropy(T, M['He3'], L0['He3']), 'extrapolated (DERIVED)'
    if T >= 2.1768:
        return its90_p(T, 'He4_high'), 'ITS-90'
    if T >= 1.25:
        return its90_p(T, 'He4_low'), 'ITS-90'
    return p_zero_entropy(T, M['He4'], L0['He4']), 'extrapolated (DERIVED)'


L0 = dict(He3=fit_L0(0.65, its90_p(0.65, 'He3'), M['He3']), He4=fit_L0(1.25, its90_p(1.25, 'He4_low'), M['He4']))


def svp_H2(T):
    """Souers et al. 1977 (UCRL-52226, osti 7110677) Table 6, normal H2 solid: ln p[Pa] = A + B/T + B' ln T (fit 4-14 K)."""
    return exp(7.570953 - 86.94152 / T + 2.860678 * log(T))


def n_ads_eq(T, theta, Eb_K=141.8, m=M['He4'], n_mono=6.37e18):
    """DERIVED: ideal 2D adsorbed gas in equilibrium with 3D gas: n3 = (theta n_mono / lambda_th) exp(-E_b/kT).
       E_b = 12.22 meV = 141.8 K (He-4 on graphite, Cole, Frankl & Goodstein RMP 53, 199 (1981));
       n_mono = sqrt(3) commensurate density 0.0637 A^-2."""
    return theta * n_mono / lam_th(T, m) * exp(-Eb_K / T)


def main():
    out = dict(stage=20)
    # ---------------- 1. requirements ----------------
    designs = [
        ('stage12_cheap_diamond', 'diamond_spin', 5.6234e-15, 3.852e-6, 2.5),
        ('stage16_Pb_meissner', 'lead_meissner', 1e-14, 2.6915348e-6, 1.0),
        ('stage19_Pino_Nb', 'niobium_meissner', 2e13 * U, 1.2302687708e-06, 0.5),
    ]
    D = {}
    for key, ck, m, dx, t in designs:
        cand = CAND[ck]
        r = s16.evaluate(m, dx, t, cand)
        T = cand['T_gas']
        R = r['R']
        G_allow = r['gamma_gas']                     # allowed gas decoherence rate (s^-1), 40% rule
        J_req = 3 * G_allow / (32 * pi ** 2 * R ** 2)  # DERIVED, species/T independent in saturated regime
        req = dict(m=m, dx=dx, t=t, R=R, T_gas=T, P_req_Pa_He=r['P_req'], n_req_He=r['n_req'], Gamma_gas_allowed=G_allow,
                   J_req=J_req, N5=r['N5'])
        # per species / temperature: n_req = 4 J_req / vbar
        req['n_req'] = {f'{sp}@{TT}K': 4 * J_req / vbar(TT, M[sp]) for sp in ('He4', 'He3', 'H2') for TT in (0.015, 0.05, 1.0, 4.0, 6.2, 300.0)}
        # regime check
        req['regime'] = {}
        for sp in ('He4', 'H2'):
            lt = lam_th(T, M[sp])
            req['regime'][sp] = dict(lambda_th=lt, dx_over_lambda=dx / lt, R_over_lambda=R / lt,
                                     saturated_fraction=s16.gas_rate(1e-15, T, R, dx, M[sp])[0] / s16.gas_rate(1e-15, T, R, dx, M[sp])[1])
        # dx needed before long-wavelength (non-saturated) regime would give a factor-2 reduction (DERIVED): u = 1
        req['dx_for_half_rate'] = s16.gas_rate(1e-15, T, R, dx)[2]
        # Rate check: stage16 gas rate reproduces J formula
        req['check_Gamma_from_J'] = (32 * pi ** 2 / 3) * R ** 2 * (r['n_req'] * vbar(T, M['He4']) / 4) / G_allow
        D[key] = req
    out['designs'] = D

    # ---------------- 3. demonstrated limits ----------------
    kHe, kH2 = langevin('He4'), langevin('H2')
    lims = []

    def lim(key, label, n, T, species, kind, src, note=''):
        J = n * vbar(T, M[species]) / 4
        lims.append(dict(key=key, label=label, n=n, T=T, species=species, kind=kind, J=J,
                         P_at_T_Pa=n * K_B * T, src=src, note=note))

    tau_base = 10.2 * YEAR
    nHe_base, nH2_base = 1 / (tau_base * kHe), 1 / (tau_base * kH2)
    lim('BASE_2017_He', 'BASE reservoir trap, 405 d pbar storage (He-only assumption)', nHe_base, 6.2, 'He4', 'upper limit (lifetime, 68% CL)',
        'https://iopscience.iop.org/article/10.1088/1367-2630/aa7e73',
        f'tau>10.2 a; Langevin k_He={kHe:.3g} m^3/s; quoted p_He<2.7e-18 mbar; DERIVED p at 6.2 K = {nHe_base*K_B*6.2/100:.2g} mbar')
    lim('BASE_2017_H2', 'BASE reservoir trap (H2-only assumption)', nH2_base, 6.2, 'H2', 'upper limit (lifetime, 68% CL)',
        'https://iopscience.iop.org/article/10.1088/1367-2630/aa7e73',
        f'k_H2={kH2:.3g}; quoted p_H<1.2e-18 mbar; DERIVED {nH2_base*K_B*6.2/100:.2g} mbar')
    tau_g = 3.4 * 30.44 * 86400
    lim('Gabrielse_1990', 'TRAP/Gabrielse 1990, 59 d pbar storage at 4.2 K', 1e8, 4.2, 'He4', 'upper limit (lifetime)',
        'https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.65.1317 ; https://inspirehep.net/literature/306801',
        f'quoted n<100 cm^-3, p<5e-17 Torr at 4.2 K; Langevin He with tau>3.4 months gives n<{1/(tau_g*kHe):.2g} m^-3 (consistent). '
        'INSPIRE abstract "more than 10^3 days" is a lost superscript of "103 days" = 3.4 months (DERIVED).')
    lim('CSR_2016', 'MPIK Cryogenic Storage Ring, ~6 K (OH- detachment on H2)', 1.4e8, 6.0, 'H2', 'upper limit (collision rate)',
        'https://doi.org/10.1063/1.4953888', 'n<140 cm^-3 at (5.5+-1) K, 35 m ring')
    lim('CryPTEx_2012', 'MPIK cryogenic Paul trap, 4 K, 30 h ion storage', 1e-13 / (K_B * 4.0), 4.0, 'H2', 'upper limit (conservative)',
        'https://doi.org/10.1063/1.4742770', 'p_H2 < ~1e-15 mbar at 4 K')
    lim('XHV_RT_plateau', 'Room-temperature XHV gauge plateau ~1e-14 Torr (Redhead 1999)', 1e-14 * 133.322 / (K_B * 300), 300.0, 'H2',
        'measured (gauge-limited)', 'https://www.jlab.org/accel/inj_group/docs/2006/XHV.pdf', 'Redhead, Vacuum 53, 137 (1999) via Stutzman 2006')
    lim('MAQRO_req', 'MAQRO requirement / early Astrium estimate 1e-13 Pa (~16-20 K)', 1e-13 / (K_B * 20), 20.0, 'H2', 'proposed',
        'https://arxiv.org/abs/1307.7021', 'outgassing directly to space behind thermal shields')
    lim('MAQRO_QPPF', 'MAQRO QPPF study: "at best 1e-11 Pa" with cover', 1e-11 / (K_B * 20), 20.0, 'H2', 'study estimate',
        'https://www.nature.com/articles/s42005-021-00656-7', 'Gasbarri et al., Commun. Phys. 4, 155 (2021)')
    lim('Pino_assumed', 'Pino 2018 assumption 1e-17 mbar at 50 mK', 1e-15 / (K_B * 0.05), 0.05, 'He4', 'assumed',
        'https://arxiv.org/abs/1603.01553', 'pressure taken from 4 K Penning-trap limits, applied at 50 mK')
    lim('Stage16_SOA', 'Stage 16 SOA 1e-15 Pa at 4 K (He)', 1e-15 / (K_B * 4.0), 4.0, 'He4', 'Stage 12/16 convention', 'HYPERFIELD Stage 16', '')
    # solar wind (not a neutral gas; a directed 1 keV proton beam): flux n v
    sw = dict(n=5e6, v=4e5, J=5e6 * 4e5, src='https://pwg.gsfc.nasa.gov/Education/wsolwind.html ; http://hdl.handle.net/2060/19750016535',
              note='directed flux n v (not n vbar/4); must be shadowed by the spacecraft/shields')
    out['limits'] = lims
    out['solar_wind'] = sw
    out['langevin'] = dict(k_He=kHe, k_H2=kH2)

    # ---------------- 4. vapour pressures and cryosorption ----------------
    vp = []
    for T in (0.015, 0.05, 0.1, 1.0, 4.0):
        row = dict(T=T)
        for iso in ('He3', 'He4'):
            p, how = svp(T, iso)
            row[iso] = dict(P_Pa=p, n=p / (K_B * T), J=p / (K_B * T) * vbar(T, M[iso]) / 4, how=how)
        pH = svp_H2(T)
        row['H2'] = dict(P_Pa=pH, n=pH / (K_B * T), J=pH / (K_B * T) * vbar(T, M['H2']) / 4,
                         how='Souers 1977 fit (4-14 K)' + (', extrapolated' if T < 4 else ''))
        # max He coverage on graphite such that J <= J_req (strictest design)
        Jmin = min(d['J_req'] for d in D.values())
        n1 = n_ads_eq(T, 1.0)
        row['He_on_graphite'] = dict(n_at_full_monolayer=n1, theta_max=min(1.0, (4 * Jmin / vbar(T, M['He4'])) / n1) if n1 > 0 else 1.0)
        vp.append(row)
    out['vapour'] = vp

    # crossover temperatures where the SVP / adsorption flux equals each design's J_req (DERIVED)
    def J_of(T, kind):
        if kind == 'He3':
            p = svp(T, 'He3')[0]; m = M['He3']
        elif kind == 'He4':
            p = svp(T, 'He4')[0]; m = M['He4']
        elif kind == 'H2':
            p = svp_H2(T); m = M['H2']
        else:  # He on graphite at coverage theta
            th = float(kind.split('_')[1])
            return n_ads_eq(T, th) * vbar(T, M['He4']) / 4
        return p / (K_B * T) * vbar(T, m) / 4

    def T_cross(J, kind, lo=0.005, hi=6.0):
        for _ in range(200):
            mid = sqrt(lo * hi)
            if J_of(mid, kind) > J:
                hi = mid
            else:
                lo = mid
        return sqrt(lo * hi)

    out['T_crossover_K'] = {dk: {k: T_cross(d['J_req'], k) for k in ('He3', 'He4', 'H2', 'graphite_1.0', 'graphite_0.001')} for dk, d in D.items()}
    # allowed leak into a cold, perfectly cryopumping chamber: each atom hits the particle region once (ASSUMED);
    # steady state J at the walls = Q / (s A); take s = 1, A = 0.01 m^2 (10 cm box) -> Q_max = J_req A
    A_ch = 0.01
    out['leak_allowance'] = {dk: dict(Q_atoms_s=d['J_req'] * A_ch, Q_mbarL_s=d['J_req'] * A_ch * K_B * 300 * 10, A_m2=A_ch,
                                      note='mbar L/s at 300 K; Pfeiffer ASM 340 hard-vacuum minimum detectable He leak <=5e-13 mbar L/s (https://heliumleakdetection.net/helium-leak-detectors/pfeiffer-adixen-leak-detectors/pfeiffer-asm-340-helium-leak-detector/)') for dk, d in D.items()}
    out['charcoal_He_capacity_4K'] = dict(theta_max=vp[-1]['He_on_graphite']['theta_max'],
                                          atoms_per_g=vp[-1]['He_on_graphite']['theta_max'] * 6.37e18 * 1000.0,
                                          note='ASSUMED 1000 m^2/g and uniform graphite binding (real charcoal has stronger pore sites)')
    out['checks'] = dict(L0_fit_K=L0, L0_fit_J_per_mol={k: v * 8.314462618 for k, v in L0.items()},
                         He4_1K_extrap_vs_DonnellyBarenghi=svp(1.0, 'He4')[0] / 15.58,
                         He4_4p2K_ITS90_Pa=its90_p(4.2, 'He4_high'), He3_0p65K_ITS90_Pa=its90_p(0.65, 'He3'), He4_1p25K_ITS90_Pa=its90_p(1.25, 'He4_low'),
                         ITS90_roundtrip=[its90_T(its90_p(T, k), k) for T, k in ((0.65, 'He3'), (1.25, 'He4_low'), (4.0, 'He4_high'))],
                         H2_svp_4p2K_Pa=svp_H2(4.2), H2_svp_3K_Pa=svp_H2(3.0), H2_svp_2K_Pa=svp_H2(2.0),
                         CERN_H2_radiation_floor_Pa='1e-8..1e-7 (condensed H2 under 300 K radiation, http://cds.cern.ch/record/2746077)')

    # ---------------- 5. gaps ----------------
    gaps = {}
    for dk, d in D.items():
        gaps[dk] = {l['key']: log10(l['J'] / d['J_req']) for l in lims}
        gaps[dk]['solar_wind_unshielded'] = log10(sw['J'] / d['J_req'])
    out['gaps_orders'] = gaps  # positive = demonstrated/assumed flux exceeds requirement (short); negative = margin

    json.dump(out, open(os.path.join(ROOT, 'results/stage20.json'), 'w'), indent=1, default=float)

    # ---------------- print ----------------
    for dk, d in D.items():
        print(f"{dk}: R={d['R']:.3g} T={d['T_gas']} P_req(He)={d['P_req_Pa_He']:.3g} n_req(He)={d['n_req_He']:.3g} J_req={d['J_req']:.3g} "
              f"chk={d['check_Gamma_from_J']:.3f} dx_half={d['dx_for_half_rate']:.3g}")
        for sp, g in d['regime'].items():
            print('   ', sp, {k: f'{v:.3g}' for k, v in g.items()})
        print('    n_req', {k: f'{v:.2g}' for k, v in d['n_req'].items() if k.split('@')[0] in ('He4', 'H2')})
    print('Langevin', kHe, kH2)
    for l in lims:
        print(f"{l['key']:18s} n={l['n']:.3g} T={l['T']} J={l['J']:.3g} P@T={l['P_at_T_Pa']:.3g} Pa  {l['note']}")
    print('solar wind J', sw['J'])
    for r in vp:
        print(r['T'], {k: (f"{v['P_Pa']:.3g} Pa n={v['n']:.3g} ({v['how']})" if isinstance(v, dict) and 'P_Pa' in v else v) for k, v in r.items() if k != 'T'})
    print('checks', out['checks'])
    print('T_cross', json.dumps(out['T_crossover_K'], indent=0))
    print('leak', out['leak_allowance']); print('charcoal', out['charcoal_He_capacity_4K'])
    for dk, g in gaps.items():
        print(dk, {k: round(v, 2) for k, v in g.items()})


if __name__ == '__main__':
    main()
