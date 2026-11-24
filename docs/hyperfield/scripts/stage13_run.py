#!/usr/bin/env python3
"""Stage 13 — measured diamond Im((ε-1)/(ε+2)) and the blackbody budget
at the Stage 12 cheapest QGEM point.

No parameters are fitted. Values are tagged CITED / DERIVED / ASSUMED.
"""
from __future__ import annotations
import json, os, sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from stage12_lib import (
    Point, evaluate, lambda_bb_scatter, lambda_bb_abs_em, gamma_bb,
    EPS_DIAMOND, RHO_DIAMOND, K_B, HBAR, C, concurrence_dephased,
    dephased_ppt_lmin,
)

# ---------------------------------------------------------------------------
# Conversion (DERIVED)
# χ = (ε-1)/(ε+2). For ε = ε' + iε'' with ε'' ≪ ε':
#   Im χ = 3 ε'' / (ε'+2)^2 ,  ε'' = ε' tanδ
#   => Im χ = 3 ε' tanδ / (ε'+2)^2
# ---------------------------------------------------------------------------
EPS_P = EPS_DIAMOND  # 5.7, CITED (standard diamond RF/mm-wave value; Millar: 5.66–5.67)


def im_chi_from_tand(tan_delta, eps=EPS_P):
    """DERIVED, low-loss limit."""
    return 3.0 * eps * tan_delta / (eps + 2.0) ** 2


def tand_from_im_chi(im_chi, eps=EPS_P):
    return im_chi * (eps + 2.0) ** 2 / (3.0 * eps)


IMCHI_PREFACTOR = 3.0 * EPS_P / (EPS_P + 2.0) ** 2  # ≈ 0.2882


# Wien peak of B_ν (frequency spectral radiance): ν_max / T ≈ 5.879e10 Hz/K
WIEN_GHZ_PER_K = 5.879e1  # CITED, Planck's law for B_ν

# Stage 12 cheapest point (DERIVED in Stage 12; frozen here)
CHEAP = dict(
    m=5.6234132519034904e-15,
    d=11.0e-6,
    dx=3.852248420036749e-6,
    t_int=2.5,
    tau_acc=0.25,
    t_tot=3.0,
    P=2.3958407731654315e-15,
    z_plate=5.0e-6,
    T_gas=4.0,
    R=7.25885450216044e-7,
    Phi=-0.5112780987191368,
)

# ---------------------------------------------------------------------------
# Literature cases. tanδ is the number we feed the converter unless im_chi is
# given directly (silica / Stage 12 assumption).
# ---------------------------------------------------------------------------
CASES = [
    dict(
        id='stage12_assumed',
        label='Stage 12 assumption',
        kind='assumed',
        im_chi=0.01,
        tan_delta=None,  # filled from converter: equivalent tanδ
        band='n/a (assumed, frequency-independent)',
        note='ASSUMED in Stage 12, not a measurement. Equivalent tanδ ≈ 0.035 at ε′=5.7, in the range of heavily defective / graphitic carbon rather than diamond.',
        source='Stage 12, chosen as a conservative 10× reduction from Romero-Isart silica Im χ = 0.1',
    ),
    dict(
        id='silica_romero',
        label='Fused silica (Romero-Isart 2011)',
        kind='cited_comparison',
        im_chi=0.1,
        tan_delta=None,
        band='thermal / IR, as used for silica nanospheres',
        note='CITED: Romero-Isart, Phys. Rev. A 84, 052121 (2011), Table III, ε_bb = 2.1 + i 0.57, which gives Im[(ε-1)/(ε+2)] = 0.1. This is silica, not diamond. Included only as the number QGEM-adjacent optomechanics inherited.',
        source='Romero-Isart PRA 84, 052121 (2011) Table III; arXiv:1110.4495',
    ),
    dict(
        id='bulk_best_cvd',
        label='Best optical-grade poly-CVD (window)',
        kind='bulk',
        tan_delta=4.5e-6,
        band='170 GHz, 300 K (record); window-grade material is similar or slightly better at 70–100 K',
        note='CITED Garin / Parshin / Derkach et al.: record tanδ ≈ 4.5×10^{-6} at 170 GHz, 300 K, still ~10× above the theoretical phonon floor. Typical optical-grade “Diamond Materials” poly-CVD is specified at tanδ ≈ 10^{-5} in the 100–200 GHz band (CITED Mazzocchi et al. 2026, quoting the manufacturer).',
        source='Garin et al., NASA/CR and Diamond Relat. Mater. (2003–2005); Mazzocchi et al., arXiv:2512.03780 (2026)',
    ),
    dict(
        id='bulk_typical_window',
        label='Typical window-grade CVD',
        kind='bulk',
        tan_delta=1.0e-5,
        band='100–200 GHz, 300 K (and the usual gyrotron-window specification)',
        note='CITED typical specification for optical-grade polycrystalline CVD used as gyrotron windows. JJAP 2025 reports 2×10^{-6}–3.5×10^{-5} over 10–330 GHz on a 2-inch poly wafer, decreasing with frequency in the 220–330 GHz band. Single-crystal CVD wafers: tanδ ≤ 2×10^{-5} at 145–170 GHz (Yamada et al., Diamond Relat. Mater. 58, 1 (2015)).',
        source='Diamond Materials spec as quoted by Mazzocchi arXiv:2512.03780; Dang et al. JJAP (2025); Yamada et al. DRM 2015',
    ),
    dict(
        id='bulk_millar_22ghz_cryo',
        label='Poly-CVD at 22 GHz, 35–100 K (Millar)',
        kind='bulk_conservative',
        tan_delta=9.5e-4,
        band='22.3 GHz, 35–100 K (this sample); 22 GHz sits near the 0.4 K Wien peak and is relevant to 1 K tails',
        note='CITED Millar, le Floch, Tobar et al., arXiv:1302.4149 (Diamond Materials poly-CVD disk). Fits: tanδ(T)=7.950×10^{-4} T^{0.045} for T<100 K (almost flat; ≈ 9.3×10^{-4} at 35 K) and 2.137×10^{-4} T^{0.329} for T>100 K. This sample is much lossier at 22 GHz than window-grade material at 170 GHz, consistent with a conductivity / defect 1/f tail. VDK 2020 cites this group (Le Floch / Hartnett / Tobar) when treating ε as “real” at O(1) K — that is fair for Re(ε)≈5.67, but Im χ is not zero.',
        source='Millar et al. arXiv:1302.4149; van de Kamp et al. PRA 102, 062807 (2020) ref. [53]',
    ),
    dict(
        id='bulk_millar_thz_300K',
        label='Same poly-CVD at 0.67–1 THz, 300 K',
        kind='bulk',
        tan_delta=6.3e-5,
        band='672–1000 GHz, 300 K (Wien peak at ~11–17 K)',
        note='CITED Millar et al.: ε* = 5.661 (1 − j 6.3×10^{-5}) around 1 THz. Loss decreases with frequency relative to their 22 GHz data. Room-temperature THz data are a useful upper envelope for 10 K thermal photons, but the cryogenic 22 GHz number is more conservative near 1 K.',
        source='Millar et al. arXiv:1302.4149, eq. (11)',
    ),
    dict(
        id='bulk_scd_subkelvin',
        label='Single-crystal CVD, sub-K microwave',
        kind='bulk',
        tan_delta=2.0e-6,
        band='microwave (superconducting microstrip), sub-kelvin; TLS analysis',
        note='CITED Mazzocchi et al., arXiv:2512.03780 (2026), Table I: lowest measured total tanδ ≈ 2.04×10^{-6} (IrSCD), 5.94×10^{-6} (HPHT-seeded CVD clone); optical-grade poly-CVD 1.04×10^{-5}. Nitrogen-related TLS floors tanδ_TLS^0 ≈ 7×10^{-6}–2.4×10^{-5}. TLS loss typically *falls* from millikelvin to a few kelvin as two-level systems saturate, so these are conservative for a 1–10 K thermal-photon budget if other channels do not take over. A tiled-clone joint reached tanδ = 1.44×10^{-4} from boundary scattering, a warning that defects dominate when they are present.',
        source='Mazzocchi, Yamada, Scherer et al., arXiv:2512.03780 (24 Aug 2026 / dated in preprint)',
    ),
    dict(
        id='nanodiamond_surface',
        label='Nanodiamond with surface sp² (dirty)',
        kind='nanodiamond',
        tan_delta=1.0e-2,
        band='microwave 2.5–9 GHz (powders); used as a conservative effective tanδ for a micron particle with a graphitic / disordered shell',
        note='CITED Cuenca et al., ACS Omega 3, 2183 (2018): microwave tanδ of sub-micron diamond powders tracks surface sp²; 5% sp² is detectable; they extrapolate tanδ ≈ 0.016 at 100% sp² (2.5 GHz). ASSUMED effective tanδ = 10^{-2} for a “dirty” NV nanodiamond (HPHT-crushed, residual graphitic shell). A DERIVED core-shell check: a 2 nm shell with order-unity effective loss on a 0.73 μm radius core gives a volume fraction 3 ΔR/R ≈ 0.8% and an effective tanδ of the same order (10^{-3}–10^{-2}) if the shell is graphitic. Detonation ND and uncleaned commercial NV powder can be in this range or worse. THz-TDS on NCD films (e.g. ACS Appl. Electron. Mater. 2020) attributes 0.1–3 THz absorption fingerprints to sp² / grain-boundary phases.',
        source='Cuenca et al., ACS Omega 3, 2183 (2018); core-shell estimate DERIVED; THz-TDS qualitative from Liu et al., ACS AEM (2020)',
    ),
    dict(
        id='nanodiamond_cleaner',
        label='Cleaner milled CVD nanodiamond',
        kind='nanodiamond',
        tan_delta=1.0e-3,
        band='effective, microwave-to-THz; ASSUMED mid-range between window-grade bulk and dirty-surface powder',
        note='ASSUMED 10^{-3}, not a direct THz measurement on a single 0.7 μm sphere. Motivated by: (i) Millar’s cryogenic 22 GHz poly-CVD already sits here even without a surface shell; (ii) milled high-purity CVD nanodiamonds show far less near-IR absorption / laser heating than commercial HPHT NV powder (Frangeskou et al., NJP 20, 043016 (2018), arXiv:1608.04724); (iii) a thinner or less graphitic shell (0.5–1 nm, moderate tanδ) gives a similar effective number. This is the cautious “realistic NV particle” case, not a claim of a measured 10^{-3}.',
        source='Scale set by Millar 2013 + Cuenca 2018 + Frangeskou 2018; the 10^{-3} value itself is ASSUMED',
    ),
]


def _fill_case(c):
    c = dict(c)
    if c.get('im_chi') is None:
        c['im_chi'] = float(im_chi_from_tand(c['tan_delta']))
    if c.get('tan_delta') is None:
        c['tan_delta'] = float(tand_from_im_chi(c['im_chi'])) if c['im_chi'] and c['id'] != 'silica_romero' else None
    # silica uses a different ε; do not invert with diamond ε'
    if c['id'] == 'silica_romero':
        c['tan_delta'] = None
        c['tan_delta_note'] = 'silica ε=2.1+i0.57, not convertible with diamond ε′=5.7'
    c['im_chi'] = float(c['im_chi'])
    return c


CASES = [_fill_case(c) for c in CASES]


def cheap_point(T_env, T_int, im_chi):
    return Point(
        name='stage12_cheapest',
        geometry='parallel', shielded=True, trapped=True,
        m=CHEAP['m'], d=CHEAP['d'], dx=CHEAP['dx'],
        t_int=CHEAP['t_int'], tau_acc=CHEAP['tau_acc'],
        P=CHEAP['P'], T_gas=CHEAP['T_gas'], z_plate=CHEAP['z_plate'],
        T_env=float(T_env), T_int_temp=float(T_int), im_chi=float(im_chi),
    )


def concurrence_at(im_chi, T_env, T_int):
    return float(evaluate(cheap_point(T_env, T_int, im_chi))['concurrence'])


def tmax(im_chi, thresh, mode='eq', T_env_fixed=1.0, lo=0.15, hi=80.0, n=50):
    """Largest T with C > thresh. mode='eq' uses T_env=T_int=T; 'int' varies T_int at fixed T_env."""
    # If already dead at lo, return lo
    def C(T):
        if mode == 'eq':
            return concurrence_at(im_chi, T, T)
        return concurrence_at(im_chi, T_env_fixed, T)
    if C(lo) <= thresh:
        return float(lo)
    if C(hi) > thresh:
        return float(hi)
    a, b = lo, hi
    for _ in range(n):
        m = 0.5 * (a + b)
        if C(m) > thresh:
            a = m
        else:
            b = m
    return float(0.5 * (a + b))


def budget_parts(T_env, T_int, R, dx, im_chi):
    """Split Γ_BB into scatter / absorb / emit (DERIVED from Stage 12 formulae)."""
    lam = HBAR * C / (K_B * max(T_env, 1e-9))
    cap = dx * dx / (1.0 + (dx / lam) ** 2)
    g_sc = lambda_bb_scatter(T_env, R) * cap if T_env > 0 else 0.0
    g_abs = lambda_bb_abs_em(T_env, R, im_chi) * cap if T_env > 0 else 0.0
    g_em = lambda_bb_abs_em(T_int, R, im_chi) * cap if T_int > 0 else 0.0
    return dict(gamma_scatter=float(g_sc), gamma_absorb=float(g_abs),
                gamma_emit=float(g_em), gamma_bb=float(g_sc + g_abs + g_em),
                lambda_th_m=float(lam))


# Debye heat capacity for a diamond sphere (CITED θ_D ≈ 2200 K; formula DERIVED)
THETA_D = 2200.0  # K, CITED standard diamond Debye temperature
M_C = 12.011 * 1.66053906660e-27  # kg per C atom


def heat_capacity(T, m=CHEAP['m']):
    """C_V = (12 π^4 / 5) N k (T/θ_D)^3  (Debye low-T, DERIVED)."""
    N = m / M_C
    return (12.0 * np.pi ** 4 / 5.0) * N * K_B * (T / THETA_D) ** 3


def radiative_timescale(T, im_chi, R=CHEAP['R']):
    """Order-of-magnitude thermal time constant from BB emission (DERIVED, rough).

    Photon number emission rate in the long-wavelength limit is related to
    Λ_em; a typical thermal photon carries ~ 2.7 kT. Then
        τ ~ C_V T / (2.7 kT × Ṙ) ~ C_V / (2.7 k Ṙ).
    This is an order-of-magnitude check, not a full spectral cooling curve.
    """
    # total photon emission rate ≈ (something) * Λ; use energy: P ≈ ∫ ħω R(ω) dω
    # A standard small-sphere result: P_em = 4 π R^3 * (4 σ T^4 / c) * 3 ω-averaged
    # We integrate the known R(k) prefactor from MAQRO / Romero-Isart:
    #   R_tot = ∫ dk  (3 V k^3 c / π) n_th(k) Imχ , V=4πR^3/3
    #   = 8 π^5 R^3 c (kT/ℏc)^4 Imχ / 15     (DERIVED, Bose integral)
    # Then <E> = (π^4/30) / (π^2/6) kT = (π^2/5) kT ≈ 1.97 kT for the k^3 n_th weight?
    # Use 2.7 kT (mean energy of a Planck photon) as a round number (ASSUMED O(1)).
    im = max(im_chi, 1e-18)
    Rtot = (8.0 * np.pi ** 5 / 15.0) * (R ** 3) * C * (K_B * T / (HBAR * C)) ** 4 * im
    Eavg = 2.7 * K_B * T  # ASSUMED O(1) Planck mean
    P = Rtot * Eavg
    Cv = heat_capacity(T)
    if P <= 0:
        return np.inf
    return float(Cv * T / P)  # τ such that C dT/dt ~ -P, so T/τ ~ P/C


def main():
    os.chdir(ROOT)
    # reference evaluate at Stage 12 conditions
    ev0 = evaluate(cheap_point(1.0, 1.0, 0.01))

    T_grid = np.unique(np.concatenate([
        np.linspace(0.2, 4, 40),
        np.linspace(4, 20, 80),
        np.linspace(20, 40, 20),
    ]))

    curves = {}
    results_cases = []
    for c in CASES:
        im = c['im_chi']
        C_eq = [concurrence_at(im, T, T) for T in T_grid]
        C_int = [concurrence_at(im, 1.0, T) for T in T_grid]
        rec = dict(
            id=c['id'], label=c['label'], kind=c['kind'],
            im_chi=im,
            tan_delta=c.get('tan_delta'),
            band=c['band'], note=c['note'], source=c['source'],
            Tmax_eq_Cgt0=tmax(im, 0.0, 'eq'),
            Tmax_eq_Cgt0p1=tmax(im, 0.1, 'eq'),
            Tmax_int_Tenv1K_Cgt0=tmax(im, 0.0, 'int', T_env_fixed=1.0),
            Tmax_int_Tenv1K_Cgt0p1=tmax(im, 0.1, 'int', T_env_fixed=1.0),
            Tmax_int_Tenv4K_Cgt0=tmax(im, 0.0, 'int', T_env_fixed=4.0),
            Tmax_int_Tenv4K_Cgt0p1=tmax(im, 0.1, 'int', T_env_fixed=4.0),
            C_at={
                '1K_eq': concurrence_at(im, 1, 1),
                '4K_eq': concurrence_at(im, 4, 4),
                '10K_eq': concurrence_at(im, 10, 10),
                'Tenv1_Tint4': concurrence_at(im, 1, 4),
                'Tenv1_Tint20': concurrence_at(im, 1, 20),
                'Tenv4_Tint4': concurrence_at(im, 4, 4),
            },
            budget_4K_eq=budget_parts(4, 4, CHEAP['R'], CHEAP['dx'], im),
            budget_1K_eq=budget_parts(1, 1, CHEAP['R'], CHEAP['dx'], im),
            tau_rad_4K_s=radiative_timescale(4.0, im),
            tau_rad_1K_s=radiative_timescale(1.0, im),
        )
        results_cases.append(rec)
        curves[c['id']] = dict(
            T=T_grid.tolist(),
            C_eq=C_eq,
            C_int_Tenv1=C_int,
        )

    # scattering-only ceiling (Im χ = 0): isolates T_env scattering
    im0 = 0.0
    scatter_ceiling = dict(
        Tmax_eq_Cgt0=tmax(im0, 0.0, 'eq'),
        Tmax_eq_Cgt0p1=tmax(im0, 0.1, 'eq'),
        note='DERIVED: Im χ = 0, so only blackbody scattering (depends on Re χ, not loss) plus air and GGN. This is the hard T_env ceiling for a perfect crystal.',
    )

    # Im χ sweep for the figure
    im_sweep = np.logspace(-8, -1, 50)
    sweep = []
    for im in im_sweep:
        sweep.append(dict(
            im_chi=float(im),
            Tmax_eq_Cgt0=tmax(im, 0.0, 'eq'),
            Tmax_eq_Cgt0p1=tmax(im, 0.1, 'eq'),
            Tmax_int_Cgt0=tmax(im, 0.0, 'int', 1.0),
            Tmax_int_Cgt0p1=tmax(im, 0.1, 'int', 1.0),
        ))

    conversion = dict(
        eps_prime=EPS_P,
        formula='Im[(ε-1)/(ε+2)] = 3 ε′ tanδ / (ε′+2)^2   (low-loss, DERIVED)',
        prefactor=float(IMCHI_PREFACTOR),
        wien_GHz_per_K=WIEN_GHZ_PER_K,
        wien_note='Frequency of the B_ν maximum is ≈ 58.8 GHz × T/K (CITED, Planck). A 0.1–10 K thermal spectrum is roughly 10 GHz–1 THz, with the 1 K peak at 59 GHz and the 4 K peak at 235 GHz.',
    )

    qgem_assumptions = dict(
        romero_isart_2011='Silica, not diamond: ε_bb=2.1+i0.57, Im[(ε-1)/(ε+2)]=0.1. CITED PRA 84, 052121 Table III. The Stage 12 default Im χ=0.01 was a 10× “diamond is cleaner than silica” guess, not a measurement.',
        yin_2013='Yin et al., Phys. Rev. A 88, 033614 (2013), arXiv:1305.1701, state that the THz emission rate of diamond is “smaller by about two orders of magnitude” than the silica sphere of Romero-Isart, implying Im χ ~ 10^{-3} for diamond at (room-temperature) THz. They do not cite a tanδ measurement. 10^{-3} is 300× larger than window-grade mm-wave Im χ, and is closer to our dirty-nanodiamond / Millar-22 GHz cases.',
        vdk_2020='van de Kamp et al., PRA 102, 062807, treat ε as real and T-independent at T~O(1) K, citing Ibarra et al. (Diamond Relat. Mater. 1997) and Le Floch / Hartnett / Tobar (the Millar 2013 data set). They note that for diamond, phonon excitations are negligible below T_int < 4 K (citing Bose 2017). They never quote Im[(ε-1)/(ε+2)]. They say air dominates and T_int can be raised to a few kelvin “without interfering with the results significantly.” That conclusion holds for window-grade Im χ; it is tight for Im χ=0.01 and for dirty nanodiamond.',
        schut_2024='Schut et al., PRR 6, 013199, cite blackbody and air as a decoherence package and take γ ~ 10^{-2} Hz at ~1 K from earlier papers. No independent Im χ.',
        bose_2017='Bose et al. PRL 2017 mention internal cooling / Tint and diamond phonons, without a loss-tangent number.',
    )

    heating = dict(
        magnetic_trap='The Stage 12 cheap point is magnetically trapped (Schut 2024 layout). There is no 1064 nm optical-trap beam during the hold, which removes the dominant heating channel of levitated-optomechanics experiments.',
        nv_laser='CITED Frangeskou et al., New J. Phys. 20, 043016 (2018), arXiv:1608.04724: commercial HPHT NV nanodiamonds in an optical dipole trap have been reported above 800 K at modest vacuum (absorption by nitrogen and amorphous/graphitic defects). Milled high-purity CVD nanodiamonds showed no measurable laser-power-dependent heating up to ~300 mW at 3–4 mbar. Any green/IR NV initialization or readout laser left on during the 2.5 s hold, in 10^{-15} Pa, would likely raise T_int far above the 15–20 K emission ceiling. Lasers must be off (or used only at the start/end, with time to cool) for the budget below to apply.',
        radiative_cooling='At 1–4 K the radiative coupling of a 0.73 μm sphere is weak, especially if Im χ is window-grade (10^{-6}). The DERIVED order-of-magnitude radiative time constant τ ~ C_V T / P_em is seconds to many hours, depending on Im χ and T (see cases). Residual gas conduction is negligible at 2.4×10^{-15} Pa. So T_int is set by whatever it was after loading / last laser pulse, not by T_env, unless a long wait or an active cooling step is built in. Internal temperature is the binding thermal variable.',
        theta_D_K=THETA_D,
        C_V_1K_J_per_K=float(heat_capacity(1.0)),
        C_V_4K_J_per_K=float(heat_capacity(4.0)),
    )

    out = dict(
        conversion=conversion,
        cheapest_point=CHEAP,
        stage12_reference=dict(
            T_env=1.0, T_int=1.0, im_chi=0.01,
            concurrence=float(ev0['concurrence']),
            gamma_air=float(ev0['gamma_air']),
            gamma_bb=float(ev0['gamma_bb']),
            gamma_t=float(ev0['gamma_t']),
            Phi=float(ev0['Phi']),
        ),
        scatter_only_ceiling=scatter_ceiling,
        cases=results_cases,
        im_chi_sweep=sweep,
        curves=curves,
        qgem_assumptions=qgem_assumptions,
        heating=heating,
        verdict_one_liner='Window-grade diamond is comfortable at a 4 K cryostat and only dies of blackbody (mostly scattering) near 12–14 K; the Stage 12 “dead at 5 K” cliff was the Im χ=0.01 assumption. Dirty nanodiamond surfaces bring the cliff back to ~4–6 K. Internal temperature, not T_env, is binding if any laser is used in vacuum.',
    )
    path = os.path.join(ROOT, 'results', 'stage13.json')
    with open(path, 'w') as f:
        json.dump(out, f, indent=2)
    print('wrote', path)
    print('scatter-only Tmax C>0 / C>0.1:', scatter_ceiling)
    for r in results_cases:
        print(f"  {r['id']:28} Imχ={r['im_chi']:.3e}  Teq C>0 {r['Tmax_eq_Cgt0']:.2f}  C>0.1 {r['Tmax_eq_Cgt0p1']:.2f}  Tint@1K C>0 {r['Tmax_int_Tenv1K_Cgt0']:.2f}  C>0.1 {r['Tmax_int_Tenv1K_Cgt0p1']:.2f}")


if __name__ == '__main__':
    main()
