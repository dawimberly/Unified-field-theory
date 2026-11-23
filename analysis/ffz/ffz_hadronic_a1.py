#!/usr/bin/env python3
"""Controlled hadronic a1: ChPT / anomaly / lattice-AdS — replace free π,N proxy.

Stacked on #35. Key: does confinement remove 8 gluon Kabat (−2/3 each)?
Is net hadronic a1 positive? Quantify uncertainty.

No tuning. DERIVED vs ASSUMED.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_hadronic_a1.py
"""

from __future__ import annotations

from fractions import Fraction as F

A1_S = F(1, 6)
A1_W = F(1, 12)
A1_V = F(-2, 3)
A1_PROXY = F(2, 3)  # #31 free π,N proxy


def confinement_removes_gluon_kabat() -> dict:
    """IR hadronic EFT contains no free massless SU(3) gluons.

    Derived: Appelquist–Carazzone + confinement ⇒ gluons are not light
    degrees of freedom below Λ_QCD; the free-field Kabat contact for 8
    massless vectors is absent from the IR free-field sum.

    Not derived: a complete UV→IR matching theorem that every bit of the
    UV gluon a1 is canceled with no leftover curvature contact. Glueballs
    and multi-hadron states can carry residual positive scalar-like a1.
    """
    delta_remove_8 = -8 * A1_V  # removing 8×(−2/3) adds +16/3
    assert delta_remove_8 == F(16, 3)
    return {
        "free_gluons_in_IR_EFT": False,
        "Kabat_8_absent_from_IR_sum": True,
        "delta_a1_from_removing_8": delta_remove_8,
        "complete_UV_IR_matching": False,
        "label_absent": "derived (IR EFT content)",
        "label_matching": "assumed incomplete — residual glueball contacts possible",
        "refs": "AC decoupling; QCD confinement; HK matching literature",
    }


def chpt_pion_sector() -> dict:
    """Pions as pseudo-Goldstones; gravitational couplings from ChPT.

    Donoghue–Leutwyler gravitational form factors / ChPT: pion–curvature
    couplings exist; the effective ξ_π is not fixed to 0 or 1/6 by SM UV
    alone (ASSUMED range). Three light pions (π±,π0) as real scalars:
      a1_π ∈ [0, 3/6] = [0, 1/2] for ξ∈[1/6, 0].
    """
    a1_xi_conf = 3 * (A1_S - A1_S)  # ξ=1/6 → 0
    a1_xi_min = 3 * A1_S  # ξ=0 → 1/2
    assert a1_xi_conf == 0
    assert a1_xi_min == F(1, 2)
    return {
        "n_pi": 3,
        "a1_xi_1_6": a1_xi_conf,
        "a1_xi_0": a1_xi_min,
        "a1_band": (a1_xi_conf, a1_xi_min),
        "xi_fixed_by_SM": False,
        "label_band": "derived arithmetic; ξ_π range assumed (ChPT/DL)",
        "refs": "Donoghue–Leutwyler; ChPT gravitational FFs",
    }


def nucleon_and_leptons_photon() -> dict:
    """Stable IR matter + photon.

    Nucleons p,n as 2 Dirac ⇒ a1 = 2×(1/6) = 1/3.
    Light leptons e,μ,νe,νμ: 6 Weyl ⇒ 6/12 = 1/2.
    Photon: −2/3 (still massless; Kabat remains).
    """
    a1_N = 2 * F(1, 6)
    a1_lep = 6 * A1_W
    a1_g = A1_V
    assert a1_N == F(1, 3)
    assert a1_lep == F(1, 2)
    assert a1_g == F(-2, 3)
    return {
        "a1_nucleons": a1_N,
        "a1_leptons": a1_lep,
        "a1_photon": a1_g,
        "sum_fixed": a1_N + a1_lep + a1_g,  # 1/3+1/2−2/3 = 1/6
        "label": "derived (IR content; nucleons as free Dirac — mild assumption)",
    }


def glueball_residual() -> dict:
    """Lightest glueballs: massive scalars (lattice/AdS/QCD).

    O(1) light glueballs (f0-like), not 8 massless vectors. Each minimal
    scalar +1/6. Band: 0–3 light glueballs ⇒ a1 ∈ [0, 1/2] ASSUMED count.
    """
    return {
        "n_glueball_band": (0, 3),
        "a1_band": (0, F(1, 2)),
        "replaces_8_gluon_Kabat": False,  # they don't restore −16/3
        "label": "assumed O(1) count from lattice/AdS/QCD spectrum lore",
        "refs": "lattice glueball spectrum; AdS/QCD reviews",
    }


def controlled_hadronic_a1() -> dict:
    """Assemble band for net hadronic a1."""
    conf = confinement_removes_gluon_kabat()
    pi = chpt_pion_sector()
    mat = nucleon_and_leptons_photon()
    gb = glueball_residual()
    fixed = mat["sum_fixed"]  # 1/6
    # low = fixed + pi_low + gb_low = 1/6 + 0 + 0
    # high = fixed + pi_high + gb_high = 1/6 + 1/2 + 1/2 = 7/6
    a1_low = fixed + pi["a1_band"][0] + gb["a1_band"][0]
    a1_high = fixed + pi["a1_band"][1] + gb["a1_band"][1]
    assert a1_low == F(1, 6)
    assert a1_high == F(7, 6)
    return {
        "a1_low": a1_low,
        "a1_high": a1_high,
        "central_proxy_compare": A1_PROXY,
        "sign_OK_entire_band": a1_low > 0,
        "gluon_Kabat_removed": conf["Kabat_8_absent_from_IR_sum"],
        "uncertainty": (
            "ξ_π, glueball multiplicity, and UV→IR matching incompleteness; "
            "band [1/6, 7/6] — all positive but width O(1)"
        ),
        "plain": (
            "With gluon Kabats absent from the IR sum (derived as EFT content) "
            "and photon+leptons+nucleons fixed, net a1 ≥ 1/6 > 0 across the "
            "ChPT/glueball uncertainty band. Positivity is robust inside this "
            "framework; absolute value is NOT precise. Complete matching "
            "remains assumed."
        ),
    }


def self_checks() -> None:
    c = confinement_removes_gluon_kabat()
    assert c["Kabat_8_absent_from_IR_sum"] is True
    assert c["complete_UV_IR_matching"] is False
    h = controlled_hadronic_a1()
    assert h["a1_low"] == F(1, 6)
    assert h["a1_high"] == F(7, 6)
    assert h["sign_OK_entire_band"] is True
    assert h["gluon_Kabat_removed"] is True


def main():
    print("=== Controlled hadronic a1 (ChPT / anomaly / lattice-AdS) ===\n")
    conf = confinement_removes_gluon_kabat()
    print(f"  Free gluons in IR EFT? {conf['free_gluons_in_IR_EFT']}")
    print(f"  8×Kabat absent from IR sum? {conf['Kabat_8_absent_from_IR_sum']}  Δ={conf['delta_a1_from_removing_8']} [{conf['label_absent']}]")
    print(f"  Complete UV→IR matching? {conf['complete_UV_IR_matching']} [{conf['label_matching']}]")
    pi = chpt_pion_sector()
    print(f"  Pion a1 band [ξ=1/6..0]: [{pi['a1_band'][0]}, {pi['a1_band'][1]}] [{pi['label_band']}]")
    mat = nucleon_and_leptons_photon()
    print(f"  N+ℓ+γ fixed sum: {mat['sum_fixed']} [{mat['label']}]")
    gb = glueball_residual()
    print(f"  Glueball a1 band: [{gb['a1_band'][0]}, {gb['a1_band'][1]}] [{gb['label']}]")
    h = controlled_hadronic_a1()
    print(f"\n  Controlled a1 band: [{h['a1_low']}, {h['a1_high']}]  vs proxy {h['central_proxy_compare']}")
    print(f"  Sign OK entire band? {h['sign_OK_entire_band']}")
    print(f"  Uncertainty: {h['uncertainty']}")
    print(f"  Plain: {h['plain']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
Confinement removes free-gluon Kabats from the IR EFT sum (derived as
content). Controlled band a1 ∈ [1/6, 7/6] — positive throughout.
Uncertainty O(1) from ξ_π and glueballs; UV→IR matching incomplete.
Proxy +2/3 sits inside the band (not a coincidence claim — consistent).
"""
    )
    print("=== HADRONIC A1 COMPLETE ===")


if __name__ == "__main__":
    main()
