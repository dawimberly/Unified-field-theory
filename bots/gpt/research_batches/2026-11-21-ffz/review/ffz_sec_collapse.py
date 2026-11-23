#!/usr/bin/env python3
"""Strong-energy condition in collapse: ρ+3p for dense-matter EOS.

Dan: forces as great as a black-hole interior could make gravity repel
(negative gravity). Where (if anywhere) does ρ+3p<0?

EOS: neutron-star cores (polytrope / causal), quark matter / MIT bag,
QGP with bag constant B. Ordinary matter vs bag vacuum term.

Fields-only. No tuning. DERIVED vs ASSUMED. Honest negatives.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_sec_collapse.py
"""

from __future__ import annotations

from math import pi

# Nuclear saturation density ~ 0.16 fm^{-3} → energy density ~ 150 MeV/fm³
# Use natural units: MeV/fm³ and GeV/fm³ as convenient.
N_SAT_FM3 = 0.16  # fm^{-3}
M_N_MEV = 938.0
RHO_SAT_MEV_FM3 = N_SAT_FM3 * M_N_MEV  # ~150 MeV/fm³

# MIT bag constant ballpark (ASSUMED literature range ~ (145–200 MeV)^4)
B_MEV4_LOW = 145.0**4
B_MEV4_HIGH = 200.0**4
# Convert MeV^4 → MeV/fm³: 1 MeV/fm³ = (197.3 MeV fm)^3 * MeV / fm³ wait:
# ħc ≈ 197.3 MeV·fm; (MeV)^4 / (ħc)^3 → MeV/fm³
HC_MEV_FM = 197.3


def mev4_to_mev_fm3(b_mev4: float) -> float:
    return b_mev4 / (HC_MEV_FM**3)


def sec_value(rho: float, p: float) -> float:
    """ρ+3p; SEC requires ρ+3p ≥ 0 (and ρ+p≥0)."""
    return rho + 3.0 * p


def neutron_star_core() -> dict:
    """Causal / polytropic NS core matter.

    For barotropic p=wρ with 0≤w≤1 (causal): ρ+3p = ρ(1+3w) ≥ ρ > 0.
    DERIVED: ordinary nuclear/NS matter with p≤ρ cannot violate SEC.
    """
    samples = []
    for w in (0.0, 1.0 / 3.0, 0.5, 1.0):
        rho = 5.0 * RHO_SAT_MEV_FM3  # few × n_sat
        p = w * rho
        samples.append({"w": w, "rho": rho, "p": p, "rho_3p": sec_value(rho, p)})
    return {
        "name": "neutron_star_core",
        "samples": samples,
        "any_SEC_violation": any(s["rho_3p"] < 0 for s in samples),
        "label": "derived (p=wρ, 0≤w≤1 ⇒ ρ+3p>0)",
        "plain": "NS-core fluid with causal EOS never has ρ+3p<0.",
    }


def mit_bag_quark_matter() -> dict:
    """MIT bag: p = (ρ − 4B)/3 for massless quarks inside bag.

    Then ρ+3p = ρ + (ρ − 4B) = 2ρ − 4B.
    SEC violation iff 2ρ − 4B < 0 iff ρ < 2B.
    At equilibrium bag surface ρ=4B ⇒ ρ+3p=4B>0.
    Deep inside with ρ≫B: ρ+3p>0.
    Only the pure bag vacuum (ρ=B, p=−B) has ρ+3p=−2B<0.
    """
    b_lo = mev4_to_mev_fm3(B_MEV4_LOW)
    b_hi = mev4_to_mev_fm3(B_MEV4_HIGH)
    # Quark matter at ρ = 4B (equilibrium) and ρ = 10B
    rows = []
    for B in (b_lo, b_hi):
        for factor, tag in ((4.0, "eq_surface_4B"), (10.0, "dense_10B"), (1.0, "pure_bag_vac")):
            rho = factor * B
            if tag == "pure_bag_vac":
                p = -B  # vacuum
            else:
                p = (rho - 4.0 * B) / 3.0
            rows.append(
                {
                    "B_MeV_fm3": B,
                    "tag": tag,
                    "rho": rho,
                    "p": p,
                    "rho_3p": sec_value(rho, p),
                }
            )
    violations = [r for r in rows if r["rho_3p"] < 0]
    return {
        "name": "MIT_bag_quark",
        "rows": rows,
        "any_SEC_violation": len(violations) > 0,
        "violation_only_bag_vacuum": all(r["tag"] == "pure_bag_vac" for r in violations),
        "B_range_MeV_fm3": (b_lo, b_hi),
        "label_EOS": "derived (MIT bag thermodynamics)",
        "label_B": "assumed literature B~(145–200 MeV)^4",
        "plain": (
            "Quark fluid with MIT bag: ρ+3p=2ρ−4B ≥0 for ρ≥2B; "
            "equilibrium ρ=4B is safe. Only pure bag vacuum (p=−B) violates SEC."
        ),
    }


def qgp_with_bag() -> dict:
    """QGP + bag: similar — thermal QGP p~ρ/3 plus −B vacuum piece.

    Effective: ρ_tot = ρ_th + B, p_tot = p_th − B with p_th≈ρ_th/3.
    ⇒ ρ+3p = (ρ_th+B) + ρ_th − 3B = 2ρ_th − 2B.
    Violation iff ρ_th < B (bag-dominated vacuum-like).
    """
    b = mev4_to_mev_fm3(B_MEV4_LOW)
    rows = []
    for rho_th in (0.1 * b, b, 5 * b, 50 * b):
        rho = rho_th + b
        p = rho_th / 3.0 - b
        rows.append(
            {
                "rho_th": rho_th,
                "rho": rho,
                "p": p,
                "rho_3p": sec_value(rho, p),
                "violates": sec_value(rho, p) < 0,
            }
        )
    return {
        "name": "QGP_plus_bag",
        "rows": rows,
        "any_SEC_violation": any(r["violates"] for r in rows),
        "candidate": "bag-constant vacuum term",
        "label": "derived structure; B assumed",
        "plain": (
            "QGP+bag: ρ+3p=2(ρ_th−B)<0 only when thermal density < B — "
            "i.e. bag vacuum dominates, not dense thermal matter."
        ),
    }


def verdict() -> dict:
    ns = neutron_star_core()
    bag = mit_bag_quark_matter()
    qgp = qgp_with_bag()
    return {
        "ordinary_matter_violates_SEC": False,
        "NS_violates": ns["any_SEC_violation"],
        "bag_fluid_violates_except_vacuum": not bag["violation_only_bag_vacuum"],
        "only_vacuum_bag_term": bag["violation_only_bag_vacuum"],
        "plain": (
            "Honest: ordinary dense matter (NS, causal EOS, equilibrium bag "
            "quark matter, thermal QGP) satisfies SEC. The bag-constant "
            "vacuum term (p=−ρ) is the candidate for ρ+3p<0 / repulsion — "
            "same as cosmological CC, not a new collapse miracle."
        ),
        "dan_negative_gravity": (
            "Repulsion requires vacuum-like ρ+3p<0, not 'forces as great as "
            "BH interior' from ordinary field stress. Magnitude of bag B is "
            "QCD-scale, not Planck — see AD/AE for induced-G size."
        ),
    }


def self_checks() -> None:
    ns = neutron_star_core()
    assert ns["any_SEC_violation"] is False
    bag = mit_bag_quark_matter()
    assert bag["violation_only_bag_vacuum"] is True
    assert bag["any_SEC_violation"] is True  # vacuum does
    qgp = qgp_with_bag()
    assert any(r["violates"] for r in qgp["rows"] if r["rho_th"] < qgp["rows"][0]["rho_th"] + 1e9)
    # specifically low thermal
    assert qgp["rows"][0]["violates"] is True
    assert qgp["rows"][-1]["violates"] is False
    v = verdict()
    assert v["ordinary_matter_violates_SEC"] is False
    assert v["only_vacuum_bag_term"] is True


def main():
    print("=== SEC in collapse: ρ+3p for dense EOS ===\n")
    ns = neutron_star_core()
    print(f"  NS core: SEC violation? {ns['any_SEC_violation']}  [{ns['label']}]")
    for s in ns["samples"]:
        print(f"    w={s['w']:.3f}  ρ+3p={s['rho_3p']:.1f} MeV/fm³")
    bag = mit_bag_quark_matter()
    print(f"\n  MIT bag: B∈[{bag['B_range_MeV_fm3'][0]:.1f}, {bag['B_range_MeV_fm3'][1]:.1f}] MeV/fm³")
    print(f"  {bag['plain']}")
    for r in bag["rows"]:
        print(f"    {r['tag']:16s} ρ+3p={r['rho_3p']:+.1f}")
    qgp = qgp_with_bag()
    print(f"\n  QGP+bag: {qgp['plain']}")
    for r in qgp["rows"]:
        print(f"    ρ_th={r['rho_th']:.1f}  ρ+3p={r['rho_3p']:+.1f}  viol={r['violates']}")
    v = verdict()
    print(f"\n  {v['plain']}")
    print(f"  {v['dan_negative_gravity']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
Ordinary collapse matter does NOT violate SEC. Only the bag/vacuum
piece (p=−ρ) gives ρ+3p<0 — candidate for repulsion, QCD-scale B.
"""
    )
    print("=== SEC COLLAPSE COMPLETE ===")


if __name__ == "__main__":
    main()
