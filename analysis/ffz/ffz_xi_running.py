#!/usr/bin/env python3
"""Does ξ run away from 1/6 in WZ+Yukawa? (new-idea rung)

Stacked on the surface/contact analysis (#20). Question: if the nonminimal
coupling ξ runs, a1∝(1/6−ξ) could feed a finite/simple contribution into the
induced-R pole and maybe explain GAP=+0.00965±0.0024.

(1) β_ξ / anomalous running of (ξ−1/6) through feasible loops, using:
      * SUSY chiral-multiplet kinetic term (superspace fixes conformal ξ)
      * Jack–Osborn local RG (δ=0 ⇒ additive b=0 in τ=0)
      * packet fixed-1/6 / local-Weyl checkpoints through O(h³)
(2) If ξ deviates, compute δξ with coefficients fixed by couplings, then
    a1∝(1/6−ξ) contribution to the h²z/K simple pole.
(3) Compare to GAP; test whether ln2/72 arises from anything derived here.
(4) Higher poles unaffected.

No tuning to 1/24. Derived vs assumed labeled.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_xi_running.py

Requires: sympy, mpmath
"""

from __future__ import annotations

from fractions import Fraction as F
from math import log, pi

import sympy as s

K = 16.0 * pi**2
C_LEDGER = 0.05131719
C_TARGET = 1.0 / 24.0
GAP = C_LEDGER - C_TARGET  # +0.00965052
GAP_ERR = 0.0024
LN2_OVER_72 = log(2) / 72.0


# ---------------------------------------------------------------------------
# (1) β_ξ and RG stability of ξ=1/6
# ---------------------------------------------------------------------------

def susy_chiral_fixes_xi():
    """In rigid N=1 on old-minimal curved superspace, the chiral kinetic term

        ∫ d⁴θ Ē Φ̄ Φ

    expands to a complex-scalar curvature coupling at the conformal value
    ξ=1/6 (together with the fermion / F-term partners). An independent
    deformation δξ R|φ|² without superpartners breaks SUSY.

    Therefore, within the SUSY WZ model used by the Newton test, ξ is not
    an independent running coupling: it is fixed to 1/6 by the superspace
    kinetic term.

    Label: DERIVED from superspace field content; ASSUMED that the benchmark
    stays in the SUSY completion (no extra non-SUSY ξ deformation).
    """
    return {
        "xi_independent_in_SUSY": False,
        "xi_fixed_value": F(1, 6),
        "label": (
            "derived (superspace kinetic ⇒ ξ=1/6); "
            "assumed (no SUSY-breaking ξ deformation)"
        ),
    }


def jack_osborn_beta_a():
    """Jack–Osborn local RG for a=1−6ξ.

    β_a = Γ a + b, with Γ the soft-mass / improvement anomalous rate.
    SUSY flat local-source prescription: Eq. (7.5) ⇒ δ=0 ⇒ additive
    curvature-source beta vanishes in the τ=0 scheme ⇒ b=0.

    Component check (packet local-Weyl checkpoint, JO 9.16 / 9.26):
      δ^(1)=0; three radial δ^(2) pieces −2λT + Th + Th cancel with
      T=tr(ȳy)=h/2, λ=h.

    Packet fixed-1/6: additive coefficients at h² and h³ vanish under the
    specified matter subtraction; a first new d_I pole at three loops in
    another continuation would first affect 4d additive beta at O(h⁴) —
    beyond the O(h²) Newton simple pole.

    Therefore β_a = Γ a  and  β_ξ = −β_a/6 = Γ (ξ − 1/6).
    Fixed point: ξ_*=1/6 is RG-stable (β_ξ=0 there).
    """
    h, u, v, c1, c2 = s.symbols("h u v c1 c2")
    # JO (5.21) / packet: γ = h/2 − h²/2 + O(h³)
    # B = 6 h γ = 3 h² − 3 h³ + O(h⁴)  (ungauged one-chiral; β_h|=4)
    gamma = h / 2 - h**2 / 2
    B = s.expand(6 * h * gamma)
    assert B == 3 * h**2 - 3 * h**3
    # Soft Γ = h ∂_h (B/h)  (Jack–Jones–Pickering / packet)
    Gamma = s.expand(h * s.diff(B / h, h))
    assert Gamma == 3 * h - 6 * h**2

    # Component δ^(2) radial cancelation (JO 9.26)
    T = h / 2
    lam = h
    delta2 = s.expand(-2 * lam * T + T * h + T * h)
    assert delta2 == 0

    # Additive b=0 ⇒ β_a = Γ a, β_ξ = Γ (ξ − 1/6)
    a, xi = s.symbols("a xi")
    b = 0
    beta_a = Gamma * a + b
    beta_xi = s.simplify(-beta_a.subs(a, 1 - 6 * xi) / 6)
    # Should equal Γ (ξ − 1/6)
    expected = s.expand(Gamma * (xi - s.Rational(1, 6)))
    assert s.expand(beta_xi - expected) == 0
    at_fixed = s.simplify(beta_xi.subs(xi, s.Rational(1, 6)))
    assert at_fixed == 0

    # Finite scheme shift a'=a+f cannot create a nonzero additive at h²
    # and only (u−v)c1+3c2 at h³ (packet); with SUSY b=0 both vanish.
    beta_h = 3 * h**2 + u * h**3
    Gamma_gen = 3 * h + v * h**2
    f = c1 * h + c2 * h**2
    shift = s.expand(beta_h * s.diff(f, h) - Gamma_gen * f)
    assert shift.coeff(h, 2) == 0

    return {
        "gamma_series": gamma,
        "B_series": B,
        "Gamma_series": Gamma,
        "b_additive": 0,
        "beta_xi": beta_xi,
        "beta_xi_at_1_6": 0,
        "delta2_radial": 0,
        "xi_star": F(1, 6),
        "RG_stable": True,
        "label": (
            "derived: β_ξ=Γ(ξ−1/6) with b=0 (JO 7.5 + component check); "
            "derived: ξ=1/6 is a fixed point; "
            "assumed: τ=0 SUSY local-source prescription of the packet"
        ),
    }


def packet_fixed_one_sixth_additive():
    """Packet fixed-1/6 checkpoint: additive a-beta coeffs vanish at h², h³."""
    h, e = s.symbols("h eps")
    B = 3 * h**2 - 3 * h**3
    H = B - 2 * e * h
    Gamma = h * s.diff(B / h, h)
    Zinv = 1 - B / (2 * e * h)
    assert s.expand(H * s.diff(Zinv, h) - Gamma * Zinv) == 0
    f = B / (6 * h)  # = h/2 − h²/2
    assert s.expand(f - (h / 2 - h**2 / 2)) == 0
    assert s.expand(B * s.diff(f, h) - Gamma * f) == 0
    # a_C conversion cancels additive bare part through h²
    Z = 1 + 3 * h / (2 * e) + h**2 * (9 / (4 * e**2) - 3 / (2 * e))
    ad = e / (3 - 2 * e)
    a = s.symbols("a")
    aC = a - ad + B / (2 * h * (3 - 2 * e))
    series = s.series(ad + Z * aC - Z * a, h, 0, 3).removeO()
    assert s.simplify(series) == 0
    return {
        "f_conversion": f,
        "additive_h2": 0,
        "additive_h3": 0,
        "label": "derived (packet fixed-1/6 runnable identities)",
    }


# ---------------------------------------------------------------------------
# (2) δξ and a1 contribution
# ---------------------------------------------------------------------------

def delta_xi_and_a1():
    """If the theory starts at ξ=1/6, β_ξ=0 ⇒ δξ=0 at all feasible orders.

    Then a1_bulk=(1/6−ξ)R=0 still, and there is no derived (1/6−ξ) feed into
    the h²z/K simple pole.

    A SUSY-breaking independent ξ deformation is OUTSIDE the model; its
    coefficient would be a free knob (tuning), not fixed by h,z,K.
    """
    susy = susy_chiral_fixes_xi()
    jo = jack_osborn_beta_a()
    # Starting at the fixed point
    xi0 = F(1, 6)
    delta_xi = 0  # β_ξ(xi0)=0 and SUSY fixes xi0
    a1_coeff = F(1, 6) - xi0  # per real scalar density /R
    assert a1_coeff == 0
    # Contribution to the three-loop Newton residue (units h²z/K):
    # a pure vacuum (1/6−ξ)/ε term has no z and is not an h²z/K piece;
    # with δξ=0 it vanishes anyway.
    delta_c_h2zK = 0.0
    return {
        "xi0": xi0,
        "delta_xi": delta_xi,
        "a1_coeff": a1_coeff,
        "delta_c": delta_c_h2zK,
        "free_knob": False,  # within SUSY; deformation would be a free knob
        "susy": susy,
        "jo": jo,
        "label": "derived: δξ=0 at ξ=1/6; δc=0",
    }


def hypothetical_nonsusy_shift_needed():
    """Honesty check: what additive b would be needed to move ξ enough
    for a1 to look like GAP? This is NOT a prediction — it shows tuning.

    At one loop, vacuum δ(1/(16πG)) ∼ (1/6−ξ)/(16π² ε) per real scalar.
    That structure is NOT the h²z/K ledger entry. To fake a match one would
    still need an extra conversion into h²z/K (free / scheme). So even a
    shifted ξ_* cannot derive GAP without further assumptions.
    """
    return {
        "can_nonsusy_xi_shift_derive_GAP": False,
        "reason": (
            "vacuum a1 pole lacks the z weight of the h²z/K test; "
            "converting it into that channel reintroduces a free scheme/"
            "matching choice. SUSY moreover forbids the shift."
        ),
        "label": "derived (operator mismatch + SUSY protection)",
    }


# ---------------------------------------------------------------------------
# (3) vs GAP; ln2/72
# ---------------------------------------------------------------------------

def compare_to_gap(delta_c: float):
    ratio = delta_c / GAP if GAP else float("nan")
    ln2_from_xi = False  # honest: nothing in β_ξ produces ln2/72
    return {
        "GAP": GAP,
        "GAP_err": GAP_ERR,
        "predicted": delta_c,
        "ratio": ratio,
        "ln2_72": LN2_OVER_72,
        "ln2_72_ratio_to_gap": LN2_OVER_72 / GAP,
        "ln2_72_from_xi_running": ln2_from_xi,
        "ln2_note": (
            "ln2/72 is NOT produced by β_ξ, JO δ=0, or a1∝(1/6−ξ). "
            "It remains an unexplained numerical proximity to GAP."
        ),
    }


def higher_poles():
    return {
        "triple_ok": True,
        "double_ok": True,
        "reason": (
            "ξ fixed at 1/6: no new (1/6−ξ) tadpole; soft dictionary and "
            "CT unchanged. Finite ξ deformation would be simple/finite only."
        ),
        "label": "derived",
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=== ξ running away from 1/6? ===\n")
    print("No tuning to 1/24. Stacked on surface/contact (#20).\n")

    print("-- (1) β_ξ / RG stability --")
    susy = susy_chiral_fixes_xi()
    print(f"  SUSY chiral kinetic fixes ξ? (independent? {susy['xi_independent_in_SUSY']})")
    print(f"  fixed value: {susy['xi_fixed_value']}  [{susy['label']}]")
    jo = jack_osborn_beta_a()
    print(f"  Γ series: {jo['Gamma_series']}")
    print(f"  additive b: {jo['b_additive']}")
    print(f"  β_ξ = Γ(ξ−1/6);  β_ξ(1/6) = {jo['beta_xi_at_1_6']}")
    print(f"  δ^(2) radial cancelation: {jo['delta2_radial']}")
    print(f"  ξ_*=1/6 RG-stable? {jo['RG_stable']}  [{jo['label']}]")
    pkt = packet_fixed_one_sixth_additive()
    print(
        f"  packet additive at h²,h³: {pkt['additive_h2']}, {pkt['additive_h3']}"
        f"  [{pkt['label']}]"
    )

    print("\n-- (2) δξ and a1 → simple pole --")
    d = delta_xi_and_a1()
    print(f"  ξ0 = {d['xi0']}; δξ = {d['delta_xi']}")
    print(f"  a1 coeff (1/6−ξ) = {d['a1_coeff']}")
    print(f"  δc in units h²z/K = {d['delta_c']:+.6e}  [{d['label']}]")
    hyp = hypothetical_nonsusy_shift_needed()
    print(f"  non-SUSY shift can derive GAP? {hyp['can_nonsusy_xi_shift_derive_GAP']}")
    print(f"  reason: {hyp['reason']}")

    print("\n-- (3) vs GAP +0.00965 ± 0.0024; ln2/72 --")
    cmp_ = compare_to_gap(d["delta_c"])
    print(f"  GAP = {cmp_['GAP']:+.8f} ± {cmp_['GAP_err']}")
    print(f"  predicted = {cmp_['predicted']:+.8f}")
    print(f"  ratio predicted/GAP = {cmp_['ratio']:.6f}")
    print(
        f"  ln2/72 = {cmp_['ln2_72']:.8f} (ratio to GAP {cmp_['ln2_72_ratio_to_gap']:.4f})"
    )
    print(f"  ln2/72 from ξ running? {cmp_['ln2_72_from_xi_running']}")
    print(f"  {cmp_['ln2_note']}")

    print("\n-- (4) Higher poles --")
    hp = higher_poles()
    print(
        f"  triple_ok={hp['triple_ok']}, double_ok={hp['double_ok']}  [{hp['label']}]"
    )
    print(f"  {hp['reason']}")

    print("\nVerdict")
    print("-------")
    print(
        "(1) ξ=1/6 is RG-stable in SUSY WZ: superspace fixes ξ, and\n"
        "    Jack–Osborn gives β_ξ=Γ(ξ−1/6) with additive b=0 (τ=0).\n"
        "(2) Starting at ξ=1/6 ⇒ δξ=0 ⇒ a1 contribution δc=0.\n"
        "(3) ratio predicted/GAP = 0. ln2/72 does NOT arise from ξ running.\n"
        "(4) Triple/double poles unaffected.\n"
        "Plain: ξ does not run away from 1/6 here; this does not explain GAP."
    )
    print("\n=== XI RUNNING COMPLETE ===")
    return cmp_


if __name__ == "__main__":
    main()
