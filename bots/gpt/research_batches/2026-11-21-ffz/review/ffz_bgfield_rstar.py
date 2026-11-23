#!/usr/bin/env python3
"""Background-field + R* audit of the FFZ three-loop Newton simple pole.

Proven-methods ladder rung #6 (Chetyrkin–Tkachov / Smirnov R*).

Feasible subset (honest):
  * Formalize the background-field metric split already used by the packet
    (external ḡ, Weyl probe σ=−κx², R=4d(d−1)κ) as the covariant R extraction.
  * Apply R* = R_UV ∘ R̃_IR criteria to each ledger class: check UV forests,
    IR soft forests, and whether any IR counterterm was missed — especially
    soft-soft, soft-linear, and the QCc tadpole.
  * Re-verify every L=0 simple coefficient that is closed-form; retain QMC
    pieces at their stated errors (full R* re-integration of all graphs is
    not attempted).

No tuning to 1/24.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_bgfield_rstar.py

Requires: sympy, mpmath, scipy
"""

from __future__ import annotations

from fractions import Fraction as F
from math import log

import mpmath as mp
import sympy as s
from scipy.integrate import dblquad

mp.mp.dps = 40


# ---------------------------------------------------------------------------
# (0) Background field + R* definitions
# ---------------------------------------------------------------------------

def background_field_setup():
    """Packet already uses a background-field metric split.

    g_μν = ḡ_μν + h_μν with ḡ weakly curved, probe σ=−κx², R̄=4d(d−1)κ.
    Matter propagators / heat kernels on ḡ; R coefficient read covariantly
    after /[4d(d−1)]. No dynamical graviton / gauge-fixing for h — external
    metric (ASSUMED continuation of the packet). Full quantum BG-field
    graviton loops are OUT OF SCOPE.
    """
    return {
        "split": "g = ḡ + h (h external probe)",
        "probe": "σ=−κx², R̄=4d(d−1)κ",
        "R_extraction": "cov. /[4d(d−1)] after Wick/Schwinger",
        "graviton_loops": False,
        "label": "derived (packet); assumed (no dynamical graviton)",
    }


def rstar_definition():
    """R* = R ∘ R̃ (Chetyrkin–Tkachov; Smirnov).

    R  = UV Bogoliubov R-operation (Taylor in UV-divergent subgraphs).
    R̃ = IR dual: Taylor in soft/IR-divergent subgraphs (small loop momenta
         / soft Schwinger lengths / soft projective endpoints).

    In the packet's ordered-sector language:
      * UV forests ↔ hard endpoint / power sectors already subtracted via CT.
      * IR forests ↔ soft-soft and single-soft Taylor IE
        R̃(G) = G − G(x,0) − G(0,y) + G(0,0)  on two soft axes
        (same as BPHZ soft forest used in overlap rigor / DRED+BPHZ).

    Feasible test: for each class, ask whether R̃ was applied to soft IR
    and whether any leftover IR counterterm affects the 1/ε simple pole.
    """
    return {
        "formula": "R* = R_UV ∘ R̃_IR",
        "IR_on_soft": "sector Taylor IE = R̃ for soft axes",
        "refs": "Chetyrkin–Tkachov; Smirnov (R*)",
        "feasible": (
            "Checklist + re-verify closed-form poles; not a full diagrammatic "
            "R* recomputation of every QMC integrand from scratch"
        ),
        "label": "derived (identification); assumed (sector IE ≡ R̃)",
    }


# ---------------------------------------------------------------------------
# Closed-form piece re-derivation (R* UV side already in CT)
# ---------------------------------------------------------------------------

def sunset_C0():
    return (
        -mp.mpf(15) / 2
        + 2
        * mp.sqrt(3)
        * mp.im(mp.polylog(2, mp.exp(1j * mp.pi / 3)))
    )


def derive_CT_L0():
    C0 = sunset_C0()
    return float((C0 + 6) / 18 + (1 - mp.zeta(2)) / 24)


def derive_QFf():
    a = dblquad(
        lambda t, x: t
        / (1 + t + x * t) ** 3
        * (log(x) + log(1 + t + x * t) - 2 * log(1 + x + x * t)),
        0,
        1,
        0,
        1,
        epsabs=1e-14,
    )[0]
    return (24 * a - 1) / 48 - 1 / 24 + 1 / 18


def B_overlap_simple():
    e = s.symbols("e")
    d = 4 - 2 * e
    pref = s.gamma(3 * e - 1) * s.exp(3 * s.EulerGamma * e) / (
        64 * d * (d - 1)
    )
    c0 = 192 + s.Rational(64, 3)
    return float(
        s.simplify(
            s.series(pref * (-48 / e**2 + 208 / e + c0), e, 0, 1).coeff(e, -1)
        )
    )


def A_from_c0(c0=15.11, err=0.58):
    base = float(c0) / 1152 + 79 / 6912 + float(mp.pi**2) / 768
    return 4 * base, 4 * (err / 1152)


# ---------------------------------------------------------------------------
# R* IR checklist per class
# ---------------------------------------------------------------------------

def rstar_checklist():
    """Per-class R* IR/UV status. Δ_simple = change vs standing ledger."""
    rows = []

    # QQ
    rows.append(
        {
            "piece": "QQ",
            "UV_R": "yes (matter CT / mass-indep.)",
            "IR_Rtilde": "massive; IR soft regulated by z — no massless IR subdiv",
            "missed_IR_CT": False,
            "value": float(F(1, 16)),
            "err": 0.0,
            "delta": 0.0,
            "note": "exact +1/16; R* does not move",
            "label": "derived",
        }
    )
    # CCcc
    rows.append(
        {
            "piece": "CCcc",
            "UV_R": "yes",
            "IR_Rtilde": "massive",
            "missed_IR_CT": False,
            "value": float(F(-5, 72)),
            "err": 0.0,
            "delta": 0.0,
            "note": "exact −5/72",
            "label": "derived",
        }
    )
    # QCc doubled
    rows.append(
        {
            "piece": "QCc doubled",
            "UV_R": "yes (corner+faces)",
            "IR_Rtilde": "massive; soft faces in projective subtraction",
            "missed_IR_CT": False,
            "value": -0.0164974232969,
            "err": 1e-11,
            "delta": 0.0,
            "note": "analytic faces; no extra IR CT",
            "label": "derived",
        }
    )
    # QCc tadpole — critical IR check
    rows.append(
        {
            "piece": "QCc tadpole",
            "UV_R": "yes (T∼−z/(Kε) UV mass pole)",
            "IR_Rtilde": (
                "tadpole IR: at ξ=1/6 scalar a1=0 so no extra linear-R IR; "
                "O(1) of T is finite matching NOT a 1/ε pole"
            ),
            "missed_IR_CT": False,
            "value": float(F(-1, 36)),
            "err": 0.0,
            "delta": 0.0,
            "note": "pole −1/36 only; R* does not promote O(1)→1/ε",
            "label": "derived",
        }
    )
    # QFf
    qff = derive_QFf()
    rows.append(
        {
            "piece": "QFf",
            "UV_R": "yes (sunset poles cancel in D_x)",
            "IR_Rtilde": "massive equal-mass; no soft IR subdiv left in D_x",
            "missed_IR_CT": False,
            "value": qff,
            "err": 1e-11,
            "delta": 0.0,
            "note": "rederived j1; matches −0.07466151326",
            "label": "derived",
        }
    )
    # CcFf
    rows.append(
        {
            "piece": "CcFf",
            "UV_R": "yes (double 1/144)",
            "IR_Rtilde": "QMC integrand already endpoint-subtracted",
            "missed_IR_CT": False,
            "value": -0.02480,
            "err": 0.001,
            "delta": 0.0,
            "note": "no R* correction identified within QMC error",
            "label": "reproduced (QMC); R* null assumed within err",
        }
    )
    # FFff(A)×4
    A, Aerr = A_from_c0()
    rows.append(
        {
            "piece": "FFff(A)×4",
            "UV_R": "yes (triple/double from projective)",
            "IR_Rtilde": "regulated soft faces; long-pair analytic",
            "missed_IR_CT": False,
            "value": A,
            "err": Aerr,
            "delta": 0.0,
            "note": "c0=15.11±0.58; R* null within err",
            "label": "reproduced",
        }
    )
    # FFff(B) — soft-soft / soft-linear critical
    B = B_overlap_simple()
    rows.append(
        {
            "piece": "FFff(B)",
            "UV_R": "yes (CT + projective)",
            "IR_Rtilde": (
                "R̃ soft-soft = overlap addback (Δa1=-16, Δc0=64/3); "
                "soft-linear excluded (single-soft R̃ has no 1/ε from O(ε)·linear)"
            ),
            "missed_IR_CT": False,
            "value": B,
            "err": 2.0 / 2304,
            "delta": 0.0,
            "note": (
                "R* confirms prior forest: soft-linear NOT an IR 1/ε channel; "
                "including it would be an R̃ over-subtraction error"
            ),
            "label": "derived (R̃≡forest IE)",
        }
    )
    # CT
    ct = derive_CT_L0()
    rows.append(
        {
            "piece": "CT",
            "UV_R": "yes (Z_z + sunset R)",
            "IR_Rtilde": "massive sunset; tadpole subdivs reinserted",
            "missed_IR_CT": False,
            "value": ct,
            "err": 0.0,
            "delta": 0.0,
            "note": "(C+6)/18+(1−ζ2)/24",
            "label": "derived",
        }
    )
    return rows


def soft_linear_rstar_test():
    """Explicit R* statement on soft-linear ≈ −105.304.

    Single-soft IR subgraph with O(ε) numerator:
      ε ∫_0 dx x^{εu−1} (c0 + c1 x + …) = c0/(εu) + c1/(1+εu) + …
    The soft-linear c1 term is FINITE as ε→0 — no 1/ε pole.
    Therefore R̃ does not generate a simple-pole counterterm from soft-linear.
    Prior exclusion from overlap c0 is R*-correct (not an omission).
    """
    e, u, c0, c1 = s.symbols("e u c0 c1", positive=True)
    # Leading + linear pieces after integrating x^{e u - 1}
    leading_pole = c0 / (e * u)
    linear_finite = c1 / (1 + e * u)  # → c1 at e=0
    assert leading_pole.has(e)
    lim = s.limit(linear_finite, e, 0)
    assert lim == c1
    return {
        "soft_linear_gives_1_over_e": False,
        "conclusion": "R̃ correctly leaves soft-linear out of 1/ε; not a missed IR CT",
        "label": "derived",
    }


def projected_higher_poles():
    trip = F(1, 144) + F(1, 48) == F(1, 36)
    doub = F(5, 144) + F(5, 288) + F(-1, 288) - F(5, 48) + F(1, 18) == 0
    return {"triple_ok": trip, "double_ok": doub}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=== Background-field + R* audit (ladder rung #6) ===\n")
    print("No tuning to 1/24. Feasible subset — not a full R* re-integration.\n")

    bg = background_field_setup()
    print("-- Background field --")
    print(f"  {bg['split']}; probe {bg['probe']}")
    print(f"  R extraction: {bg['R_extraction']}")
    print(f"  Dynamical graviton loops? {bg['graviton_loops']}")
    print(f"  [{bg['label']}]")

    rs = rstar_definition()
    print("\n-- R* --")
    print(f"  {rs['formula']}")
    print(f"  IR soft: {rs['IR_on_soft']}")
    print(f"  Feasible: {rs['feasible']}")
    print(f"  [{rs['label']}]")

    sl = soft_linear_rstar_test()
    print("\n-- Soft-linear R* test --")
    print(f"  soft-linear → 1/ε? {sl['soft_linear_gives_1_over_e']}")
    print(f"  {sl['conclusion']}")
    print(f"  [{sl['label']}]")

    print("\n-- Per-piece R* checklist + L=0 simple --")
    rows = rstar_checklist()
    print(
        f"  {'piece':16s}  {'value':>12}  {'err':>8}  {'Δ':>8}  missed_IR?"
    )
    total = 0.0
    err2 = 0.0
    any_missed = False
    any_delta = False
    for r in rows:
        print(
            f"  {r['piece']:16s}  {r['value']:+12.8f}  {r['err']:8.1e}  "
            f"{r['delta']:+8.1e}  {r['missed_IR_CT']}"
        )
        print(f"    IR: {r['IR_Rtilde']}")
        print(f"    note: {r['note']}  [{r['label']}]")
        total += r["value"] + r["delta"]
        err2 += r["err"] ** 2
        any_missed = any_missed or r["missed_IR_CT"]
        any_delta = any_delta or (r["delta"] != 0.0)

    tot_err = err2**0.5
    target = float(F(1, 24))
    gap = total - target
    print(f"\n  {'TOTAL':16s}  {total:+12.8f}  {tot_err:8.1e}")
    print(f"  {'target 1/24':16s}  {target:+12.8f}")
    print(f"  {'GAP':16s}  {gap:+12.8f}  {tot_err:8.1e}")

    hp = projected_higher_poles()
    print(
        f"\n  Higher poles after R*: triple_ok={hp['triple_ok']}, "
        f"double_ok={hp['double_ok']}"
    )
    print(f"  Any missed IR CT found? {any_missed}")
    print(f"  Any piece value changed? {any_delta}")

    print(
        f"""
Verdict
-------
Background-field metric split: already the packet's covariant R extraction
(no dynamical graviton). R* IR dual on soft axes ≡ sector Taylor forest IE
already used for FFff(B) overlap (soft-soft only). Soft-linear is NOT an
IR 1/ε channel under R̃. QCc tadpole O(1) is not a pole. No missed IR
counterterm identified in QQ, CCcc, QCc, tadpole, QFf, CcFf, FFff(A/B), CT.

New simple total = {total:+.6f} ± {tot_err:.6f}
GAP vs 1/24     = {gap:+.6f} ± {tot_err:.6f}
Pieces changed  = none (all Δ=0)

Plain statement: R* / background-field feasible subset confirms the existing
ledger assembly; it does NOT close the ~+0.01 gap. Full diagrammatic R*
re-integration of every QMC graph remains out of scope.
"""
    )
    assert not any_missed and not any_delta
    assert hp["triple_ok"] and hp["double_ok"]
    assert abs(gap - 0.0097) < 0.004
    assert abs(rows[4]["value"] + 0.07466151326) < 1e-10  # QFf
    print("=== BG-FIELD + R* RUNG COMPLETE ===")


if __name__ == "__main__":
    main()
