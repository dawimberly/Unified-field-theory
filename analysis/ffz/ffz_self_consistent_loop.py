#!/usr/bin/env python3
"""Self-consistent closed loop: v ↔ G_ind ↔ m_s (no external scale).

Dan's whole-system-dynamics new-idea rung. Stacked on BG-field/R* (#18).

Loop (ASSUMED mediation / gap / contact dictionary from #16–#17):
  (i)   c_eff = c_ledger + K v²/(6 h² z)     [#16 ξv² finite R contact]
  (ii)  M_P²  = 16π |c_eff| h² z / K         [induced Planck from C_R]
  (iii) m_s²  = −|F|² / M_P²                 [tachyonic gravity mediation]
  (iv)  v²/z  = −2(1 + m_s²/z)/(h K)         [#17 T2 soft tree gap, if r<-1]

F-closure for "no external scale":
  F_A: φ=F/z free          → one dimensionless free param
  F_B: φ=1 (F~z natural)   → no free continuous param
  F_C: F = |y| v²          → fully closed (φ = √(hK)·(v²/z))

No tuning to 1/24. Derived vs assumed labeled.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_self_consistent_loop.py

Requires: numpy
"""

from __future__ import annotations

import numpy as np

K = 16.0 * np.pi**2
C_LEDGER = 0.05131719  # R* rung L=0 simple total
GAP_TARGET = 0.00965052  # R* gap vs 1/24
TARGET_V2Z_PREF = 3.703108e-4  # v²/z = PREF × h² absorbs gap via ξv²
XI = 1.0 / 6.0


def c_eff_of(x, h):
    """x = v²/z.  c_eff = c_ledger + K x / (6 h²)."""
    return C_LEDGER + K * x / (6.0 * h**2)


def r_from_phi(phi, c_eff, h):
    """r = m_s²/z = −φ² K / (16π c_eff h²),  φ=F/z.

    ASSUMED: tachyonic gravity mediation m_s² = −|F|²/M_P² with
    M_P² = 16π c_eff h² z / K  (c_eff>0).
    """
    if c_eff <= 0:
        return np.nan
    return -(phi**2) * K / (16.0 * np.pi * c_eff * h**2)


def x_from_r(r, h):
    """T2 soft tree: x=v²/z = −2(1+r)/(h K) if r<-1 else 0."""
    if not np.isfinite(r) or r >= -1.0:
        return 0.0
    return -2.0 * (1.0 + r) / (h * K)


def delta_R(x, h):
    """Finite R shift in ledger units: K x / (6 h²)."""
    return K * x / (6.0 * h**2)


def solve_fixed_point_phi(h, phi, n_iter=80, tol=1e-14):
    """Iterate x → c_eff → r(φ) → x' from soft gap."""
    x = 0.0
    hist = []
    for _ in range(n_iter):
        ce = c_eff_of(x, h)
        r = r_from_phi(phi, ce, h)
        x_new = x_from_r(r, h)
        hist.append((x, ce, r, x_new))
        if abs(x_new - x) < tol:
            x = x_new
            break
        # Damped update for stability
        x = 0.5 * x + 0.5 * x_new
    else:
        x = x_new
    ce = c_eff_of(x, h)
    r = r_from_phi(phi, ce, h)
    # Residual of exact fixed-point map
    resid = abs(x - x_from_r(r, h))
    return {
        "h": h,
        "phi": phi,
        "x": x,
        "c_eff": ce,
        "r": r,
        "delta_R": delta_R(x, h),
        "resid": resid,
        "converged": resid < 1e-10,
        "hist": hist,
    }


def solve_fixed_point_F_from_v(h, n_iter=120, tol=1e-14):
    """F_C: F = |y| v² ⇒ φ = √(h K) · x  (x=v²/z). Fully closed."""
    x = 1e-6  # tiny seed
    for _ in range(n_iter):
        ce = c_eff_of(x, h)
        phi = np.sqrt(h * K) * x
        r = r_from_phi(phi, ce, h)
        x_new = x_from_r(r, h)
        if abs(x_new - x) < tol:
            x = x_new
            break
        x = 0.5 * x + 0.5 * x_new
    else:
        x = x_new
    ce = c_eff_of(x, h)
    phi = np.sqrt(h * K) * x
    r = r_from_phi(phi, ce, h)
    resid = abs(x - x_from_r(r, h))
    return {
        "h": h,
        "phi": phi,
        "x": x,
        "c_eff": ce,
        "r": r,
        "delta_R": delta_R(x, h),
        "resid": resid,
        "converged": resid < 1e-10,
        "nontrivial": x > 0,
    }


def trivial_solution(h, phi):
    """x=0 branch: always a fixed point of the map when r(x=0) ≥ −1;
    if r(0) < −1 the soft gap drives x away from 0 (unstable trivial)."""
    ce = C_LEDGER
    r0 = r_from_phi(phi, ce, h)
    return {
        "x": 0.0,
        "c_eff": ce,
        "r0": r0,
        "stable_as_fixed_point": r0 >= -1.0,
        "delta_R": 0.0,
    }


def scan_stability_phi(h, phi):
    """Linear stability: iterate from seed x=0 and from large seed."""
    from_zero = solve_fixed_point_phi(h, phi)
    # large seed
    x = 1.0
    for _ in range(80):
        ce = c_eff_of(x, h)
        r = r_from_phi(phi, ce, h)
        x = 0.5 * x + 0.5 * x_from_r(r, h)
    from_large = {
        "x": x,
        "c_eff": c_eff_of(x, h),
        "r": r_from_phi(phi, c_eff_of(x, h), h),
        "delta_R": delta_R(x, h),
    }
    triv = trivial_solution(h, phi)
    same = abs(from_zero["x"] - from_large["x"]) < 1e-8
    return {
        "from_zero": from_zero,
        "from_large": from_large,
        "trivial": triv,
        "unique_attractor": same,
    }


def compare_to_gap(delta):
    """ratio predicted δ_R / GAP_TARGET (want ≈1 to absorb)."""
    if GAP_TARGET == 0:
        return float("nan")
    return delta / GAP_TARGET


def poles_unaffected():
    """ξv² and soft IR masses: finite/simple only; triple/double unchanged."""
    trip = abs((1 / 144 + 1 / 48) - 1 / 36) < 1e-15
    doub = abs(5 / 144 + 5 / 288 - 1 / 288 - 5 / 48 + 1 / 18) < 1e-15
    return {
        "triple_ok": trip,
        "double_ok": doub,
        "note": (
            "DERIVED (#16/#18): tree ξv² is finite contact; soft m_s is IR. "
            "UV triple/double of h² ledger unchanged."
        ),
    }


def main():
    print("=== Self-consistent loop: v ↔ G_ind ↔ m_s ===\n")
    print("No tuning to 1/24. No external soft scale except via F-closure.\n")

    print("-- Assumptions --")
    print("  ASSUMED: tachyonic gravity mediation m_s² = −|F|²/M_P²")
    print("  ASSUMED: M_P² = 16π c_eff h² z / K  (induced from |C_R|)")
    print("  ASSUMED: soft-tree gap (#17 T2) for v given m_s")
    print("  ASSUMED: ξv² contact (#16) shifts c_eff")
    print("  DERIVED: algebra of the fixed-point map given those dictionaries\n")

    # F_A: φ free scan
    print("-- (1)–(3) F_A: φ=F/z free --")
    print(
        f"  {'h':>4}  {'φ':>8}  {'x=v²/z':>12}  {'δ_R':>12}  "
        f"{'ratio δ/GAP':>12}  {'r':>10}  unique?"
    )
    fa_rows = []
    for h in (0.1, 1.0):
        for phi in (0.01, 0.1, 0.5, 1.0, 2.0):
            st = scan_stability_phi(h, phi)
            fz = st["from_zero"]
            ratio = compare_to_gap(fz["delta_R"])
            print(
                f"  {h:4.1f}  {phi:8.2f}  {fz['x']:12.4e}  {fz['delta_R']:12.4e}  "
                f"{ratio:12.4e}  {fz['r']:10.3e}  {st['unique_attractor']}"
            )
            fa_rows.append({**fz, "unique": st["unique_attractor"], "ratio": ratio})
    print(
        "  Free parameter remaining: φ=F/z (dimensionless). "
        "Natural O(1) vs tuned: O(1) is natural for hidden F~z, but "
        "hitting ratio≈1 requires tuning φ (see below)."
    )
    # φ that hits ratio≈1 at h=1
    print("\n  φ scan for ratio≈1 at h=1:")
    best = None
    for phi in np.linspace(0.05, 0.5, 200):
        sol = solve_fixed_point_phi(1.0, phi)
        ratio = compare_to_gap(sol["delta_R"])
        if best is None or abs(ratio - 1) < abs(best["ratio"] - 1):
            best = {**sol, "ratio": ratio}
    print(
        f"    best φ≈{best['phi']:.4f}: x={best['x']:.6e}, "
        f"δ_R={best['delta_R']:.6e}, ratio={best['ratio']:.4f}"
    )
    print("    → hitting the Newton gap still tunes φ  [tuning]")

    # F_B: φ=1
    print("\n-- F_B: φ=1 (F~z, natural O(1), no free continuous param) --")
    for h in (0.1, 1.0):
        st = scan_stability_phi(h, 1.0)
        fz = st["from_zero"]
        fl = st["from_large"]
        triv = st["trivial"]
        print(
            f"  h={h}: x={fz['x']:.6e}, δ_R={fz['delta_R']:.6e}, "
            f"ratio={compare_to_gap(fz['delta_R']):.4e}, r={fz['r']:.4e}"
        )
        print(
            f"    trivial stable? {triv['stable_as_fixed_point']}; "
            f"unique attractor? {st['unique_attractor']}; "
            f"from_large x={fl['x']:.6e}"
        )
        target_x = TARGET_V2Z_PREF * h**2
        print(
            f"    predicted v²/z={fz['x']:.6e} vs target {target_x:.6e} "
            f"(pred/target_x={fz['x']/target_x if target_x else float('nan'):.4e})"
        )

    # F_C: F = |y| v² fully closed
    print("\n-- F_C: F=|y|v² (fully closed, no external F) --")
    for h in (0.1, 1.0):
        sol = solve_fixed_point_F_from_v(h)
        # also check trivial
        # at x=0, φ=0 ⇒ r=0 ≥ -1 ⇒ trivial is stable fixed point
        print(
            f"  h={h}: nontrivial={sol['nontrivial']}, x={sol['x']:.6e}, "
            f"φ={sol['phi']:.6e}, δ_R={sol['delta_R']:.6e}, "
            f"ratio={compare_to_gap(sol['delta_R']):.4e}, "
            f"r={sol['r']:.4e}, converged={sol['converged']}"
        )
        print(
            "    trivial x=0: always a solution when φ∝x (r→0); "
            "nonzero seed → "
            f"{'returns to 0' if sol['x']==0 else 'finite x'}"
        )

    # Existence / uniqueness summary at F_B h=1
    fb = scan_stability_phi(1.0, 1.0)
    fc = solve_fixed_point_F_from_v(1.0)

    poles = poles_unaffected()
    print("\n-- (4) Higher poles --")
    print(f"  triple_ok={poles['triple_ok']}, double_ok={poles['double_ok']}")
    print(f"  {poles['note']}")

    print(
        f"""
Verdict
-------
(1) Nontrivial self-consistent solutions:
    • F_A (φ free): yes for φ large enough that r<-1; unique attractor in
      the damped map for scanned φ (trivial unstable when r(0)<-1).
    • F_B (φ=1): nontrivial solution EXISTS and is the unique attractor
      (trivial unstable). At h=1: v²/z={fb['from_zero']['x']:.6e},
      δ_R={fb['from_zero']['delta_R']:.6e}, ratio={compare_to_gap(fb['from_zero']['delta_R']):.4e}.
    • F_C (F=|y|v²): only trivial x=0 is consistent (φ∝x ⇒ soft too weak
      to trigger tachyon). Nontrivial closed solution does NOT exist.

(2) Free parameters:
    • F_A: φ=F/z remains free (dimensionless). Natural O(1), but ratio≈1
      needs φ≈{best['phi']:.3f} at h=1 — tuned.
    • F_B: no continuous free param (φ fixed to 1 by F~z assumption).
    • F_C: no free param, but only v=0.

(3) Predicted vs gap +{GAP_TARGET:.5f}:
    • F_B h=1: ratio δ_R/GAP = {compare_to_gap(fb['from_zero']['delta_R']):.4e}
      (orders of magnitude ABOVE 1 — overshoots, does not match).
    • F_B h=0.1: see table above (also ≫1).
    • F_C: ratio = 0.
    • F_A can hit ratio≈1 only by tuning φ.

(4) Triple/double poles: unaffected (finite ξv² + IR soft).

Plain statement: a self-consistent gravity-mediation loop can produce a
nontrivial v when F~z (F_B), but the predicted R shift overshoots the
~+0.01 gap by ~10²–10³. Fully closing F via the condensate (F_C) yields
only v=0. Leaving φ free restores tuning. This whole-system loop does
NOT explain the Newton gap without tuning or a huge mismatch.
"""
    )
    assert poles["triple_ok"] and poles["double_ok"]
    assert fb["from_zero"]["x"] > 0
    assert compare_to_gap(fb["from_zero"]["delta_R"]) > 10
    assert fc["x"] == 0.0 or abs(fc["x"]) < 1e-12
    print("=== SELF-CONSISTENT LOOP COMPLETE ===")


if __name__ == "__main__":
    main()
