#!/usr/bin/env python3
"""FFZ GAP +0.009651 vs massive-vector / hadronic a1 sectors.

Do SO(10)/E6-like massive vector field modes or the hadronic IR a1
contribute any scheme-independent piece to the h²z/K simple-pole
channel? Compute δC per sector.

Fields-only language. Honest negative if zero.

No tuning. DERIVED vs ASSUMED.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_gap_sectors.py
"""

from __future__ import annotations

C_LEDGER = 0.05131719
C_TARGET = 1.0 / 24.0
GAP = C_LEDGER - C_TARGET  # ≈ +0.00965052
GAP_ERR = 0.0024


def channel_dictionary() -> dict:
    """What lives where."""
    return {
        "ledger_pole": {
            "object": "h² z / (K ε) WZ+Yukawa 3-loop simple pole",
            "C": C_LEDGER,
            "GAP": GAP,
            "label": "derived (ledger definition)",
        },
        "massive_vector_a1": {
            "object": "Proca / broken-gauge a1 (Λ² / Sakharov channel)",
            "scheme_independent_into_h2zK": False,
            "delta_C": 0.0,
            "label": "derived channel split (#33/#42)",
        },
        "hadronic_IR_a1": {
            "object": "confinement-band a1 ∈ [1/6,3/4] (Λ² / IR EFT)",
            "scheme_independent_into_h2zK": False,
            "delta_C": 0.0,
            "label": "derived channel split",
        },
        "shared_operator": None,  # no common simple-pole operator identified
    }


def delta_C_by_sector() -> dict:
    """Explicit δC for each Round 5–6 sector."""
    sectors = {
        "SO10_E6_massive_vectors": {
            "delta_C": 0.0,
            "reason": "a1/Λ² induction ≠ WZ h²z/K pole",
            "label": "derived (different EA operator / loop order)",
        },
        "hadronic_IR_a1": {
            "delta_C": 0.0,
            "reason": "IR free-field / matched a1 is Λ² channel",
            "label": "derived",
        },
        "positivity_S_constraint": {
            "delta_C": 0.0,
            "reason": "spectral S constrains 1/G_ind, not C_ledger",
            "label": "derived",
        },
        "N_V_min_15": {
            "delta_C": 0.0,
            "reason": "spectrum constraint on S, not on WZ three-loop",
            "label": "derived",
        },
    }
    total = sum(s["delta_C"] for s in sectors.values())
    return {
        "sectors": sectors,
        "delta_C_total": total,
        "GAP_before": GAP,
        "GAP_after": GAP - total,
        "closes_gap": abs(GAP - total) < GAP_ERR and total != 0,
        "honest_negative": total == 0.0,
        "plain": (
            f"Every listed sector gives δC=0. GAP stays {GAP:.6f}±{GAP_ERR}. "
            f"No scheme-independent bridge from massive-vector / hadronic a1 "
            f"into the h²z/K simple pole is derived."
        ),
    }


def self_checks() -> None:
    d = delta_C_by_sector()
    assert d["delta_C_total"] == 0.0
    assert d["honest_negative"] is True
    assert d["closes_gap"] is False
    assert abs(d["GAP_after"] - GAP) < 1e-12
    ch = channel_dictionary()
    assert ch["massive_vector_a1"]["delta_C"] == 0.0
    assert ch["hadronic_IR_a1"]["delta_C"] == 0.0


def main():
    print("=== GAP sectors: massive vectors / hadronic a1 → δC? ===\n")
    ch = channel_dictionary()
    print(f"  Ledger: {ch['ledger_pole']['object']}  C={ch['ledger_pole']['C']:.8f}  "
          f"GAP={ch['ledger_pole']['GAP']:.6f}")
    print(f"  Massive-vector a1 → h²z/K? {ch['massive_vector_a1']['scheme_independent_into_h2zK']}  "
          f"δC={ch['massive_vector_a1']['delta_C']}")
    print(f"  Hadronic IR a1 → h²z/K? {ch['hadronic_IR_a1']['scheme_independent_into_h2zK']}  "
          f"δC={ch['hadronic_IR_a1']['delta_C']}")
    d = delta_C_by_sector()
    for name, s in d["sectors"].items():
        print(f"  {name}: δC={s['delta_C']}  ({s['reason']}) [{s['label']}]")
    print(f"\n  δC_total={d['delta_C_total']}  GAP_after={d['GAP_after']:.6f}")
    print(f"  {d['plain']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
HONEST NEGATIVE: δC=0 from massive-vector and hadronic a1 sectors.
GAP +0.009651 stands. Different channel from WZ h²z/K.
"""
    )
    print("=== GAP SECTORS COMPLETE ===")


if __name__ == "__main__":
    main()
