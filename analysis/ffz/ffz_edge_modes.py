#!/usr/bin/env python3
"""AH: Donnelly–Wall-type edge modes of the *original* FFZ fields.

After AG (surface/TD → δC=0), ask whether edge-mode / contact contributions
of the WZ field content itself (complex scalar ξ=1/6 + Weyl fermion; no
gauge) can shift the h²z/K simple pole by GAP=+0.009651.

Compare any nonzero δC to GAP and to scheme κ=−GAP from #23.
Flag ln2/72 numerical coincidence.

No tuning. DERIVED vs ASSUMED. Fields-only. Honest negatives.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_edge_modes.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import log

C_LEDGER = 0.05131719
C_TARGET = 1.0 / 24.0
GAP = C_LEDGER - C_TARGET  # ≈ +0.00965052
GAP_ERR = 0.0024
KAPPA_SCHEME = -GAP  # #23
LN2_OVER_72 = log(2) / 72.0


def original_fields() -> dict:
    return {
        "complex_scalar_xi": F(1, 6),
        "weyl_fermion": True,
        "massless_gauge": False,
        "label": "derived (FFZ Newton test content)",
    }


def scalar_edge_contact() -> dict:
    """Solodukhin / EE for nonminimal scalar.

    Area-law quad ∝ (1/6 − ξ). At ξ=1/6 → 0.
    No Kabat-type gauge contact for scalars (Kabat is spin-1).
    DW edge-mode story for Maxwell does not apply to WZ scalar.
    """
    xi = F(1, 6)
    area_quad = F(1, 6) - xi
    assert area_quad == 0
    return {
        "field": "complex scalar (WZ)",
        "a1_bulk_xi_piece": area_quad,
        "kabat_contact": 0,
        "DW_edge_applicable": False,
        "delta_C_h2zK": 0.0,
        "label": "derived (Solodukhin EE ∝ 1/6−ξ; no scalar Kabat)",
        "dw_reading": "assumed N/A — DW is gauge edge flux",
    }


def fermion_edge_contact() -> dict:
    """Weyl / Majorana: Solodukhin — spin-1/2 has no Kabat contact."""
    return {
        "field": "Weyl fermion (WZ)",
        "kabat_contact": 0,
        "DW_edge_applicable": False,
        "delta_C_h2zK": 0.0,
        "label": "derived (no spin-1/2 contact)",
    }


def channel_split() -> dict:
    """Even a nonzero free-field contact lives in Λ² / a1, not h²z/K."""
    return {
        "edge_a1_channel": "Λ² / Sakharov 1/G",
        "ledger_channel": "h² z / K three-loop simple pole",
        "shared_operator": None,
        "delta_C_from_channel_map": 0.0,
        "label": "derived channel split (#33/#42/#54)",
    }


def compare_to_gap_and_kappa(delta_C: float) -> dict:
    coin = {
        "name": "ln2/72",
        "value": LN2_OVER_72,
        "diff_to_GAP": LN2_OVER_72 - GAP,
        "within_1sigma": abs(LN2_OVER_72 - GAP) <= GAP_ERR,
        "label": "NOT derived — numerical coincidence only (#20/#22)",
    }
    return {
        "delta_C": delta_C,
        "GAP": GAP,
        "kappa_scheme_23": KAPPA_SCHEME,
        "delta_C_equals_GAP": abs(delta_C - GAP) <= 1e-15 and delta_C != 0,
        "delta_C_equals_minus_kappa": abs(delta_C - (-KAPPA_SCHEME)) <= 1e-15
        and delta_C != 0,
        "closes_gap": abs(GAP - delta_C) <= GAP_ERR and delta_C != 0,
        "coincidence_ln2_72": coin,
        "plain": (
            f"δC={delta_C}. GAP={GAP:+.8f}. κ(#23)={KAPPA_SCHEME:+.8f}. "
            f"ln2/72={LN2_OVER_72:.8f} within 1σ of GAP — coincidence only."
        ),
    }


def edge_mode_result() -> dict:
    s = scalar_edge_contact()
    f = fermion_edge_contact()
    ch = channel_split()
    total = s["delta_C_h2zK"] + f["delta_C_h2zK"] + ch["delta_C_from_channel_map"]
    cmp_ = compare_to_gap_and_kappa(total)
    return {
        "scalar": s,
        "fermion": f,
        "channel": ch,
        "delta_C_total": total,
        "compare": cmp_,
        "honest_negative": total == 0.0,
        "plain": (
            "Original FFZ fields have no DW/Kabat edge contact. "
            "ξ=1/6 kills scalar area quad. Channel split ⇒ δC=0 into h²z/K. "
            f"GAP {GAP:.6f} stands. κ=−GAP remains scheme (#23)."
        ),
    }


def self_checks() -> None:
    r = edge_mode_result()
    assert r["delta_C_total"] == 0.0
    assert r["honest_negative"] is True
    assert r["compare"]["closes_gap"] is False
    assert scalar_edge_contact()["a1_bulk_xi_piece"] == 0
    assert scalar_edge_contact()["DW_edge_applicable"] is False
    assert fermion_edge_contact()["kabat_contact"] == 0
    assert abs(KAPPA_SCHEME + GAP) < 1e-15
    coin = r["compare"]["coincidence_ln2_72"]
    assert coin["within_1sigma"] is True  # known coincidence


def main():
    print("=== AH: Edge modes of original FFZ fields vs GAP ===\n")
    fc = original_fields()
    print(f"  fields: scalar ξ={fc['complex_scalar_xi']}, Weyl, gauge={fc['massless_gauge']}")

    r = edge_mode_result()
    print("\n-- Scalar --")
    s = r["scalar"]
    print(f"  a1 ξ-piece={s['a1_bulk_xi_piece']}  Kabat={s['kabat_contact']}")
    print(f"  DW applicable? {s['DW_edge_applicable']}  [{s['label']}]")

    print("\n-- Fermion --")
    f = r["fermion"]
    print(f"  Kabat={f['kabat_contact']}  [{f['label']}]")

    print("\n-- Channel --")
    ch = r["channel"]
    print(f"  edge channel: {ch['edge_a1_channel']}")
    print(f"  ledger: {ch['ledger_channel']}  δC_map={ch['delta_C_from_channel_map']}")

    print("\n-- vs GAP / κ --")
    c = r["compare"]
    print(f"  δC_total={c['delta_C']}")
    print(f"  GAP={c['GAP']:+.8f}  κ(#23)={c['kappa_scheme_23']:+.8f}")
    print(f"  closes gap? {c['closes_gap']}")
    coin = c["coincidence_ln2_72"]
    print(
        f"  coincidence {coin['name']}={coin['value']:.8f}  "
        f"diff_GAP={coin['diff_to_GAP']:+.2e}  within_1σ={coin['within_1sigma']}"
    )
    print(f"\n  {r['plain']}")

    print("\n=== EDGE MODES COMPLETE ===")
    return r


if __name__ == "__main__":
    main()
