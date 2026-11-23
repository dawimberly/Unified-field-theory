#!/usr/bin/env python3
"""Scale-dependent 1/G_ind(μ): sign crossover across SM → hadronic scales.

Uses free-field / controlled a1 at each threshold. Flag coincidences
with Λ_QCD and m_N.

No tuning. DERIVED vs ASSUMED.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_running_gind.py
"""

from __future__ import annotations

from fractions import Fraction as F

# a1 proxies as "sign of 1/G_ind" (overall positive prefactor ASSUMED Sakharov)
LAMBDA_QCD_GEV = 0.3
M_N_GEV = 0.938


def running_table() -> list[dict]:
    """a1(μ) waypoints (from #31/#G)."""
    return [
        {"mu": 1e16, "name": "UV_SM", "a1": F(-43, 12), "label": "derived"},
        {"mu": 50.0, "name": "EW_Proca", "a1": F(-1, 12), "label": "derived"},
        {"mu": 2.0, "name": "AC_above_QCD", "a1": F(-13, 4), "label": "derived"},
        {
            "mu": LAMBDA_QCD_GEV,
            "name": "confinement",
            "a1": F(1, 6),  # low end of controlled band
            "a1_high": F(7, 6),
            "label": "controlled band (#G); matching assumed",
        },
        {"mu": M_N_GEV, "name": "nucleon", "a1": F(1, 6), "label": "same IR band"},
        {"mu": 0.001, "name": "QED_only", "a1": F(-1, 6), "label": "derived; if hadrons frozen out"},
    ]


def crossover() -> dict:
    rows = running_table()
    # First scale where a1>0 when scanning down
    cross = None
    for r in rows:
        a1 = r["a1"]
        if a1 > 0 and cross is None:
            cross = r
    return {
        "crossover": cross,
        "at_Lambda_QCD": cross is not None and abs(cross["mu"] - LAMBDA_QCD_GEV) < 1e-9,
        "coincidence_flags": [
            "Crossover tied to confinement by construction of the IR band — not an independent prediction",
            "Proximity to m_N is the same IR regime, not a new coincidence claim",
        ],
        "QED_only_negative": True,
        "plain": (
            "Sign becomes positive at/below confinement in the controlled "
            "hadronic counting; negative again in a pure QED+lepton theory "
            "without light hadrons. Crossover = Λ_QCD by construction of "
            "when gluon Kabats drop — flag as structural, not miraculous."
        ),
    }


def self_checks() -> None:
    rows = running_table()
    assert rows[0]["a1"] < 0 and rows[2]["a1"] < 0
    assert rows[3]["a1"] > 0
    c = crossover()
    assert c["crossover"]["name"] == "confinement"
    assert c["at_Lambda_QCD"] is True


def main():
    print("=== Running 1/G_ind(μ) sign ===\n")
    for r in running_table():
        hi = r.get("a1_high")
        extra = f" .. {hi}" if hi is not None else ""
        print(f"  μ={r['mu']:.3g} GeV  {r['name']:16s}  a1={r['a1']}{extra}  [{r['label']}]")
    c = crossover()
    print(f"\n  Crossover: {c['crossover']['name']} at μ={c['crossover']['mu']} GeV")
    print(f"  Equals Λ_QCD by construction? {c['at_Lambda_QCD']}")
    for f in c["coincidence_flags"]:
        print(f"  coincidence: {f}")
    print(f"  {c['plain']}")
    self_checks()
    print("  self-checks PASS")
    print("=== RUNNING G_IND COMPLETE ===")


if __name__ == "__main__":
    main()
