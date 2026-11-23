#!/usr/bin/env python3
"""Low-energy IR roster with Appelquist–Carazzone decoupling.

Stacked on shared-boundary (#30). Track active fields by scale; net a1.
Tests Dan's 'mass switches on gravity'.

No tuning. DERIVED vs ASSUMED.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_ir_roster.py
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F

A1_S = F(1, 6)
A1_W = F(1, 12)
A1_V = F(-2, 3)
A1_P = F(1, 2)


@dataclass(frozen=True)
class ScaleRoster:
    name: str
    mu_GeV: float
    n_scalar: int
    n_weyl: int
    n_massless_V: int
    n_proca: int
    note: str
    label: str

    @property
    def a1(self) -> F:
        return (
            self.n_scalar * A1_S
            + self.n_weyl * A1_W
            + self.n_massless_V * A1_V
            + self.n_proca * A1_P
        )


def corrected_ladder() -> list[ScaleRoster]:
    """Scale-by-scale free-field a1 with AC decoupling.

    Hadronic replacement of quarks+gluons by π,N is ASSUMED (not a
    controlled QCD heat-kernel theorem).
    """
    return [
        ScaleRoster(
            "UV_SM", 1e16, 4, 45, 12, 0,
            "full SM massless gauge", "derived",
        ),
        ScaleRoster(
            "EW_Proca", 50.0, 4, 45, 9, 3,
            "W±Z as Proca (#29)", "derived arithmetic",
        ),
        ScaleRoster(
            "AC_above_QCD", 2.0, 0, 33, 9, 0,
            "AC remove H,W,Z,t; γ+8g free; 33 Weyl remain",
            "derived count; threshold order assumed",
        ),
        ScaleRoster(
            "hadronic_proxy", 0.1, 3, 10, 1, 0,
            "γ + light leptons + π(3) + N(p,n); gluons confined",
            "ASSUMED free-field hadrons",
        ),
        ScaleRoster(
            "QED_leptons", 0.001, 0, 6, 1, 0,
            "γ + e,μ,νe,νμ only", "derived QED content",
        ),
    ]


def mass_switches_gravity_test() -> dict:
    rows = {r.name: r for r in corrected_ladder()}
    signs = {n: (r.a1 > 0) for n, r in rows.items()}
    return {
        "rows": rows,
        "signs": signs,
        "any_positive": any(signs.values()),
        "hadronic_positive": signs["hadronic_proxy"],
        "QED_positive": signs["QED_leptons"],
        "derived_positive": False,  # only assumed hadronic proxy is +
        "plain": (
            "Hadronic proxy a1=+2/3 is positive but ASSUMED (π,N as free HK). "
            "QED+leptons a1=−1/6. Mass-on does not derive positive induced G."
        ),
    }


def self_checks() -> None:
    rows = corrected_ladder()
    assert rows[0].a1 == F(-43, 12)
    assert rows[1].a1 == F(-1, 12)
    assert rows[2].a1 == F(33, 12) - 6 == F(-13, 4)
    assert rows[3].a1 == F(2, 3)
    assert rows[3].a1 > 0
    assert rows[4].a1 == F(-1, 6)
    t = mass_switches_gravity_test()
    assert t["hadronic_positive"] is True
    assert t["QED_positive"] is False
    assert t["derived_positive"] is False


def main():
    print("=== IR roster / Appelquist–Carazzone a1 ===\n")
    print("No tuning. Stacked on #30.\n")
    print(f"  {'scale':16s} {'μ/GeV':>10} {'a1':>10} sign  label")
    for r in corrected_ladder():
        print(
            f"  {r.name:16s} {r.mu_GeV:10.3g} {str(r.a1):>10} "
            f"{'OK' if r.a1 > 0 else 'NO':>4}  [{r.label}]"
        )
        print(f"    {r.note}")
    t = mass_switches_gravity_test()
    print(f"\n  Hadronic +? {t['hadronic_positive']}  QED +? {t['QED_positive']}  derived +? {t['derived_positive']}")
    print(f"  {t['plain']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
UV −43/12, EW-Proca −1/12, AC-above-QCD −13/4: all wrong-sign (derived).
Hadronic proxy +2/3 positive only under ASSUMED free π,N replacement —
flag as coincidence risk, not a derivation that mass switches gravity on.
QED+leptons −1/6 (photon Kabat remains).
"""
    )
    print("=== IR ROSTER COMPLETE ===")


if __name__ == "__main__":
    main()
