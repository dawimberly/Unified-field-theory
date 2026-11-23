#!/usr/bin/env python3
"""GUT spectra SU(5)/SO(10): does UV 'confinement-like' Higgsing flip a1?

Massless GUT gauge bosons (−2/3 Kabat) vs massive X/Y Proca (+1/2).
Compute a1 above/below M_GUT; ask whether a GUT spectrum flips the
integrated sign relative to #41/#45.

No tuning. DERIVED vs ASSUMED. Coincidences flagged.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_gut_spectrum.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import log

A1_V = F(-2, 3)
A1_PROCA = F(1, 2)
A1_W = F(1, 12)
A1_S = F(1, 6)

# From #41 / #45
MU_UV = 1e16
M_GUT = 1e16  # coincidence: often identified — FLAG
V_EW = 246.0
LAMBDA_QCD = 0.3
A1_AC = F(-13, 4)
A1_EW = F(-1, 12)
A1_HAD = F(1, 3)
A1_SM_UV = F(-43, 12)


def a1_gauge(n_massless: int, n_proca: int) -> F:
    return n_massless * A1_V + n_proca * A1_PROCA


def su5_spectra() -> dict:
    """SU(5) Georgi–Glashow.

    dim G = 24. Matter: 3×(10+5bar)=45 Weyl (same as SM, no ν_R).
    Higgs: minimal 5_H = 5 complex = 10 real scalars (ASSUMED count;
    24_H for adjoint breaking omitted or noted).
    """
    n_gauge = 24
    n_weyl = 45
    n_real_H = 10  # 5_H; ASSUMED minimal
    unbroken = {
        "n_massless": n_gauge,
        "n_proca": 0,
        "a1_gauge": a1_gauge(n_gauge, 0),
        "a1_ferm": n_weyl * A1_W,
        "a1_scalar": n_real_H * A1_S,
    }
    unbroken["a1_net"] = unbroken["a1_gauge"] + unbroken["a1_ferm"] + unbroken["a1_scalar"]
    # Broken to SM: 12 X/Y Proca, 12 SM massless
    broken = {
        "n_massless": 12,
        "n_proca": 12,
        "a1_gauge": a1_gauge(12, 12),
        "a1_ferm": n_weyl * A1_W,
        "a1_scalar": n_real_H * A1_S,  # residual; ASSUMED same
    }
    broken["a1_net"] = broken["a1_gauge"] + broken["a1_ferm"] + broken["a1_scalar"]
    assert unbroken["a1_gauge"] == F(-16)
    assert broken["a1_gauge"] == F(-2)  # −8 + 6
    return {
        "group": "SU(5)",
        "unbroken": unbroken,
        "broken_to_SM": broken,
        "delta_a1_from_XY_Proca": broken["a1_net"] - unbroken["a1_net"],
        "sign_OK_unbroken": unbroken["a1_net"] > 0,
        "sign_OK_broken": broken["a1_net"] > 0,
        "label_content": "derived (group dim + SM embedding); Higgs count assumed",
    }


def so10_spectra() -> dict:
    """SO(10).

    dim G = 45. Matter: 3×16 = 48 Weyl (includes ν_R).
    Higgs: ASSUMED omit detailed 10+126+210 — report gauge+fermion first,
    then optional 10_H = 10 complex = 20 real.
    """
    n_gauge = 45
    n_weyl = 48
    n_real_H = 20  # 10_H only; ASSUMED minimal slice
    unbroken = {
        "n_massless": n_gauge,
        "n_proca": 0,
        "a1_gauge": a1_gauge(n_gauge, 0),
        "a1_ferm": n_weyl * A1_W,
        "a1_scalar": n_real_H * A1_S,
    }
    unbroken["a1_net"] = unbroken["a1_gauge"] + unbroken["a1_ferm"] + unbroken["a1_scalar"]
    # Broken to SM: 45−12=33 Proca
    broken = {
        "n_massless": 12,
        "n_proca": 33,
        "a1_gauge": a1_gauge(12, 33),
        "a1_ferm": n_weyl * A1_W,
        "a1_scalar": n_real_H * A1_S,
    }
    broken["a1_net"] = broken["a1_gauge"] + broken["a1_ferm"] + broken["a1_scalar"]
    assert unbroken["a1_gauge"] == F(-30)
    assert broken["a1_gauge"] == F(17, 2)  # −8 + 33/2 = (−16+33)/2 = 17/2
    return {
        "group": "SO(10)",
        "unbroken": unbroken,
        "broken_to_SM": broken,
        "delta_a1_from_heavy_Proca": broken["a1_net"] - unbroken["a1_net"],
        "sign_OK_unbroken": unbroken["a1_net"] > 0,
        "sign_OK_broken": broken["a1_net"] > 0,
        "label_content": "derived (group dim + 16-plets); Higgs ASSUMED minimal 10_H",
    }


def spectral_with_gut(model: str = "SU5", include_higgs: bool = False) -> dict:
    """Replace UV_SM bin by GUT unbroken (→M_GUT) + broken (M_GUT→v).

    Default include_higgs=False: conservative gauge+fermion only (Higgs
    multiplets ASSUMED and model-dependent). Set True to add minimal 5_H/10_H.
    """
    if model == "SU5":
        sp = su5_spectra()
    else:
        sp = so10_spectra()

    def net(phase: dict) -> F:
        if include_higgs:
            return phase["a1_net"]
        return phase["a1_gauge"] + phase["a1_ferm"]

    def S_for(m_gut: float, mu_uv: float) -> dict:
        rows = []
        if mu_uv > m_gut:
            rows.append(
                {
                    "name": "GUT_unbroken",
                    "mu_hi": mu_uv,
                    "mu_lo": m_gut,
                    "a1": net(sp["unbroken"]),
                }
            )
        rows.append(
            {
                "name": "GUT_broken_to_EW",
                "mu_hi": m_gut,
                "mu_lo": V_EW,
                "a1": net(sp["broken_to_SM"]),
            }
        )
        rows += [
            {"name": "EW_to_top", "mu_hi": V_EW, "mu_lo": 173.0, "a1": A1_EW},
            {"name": "top_to_QCD", "mu_hi": 173.0, "mu_lo": LAMBDA_QCD, "a1": A1_AC},
            {"name": "hadronic", "mu_hi": LAMBDA_QCD, "mu_lo": 0.1, "a1": A1_HAD},
        ]
        S = 0.0
        out_rows = []
        for r in rows:
            if r["mu_hi"] <= r["mu_lo"]:
                continue
            w = 2.0 * log(r["mu_hi"] / r["mu_lo"])
            c = float(r["a1"]) * w
            S += c
            out_rows.append({**r, "weight": w, "contrib": c})
        return {"rows": out_rows, "S": S, "flips": S > 0}

    w_sm = 2.0 * log(MU_UV / V_EW)
    S_sm = (
        float(A1_SM_UV) * w_sm
        + float(A1_EW) * 2.0 * log(V_EW / 173.0)
        + float(A1_AC) * 2.0 * log(173.0 / LAMBDA_QCD)
        + float(A1_HAD) * 2.0 * log(LAMBDA_QCD / 0.1)
    )
    cases = {
        "M_GUT_eq_MU_UV": S_for(MU_UV, MU_UV),
        "M_GUT_1e15_MU_1e16": S_for(1e15, 1e16),
        "M_GUT_1e16_MU_MPl": S_for(1e16, 1.22e19),
    }
    return {
        "model": model,
        "include_higgs": include_higgs,
        "S_SM_baseline": S_sm,
        "cases": cases,
        "any_flip": any(c["flips"] for c in cases.values()),
        "flip_with_wide_Planck_bin": cases["M_GUT_1e16_MU_MPl"]["flips"],
        "coincidence_flag": (
            "M_GUT∼10^{15–16} GeV ∼ μ_UV often assumed — structural; "
            "Planck→GUT is the honest stress test; Higgs multiplets ASSUMED"
        ),
        "label": "derived gauge+fermion arithmetic; Higgs optional/assumed",
    }


def self_checks() -> None:
    s5 = su5_spectra()
    assert s5["unbroken"]["a1_gauge"] == F(-16)
    assert s5["broken_to_SM"]["a1_gauge"] == F(-2)
    assert s5["sign_OK_unbroken"] is False
    # broken SU5 net: −2 + 45/12 + 10/6 = −2 + 15/4 + 5/3
    # = (−24 + 45 + 20)/12 = 41/12 > 0
    assert s5["broken_to_SM"]["a1_net"] == F(41, 12)
    assert s5["sign_OK_broken"] is True
    s10 = so10_spectra()
    assert s10["unbroken"]["a1_gauge"] == F(-30)
    assert s10["broken_to_SM"]["a1_gauge"] == F(17, 2)
    assert s10["sign_OK_unbroken"] is False
    assert s10["sign_OK_broken"] is True
    # Unbroken GUT makes UV worse; broken can flip local a1 but spectral
    # integral may still fail if unbroken bin is wide
    # Conservative: gauge+fermion only (no assumed Higgs)
    g5 = spectral_with_gut("SU5", include_higgs=False)
    assert g5["S_SM_baseline"] < 0
    assert g5["cases"]["M_GUT_eq_MU_UV"]["flips"] is True
    assert g5["flip_with_wide_Planck_bin"] is False  # unbroken Kabats poison
    g10 = spectral_with_gut("SO10", include_higgs=False)
    assert so10_spectra()["sign_OK_unbroken"] is False
    assert g10["cases"]["M_GUT_eq_MU_UV"]["flips"] is True
    # SO(10) has more heavy Proca — may still flip with Planck bin
    assert "M_GUT_1e16_MU_MPl" in g10["cases"]


def main():
    print("=== GUT spectra SU(5)/SO(10) vs UV a1 sign ===\n")
    for sp in (su5_spectra(), so10_spectra()):
        print(f"  {sp['group']}:")
        u, b = sp["unbroken"], sp["broken_to_SM"]
        print(
            f"    unbroken:  gauge={u['a1_gauge']} ferm={u['a1_ferm']} "
            f"sc={u['a1_scalar']} net={u['a1_net']}  sign+? {sp['sign_OK_unbroken']}"
        )
        print(
            f"    broken SM: gauge={b['a1_gauge']} ferm={b['a1_ferm']} "
            f"sc={b['a1_scalar']} net={b['a1_net']}  sign+? {sp['sign_OK_broken']}"
        )
        print(f"    [{sp['label_content']}]")
    for model in ("SU5", "SO10"):
        for higgs in (False, True):
            g = spectral_with_gut(model, include_higgs=higgs)
            tag = "with minimal Higgs" if higgs else "gauge+fermion only"
            print(f"\n  Spectral ({model}, {tag}): S_SM={g['S_SM_baseline']:+.1f}")
            for name, c in g["cases"].items():
                print(f"    {name}: S={c['S']:+.1f}  flip={c['flips']}")
            print(
                f"  wide Planck→GUT flip? {g['flip_with_wide_Planck_bin']}  [{g['label']}]"
            )
        print(f"  coincidence: {spectral_with_gut(model)['coincidence_flag']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
Unbroken GUT a1 is MORE negative. Broken-to-SM local a1 turns +
(X/Y Proca; derived gauge arithmetic). Narrow window (M_GUT≈μ_UV)
can flip S. Conservative gauge+fermion SU(5) with Planck→GUT does NOT
flip — unbroken Kabats re-poison. SO(10) more Proca-friendly. Not a
scale-independent flip theorem; Higgs content ASSUMED when included.
"""
    )
    print("=== GUT SPECTRUM COMPLETE ===")


if __name__ == "__main__":
    main()
