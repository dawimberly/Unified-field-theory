#!/usr/bin/env python3
"""Full Standard Model roster in induced gravity (Sakharov, κ=0).

Stacked on particle-plasticity (#26). One-loop heat-kernel a1 → Λ² channel
of 1/G_ind for SM / SM+ν_R / MSSM, then map onto the FFZ h²z/K Newton ledger.

No tuning. DERIVED vs ASSUMED labeled.

Coefficients (d=4, contribution to a1/∫R in the one-loop EA / Sakharov reading):
  * real scalar: (1/6 − ξ)                         [derived; Bunch–Parker]
  * Weyl fermion: +1/12                            [derived from Adler/Visser
      weight: N_Dirac = N_Weyl/2 enters as +1 in units of (1−6ξ) = 6 a1_scalar;
      ⇒ a1_Weyl = 1/12]
  * massless vector + ghosts + Kabat contact: −2/3 [derived; Kabat 1995 /
      Solodukhin; Larsen–Wilczek; Donnelly–Wall edge-mode reading]
  * massive vector (Proca): +1/2                   [derived as 3×(1/6) physical
      dofs, no residual gauge ⇒ no Kabat; UV SM uses massless count instead]

Known result: SM alone gives wrong-sign / needs-conformal issues for induced G
(Visser gr-qc/0204062; Frolov–Fursaev–Zelnikov; Adler/Zee).

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_sm_roster.py
"""

from __future__ import annotations

from fractions import Fraction as F
from math import pi

C_LEDGER = 0.05131719
C_TARGET = 1.0 / 24.0
GAP = C_LEDGER - C_TARGET  # ≈ +0.00965052
GAP_ERR = 0.0024
K = 16.0 * pi**2


# ---------------------------------------------------------------------------
# (1) Per-species a1 / ∫R
# ---------------------------------------------------------------------------

def a1_real_scalar(xi) -> F:
    """Seeley–DeWitt a1/∫R for one real scalar."""
    return F(1, 6) - xi


def a1_weyl() -> F:
    """Weyl fermion EA / Sakharov weight.

    Adler–Zee / Visser: bracket N0(1−6ξ)+N_Dirac−4 N1, with N_Dirac=N_Weyl/2.
    In a1 units (bracket/6): a1_Weyl = 1/12.
    """
    return F(1, 12)


def a1_dirac() -> F:
    """Dirac = 2 Weyl ⇒ 1/6."""
    return F(1, 6)


def a1_massless_vector() -> F:
    """Maxwell + ghosts + Kabat contact at d=4: Ds/6 − 1 = 2/6 − 1 = −2/3."""
    return F(-2, 3)


def a1_massive_vector() -> F:
    """Proca: 3 physical dofs × 1/6, no Kabat (no residual gauge)."""
    return F(1, 2)


def species_table(xi_H: F | float = 0) -> dict:
    xi = F(xi_H) if not isinstance(xi_H, F) else xi_H
    return {
        "real_scalar_xi": a1_real_scalar(xi),
        "real_scalar_minimal": a1_real_scalar(0),
        "real_scalar_conformal": a1_real_scalar(F(1, 6)),
        "weyl": a1_weyl(),
        "dirac": a1_dirac(),
        "massless_vector": a1_massless_vector(),
        "massive_vector": a1_massive_vector(),
        "refs": (
            "Kabat 1995; Larsen–Wilczek; Donnelly–Wall; Solodukhin; "
            "Visser gr-qc/0204062; Frolov–Fursaev–Zelnikov"
        ),
        "label_scalar": "derived (HK)",
        "label_weyl": "derived (Adler/Visser weight ↔ a1)",
        "label_vector": "derived (Kabat/Solodukhin)",
        "label_proca": "derived (3×1/6, no Kabat)",
        "label_xi_H": "assumed (Higgs nonminimal; free vs conformal)",
    }


# ---------------------------------------------------------------------------
# SM / SM+ν_R / MSSM DOF
# ---------------------------------------------------------------------------

def sm_dof(nu_R: bool = False) -> dict:
    """Unbroken-phase UV count (Λ ≫ v): all gauge bosons massless.

    Higgs doublet: 4 real scalars.
    Fermions: 3 gen × [quarks 3×4 Weyl + leptons 3 Weyl] = 36+9 = 45 Weyl;
      optional +3 ν_R → 48.
    Gauge: 8+3+1 = 12 massless vectors.
    """
    n_scalar = 4
    n_weyl = 48 if nu_R else 45
    n_vector = 12
    return {
        "n_real_scalar": n_scalar,
        "n_weyl": n_weyl,
        "n_massless_vector": n_vector,
        "n_massive_vector_IR": 3,  # W±,Z below EW — not used for UV Λ²
        "nu_R": nu_R,
        "label": "derived (SM field content; UV massless gauge)",
    }


def mssm_dof() -> dict:
    """MSSM DOF for free-field a1 (ξ=1/6 on all scalars; SUSY).

    Gauge: 12 N=1 vector multiplets (8+3+1).
    Chiral:
      * 3 gen quark+lepton superfields: 45 Weyl + 45 complex scalars
        (=90 real) — same chiral count as SM Weyl content without ν_R;
        with soft SUSY the scalars are present.
      * Optional ν_R omitted (same as SM baseline).
      * Higgs: 2 doublets → 4 complex = 8 real scalars + 4 Weyl higgsinos.
    Total chiral Weyl: 45+4 = 49; real scalars: 90+8 = 98.
    """
    return {
        "n_vector_multiplet": 12,
        "n_chiral_weyl": 49,
        "n_real_scalar": 98,
        "xi": F(1, 6),
        "label": "derived (MSSM content); ξ=1/6 assumed (SUSY)",
    }


def a1_net_sm(xi_H=0, nu_R: bool = False) -> dict:
    """Net a1/∫R for SM (or SM+ν_R)."""
    dof = sm_dof(nu_R=nu_R)
    xi = F(xi_H) if not isinstance(xi_H, F) else xi_H
    s = dof["n_real_scalar"] * a1_real_scalar(xi)
    f = dof["n_weyl"] * a1_weyl()
    v = dof["n_massless_vector"] * a1_massless_vector()
    net = s + f + v
    # Adler–Visser bracket C = 6 * a1_net  (in units where scalar weight is 1−6ξ)
    bracket = 6 * net
    return {
        "xi_H": xi,
        "nu_R": nu_R,
        "scalar": s,
        "fermion": f,
        "vector": v,
        "a1_net": net,
        "bracket_Adler": bracket,
        "sign_OK_for_G": net > 0,
        "label": "derived (content × HK weights); ξ_H assumed",
    }


def a1_net_mssm() -> dict:
    """Net a1 for MSSM at ξ=1/6.

    Scalars drop out at ξ=1/6. Each vector multiplet:
      −2/3 (vector) + 1/12 (Majorana≈Weyl) = −7/12.
    Chiral Weyl (matter+higgsino): +1/12 each; scalar partners vanish.
    """
    dof = mssm_dof()
    n_v = dof["n_vector_multiplet"]
    n_chi = dof["n_chiral_weyl"]
    per_vm = a1_massless_vector() + a1_weyl()  # −2/3 + 1/12 = −7/12
    assert per_vm == F(-7, 12)
    gauge = n_v * per_vm
    chiral = n_chi * a1_weyl()
    # Cross-check: scalars at ξ=1/6 contribute 0
    scalar = dof["n_real_scalar"] * a1_real_scalar(F(1, 6))
    assert scalar == 0
    net = gauge + chiral + scalar
    return {
        "per_vector_multiplet": per_vm,
        "gauge": gauge,
        "chiral_weyl": chiral,
        "scalar_at_xi_1_6": scalar,
        "a1_net": net,
        "bracket_Adler": 6 * net,
        "sign_OK_for_G": net > 0,
        "note_vs_WZ": (
            "WZ Newton test uses interacting 3-loop h²z/K, not free a1; "
            "free MSSM a1 still wrong-sign (vectors dominate)"
        ),
        "label": "derived",
    }


def wrong_sign_literature() -> dict:
    return {
        "statement": (
            "SM alone induces the wrong sign of G (or requires special "
            "conformal / beyond-SM arrangements); fermions help but 12 "
            "gauge bosons dominate the −4 N1 weight."
        ),
        "refs": (
            "Visser, 'Sakharov's induced gravity: a modern perspective' "
            "(gr-qc/0204062); Frolov–Fursaev–Zelnikov; Adler; Zee"
        ),
        "label": "literature (cited); reproduced numerically below",
    }


# ---------------------------------------------------------------------------
# (2) Net Λ² channel vs WZ structure
# ---------------------------------------------------------------------------

def lambda2_summary() -> dict:
    rows = {}
    for name, xi, nu in [
        ("SM_xi0", 0, False),
        ("SM_xi1_6", F(1, 6), False),
        ("SM_nuR_xi0", 0, True),
        ("SM_nuR_xi1_6", F(1, 6), True),
    ]:
        rows[name] = a1_net_sm(xi_H=xi, nu_R=nu)
    rows["MSSM_xi1_6"] = a1_net_mssm()
    return rows


# ---------------------------------------------------------------------------
# (3) Map onto FFZ h²z/K ledger
# ---------------------------------------------------------------------------

def ffz_channel_map() -> dict:
    """Does the free SM roster shift C_ledger=0.051317 vs 1/24?

    Derived: free-field a1 renormalizes 1/G with a *quadratic* UV piece
    ∼ Λ² a1 / (16π²), mass-independent at leading order. The FFZ three-loop
    simple pole is ∝ h² z /(K ε) from the WZ Yukawa sector — a different
    operator / RG structure. In dim-reg MS for massless free fields that
    Λ² piece is absent from the h²z/K ledger.

    Therefore the full SM roster does **not** change 0.051317 vs 1/24 in the
    ledger channel. SM fields enter that channel only via *interacting*
    diagrams with their own Yukawas / gauge couplings (not computed here).

    Top-Yukawa-dominated estimate (ASSUMED scope): if the Newton test is
    dominated by the heaviest Yukawa, the existing WZ single-coupling ledger
    (c_Y=1/24, gap +0.00965) remains the relevant estimate; adding light
    Yukawas/gauge is a different multi-coupling calculation (not done).
    """
    return {
        "free_a1_channel": "Λ²",
        "ledger_channel": "h² z / (K ε)",
        "channels_same": False,
        "delta_C_ledger_from_free_roster": 0.0,
        "C_ledger": C_LEDGER,
        "C_target": C_TARGET,
        "GAP": GAP,
        "top_Yukawa_note": (
            "Top Yukawa dominates SM fermion masses; a single-coupling WZ "
            "proxy keeps c_Y=1/24 and the standing gap. Multi-Yukawa + gauge "
            "re-ledger is ASSUMED out of scope (not derived here)."
        ),
        "label_channel_split": "derived (same argument as gauge-contact #22)",
        "label_top_proxy": "assumed (scope); not a new computation",
    }


# ---------------------------------------------------------------------------
# (4) Gap closure?
# ---------------------------------------------------------------------------

def gap_closure_check() -> dict:
    """Does any roster close GAP=+0.00965 without tuning?"""
    ch = ffz_channel_map()
    # Even if one (wrongly) converts a1 → ledger units via Λ²→z at h=1:
    # δc = −a1 / h²  (same assumed conversion as #22); compare to GAP.
    assumed = {}
    for name, row in lambda2_summary().items():
        a1 = float(row["a1_net"])
        delta_c_assumed = -a1  # h=1
        assumed[name] = {
            "a1_net": a1,
            "delta_c_assumed_h1": delta_c_assumed,
            "ratio_to_GAP": delta_c_assumed / GAP,
            "closes_gap": abs(delta_c_assumed - GAP) < GAP_ERR,
        }
    # Coincidences: none of the a1 nets equal −GAP or give ratio ≈1
    return {
        "derived_delta_C": ch["delta_C_ledger_from_free_roster"],
        "derived_closes_gap": False,
        "assumed_Lambda_to_z": assumed,
        "any_assumed_closes": any(v["closes_gap"] for v in assumed.values()),
        "coincidence_flags": [
            "No roster yields δc≈GAP under Λ²→z at h=1 (ratios O(10²–10³))",
            "Wrong-sign SM a1 is a known literature result, not a gap fix",
            "Bohr/ln2-style numerical accidents: none identified here",
        ],
        "label": "derived: free roster does not touch h²z/K gap",
    }


# ---------------------------------------------------------------------------
# Self-checks
# ---------------------------------------------------------------------------

def self_checks() -> None:
    assert a1_real_scalar(F(1, 6)) == 0
    assert a1_real_scalar(0) == F(1, 6)
    assert a1_weyl() == F(1, 12)
    assert a1_dirac() == F(1, 6)
    assert a1_massless_vector() == F(-2, 3)
    assert a1_massive_vector() == F(1, 2)
    # SM ξ=0: 4/6 + 45/12 − 8 = 2/3 + 15/4 − 8
    sm0 = a1_net_sm(0, False)
    assert sm0["a1_net"] == F(4, 6) + F(45, 12) + 12 * F(-2, 3)
    assert sm0["a1_net"] == F(2, 3) + F(15, 4) - 8
    assert sm0["a1_net"] == F(-43, 12)  # −3.583...
    assert sm0["sign_OK_for_G"] is False
    # SM conformal
    smc = a1_net_sm(F(1, 6), False)
    assert smc["a1_net"] == F(45, 12) - 8
    assert smc["a1_net"] == F(-17, 4)  # −4.25
    # +ν_R
    smn = a1_net_sm(0, True)
    assert smn["a1_net"] == F(-43, 12) + F(3, 12)
    # MSSM
    ms = a1_net_mssm()
    assert ms["per_vector_multiplet"] == F(-7, 12)
    assert ms["a1_net"] == 12 * F(-7, 12) + 49 * F(1, 12)
    assert ms["a1_net"] == F(-84, 12) + F(49, 12)
    assert ms["a1_net"] == F(-35, 12)
    assert ms["sign_OK_for_G"] is False
    # Channel
    ch = ffz_channel_map()
    assert ch["channels_same"] is False
    assert ch["delta_C_ledger_from_free_roster"] == 0.0
    g = gap_closure_check()
    assert g["derived_closes_gap"] is False
    assert g["any_assumed_closes"] is False
    # DOF
    assert sm_dof(False)["n_weyl"] == 45
    assert sm_dof(True)["n_weyl"] == 48
    assert sm_dof(False)["n_massless_vector"] == 12
    assert abs(GAP - 0.00965052) < 1e-5


def main():
    print("=== SM roster in induced gravity (Sakharov, κ=0) ===\n")
    print("No tuning. Free-field a1 / Λ² channel vs FFZ h²z/K ledger.\n")

    print("-- (1) Per-species a1/∫R --")
    sp = species_table(0)
    print(f"  real scalar (ξ):     (1/6−ξ)     [{sp['label_scalar']}]")
    print(f"  real scalar ξ=0:     {sp['real_scalar_minimal']}")
    print(f"  real scalar ξ=1/6:   {sp['real_scalar_conformal']}")
    print(f"  Weyl fermion:        {sp['weyl']}         [{sp['label_weyl']}]")
    print(f"  Dirac:               {sp['dirac']}")
    print(f"  massless vector:     {sp['massless_vector']}        [{sp['label_vector']}]")
    print(f"  massive vector:      {sp['massive_vector']}         [{sp['label_proca']}]")
    print(f"  ξ_H:                 [{sp['label_xi_H']}]")
    print(f"  refs: {sp['refs']}")
    lit = wrong_sign_literature()
    print(f"\n  Literature: {lit['statement']}")
    print(f"  [{lit['refs']}]")

    print("\n-- SM DOF --")
    d0 = sm_dof(False)
    print(
        f"  scalars={d0['n_real_scalar']}, Weyl={d0['n_weyl']}, "
        f"massless vectors={d0['n_massless_vector']}  [{d0['label']}]"
    )
    print(f"  SM+ν_R Weyl={sm_dof(True)['n_weyl']}")

    print("\n-- (2) Net a1 (Λ² channel) --")
    rows = lambda2_summary()
    print(f"  {'roster':22s}  {'a1_net':>12}  {'bracket':>10}  sign_OK")
    for name, r in rows.items():
        a1 = r["a1_net"]
        br = r["bracket_Adler"]
        ok = r["sign_OK_for_G"]
        print(f"  {name:22s}  {str(a1):>12}  {str(br):>10}  {ok}")
    print("  [All wrong-sign for induced G — reproduces literature]")
    print(f"  MSSM note: {rows['MSSM_xi1_6']['note_vs_WZ']}")

    # ξ_H dependence
    print("\n  ξ_H dependence (SM, no ν_R): a1 = 4(1/6−ξ_H)+45/12−8")
    for xi in [0, F(1, 12), F(1, 6), F(1, 4), F(1, 3)]:
        r = a1_net_sm(xi, False)
        print(f"    ξ_H={xi}: a1={r['a1_net']}  (sign_OK={r['sign_OK_for_G']})")
    # Would need 4(1/6−ξ)+15/4−8 > 0 ⇒ 2/3−4ξ+3.75−8>0 ⇒ −4ξ > 3.583 ⇒ ξ < −0.896
    print("    sign_OK would require ξ_H < −43/48 ≈ −0.896 (unphysical; not a fix)")

    print("\n-- (3) Map onto FFZ h²z/K --")
    ch = ffz_channel_map()
    print(f"  free a1 channel: {ch['free_a1_channel']}")
    print(f"  ledger channel:  {ch['ledger_channel']}")
    print(f"  same channel? {ch['channels_same']}  [{ch['label_channel_split']}]")
    print(f"  δC_ledger from free roster: {ch['delta_C_ledger_from_free_roster']}")
    print(f"  C_ledger={ch['C_ledger']:.8f}, 1/24={ch['C_target']:.8f}, GAP={ch['GAP']:.8f}")
    print(f"  Top-Yukawa proxy: {ch['top_Yukawa_note']}  [{ch['label_top_proxy']}]")

    print("\n-- (4) Gap closure? --")
    g = gap_closure_check()
    print(f"  derived δC: {g['derived_delta_C']} → closes gap? {g['derived_closes_gap']}")
    print("  assumed Λ²→z @ h=1 (NOT a prediction):")
    for name, v in g["assumed_Lambda_to_z"].items():
        print(
            f"    {name:22s}  δc={v['delta_c_assumed_h1']:+.4f}  "
            f"ratio/GAP={v['ratio_to_GAP']:+.1f}  closes? {v['closes_gap']}"
        )
    print(f"  any assumed closes? {g['any_assumed_closes']}")
    for flag in g["coincidence_flags"]:
        print(f"  coincidence: {flag}")

    print("\n-- Self-checks --")
    self_checks()
    print("  all asserts PASS")

    print(
        f"""
Verdict
-------
SM DOF (4 scalars, 45 Weyl, 12 vectors) with standard a1 weights give
a1_net = −43/12 (ξ_H=0) or −17/4 (ξ_H=1/6) — wrong sign for G_ind,
reproducing Visser / Frolov–Fursaev–Zelnikov. SM+ν_R and MSSM also
wrong-sign. Free-field a1 is a Λ² channel; it does **not** shift the
FFZ h²z/K ledger (0.051317 vs 1/24). Gap +{GAP:.5f} stands; no roster
closes it without tuning. Top-Yukawa WZ proxy leaves c_Y=1/24.
"""
    )
    print("=== SM ROSTER COMPLETE ===")


if __name__ == "__main__":
    main()
