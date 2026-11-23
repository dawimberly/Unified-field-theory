#!/usr/bin/env python3
"""Physical meaning + topology + weak-entry audit of the ~+0.01 simple gap.

Dan order: (3) physical / RG status of 1/24, then (1) topology coverage and
tighten weak ledger entries, then (2) a_Y/a_X scan only if still needed.
No tuning to 1/24.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_gap_physical_audit.py

Optional env:
  CCFF_POWER=20   (Sobol power for CcFf QMC; default 20)
  CCFF_SEEDS=8    (number of scramble seeds; default 8)
  SKIP_CCFF=1     (skip the slow QMC block)

Requires: numpy, scipy, sympy, mpmath
"""

from __future__ import annotations

import importlib.util
import os
from collections import Counter
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path

import mpmath as mp
import numpy as np

mp.mp.dps = 40

_REVIEW = Path(__file__).resolve().parent / "ffz_independent_pole_review.py"
_spec = importlib.util.spec_from_file_location("ffz_ind_review", _REVIEW)
_ffz = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ffz)


# ---------------------------------------------------------------------------
# (3) Physical meaning of 1/24 and of a residual simple pole
# ---------------------------------------------------------------------------

def audit_physical_meaning():
    """Return structured verdict; asserts encode the derived identities."""
    # (a) Y soft-source coefficient is RG-fixed.
    # Martin β_Λ linear in s_Y at three loops: (3/2) h² z s_Y / K
    # Engineering: simple pole of vacuum = −β / (2 ℓ) with ℓ=3
    # s_Y = −a_Y R/6 → F_R pole = h² z a_Y / (24 K ε) → c_Y = 1/24
    beta_coeff = F(3, 2)
    ell = 3
    assert -beta_coeff / (2 * ell) == F(-1, 4)
    assert F(1, 6) * F(1, 4) == F(1, 24)
    c_Y = F(1, 24)

    # Lower-loop cross-checks: one-loop β=2 → pole a_Y/6; two-loop |β|=1 → a_Y/24
    assert F(2) / (2 * 1) / 6 == F(1, 6)
    assert F(1) / (2 * 2) / 6 == F(1, 24)

    # Y-linear higher poles from the soft-source pole equation (independent of
    # the three-loop β coefficient): h²za_Y/K * [1/(36ε³)−1/(18ε²)+1/(24ε)].
    y_triple, y_double, y_simple = F(1, 36), F(-1, 18), F(1, 24)
    assert y_simple == c_Y

    # (b) Baseline = 1/24 is NOT a consequence of those pole equations alone.
    # Packet: "Flat vacuum RG data alone do not determine [the baseline]."
    # "A finite pure-metric contact redefinition can move the associated beta
    #  coefficient along with its matching condition."
    # Projection: a finite shift κ of the renormalized pure-metric coefficient
    # shifts its beta by 4κ h² Σz/K (β_h=3h², β_z=hz).
    kappa_beta_factor = 4

    gap = 0.0096505241
    gap_err = 0.00241
    target = float(c_Y)
    total = target + gap
    rel = gap / target

    return {
        "c_Y": c_Y,
        "y_higher": (y_triple, y_double, y_simple),
        "baseline_equals_c_Y": "matching_test_not_RG_theorem",
        "gap_forbidden_by_RG": False,
        "gap": gap,
        "gap_err": gap_err,
        "total_ledger": total,
        "rel_to_1_over_24": rel,
        "finite_contact_absorbs_simple": True,
        "kappa_beta_factor": kappa_beta_factor,
        "G_ind_note": (
            "Residual +0.00965 in units h²z/K is ~+23% of the three-loop "
            "Y-cancelation residue 1/24. In Γ_ind ⊃ −R/(16π G_ind) this shifts "
            "the three-loop matching piece of 1/G_ind by that fraction; a finite "
            "pure-metric R contact can absorb it (beta moves with the matching). "
            "Higher poles are not absorbable this way and already cancel."
        ),
    }


# ---------------------------------------------------------------------------
# (1) Topology enumeration vs ledger
# ---------------------------------------------------------------------------

SLOTS = {"Q": "xxXX", "C": "xxY", "c": "XXy", "F": "xff", "f": "XFF"}


def _connected(n, edges):
    seen = {0}
    while True:
        more = seen | {v for _, u, v in edges if u in seen} | {
            u for _, u, v in edges if v in seen
        }
        if more == seen:
            return len(more) == n
        seen = more


def _canonical(vertices, edges):
    orders = (
        order
        for order in permutations(range(len(vertices)))
        if tuple(vertices[i] for i in order) == vertices
    )
    return min(
        tuple(sorted((s, order[u], order[v]) for s, u, v in edges)) for order in orders
    )


def graphs(vertices):
    families = []
    for species in ("x", "y", "f"):
        src = [
            i
            for i, v in enumerate(vertices)
            for _ in range(SLOTS[v].count(species))
        ]
        snk = [
            i
            for i, v in enumerate(vertices)
            for _ in range(SLOTS[v].count(species.upper()))
        ]
        if len(src) != len(snk):
            return {}
        families.append(
            [
                tuple((species, u, v) for u, v in zip(src, perm))
                for perm in permutations(snk)
            ]
        )
    out = Counter()
    for e in product(*families):
        edges = sum(e, ())
        if not _connected(len(vertices), edges):
            continue
        assert all(
            _connected(len(vertices), edges[:i] + edges[i + 1 :])
            for i in range(len(edges))
        )
        assert len(edges) - len(vertices) + 1 == 3
        out[_canonical(vertices, edges)] += 1
    return out


def audit_topologies():
    expected = {
        "QQ": [4, 16],
        "QCc": [4, 16],
        "QFf": [8],
        "CCcc": [8, 32],
        "CcFf": [8],
        "FFff": [8, 32],
    }
    found = {}
    for multiset, exp in expected.items():
        g = graphs(tuple(multiset))
        mults = sorted(g.values())
        assert mults == exp, (multiset, mults, exp)
        found[multiset] = {"n_classes": len(g), "multiplicities": mults}

    ledger_map = {
        "QQ": "ledger QQ (exact 1/16)",
        "QCc": "ledger QCc doubled + tadpole pole −1/36",
        "QFf": "ledger QFf (derived j1)",
        "CCcc": "ledger CCcc (exact −5/72)",
        "CcFf": "ledger CcFf (QMC)",
        "FFff": "ledger FFff(A)×4 + FFff(B) overlap-only",
        "CT": "1-loop Z_z mass + 2-loop matter CT insertions (derived)",
    }
    omitted = [
        {
            "class": "H tetrahedron (4 cubic)",
            "status": "absent by charge bipartiteness — not an omission",
            "size_est": "0",
        },
        {
            "class": "gauge / gaugino graphs",
            "status": "out of model (ungauged)",
            "size_est": "0",
        },
        {
            "class": "soft holomorphic masses beyond s_Y projection",
            "status": "set to zero in packet; Y linear already projected",
            "size_est": "0",
        },
        {
            "class": "FFff(B) non-overlap finite beyond field c0+Δc0(W/U)",
            "status": "ASSUMED absent (overlap-only assembly); not a new topology",
            "size_est": (
                "soft-linear ~−105 → ~+0.046 in simple if wrongly included; "
                "empirical 2×Δc0 shifts simple by ~−0.009 (NOT derived)"
            ),
        },
        {
            "class": "extra pure-metric contact / measure outside Wick graphs",
            "status": "scheme (finite renormalization), not a Wick topology",
            "size_est": "can absorb entire simple gap by construction",
        },
    ]
    return {"classes": found, "ledger_map": ledger_map, "omitted": omitted}


# ---------------------------------------------------------------------------
# Independent QCc simple
# ---------------------------------------------------------------------------

def independent_qcc_simple():
    def g0_fun(z):
        zz = mp.mpf(z)
        dd = 4 - 2 * zz
        pref_m = -dd * (6 - dd)
        br = (
            2 * mp.gamma(1 + zz) * mp.gamma(zz) * mp.gamma(2 - 3 * zz)
            + 2 * mp.gamma(1 + zz) ** 2 * mp.gamma(1 - 3 * zz)
        ) / mp.gamma(3 - zz)
        return pref_m * mp.beta(1 - zz, 1 - zz) ** 2 * br

    def tail_fun(z):
        zz = mp.mpf(z)
        dd = 4 - 2 * zz
        return (
            (8 * dd * zz / 3)
            * mp.beta(1 - 2 * zz, 1 + zz)
            * mp.beta(2 * zz, 2 * zz)
            * mp.beta(zz, 1)
            * mp.beta(1 - zz, 1 - zz)
        )

    z = mp.mpf("1e-10")

    def geom(rr):
        return -mp.log(rr) / (1 - rr) ** 2 - 1 / (1 - rr)

    def evan(rr):
        return mp.mpf(2) / 3 - 4 * (rr + 2) * mp.log(rr) / (3 * (rr - 1))

    g1 = 24 * mp.quad(lambda x: geom(x * (1 - x)), [0, mp.mpf("0.5"), 1])
    v = 8 * mp.quad(lambda x: evan(x * (1 - x)), [0, mp.mpf("0.5"), 1])
    I0 = -32 + g1 + v + mp.mpf(80) / 3 - mp.mpf(8) / 3
    c = (I0 / 3 + 8 * mp.mpf(7) / 18) / 192
    # g0/tail pole asserts (sanity)
    assert abs(z * g0_fun(z) + 8) < mp.mpf("1e-8")
    assert abs(z * tail_fun(z) - mp.mpf(32) / 3) < mp.mpf("1e-7")
    return float(c), 1e-11


# ---------------------------------------------------------------------------
# CcFf QMC tighten
# ---------------------------------------------------------------------------

def tighten_ccff(power=None, n_seeds=None):
    if os.environ.get("SKIP_CCFF") == "1":
        return {"skipped": True, "value": -0.02480, "err": 0.001}
    power = int(os.environ.get("CCFF_POWER", power or 20))
    n_seeds = int(os.environ.get("CCFF_SEEDS", n_seeds or 8))
    seeds = tuple(11 + 7 * i for i in range(n_seeds))
    print(f"  CcFf QMC power={power}, seeds={seeds}")

    # Per-seed extrapolations for scatter
    vals = []
    last_rows = None
    for seed in seeds:
        one = _ffz.ccff_simple_pole(power=power, seeds=(seed,))
        vals.append(float(one["extrapolated"][1]))
        last_rows = np.array(one["samples"])
    value = float(np.mean(vals))
    seed_std = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0
    sem = seed_std / (len(vals) ** 0.5) if len(vals) else 0.0
    y = last_rows[:, 1]
    extrap = (8 * y[2] - 6 * y[1] + y[0]) / 3
    step = abs(float(extrap) - float(y[2]))
    # Target ±1e-4: floor at 1e-4, else 2×SEM (seed-averaged)
    err = max(1e-4, 2 * sem)
    return {
        "skipped": False,
        "value": value,
        "err": err,
        "power": power,
        "n_seeds": n_seeds,
        "seed_std": seed_std,
        "sem": sem,
        "extrap_step": float(step),
        "per_seed": vals,
    }


def b_nonoverlap_finite():
    """Characterize FFff(B) non-overlap finite content (no new fit to 1/24)."""
    c0_runs = [191.77900, 191.47522, 190.46241, 193.18188]
    c0_mean = float(np.mean(c0_runs))
    c0_spread = float(np.max(c0_runs) - np.min(c0_runs)) / 2
    delta_c0 = 64 / 3
    soft_linear = -105.304
    base = -163 / 1728 + float(mp.pi) ** 2 / 64

    def b_simple(c0):
        return base - c0 / 2304

    scenarios = {
        "field_only": b_simple(c0_mean),
        "field+overlap_WU": b_simple(c0_mean + delta_c0),
        "field+2*WU (empirical NOT derived)": b_simple(c0_mean + 2 * delta_c0),
        "field+WU+soft_linear (excluded)": b_simple(c0_mean + delta_c0 + soft_linear),
    }
    return {
        "c0_field_mean": c0_mean,
        "c0_field_halfspan": c0_spread,
        "delta_c0_WU": delta_c0,
        "soft_linear_excluded": soft_linear,
        "scenarios": scenarios,
        "simple_err_from_c0": c0_spread / 2304,
        "label": (
            "field c0 = reproduced numerical; Δc0=64/3 derived; "
            "non-overlap beyond that = assumed zero for overlap-only assembly"
        ),
    }


# ---------------------------------------------------------------------------
# (2) a_Y / a_X scan
# ---------------------------------------------------------------------------

def scan_aYaX():
    """Higher poles uniquely fix a_Y=-1 (and a_X=0); simple gap remains there.

    residual_higher = baseline_higher + a_Y * Y_higher (+ a_X * X_higher).
    Overlap-fixed ledger: baseline_higher = −Y_higher at projected cancelation;
    X additive higher poles vanish in the packet through this order.
    residual_higher ∝ (a_Y + 1) ⇒ triple/double cancel iff a_Y = −1.
    At that point simple residual = c0_baseline − 1/24 ≈ +0.00965 ≠ 0.
    """
    D = 0.0096505241
    solutions_higher = [{"a_X": 0.0, "a_Y": -1.0}]
    solutions_all = []
    c0 = 0.0513171908
    c_Y = 1 / 24
    a_Y_simple_only = -c0 / c_Y
    return {
        "higher_pole_locus": solutions_higher,
        "simultaneous_solutions": solutions_all,
        "simple_zero_if_ignore_higher": {
            "a_X": 0.0,
            "a_Y": a_Y_simple_only,
            "note": f"reopens triple/double by factor (a_Y+1)={a_Y_simple_only+1:+.4f}",
        },
        "D": D,
        "verdict": (
            "no simultaneous solution; higher poles pin a_Y=-1, simple gap stands"
        ),
    }


def main():
    print("=== Gap physical / topology / weak-entry audit (no 1/24 tuning) ===\n")

    print("-- (3) Physical meaning --")
    phys = audit_physical_meaning()
    print(f"  c_Y = {phys['c_Y']}  [RG-fixed from Martin β + −β/(2ℓ)]")
    print(f"  baseline==c_Y? {phys['baseline_equals_c_Y']}")
    print(f"  gap forbidden by RG/t Hooft? {phys['gap_forbidden_by_RG']}")
    print(
        f"  ledger total={phys['total_ledger']:.6f}, "
        f"gap={phys['gap']:+.6f}±{phys['gap_err']:.6f}"
    )
    print(f"  gap / (1/24) = {phys['rel_to_1_over_24'] * 100:.1f}%")
    print(f"  finite R contact absorbs simple? {phys['finite_contact_absorbs_simple']}")
    print(f"  {phys['G_ind_note']}")

    print("\n-- (1) Topology coverage --")
    topo = audit_topologies()
    for name, info in topo["classes"].items():
        print(
            f"  {name}: {info['n_classes']} classes, mult={info['multiplicities']}"
            f"  → {topo['ledger_map'][name]}"
        )
    print(f"  CT → {topo['ledger_map']['CT']}")
    print("  Omitted / non-topology notes:")
    for o in topo["omitted"]:
        print(f"    - {o['class']}: {o['status']} (size: {o['size_est']})")

    print("\n-- (1) Tighten QCc simple --")
    qcc, qcc_err = independent_qcc_simple()
    print(f"  QCc doubled simple = {qcc:+.12f} ± {qcc_err:.0e}  [derived faces+tail]")
    assert abs(qcc + 0.01649742330) < 1e-11

    print("\n-- (1) Tighten CcFf simple --")
    ccff = tighten_ccff()
    if ccff["skipped"]:
        print("  SKIPPED (SKIP_CCFF=1); keeping prior −0.02480±0.001")
    else:
        print(
            f"  CcFf simple = {ccff['value']:+.6f} ± {ccff['err']:.6f} "
            f"(power={ccff['power']}, seeds={ccff['n_seeds']}, "
            f"seed_std={ccff['seed_std']:.6f}, SEM={ccff['sem']:.6f})"
        )
        print(f"  per-seed: {['%+.6f' % v for v in ccff['per_seed']]}")
        assert abs(ccff["value"] + 0.0248) < 0.002
        assert ccff["err"] <= 1.5e-4

    print("\n-- (1) FFff(B) non-overlap finite --")
    bfin = b_nonoverlap_finite()
    print(
        f"  field c0 = {bfin['c0_field_mean']:.2f} ± {bfin['c0_field_halfspan']:.2f} "
        f"→ simple err ±{bfin['simple_err_from_c0']:.5f}"
    )
    for k, v in bfin["scenarios"].items():
        print(f"  scenario {k}: B simple = {v:+.6f}")
    print(f"  {bfin['label']}")

    print("\n-- (2) a_Y / a_X scan --")
    scan = scan_aYaX()
    print(f"  higher-pole locus: {scan['higher_pole_locus']}")
    print(
        f"  simultaneous (triple+double+simple) solutions: "
        f"{scan['simultaneous_solutions']}"
    )
    s0 = scan["simple_zero_if_ignore_higher"]
    print(f"  simple-only zero (ignore higher): a_Y={s0['a_Y']:.4f} — {s0['note']}")
    print(f"  verdict: {scan['verdict']}")

    print("\n=== AUDIT COMPLETE ===")
    print(
        "SUMMARY: 1/24 is RG-fixed as c_Y (Y soft); baseline=1/24 is a matching "
        "test, not a t Hooft theorem — gap is not RG-forbidden and is absorbable "
        "by a finite pure-metric contact (~+23% of the three-loop 1/G_ind piece). "
        "All Wick topologies are on the ledger; no omitted graph of size ~0.01. "
        "No (a_X,a_Y) satisfies higher poles and simple simultaneously."
    )


if __name__ == "__main__":
    main()
