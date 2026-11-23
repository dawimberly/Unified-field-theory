#!/usr/bin/env python3
"""Positivity of induced 1/G as a selection principle.

Treat S>0 (and rough M_Pl² magnitude) as constraints; map allowed
region in (N_V, M_*, Λ_UV). Quantify remaining freedom / tuning.
Flag coincidences. Fields-only language.

No tuning dials — report how much parameter volume survives.

Run:
  OPENBLAS_NUM_THREADS=1 python3 ffz_positivity_region.py
"""

from __future__ import annotations

from math import log, log10

# Reuse spectral from minimal spectrum
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_minimal_spectrum import (  # noqa: E402
    M_PL,
    spectral_S,
)

# Rough Sakharov: 1/(16πG) ~ |a1| Λ_UV² / (12π)  → need Λ_UV ~ M_Pl for magnitude
# We flag magnitude as ASSUMED cutoff coincidence, not a derived win.


def scan_region() -> dict:
    """Grid scan: N_V, M_*, Λ_UV."""
    n_vals = list(range(0, 70, 1))
    m_star_vals = [1e14, 3e14, 1e15, 3e15, 1e16, 3e16, 1e17]
    mu_uv_vals = [1e17, 1e18, 1e19, M_PL, 5e19]

    allowed = []
    total = 0
    for n in n_vals:
        for ms in m_star_vals:
            for mu in mu_uv_vals:
                if mu <= ms:
                    continue
                total += 1
                s = spectral_S(n, ms, mu)
                if s > 0:
                    allowed.append(
                        {
                            "n_v": n,
                            "m_star": ms,
                            "mu_uv": mu,
                            "S": s,
                        }
                    )

    # Volume fraction
    frac = len(allowed) / total if total else 0.0
    # N_V range among allowed at fiducial (M_*=1e16, mu=M_Pl)
    fid = [a for a in allowed if a["m_star"] == 1e16 and abs(a["mu_uv"] - M_PL) < 1e10]
    n_fid = sorted({a["n_v"] for a in fid})
    return {
        "n_grid": len(n_vals),
        "m_grid": len(m_star_vals),
        "mu_grid": len(mu_uv_vals),
        "total_points": total,
        "n_allowed": len(allowed),
        "fraction_allowed": frac,
        "n_v_fiducial_allowed": (n_fid[0], n_fid[-1]) if n_fid else None,
        "n_v_span_fiducial": (n_fid[-1] - n_fid[0] + 1) if n_fid else 0,
        "sample_min_S_allowed": min((a["S"] for a in allowed), default=None),
        "label": "derived scan on ASSUMED free-field spectral model",
    }


def magnitude_constraint() -> dict:
    """1/G ~ M_Pl² needs Λ_UV ~ M_Pl (Sakharov). Freedom?"""
    return {
        "requires_Lambda_UV_near_MPl": True,
        "coincidence_flag": (
            "Λ_UV∼M_Pl is the usual Sakharov coincidence — selecting it "
            "for magnitude is ASSUMED, not derived from S>0 alone"
        ),
        "tuning_assessment": (
            "Sign constraint S>0 cuts a large fraction of (N_V,M_*,Λ_UV) "
            "but leaves an O(1)–decade band in N_V and M_* — not a unique "
            "prediction. Magnitude forces Λ_UV∼M_Pl by hand."
        ),
        "label": "assumed magnitude criterion",
    }


def freedom_report(scan: dict) -> dict:
    """Is the surviving region 'tuning'?"""
    frac = scan["fraction_allowed"]
    span = scan["n_v_span_fiducial"]
    # Heuristic: if >20% of grid allowed and N_V span >10, call "loose"
    loose = frac > 0.2 and span > 10
    return {
        "fraction_allowed": frac,
        "n_v_span_at_fiducial": span,
        "characterized_as": "loose selection / residual freedom" if loose else "moderately selective",
        "is_fine_tuning_claim": False,  # not a % tuning of a continuous dial to 1e-N
        "plain": (
            f"Allowed fraction of scanned grid ≈ {100*frac:.1f}%. "
            f"At fiducial scales N_V spans {span} integers. "
            f"This is a coarse selection principle, not a derivation of "
            f"unique (N_V,M_*,Λ_UV). Calling it 'tuning' would overstate — "
            f"calling it 'predicted' would also overstate."
        ),
    }


def self_checks() -> None:
    sc = scan_region()
    assert sc["n_allowed"] > 0
    assert sc["fraction_allowed"] > 0
    assert sc["n_v_fiducial_allowed"] is not None
    lo, hi = sc["n_v_fiducial_allowed"]
    assert lo == 15  # matches Y at M_*=1e16, M_Pl
    assert hi >= 33
    fr = freedom_report(sc)
    assert fr["is_fine_tuning_claim"] is False
    mag = magnitude_constraint()
    assert mag["requires_Lambda_UV_near_MPl"] is True


def main():
    print("=== Positivity selection region (N_V, M_*, Λ_UV) ===\n")
    sc = scan_region()
    print(
        f"  Grid points={sc['total_points']}  allowed={sc['n_allowed']}  "
        f"fraction={100*sc['fraction_allowed']:.1f}%  [{sc['label']}]"
    )
    print(f"  Fiducial (M_*=1e16, Λ_UV=M_Pl) N_V ∈ {sc['n_v_fiducial_allowed']} "
          f"(span {sc['n_v_span_fiducial']})")
    mag = magnitude_constraint()
    print(f"  Magnitude: {mag['coincidence_flag']}")
    print(f"  {mag['tuning_assessment']}")
    fr = freedom_report(sc)
    print(f"\n  Freedom: {fr['plain']}")
    print(f"  Characterized as: {fr['characterized_as']}")
    self_checks()
    print("\n  self-checks PASS")
    print(
        """
Verdict
-------
S>0 is a useful cut but leaves substantial freedom in (N_V, M_*).
Λ_UV∼M_Pl for Planck magnitude remains an ASSUMED coincidence.
Not unique selection; not fine-tuning drama either.
"""
    )
    print("=== POSITIVITY REGION COMPLETE ===")


if __name__ == "__main__":
    main()
