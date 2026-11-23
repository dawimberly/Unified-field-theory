#!/usr/bin/env python3
"""Unit test for ffz_distortion_phase (CI hook)."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_distortion_phase import (  # noqa: E402
    A0_BOHR,
    GRAN_SASSO_R0_MIN,
    bec_healing_from_case,
    candidate_widths,
    modified_taus,
    self_checks,
    supplies_r0_without_tuning,
    testable_prediction,
)


def test_self_checks():
    self_checks()


def test_no_universal_r0():
    assert supplies_r0_without_tuning(candidate_widths())[
        "any_derived_universal_R0_from_C_ledger"
    ] is False


def test_a0_coincidence_proximity():
    # Flag proximity; do NOT claim derivation
    assert abs(A0_BOHR / GRAN_SASSO_R0_MIN - 1.0) < 0.05


def test_bec_healing_above_bound():
    assert bec_healing_from_case()["xi_m"] > GRAN_SASSO_R0_MIN


def test_modified_taus_bec_not_dominant():
    bec = modified_taus()[0]
    assert bec.R0_choice.startswith("GP")
    assert bec.gravity_dominates is False


def test_prediction_flags():
    p = testable_prediction()
    assert p["distinct_from_plain_DP"] is True
    assert p["fixes_R0"] is False


if __name__ == "__main__":
    test_self_checks()
    test_no_universal_r0()
    test_a0_coincidence_proximity()
    test_bec_healing_above_bound()
    test_modified_taus_bec_not_dominant()
    test_prediction_flags()
    print("test_ffz_distortion_phase: ALL PASS")
