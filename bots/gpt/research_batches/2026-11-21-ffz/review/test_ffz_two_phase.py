#!/usr/bin/env python3
"""Unit test for ffz_two_phase (CI hook)."""

from __future__ import annotations

import pathlib
import sys
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_two_phase import (  # noqa: E402
    A1_SM_XI0,
    MICROSCOPE_ETA,
    a1_ir_proca_replacement,
    a1_uv_sm,
    cassini_safe_version,
    compare_two_vs_three,
    fifth_force_bounds,
    self_checks,
    two_phase_formalization,
)


def test_self_checks():
    self_checks()


def test_uv_ir_a1():
    assert a1_uv_sm() == A1_SM_XI0 == F(-43, 12)
    ir = a1_ir_proca_replacement()
    assert ir["a1_IR_Proca"] == F(-1, 12)
    assert ir["sign_OK_IR"] is False


def test_formalization():
    tp = two_phase_formalization()
    assert tp["mass_gravity_locked_pair_derived"] is False
    assert tp["trace"]["T_trace_EM"] == 0


def test_cassini_safe_and_bounds():
    c = cassini_safe_version()
    assert c["light_bending_OK"] is True
    b = fifth_force_bounds()
    assert b["MICROSCOPE_eta"] == MICROSCOPE_ETA
    assert b["closes_Newton_gap"] is False


def test_vs_three_phase():
    cmp_ = compare_two_vs_three()
    assert cmp_["three_phase_28"]["derived"] is False
    assert "avoided" in cmp_["two_phase_this"]["photon_issue"]


if __name__ == "__main__":
    test_self_checks()
    test_uv_ir_a1()
    test_formalization()
    test_cassini_safe_and_bounds()
    test_vs_three_phase()
    print("test_ffz_two_phase: ALL PASS")
