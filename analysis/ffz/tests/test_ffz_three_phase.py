#!/usr/bin/env python3
"""Unit test for ffz_three_phase (CI hook)."""

from __future__ import annotations

import pathlib
import sys
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_three_phase import (  # noqa: E402
    A1_SM_XI0,
    CASSINI_GAMMA_MINUS_1,
    experimental_bounds,
    majorana_a1_reconciliation,
    self_checks,
    three_phase_formalization,
    z3_boundary_sign_test,
)


def test_self_checks():
    self_checks()


def test_majorana_reconciliation():
    r = majorana_a1_reconciliation()
    assert r["net_N1_EE_PR22"] == F(-1, 3)
    assert r["net_N1_Sakharov_PR27"] == F(-7, 12)
    assert r["for_induced_G_use"] == F(-7, 12)


def test_three_phase_not_derived():
    tp = three_phase_formalization()
    assert tp["E_B_G_120_derived"] is False
    assert tp["trace"]["T_trace_EM"] == 0
    assert tp["trace"]["clash_light_bending"] is True


def test_z3_no_free_flip():
    z = z3_boundary_sign_test()
    assert z["a1_SM_baseline"] == A1_SM_XI0 == F(-43, 12)
    assert z["derived_without_assumption_flips_sign"] is False
    assert z["option_A_shared_Kabat_a1"] == F(53, 12)


def test_cassini():
    ex = experimental_bounds()
    assert ex["Cassini_|gamma-1|_max"] == CASSINI_GAMMA_MINUS_1 == 2.3e-5
    assert ex["photon_decoupling_ruled_out"] is True


if __name__ == "__main__":
    test_self_checks()
    test_majorana_reconciliation()
    test_three_phase_not_derived()
    test_z3_no_free_flip()
    test_cassini()
    print("test_ffz_three_phase: ALL PASS")
