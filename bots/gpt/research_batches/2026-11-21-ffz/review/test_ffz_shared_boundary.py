#!/usr/bin/env python3
"""Unit test for ffz_shared_boundary (CI hook)."""

from __future__ import annotations

import pathlib
import sys
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_shared_boundary import (  # noqa: E402
    A1_SM_IR_PROCA,
    derived_a1_with_edge_counting,
    donnelly_wall_edge_modes,
    kabat_contact_per_sector,
    self_checks,
    shared_kabat_53_12_check,
)


def test_self_checks():
    self_checks()


def test_no_sharing():
    dw = donnelly_wall_edge_modes()
    assert dw["EM_shares_with_SU3"] is False
    k = kabat_contact_per_sector()
    assert k["UV_contacts"]["total"] == 12
    assert k["IR_contacts"]["total"] == 9


def test_53_12_not_derived():
    assert shared_kabat_53_12_check()["derivable"] is False


def test_ir_no_flip():
    a = derived_a1_with_edge_counting()
    assert a["a1_IR_derived"] == A1_SM_IR_PROCA == F(-1, 12)
    assert a["flip_derived"] is False


if __name__ == "__main__":
    test_self_checks()
    test_no_sharing()
    test_53_12_not_derived()
    test_ir_no_flip()
    print("test_ffz_shared_boundary: ALL PASS")
