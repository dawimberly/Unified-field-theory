#!/usr/bin/env python3
"""Unit test for ffz_hadronic_a1."""

from __future__ import annotations

import pathlib
import sys
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_hadronic_a1 import (  # noqa: E402
    confinement_removes_gluon_kabat,
    controlled_hadronic_a1,
    self_checks,
)


def test_self_checks():
    self_checks()


def test_band_positive():
    h = controlled_hadronic_a1()
    assert h["a1_low"] == F(1, 6)
    assert h["a1_high"] == F(7, 6)
    assert h["sign_OK_entire_band"] is True
    assert confinement_removes_gluon_kabat()["Kabat_8_absent_from_IR_sum"] is True
    assert confinement_removes_gluon_kabat()["complete_UV_IR_matching"] is False


if __name__ == "__main__":
    test_self_checks()
    test_band_positive()
    print("test_ffz_hadronic_a1: ALL PASS")
