#!/usr/bin/env python3
"""Unit test for ffz_curvature_match."""

from __future__ import annotations

import pathlib
import sys
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_curvature_match import (  # noqa: E402
    matched_band,
    self_checks,
    soft_wall_glueball,
)


def test_self_checks():
    self_checks()


def test_narrowed_positive():
    m = matched_band()
    assert m["a1_central"] == F(1, 3)
    assert m["a1_low"] == F(1, 6)
    assert m["a1_high"] == F(3, 4)
    assert m["sign_OK_entire_band"] is True
    assert m["width_matched"] < m["width_36"]
    assert soft_wall_glueball()["restores_8_Kabat"] is False


if __name__ == "__main__":
    test_self_checks()
    test_narrowed_positive()
    print("test_ffz_curvature_match: ALL PASS")
