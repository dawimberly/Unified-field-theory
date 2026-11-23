#!/usr/bin/env python3
"""Unit test for ffz_planck_stress."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_planck_stress import (  # noqa: E402
    all_stress_tests,
    self_checks,
    stress_complete_higgsing_short,
    stress_complete_higgsing_to_EW,
)


def test_self_checks():
    self_checks()


def test_short_fails_long_proca_flips():
    assert stress_complete_higgsing_short()["sign_positive"] is False
    assert stress_complete_higgsing_to_EW()["sign_positive"] is True
    a = all_stress_tests()
    assert "SO10_heavy_Proca" in a["flip_names"]
    assert "E6_heavy_Proca" in a["flip_names"]


if __name__ == "__main__":
    test_self_checks()
    test_short_fails_long_proca_flips()
    print("test_ffz_planck_stress: ALL PASS")
