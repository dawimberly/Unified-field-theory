#!/usr/bin/env python3
"""Unit test for ffz_ir_roster (CI hook)."""

from __future__ import annotations

import pathlib
import sys
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_ir_roster import corrected_ladder, mass_switches_gravity_test, self_checks  # noqa: E402


def test_self_checks():
    self_checks()


def test_signs():
    rows = corrected_ladder()
    assert rows[0].a1 == F(-43, 12)
    assert rows[3].a1 == F(2, 3)
    assert rows[4].a1 == F(-1, 6)
    t = mass_switches_gravity_test()
    assert t["derived_positive"] is False
    assert t["hadronic_positive"] is True


if __name__ == "__main__":
    test_self_checks()
    test_signs()
    print("test_ffz_ir_roster: ALL PASS")
