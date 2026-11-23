#!/usr/bin/env python3
"""Unit test for ffz_uv_higgsing."""

from __future__ import annotations

import pathlib
import sys
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_uv_higgsing import (  # noqa: E402
    DELTA_PER_VECTOR,
    mechanism_stueckelberg,
    self_checks,
    survey,
)


def test_self_checks():
    self_checks()


def test_survivors_and_stueckelberg_fail():
    assert DELTA_PER_VECTOR == F(7, 6)
    s = survey()
    assert s["n_survivors"] == 4
    assert mechanism_stueckelberg()["survives_U"] is False
    assert "asymp_safety_or_condensation" not in s["survivors"]


if __name__ == "__main__":
    test_self_checks()
    test_survivors_and_stueckelberg_fail()
    print("test_ffz_uv_higgsing: ALL PASS")
