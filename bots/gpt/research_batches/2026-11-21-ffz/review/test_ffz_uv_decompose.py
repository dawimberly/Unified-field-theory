#!/usr/bin/env python3
"""Unit test for ffz_uv_decompose."""

from __future__ import annotations

import pathlib
import sys
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_uv_decompose import decompose, minimal_flip, self_checks  # noqa: E402


def test_self_checks():
    self_checks()


def test_uv_dominates_proca_flips():
    d = decompose()
    assert d["dominant_negative"] == "UV_SM"
    assert d["uv_gauge_a1"] == F(-8)
    m = minimal_flip()
    assert m["option_A_all_Proca_12"]["flips"] is True
    assert m["option_B_conformal_Higgs"]["flips"] is False


if __name__ == "__main__":
    test_self_checks()
    test_uv_dominates_proca_flips()
    print("test_ffz_uv_decompose: ALL PASS")
