#!/usr/bin/env python3
"""Unit test for ffz_uv_consistency."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_uv_consistency import (  # noqa: E402
    ir_massless_crossover,
    self_checks,
    survivor_scorecard,
)


def test_self_checks():
    self_checks()


def test_proca_ew_ruled_out_so10_ok():
    ir = ir_massless_crossover()
    assert ir["complete_Proca_to_EW_OK_for_QCD"] is False
    s = survivor_scorecard()
    assert s["SO10_E6"]["overall_viable_as_assumption"] is True
    assert s["complete_Proca_to_EW"]["overall_viable_as_assumption"] is False


if __name__ == "__main__":
    test_self_checks()
    test_proca_ew_ruled_out_so10_ok()
    print("test_ffz_uv_consistency: ALL PASS")
