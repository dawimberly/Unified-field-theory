#!/usr/bin/env python3
"""Unit test for ffz_ledger_confinement."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_ledger_confinement import (  # noqa: E402
    gap_shift_check,
    self_checks,
    strong_wz_analogue,
)


def test_self_checks():
    self_checks()


def test_ruled_out():
    s = strong_wz_analogue()
    assert s["derived_pole_shift_from_confinement"] is False
    g = gap_shift_check()
    assert g["delta_GAP"] == 0.0
    assert g["closes_gap"] is False


if __name__ == "__main__":
    test_self_checks()
    test_ruled_out()
    print("test_ffz_ledger_confinement: ALL PASS")
