#!/usr/bin/env python3
"""Unit test for ffz_gap_sectors."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_gap_sectors import delta_C_by_sector, self_checks  # noqa: E402


def test_self_checks():
    self_checks()


def test_honest_negative():
    d = delta_C_by_sector()
    assert d["honest_negative"] is True
    assert d["delta_C_total"] == 0.0
    assert d["closes_gap"] is False


if __name__ == "__main__":
    test_self_checks()
    test_honest_negative()
    print("test_ffz_gap_sectors: ALL PASS")
