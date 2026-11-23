#!/usr/bin/env python3
"""Unit test for ffz_phase_lock."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_phase_lock import (  # noqa: E402
    GAP,
    phase_lock_scan,
    self_checks,
)


def test_self_checks():
    self_checks()


def test_tuning_required():
    r = phase_lock_scan()
    assert r["tuning_required"] is True
    assert r["derived_closes_gap"] is False
    assert r["honest_negative"] is True
    assert abs(r["contact_tuned"]["alpha"] - GAP) < 1e-15


if __name__ == "__main__":
    test_self_checks()
    test_tuning_required()
    print("test_ffz_phase_lock: ALL PASS")
