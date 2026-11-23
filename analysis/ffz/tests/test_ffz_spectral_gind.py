#!/usr/bin/env python3
"""Unit test for ffz_spectral_gind."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_spectral_gind import (  # noqa: E402
    integrated_sign,
    planck_cutoff_needed,
    self_checks,
)


def test_self_checks():
    self_checks()


def test_total_negative_needs_cutoff():
    s = integrated_sign()
    assert s["total_sign_positive_central"] is False
    assert s["S_central"] < 0
    p = planck_cutoff_needed()
    assert p["cutoff_assumed"] is True
    assert p["near_MPl_without_cutoff"] is False


if __name__ == "__main__":
    test_self_checks()
    test_total_negative_needs_cutoff()
    print("test_ffz_spectral_gind: ALL PASS")
