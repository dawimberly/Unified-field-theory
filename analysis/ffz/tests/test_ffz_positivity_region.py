#!/usr/bin/env python3
"""Unit test for ffz_positivity_region."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_positivity_region import (  # noqa: E402
    freedom_report,
    scan_region,
    self_checks,
)


def test_self_checks():
    self_checks()


def test_region_not_unique():
    sc = scan_region()
    lo, hi = sc["n_v_fiducial_allowed"]
    assert lo == 15
    assert hi - lo >= 10
    fr = freedom_report(sc)
    assert fr["is_fine_tuning_claim"] is False


if __name__ == "__main__":
    test_self_checks()
    test_region_not_unique()
    print("test_ffz_positivity_region: ALL PASS")
