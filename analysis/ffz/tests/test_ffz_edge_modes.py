#!/usr/bin/env python3
"""Unit test for ffz_edge_modes."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_edge_modes import (  # noqa: E402
    GAP,
    KAPPA_SCHEME,
    edge_mode_result,
    self_checks,
)


def test_self_checks():
    self_checks()


def test_delta_c_zero():
    r = edge_mode_result()
    assert r["delta_C_total"] == 0.0
    assert r["honest_negative"] is True
    assert abs(KAPPA_SCHEME + GAP) < 1e-15


if __name__ == "__main__":
    test_self_checks()
    test_delta_c_zero()
    print("test_ffz_edge_modes: ALL PASS")
