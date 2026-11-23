#!/usr/bin/env python3
"""Unit test for ffz_trace_anomaly."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_trace_anomaly import induced_G_from_anomaly, self_checks  # noqa: E402


def test_self_checks():
    self_checks()


def test_no_positive_1G():
    assert induced_G_from_anomaly()["positive_1G_derived"] is False


if __name__ == "__main__":
    test_self_checks()
    test_no_positive_1G()
    print("test_ffz_trace_anomaly: ALL PASS")
