#!/usr/bin/env python3
"""Unit test for ffz_surface_audit."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_surface_audit import (  # noqa: E402
    GAP,
    audit_table,
    self_checks,
    term_euler,
    three_loop_setup,
)


def test_self_checks():
    self_checks()


def test_no_derived_gap_close():
    a = audit_table()
    assert a["derived_surface_TD_sum"] == 0.0
    assert a["honest_negative"] is True
    assert abs(a["vacuum_kappa"]["kappa"] + GAP) < 1e-12
    assert three_loop_setup()["geometric_boundary"] is False
    assert term_euler()["varies_to_Einstein"] is False


if __name__ == "__main__":
    test_self_checks()
    test_no_derived_gap_close()
    print("test_ffz_surface_audit: ALL PASS")
