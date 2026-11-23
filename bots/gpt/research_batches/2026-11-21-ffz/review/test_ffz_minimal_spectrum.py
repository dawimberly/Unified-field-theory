#!/usr/bin/env python3
"""Unit test for ffz_minimal_spectrum."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_minimal_spectrum import (  # noqa: E402
    minimal_n_v,
    self_checks,
    so10_e6_status,
    two_loop_unification_sketch,
)


def test_self_checks():
    self_checks()


def test_nmin_and_not_unique():
    m = minimal_n_v()
    assert m["n_v_min"] == 15
    assert m["S_at_min"] > 0
    st = so10_e6_status(m["n_v_min"])
    assert st["selected_uniquely"] is False
    assert st["SO10_allowed"] is True
    assert two_loop_unification_sketch()["SM_unifies"] is False


if __name__ == "__main__":
    test_self_checks()
    test_nmin_and_not_unique()
    print("test_ffz_minimal_spectrum: ALL PASS")
