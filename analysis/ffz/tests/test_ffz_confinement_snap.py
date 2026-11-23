#!/usr/bin/env python3
"""Unit test for ffz_confinement_snap."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_confinement_snap import assessment, self_checks  # noqa: E402


def test_self_checks():
    self_checks()


def test_not_derived_map():
    assert assessment()["maps_to_FFZ_collapse_derived"] is False


if __name__ == "__main__":
    test_self_checks()
    test_not_derived_map()
    print("test_ffz_confinement_snap: ALL PASS")
