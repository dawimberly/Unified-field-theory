#!/usr/bin/env python3
"""Unit test for ffz_ee_confinement."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_ee_confinement import self_checks, supports_kabat_removal  # noqa: E402


def test_self_checks():
    self_checks()


def test_support_not_proof():
    s = supports_kabat_removal()
    assert s["supports_gluon_Kabat_absence"] is True
    assert s["proves_a1_matching"] is False


if __name__ == "__main__":
    test_self_checks()
    test_support_not_proof()
    print("test_ffz_ee_confinement: ALL PASS")
