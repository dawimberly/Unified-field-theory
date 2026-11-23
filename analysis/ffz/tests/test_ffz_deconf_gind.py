#!/usr/bin/env python3
"""Unit test for ffz_deconf_gind."""

from __future__ import annotations

import pathlib
import sys
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_deconf_gind import (  # noqa: E402
    a1_across_deconfinement,
    induced_g_shift,
    self_checks,
)


def test_self_checks():
    self_checks()


def test_local_flip_tiny_G():
    a = a1_across_deconfinement()
    assert a["a1_deconfined"] == F(-5)
    assert a["sign_flip_local_a1"] is True
    g = induced_g_shift()
    assert g["G_ind_Planck_sign_flips"] is False
    assert g["frac_shift_vs_Planck_G"] < 1e-35


if __name__ == "__main__":
    test_self_checks()
    test_local_flip_tiny_G()
    print("test_ffz_deconf_gind: ALL PASS")
