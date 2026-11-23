#!/usr/bin/env python3
"""Unit test for ffz_gut_spectrum."""

from __future__ import annotations

import pathlib
import sys
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_gut_spectrum import (  # noqa: E402
    self_checks,
    so10_spectra,
    spectral_with_gut,
    su5_spectra,
)


def test_self_checks():
    self_checks()


def test_broken_positive_unbroken_negative():
    s5 = su5_spectra()
    assert s5["sign_OK_unbroken"] is False
    assert s5["sign_OK_broken"] is True
    assert s5["broken_to_SM"]["a1_net"] == F(41, 12)
    s10 = so10_spectra()
    assert s10["sign_OK_unbroken"] is False
    assert s10["sign_OK_broken"] is True
    g = spectral_with_gut("SU5", include_higgs=False)
    assert g["cases"]["M_GUT_eq_MU_UV"]["flips"] is True
    assert g["flip_with_wide_Planck_bin"] is False


if __name__ == "__main__":
    test_self_checks()
    test_broken_positive_unbroken_negative()
    print("test_ffz_gut_spectrum: ALL PASS")
