#!/usr/bin/env python3
"""Unit test for ffz_lambda_pressure."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_lambda_pressure import dark_energy_without_tuning, lambda_ind_sign, self_checks  # noqa: E402


def test_self_checks():
    self_checks()


def test_no_DE():
    assert lambda_ind_sign()["matches_observed_DE_positive"] is False
    assert dark_energy_without_tuning()["DE_from_roster_without_tuning"] is False


if __name__ == "__main__":
    test_self_checks()
    test_no_DE()
    print("test_ffz_lambda_pressure: ALL PASS")
