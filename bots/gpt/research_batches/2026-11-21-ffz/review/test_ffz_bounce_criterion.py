#!/usr/bin/env python3
"""Unit test for ffz_bounce_criterion."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_bounce_criterion import (  # noqa: E402
    ac_ad_bounce,
    observables,
    self_checks,
    solar_mass_collapse,
)


def test_self_checks():
    self_checks()


def test_no_bounce_no_observables():
    b = ac_ad_bounce()
    assert b["AC_produces_bounce"] is False
    assert b["AD_produces_bounce"] is False
    assert observables()["GRB_FRB_echoes_derived_here"] is False
    c = solar_mass_collapse(10.0)
    assert 25 < c["R_s_km"] < 35


if __name__ == "__main__":
    test_self_checks()
    test_no_bounce_no_observables()
    print("test_ffz_bounce_criterion: ALL PASS")
