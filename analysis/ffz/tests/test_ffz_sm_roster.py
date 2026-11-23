#!/usr/bin/env python3
"""Unit test for ffz_sm_roster (CI hook)."""

from __future__ import annotations

import pathlib
import sys
from fractions import Fraction as F

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_sm_roster import (  # noqa: E402
    GAP,
    a1_massless_vector,
    a1_net_mssm,
    a1_net_sm,
    a1_real_scalar,
    a1_weyl,
    ffz_channel_map,
    gap_closure_check,
    self_checks,
    sm_dof,
)


def test_self_checks():
    self_checks()


def test_species_weights():
    assert a1_real_scalar(F(1, 6)) == 0
    assert a1_weyl() == F(1, 12)
    assert a1_massless_vector() == F(-2, 3)


def test_sm_wrong_sign():
    assert a1_net_sm(0, False)["a1_net"] == F(-43, 12)
    assert a1_net_sm(F(1, 6), False)["sign_OK_for_G"] is False
    assert a1_net_mssm()["a1_net"] == F(-35, 12)


def test_dof_counts():
    assert sm_dof(False)["n_weyl"] == 45
    assert sm_dof(True)["n_weyl"] == 48
    assert sm_dof(False)["n_massless_vector"] == 12
    assert sm_dof(False)["n_real_scalar"] == 4


def test_channel_and_gap():
    ch = ffz_channel_map()
    assert ch["channels_same"] is False
    assert ch["delta_C_ledger_from_free_roster"] == 0.0
    g = gap_closure_check()
    assert g["derived_closes_gap"] is False
    assert abs(GAP - 0.00965052) < 1e-5


if __name__ == "__main__":
    test_self_checks()
    test_species_weights()
    test_sm_wrong_sign()
    test_dof_counts()
    test_channel_and_gap()
    print("test_ffz_sm_roster: ALL PASS")
