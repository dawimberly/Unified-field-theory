#!/usr/bin/env python3
"""Unit test for ffz_particle_plasticity (CI hook)."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
FFZ = HERE.parent
for path in (FFZ, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from ffz_particle_plasticity import (  # noqa: E402
    GRAN_SASSO_R0_MIN,
    herschel_bulkley_tau,
    hysteresis_scan,
    kurkijarvi_switching,
    rate_dependence_analysis,
    r0_from_yield_check,
    self_checks,
    yield_thresholds,
)


def test_self_checks():
    self_checks()


def test_yield_robust_O1():
    yt = yield_thresholds()
    assert yt["robust_O1"] is True
    for v in yt["F_c_over_V0"].values():
        assert 0.1 < v < 2.0


def test_hysteresis_underdamped():
    rows = {h["beta_c"]: h for h in hysteresis_scan()}
    assert rows[0.1]["hysteretic"] is False
    assert rows[100.0]["hysteretic"] is True


def test_no_untuned_r0():
    r0 = r0_from_yield_check()
    assert r0["universal_from_C_ledger"] is False
    assert r0["a0_coincidence_only"] is True
    assert abs(r0["a0_ratio_to_bound"] - 1.0) < 0.05


def test_rate_thickening_and_kfd():
    rd = rate_dependence_analysis()
    assert rd["HB_thickening_n1_5"][-1] > rd["HB_thickening_n1_5"][0]
    assert rd["KFD_Isw"][-1] > rd["KFD_Isw"][1]
    assert herschel_bulkley_tau(1.0, 0.2, 1.5, 10.0) > herschel_bulkley_tau(
        1.0, 0.2, 1.5, 0.0
    )
    assert kurkijarvi_switching(1.0, 100.0, 1.0, 0.3) > kurkijarvi_switching(
        1.0, 0.1, 1.0, 0.3
    )


def test_gran_sasso_constant():
    assert GRAN_SASSO_R0_MIN == 5.4e-11


if __name__ == "__main__":
    test_self_checks()
    test_yield_robust_O1()
    test_hysteresis_underdamped()
    test_no_untuned_r0()
    test_rate_thickening_and_kfd()
    test_gran_sasso_constant()
    print("test_ffz_particle_plasticity: ALL PASS")
