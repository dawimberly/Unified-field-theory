#!/usr/bin/env python3
"""Unit test for ffz_induced_collapse (CI hook)."""

from __future__ import annotations

import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from ffz_induced_collapse import (  # noqa: E402
    C_LEDGER,
    C_TARGET,
    case_bec_rb87,
    case_flux_qubit,
    case_nanosphere,
    g_ind_scaling,
    gran_sasso_2020_note,
    induced_R0_smearing,
    induced_vs_target,
    phase_framing,
    self_checks,
    tau_DP,
    tau_far,
)


def test_self_checks():
    self_checks()


def test_tau_ratio_direction():
    """Larger C ⇒ smaller G_ind ⇒ longer τ."""
    sc = g_ind_scaling(C_LEDGER, C_TARGET)
    assert C_LEDGER > C_TARGET
    assert sc["G_ind_over_G_ref"] < 1.0
    assert sc["tau_over_tau_ref"] > 1.0
    assert abs(sc["tau_over_tau_ref"] - C_LEDGER / C_TARGET) < 1e-12


def test_induced_vs_target_number():
    inv = induced_vs_target()
    assert abs(inv["tau_induced_over_tau_1_24"] - 0.05131719 * 24.0) < 1e-8


def test_gran_sasso_bound_present():
    gs = gran_sasso_2020_note()
    assert gs["R0_lower_bound_m"] == 5.4e-11


def test_r0_not_natural():
    assert induced_R0_smearing()["naturally_supplies_viable_R0"] is False


def test_phase_is_rewrite():
    assert "rewrite" in phase_framing()["rewrite_or_content"]


def test_cases_gravity_vs_coherence():
    bec = case_bec_rb87()
    flux = case_flux_qubit()
    nano = case_nanosphere()
    # BEC / flux: gravity much slower than coherence (with stated assumptions)
    assert bec.tau_induced > bec.thermal_or_coherence_s
    assert flux.tau_induced > flux.thermal_or_coherence_s
    # Nanosphere: τ in a regime that can compete (order ms–s depending on d)
    assert nano.tau_induced > 0
    assert tau_DP(nano.M, nano.R0, nano.d) == nano.tau_G_N


def test_far_limit_matches():
    M, R0 = 1e-14, 1e-6
    # λ = d/(2R0) = 5e5 ⇒ f = 1 − 5/(12λ) ≈ 1 − 8e-7; allow that correction
    assert abs(tau_DP(M, R0, d=1.0) / tau_far(M, R0) - 1.0) < 1e-5


if __name__ == "__main__":
    test_self_checks()
    test_tau_ratio_direction()
    test_induced_vs_target_number()
    test_gran_sasso_bound_present()
    test_r0_not_natural()
    test_phase_is_rewrite()
    test_cases_gravity_vs_coherence()
    test_far_limit_matches()
    print("test_ffz_induced_collapse: ALL PASS")
