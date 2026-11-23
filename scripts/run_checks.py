"""Run research reproductions sequentially; fail on any failed assertion."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
repo = root.parent
ffz_tests = repo / "analysis" / "ffz" / "tests"

for name in [
    "scalar_full_response_scan",
    "scalar_source_projection_check",
    "metric_tensor_ward_check",
    "metric_spectral_check",
]:
    subprocess.run([sys.executable, str(root / (name + ".py"))], check=True)
    print(name + " PASS", flush=True)

# FFZ overnight / ladder unit tests (analysis/ffz/tests)
ffz_test_names = [
    "test_ffz_induced_collapse",
    "test_ffz_distortion_phase",
    "test_ffz_particle_plasticity",
    "test_ffz_sm_roster",
    "test_ffz_three_phase",
    "test_ffz_two_phase",
    "test_ffz_shared_boundary",
    "test_ffz_ir_roster",
    "test_ffz_trace_anomaly",
    "test_ffz_ledger_reconcile",
    "test_ffz_lambda_pressure",
    "test_ffz_string_kk_tower",
    "test_ffz_hadronic_a1",
    "test_ffz_ee_confinement",
    "test_ffz_running_gind",
    "test_ffz_confinement_snap",
    "test_ffz_curvature_match",
    "test_ffz_spectral_gind",
    "test_ffz_ledger_confinement",
    "test_ffz_casimir_lesage",
    "test_ffz_uv_decompose",
    "test_ffz_gut_spectrum",
    "test_ffz_uv_higgsing",
    "test_ffz_planck_stress",
    "test_ffz_uv_consistency",
    "test_ffz_minimal_spectrum",
    "test_ffz_positivity_region",
    "test_ffz_gap_sectors",
    "test_ffz_sec_collapse",
    "test_ffz_deconf_gind",
    "test_ffz_bounce_criterion",
    "test_ffz_surface_audit",
    "test_ffz_edge_modes",
    "test_ffz_phase_lock",
]
for name in ffz_test_names:
    path = ffz_tests / (name + ".py")
    subprocess.run([sys.executable, str(path)], check=True)
    print(name.replace("test_", "") + " PASS", flush=True)
