"""Run research reproductions sequentially; fail on any failed assertion."""
from pathlib import Path
import subprocess
import sys
root=Path(__file__).resolve().parent
repo=root.parent
for name in ['scalar_full_response_scan','scalar_source_projection_check','metric_tensor_ward_check','metric_spectral_check']:
    subprocess.run([sys.executable,str(root/(name+'.py'))],check=True)
    print(name+' PASS',flush=True)
# Penrose–Diósi / induced-G collapse checks
collapse_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_induced_collapse.py'
subprocess.run([sys.executable,str(collapse_test)],check=True)
print('ffz_induced_collapse PASS',flush=True)
distortion_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_distortion_phase.py'
subprocess.run([sys.executable,str(distortion_test)],check=True)
print('ffz_distortion_phase PASS',flush=True)
plasticity_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_particle_plasticity.py'
subprocess.run([sys.executable,str(plasticity_test)],check=True)
print('ffz_particle_plasticity PASS',flush=True)
sm_roster_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_sm_roster.py'
subprocess.run([sys.executable,str(sm_roster_test)],check=True)
print('ffz_sm_roster PASS',flush=True)
three_phase_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_three_phase.py'
subprocess.run([sys.executable,str(three_phase_test)],check=True)
print('ffz_three_phase PASS',flush=True)
two_phase_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_two_phase.py'
subprocess.run([sys.executable,str(two_phase_test)],check=True)
print('ffz_two_phase PASS',flush=True)
shared_boundary_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_shared_boundary.py'
subprocess.run([sys.executable,str(shared_boundary_test)],check=True)
print('ffz_shared_boundary PASS',flush=True)
ir_roster_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_ir_roster.py'
subprocess.run([sys.executable,str(ir_roster_test)],check=True)
print('ffz_ir_roster PASS',flush=True)
trace_anomaly_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_trace_anomaly.py'
subprocess.run([sys.executable,str(trace_anomaly_test)],check=True)
print('ffz_trace_anomaly PASS',flush=True)
ledger_reconcile_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_ledger_reconcile.py'
subprocess.run([sys.executable,str(ledger_reconcile_test)],check=True)
print('ffz_ledger_reconcile PASS',flush=True)
lambda_pressure_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_lambda_pressure.py'
subprocess.run([sys.executable,str(lambda_pressure_test)],check=True)
print('ffz_lambda_pressure PASS',flush=True)
string_kk_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_string_kk_tower.py'
subprocess.run([sys.executable,str(string_kk_test)],check=True)
print('ffz_string_kk_tower PASS',flush=True)
hadronic_a1_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_hadronic_a1.py'
subprocess.run([sys.executable,str(hadronic_a1_test)],check=True)
print('ffz_hadronic_a1 PASS',flush=True)
ee_confinement_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_ee_confinement.py'
subprocess.run([sys.executable,str(ee_confinement_test)],check=True)
print('ffz_ee_confinement PASS',flush=True)
running_gind_test=repo/'bots'/'gpt'/'research_batches'/'2026-11-21-ffz'/'review'/'test_ffz_running_gind.py'
subprocess.run([sys.executable,str(running_gind_test)],check=True)
print('ffz_running_gind PASS',flush=True)
