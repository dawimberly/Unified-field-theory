"""Run the four research reproductions sequentially; fail on any failed assertion."""
from pathlib import Path
import subprocess
import sys
root=Path(__file__).resolve().parent
for name in ['scalar_full_response_scan','scalar_source_projection_check','metric_tensor_ward_check','metric_spectral_check']:
    subprocess.run([sys.executable,str(root/(name+'.py'))],check=True)
    print(name+' PASS',flush=True)
