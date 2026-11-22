import json,numpy as np
from smooth_hopf_variational_2026_11_01 import Model,OUT
z=np.load(OUT/'extended_final.npz');p=z['parameters'];H=z['hessian'];ev,V=np.linalg.eigh((H+H.T)/2);d=V[:,0];m=Model(56,24,48,3,2);E=float(m.energy(p));rows=[]
for ep in [1e-3,5e-4]:
 Ep=float(m.energy(p+ep*d));Em=float(m.energy(p-ep*d));rows.append(dict(step=ep,symmetric_curvature=(Ep+Em-2*E)/ep**2))
res=dict(smallest_restricted_hessian_eigenvalue=float(ev[0]),direct_energy_curvature_checks=rows)
assert all(abs(r['symmetric_curvature']-ev[0])<2e-5 for r in rows)
(OUT/'smallest_mode_direct_check.json').write_text(json.dumps(res,indent=2));print(json.dumps(res))
