"""Test omitted smooth modes at the refined 76-parameter endpoint."""
import json,numpy as np
from smooth_hopf_variational_2026_11_01 import Model,embed,OUT
r=json.loads((OUT/'extended_order3_refined.json').read_text());p=np.array(r['parameters']);m=Model(56,24,48,3,2);large=Model(56,24,48,4,3);q=embed(p,m,large);g=large.jac(q)
newamp=[6+large.amp_exps.index(e) for e in large.amp_exps if e not in m.amp_exps];newrot=[large.offset+3*i+j for i,e in enumerate(large.exps) if e not in m.exps for j in range(3)]
ids=newamp+newrot;j=max(ids,key=lambda i:abs(g[i]));sign=-np.sign(g[j]);rows=[]
for grid in [(40,18,36),(56,24,48),(72,30,60)]:
 test=Model(*grid,order=4,amp_order=3);E=float(test.energy(q));gg=test.jac(q)
 for step in [1e-2,1e-3,1e-4]:
  pert=q.copy();pert[j]+=sign*step;Ep=float(test.energy(pert));rows.append(dict(grid=list(grid),step=float(sign*step),energy=E,derivative=float(gg[j]),perturbed_energy=Ep,energy_change=Ep-E))
res=dict(parameter_index=j,sector='amplitude' if j<large.offset else 'orientation',max_new_amplitude_gradient=float(abs(g[newamp]).max()),max_new_orientation_gradient=float(abs(g[newrot]).max()),rows=rows,limitation='New smooth modes outside the minimized family; energies are numerical quadratures, not interval-certified values.')
(OUT/'extended_omitted_mode_countercheck.json').write_text(json.dumps(res,indent=2));print(json.dumps(res),flush=True)
