"""Resolution-refined countercheck of the old radius-only equilibrium claim."""
import json
import numpy as np
from smooth_hopf_variational_2026_11_01 import Model,OUT
base=json.loads((OUT/'radial_control.json').read_text());p=np.array(base['parameters']);m=Model(order=1);g=m.jac(p);j=6+int(np.argmax(abs(g[6:])));sign=-np.sign(g[j]);rows=[]
for grid in [(28,12,24),(40,18,36),(56,24,48)]:
 mf=Model(*grid,order=1);z=p.astype(complex);z[j]+=1e-25j;der=float(mf.energy(z).imag/1e-25)
 E=float(mf.energy(p));q=p.copy();q[j]+=sign*.01;Epert=float(mf.energy(q))
 rows.append(dict(grid=grid,base_energy=E,orientation_derivative=der,perturbed_energy=Epert,energy_change=Epert-E))
assert all(r['energy_change']<-.1 for r in rows)
res=dict(parameter_index=j,monomial_exponents=m.exps[(j-6)//3],rotation_component=(j-6)%3,parameter_step=float(sign*.01),rows=rows,interpretation='A smooth topology-preserving direction lowers energy. The radius/amplitude-only optimum is not stationary in the enlarged family. Numerical quadrature refinement, not an interval-certified continuum proof.')
(OUT/'radial_nonstationarity_countercheck.json').write_text(json.dumps(res,indent=2));print(json.dumps(res))
