"""Enlarge the topology-preserving family and check integration/basis convergence."""
import json,time
import numpy as np
from smooth_hopf_variational_2026_11_01 import Model,embed,optimize,OUT
old=Model(order=2);r=json.loads((OUT/'order2.json').read_text());p=np.array(r['parameters'])
m=Model(40,18,36,3,2);q=embed(p,old,m)
g=m.jac(q);print(json.dumps(dict(stage='extended_start',parameters=m.size,energy=float(m.energy(q)),max_gradient=float(abs(g).max()),max_new_amplitude_gradient=float(abs(g[6:m.offset]).max()))),flush=True)
q,res=optimize(m,q,400,'extended_order3')
checks=[]
for grid in [(40,18,36),(56,24,48),(72,30,60)]:
 mf=Model(*grid,order=3,amp_order=2);checks.append(dict(grid=list(grid),energy=float(mf.energy(q)),max_gradient=float(abs(mf.jac(q)).max())))
print(json.dumps(dict(stage='refinement',rows=checks)),flush=True)
# Reoptimize if the finer integration changes the stationarity result.
if checks[1]['max_gradient']>1e-3:
 mf=Model(56,24,48,3,2);q,res=optimize(mf,q,300,'extended_order3_refined');m=mf
 fine=Model(72,30,60,3,2);checks.append(dict(grid=fine.grid,energy=float(fine.energy(q)),max_gradient=float(abs(fine.jac(q)).max()),after_refinement=True))
# A larger family provides a first-variation test outside the optimized space.
large=Model(56,24,48,4,3);ql=embed(q,m,large);gl=large.jac(ql)
newamps=[6+large.amp_exps.index(e) for e in large.amp_exps if e not in m.amp_exps]
newrots=[large.offset+3*i+j for i,e in enumerate(large.exps) if e not in m.exps for j in range(3)]
# Restricted Hessian. It is not a full field Hessian or a physical frequency spectrum.
ep=1e-4;H=np.empty((m.size,m.size));start=time.time()
for j in range(m.size):
 v=np.zeros(m.size);v[j]=ep;H[:,j]=(m.jac(q+v)-m.jac(q-v))/(2*ep)
 if j and j%20==0:print(json.dumps(dict(stage='restricted_hessian',column=j,columns=m.size,seconds=time.time()-start)),flush=True)
asym=float(np.linalg.norm(H-H.T)/max(1,np.linalg.norm(H)));ev=np.linalg.eigvalsh((H+H.T)/2)
summary=dict(final_label=res['label'],parameter_count=m.size,quadrature_grid=m.grid,energy=float(m.energy(q)),max_parameter_gradient=float(abs(m.jac(q)).max()),bound_active=res['bound_active'],quadrature_refinement=checks,new_amplitude_gradient=float(abs(gl[newamps]).max()) if newamps else 0,new_orientation_gradient=float(abs(gl[newrots]).max()) if newrots else 0,larger_family_energy=float(large.energy(ql)),hessian_symmetry_relative_error=asym,restricted_hessian_eigenvalues=ev.tolist(),negative_restricted_eigenvalues_below_minus_1e_5=int(np.sum(ev < -1e-5)),zero_like_restricted_eigenvalues=int(np.sum(abs(ev)<1e-5)),hessian_difference_step=ep,topology='Charge one by continuous smooth homotopy; no node-wise unwinding.',limitation='Positive restricted Hessian does not establish full PDE stability. New-mode gradients assess remaining omitted directions; quadrature-refinement error is not a rigorous bound.')
np.savez_compressed(OUT/'extended_final.npz',parameters=q,hessian=H)
(OUT/'extended_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
