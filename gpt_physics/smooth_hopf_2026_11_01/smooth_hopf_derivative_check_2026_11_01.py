"""Independent spatial-derivative control for the smooth Hopf family."""
import json
import numpy as np
from smooth_hopf_variational_2026_11_01 import Model,OUT
m=Model(12,6,12,2);rng=np.random.default_rng(17);p=np.zeros(m.size);p[1:6]=[.1,.03,.01,.005,0];p[6:]=.05*rng.normal(size=m.size-6)
def direct(y):
 r2=np.sum(y*y,axis=1);den=1+r2;z1=2*(y[:,0]+1j*y[:,1])/den;z2=(2*y[:,2]+1j*(r2-1))/den;z=z1*np.conj(z2)
 n0=np.stack([2*z.real,2*z.imag,abs(z1)**2-abs(z2)**2],axis=-1)
 B=np.stack([np.prod(y**np.array(e),axis=1)*np.exp(-r2/4) for e in m.exps],axis=1);v=B@p[6:].reshape(-1,3)
 # Independent quaternion rotation expression.
 q0=1/np.sqrt(1+np.sum(v*v,axis=1));q=v*q0[:,None]
 nt=n0+2*np.cross(q,np.cross(q,n0)+q0[:,None]*n0)
 H=np.exp(-r2[:,None]/m.widths[None,:]**2);rho=np.exp(-H@p[1:6])
 return rho,nt
rho,dr,n,G=m.fields(p);eligible=np.where(np.linalg.norm(m.y,axis=1)<3)[0];ids=rng.choice(eligible,20,replace=False)
errs=[]
for j in range(3):
 v=np.zeros(3);v[j]=1e-5
 rp,np_=direct(m.y[ids]+v);rm,nm=direct(m.y[ids]-v)
 errs.append(dict(axis=j,amplitude_max_error=float(np.max(abs((rp-rm)/2e-5-dr[ids,j]))),orientation_max_error=float(np.max(abs((np_-nm)/2e-5-G[ids,j])))))
assert max(max(e['amplitude_max_error'],e['orientation_max_error']) for e in errs)<2e-8
(OUT/'independent_spatial_derivative_checks.json').write_text(json.dumps(errs,indent=2));print(json.dumps(errs))
