"""Audit explicit Lref=10 polynomial matching; not an adopted UV action.
Run python experiments/matched_spin2_2026_10_30/check.py from repository root.
"""
from pathlib import Path
import importlib.util,json,platform
import numpy as np
import scipy
from scipy.optimize import root,brentq,minimize_scalar
p=Path(__file__).resolve().parents[1]/'spin2_stieltjes_2026_10_30/check.py'
spec=importlib.util.spec_from_file_location('spin2_reference',p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def coeff(L):
 z,c,_=m.spectrum(L)
 return m.C*np.dot(c,z*np.log(z))/24,m.C*np.dot(c,np.log(z))/240
ref1,ref2=coeff(10)
def shifts(L):
 a,b=coeff(L);return ref1-a,ref2-b
def F(q,L):
 a,b=shifts(L);return m.kernel(q,L)/q+a+b*q

def contour(L,vertices,n):
 q=np.concatenate([np.linspace(a,b,n,endpoint=False) for a,b in zip(vertices[:-1],vertices[1:])]+[np.array(vertices[-1:])])
 for iteration in range(24):
  f=F(q,L);mid=(q[:-1]+q[1:])/2;fm=F(mid,L)
  p1=np.angle(fm/f[:-1]);p2=np.angle(f[1:]/fm)
  bad=(abs(p1)>.10)|(abs(p2)>.10)|(abs(p1-p2)>.025)
  bad|=abs(f[1:]-f[:-1])>.25*np.minimum(abs(f[1:]),abs(f[:-1]))
  if not np.any(bad):break
  new=[]
  for j in range(len(mid)):
   new.append(q[j])
   if bad[j]:new.append(mid[j])
  new.append(q[-1]);q=np.array(new)
 else:raise RuntimeError('unresolved contour')
 f=F(q,L);steps=np.angle(f[1:]/f[:-1]);w=float(np.sum(steps)/(2*m.PI))
 return dict(initial_points_per_edge=n,final_points=len(q),winding=w,max_phase=float(max(abs(steps))),minimum_boundary_modulus=float(min(abs(f))))
def main():
 rows=[];counts=[];minima=[];bounds=[]
 for L in [3,5,30,100]:
  a,b=shifts(L)
  if L>10:
   k1,k2=coeff(L); aa=-b*k2;bb=b*k1-a*k2;cc=k1*ref1
   bounds.append(dict(L=L,quadratic=float(aa),linear=float(bb),constant=float(cc),discriminant=float(bb*bb-4*aa*cc)))
  if L<10:q=complex(brentq(lambda x:float(F(x,L).real),1,1000,xtol=1e-12))
  else:
   rr=root(lambda v:[float(F(complex(*v),L).real),float(F(complex(*v),L).imag)],[10,80],tol=1e-11)
   assert rr.success
   q=complex(*rr.x)
   opt=minimize_scalar(lambda x:float(F(x,L).real),bounds=(1,200),method='bounded',options={'xatol':1e-8})
   minima.append(dict(L=L,q=float(opt.x),F_min=float(opt.fun)))
  delta=a*q+b*q*q
  direct=[complex(m.direct(q,L,n)+delta) for n in [256,512,1024]]
  derivatives=[complex((F(q+h,L)-F(q-h,L))/(2*h)) for h in [.001,.0001,.00001]]
  rows.append(dict(L=L,q=[q.real,q.imag],delta_k1=float(a),delta_k2=float(b),F_residual=float(abs(F(q,L))),independent_integral_nodes=[256,512,1024],independent_K_residuals=[abs(x) for x in direct],integral_convergence=abs(direct[-1]-direct[-2]),K_prime_values=[[complex(q*d).real,complex(q*d).imag] for d in derivatives],inverse_K_residue=[complex(1/(q*derivatives[-1])).real,complex(1/(q*derivatives[-1])).imag],conjugation_error=float(abs(F(q.conjugate(),L)-F(q,L).conjugate())),expansion_parameter_abs_q_over_4=abs(q)/4,local_quadratic_ratio=float(abs(ref2*q/ref1))))
  radius=1
  v=[q-radius-1j*radius,q+radius-1j*radius,q+radius+1j*radius,q-radius+1j*radius,q-radius-1j*radius]
  for n in [128,256]:counts.append(dict(L=L,kind='isolating',**contour(L,v,n)))
 for L in [3,5,10,30,100]:
  for R in [200,1000,10000]:
   for n in [256,512]:
    h=.001;v=[-R+1j*h,R+1j*h,R+1j*R,-R+1j*R,-R+1j*h]
    counts.append(dict(L=L,kind='upper_rectangle',R=R,h=h,**contour(L,v,n)))
 result=dict(scope='Explicitly polynomial matched determinant; Lref=10. Not a physical UV completion or an EFT validity claim.',versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,longdouble_bits=np.finfo(m.LD).nmant),roots=rows,positive_real_minima=minima,cauchy_schwarz_lower_bound=bounds,contours=counts)
 (Path(__file__).resolve().parent/'results.json').write_text(json.dumps(result,indent=2))
 assert max(r['F_residual'] for r in rows)<1e-13
 assert max(max(r['independent_K_residuals']) for r in rows)<1e-9
 assert max(r['integral_convergence'] for r in rows)<1e-10
 assert all(abs(c['winding']-(1 if c['kind']=='isolating' or c['L']>10 else 0))<1e-8 for c in counts)
 assert all(c['max_phase']<.21 for c in counts)
 assert all(v['F_min']>0 for v in minima)
 assert all(v['quadratic']>0 and v['discriminant']<0 for v in bounds)
 print(json.dumps(dict(roots=rows,minima=minima,contour_checks=len(counts),max_independent_residual=max(max(r['independent_K_residuals']) for r in rows)),indent=2))
if __name__=='__main__':main()
