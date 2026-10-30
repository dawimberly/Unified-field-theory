"""Independent numerical diagnostics for the analytic/exact cut count.
These floating-point paths are cross-checks, not the global certificate.
"""
from pathlib import Path
import importlib.util,json,warnings
import numpy as np
from scipy.integrate import quad
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def load(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
m=load(ROOT/'experiments/matched_inverse_reconstruction_2026_10_30/check.py','response')
LD=m.LD;CD=m.CD;PI=m.PI
cert=json.loads((HERE/'results.json').read_text())
WARNINGS=[]
RAW_R={3:1024,5:1664,10:3264,30:38656,100:128256}
MATCH_R={3:256,5:416,10:3264,30:2416,100:8016}
def raw_lower(model,s):
 if s>=model.switch:return model.raw_uv(CD(-s),logq=np.asarray(np.log(LD(s))-1j*PI,dtype=CD))
 return m.base.boundary(s,model.L)
def lower(model,s,matched):
 a=raw_lower(model,s)
 return a-model.d1*s+model.d2*s*s if matched else a

def winding(model,matched,R,n):
 theta=np.linspace(-PI,PI,n+1,dtype=LD);q=R*np.exp(1j*theta)
 a=model.A(q) if matched else model.raw_uv(q)
 a[0]=lower(model,R,matched);a[-1]=a[0].conjugate()
 edges=list(map(LD,4*model.z))+[LD(R)];ss=[]
 for lo,hi in zip(edges[:-1],edges[1:]):
  t=np.linspace(0,1,n//8+1,dtype=LD)
  ss.extend(lo+(hi-lo)*np.sin(PI*t/2)**2)
 lip=np.array([lower(model,s,matched) for s in ss],dtype=CD)
 path=np.concatenate([a,lip[::-1].conjugate(),lip]);path=np.append(path,path[0])
 steps=np.angle(path[1:]/path[:-1]);number=np.sum(steps)/(2*PI)
 return dict(nodes=n,R=R,winding=float(number),maximum_phase_step=float(np.max(abs(steps))),minimum_sampled_modulus=float(np.min(abs(path))))

def parameter_real(model,s):
 zz=[];errors=[]
 for z in model.z:
  z=float(z)
  if s>4*z:
   b=np.sqrt(1-4*z/s);pts=[(1-b)/2,(1+b)/2]
  elif s==4*z:pts=[.5]
  else:pts=[]
  if s>=4*z:
   beta=np.sqrt(1-4*z/s);r1=(1-beta)/2;r2=(1+beta)/2
   integrand=lambda x:np.log(s)+np.log(abs(x-r1))+np.log(abs(x-r2))
  else:integrand=lambda x:np.log(z-s*x*(1-x))
  with warnings.catch_warnings(record=True) as caught:
   warnings.simplefilter('always')
   ans,err=quad(integrand,0,1,points=pts,epsabs=2e-11,epsrel=2e-11,limit=400)
  for warning in caught:WARNINGS.append(dict(L=model.L,s=s,z=z,message=str(warning.message)))
  zz.append(ans);errors.append(err)
 U=model.z+LD(s)/2
 value=np.dot(model.c,U*U*np.array(zz,dtype=LD)-U*model.z*np.log(model.z))/4
 error=np.dot(abs(model.c),U*U*np.array(errors,dtype=LD))/4
 return value,error

def main():
 contours=[];integrals=[]
 for L in [3,5,10,30,100]:
  model=m.Model(L)
  for matched in [False,True]:
   expected=next(x['total_A_zeros_counting_multiplicity'] for x in cert['global_counts'] if x['L']==L and x['matched_to_L10']==matched)
   for n in [512,2048,8192]:
    R=2*(MATCH_R[L] if matched else RAW_R[L]);r=winding(model,matched,R,n)
    assert abs(r['winding']-expected)<1e-10 and r['maximum_phase_step']<PI/2
    contours.append(dict(L=L,matched=matched,expected=expected,**r))
  inputs=[(float(x['s']),x['real_A_over_C']['False']) for x in cert['thresholds'] if x['L']==L]
  inputs += [(sum(v['s']['decimal'])/2,v['real_A_over_C']['False']) for x in cert['candidate_segments'] if x['L']==L for v in x['candidates']]
  for s,box in inputs:
   value,error=parameter_real(model,s);mid=sum(box['decimal'])/2
   rel=float(abs(value-mid)/max(1,abs(mid)))
   assert rel<1e-8
   integrals.append(dict(L=L,s=s,parameter_real_A_over_C=float(value),interval_midpoint=mid,relative_discrepancy=rel,quadrature_error_estimate=float(error)))
 out=dict(scope='Diagnostics only; exact/analytic proof is in check.py and REPORT.md.',numpy=np.__version__,longdouble_bits=np.finfo(LD).nmant,quadrature_warnings=WARNINGS,slit_windings=contours,independent_real_parameter_integrals=integrals)
 (HERE/'numerical_results.json').write_text(json.dumps(out,indent=2))
 print(json.dumps(dict(windings=len(contours),integrals=len(integrals),maximum_phase_step=max(x['maximum_phase_step'] for x in contours),max_integral_relative_error=max(x['relative_discrepancy'] for x in integrals))))
if __name__=='__main__':main()
