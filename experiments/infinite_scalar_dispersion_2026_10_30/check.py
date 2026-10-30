"""Infinite-cut scalar dispersion and matched inverse moments, with exact UV bounds.
Run python experiments/infinite_scalar_dispersion_2026_10_30/check.py.
Numerical quadrature/residues are diagnostics; tail bounds use exact inequalities.
"""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,platform,warnings
import numpy as np
import scipy
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent

def load(path,name):
 s=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
m=load('experiments/matched_inverse_reconstruction_2026_10_30/check.py','model')
iv=load('experiments/matched_inverse_reconstruction_2026_10_30/local_pole_certificate.py','interval')
uv=m.uv;LD=m.LD;CD=m.CD;PI=m.PI;C=m.C
WARNINGS=[]
def ratstr(x):return str(x)
def interval_sum(terms):
 ans=(F(0),F(0))
 for x in terms:ans=uv.add(ans,x)
 return ans
def coefficients(L):
 z,c=uv.m.spec(L)
 t0=interval_sum(uv.scale(uv.log_bounds(a),b) for a,b in zip(z,c))
 t1=interval_sum(uv.scale(uv.log_bounds(a),b*a) for a,b in zip(z,c))
 return uv.scale(t1,F(-1,8)),uv.scale(t0,F(1,16))
def bounds(L,matched,v):
 R=F(v['R']);a1,a2=coefficients(L);r1,r2=coefficients(10)
 d1=uv.add(r1,uv.scale(a1,-1));d2=uv.add(r2,uv.scale(a2,-1))
 # |Im(A/C)/pi| <= B/s for every s>=R. Each subleading term*s decreases.
 lead=abs(v['leading_cut_over_pi'])*R
 B=2*lead-R*v['cut_sign_margin'];assert B>0 and v['cut_sign_margin']>0
 if matched and L!=10:
  power=5
  margin=uv.minimum(d2)-uv.absolute(d1)/R-(uv.absolute(v['constant'])+v['raw_difference_bound'])/R**2
 else:
  power=1; margin=v['raw_no_zero_margin']
 assert margin>0
 lower=iv.iv(margin)[0];upper=iv.iv(B)[1];assert lower>0
 invCupper=16*iv.PI[1]**2
 K=iv.iv(invCupper*upper/lower**2)[1]
 return dict(L=L,matched_to_L10=matched,R=int(R),rho_decay_power=power,scaled_kernel_lower_coefficient=str(lower),cut_imag_B_upper=str(upper),rho_K_upper=str(K),rho_K_decimal=float(K),interpretation='For all s>=R, 0<rho(s)<=K/s^power. Exact Fractions, pi enclosure, and monotonic UV bounds; not a quadrature error certificate.')

class Response:
 def __init__(self,model,matched):self.model=model;self.matched=matched;self.L=model.L
 def A(self,q):
  q=CD(q)
  if self.matched:return self.model.A(q)
  return self.model.raw_uv(q) if abs(q)>=self.model.switch else m.base.A(q,self.L)
 def derivative(self,q):
  v=self.model.derivative(q)
  return v if self.matched else v-self.model.d1-2*self.model.d2*q
 def lower(self,s):
  s=LD(s)
  if self.matched:return self.model.lower(s)
  if s>=self.model.switch:return self.model.raw_uv(CD(-s),logq=np.asarray(np.log(s)-1j*PI,dtype=CD))
  return m.base.boundary(s,self.L)
 def rho(self,s):return (1/self.lower(s)).imag/PI
 def poles(self):
  old=json.loads((ROOT/'experiments/matched_inverse_reconstruction_2026_10_30/results.json').read_text())
  if self.matched:
   seeds=[CD(complex(*p['q'])) for row in old['pole_residues'] if row['L']==self.L for p in row['isolated']]
   slope=self.model.ref[0]
  else:
   raw=json.loads((ROOT/'experiments/inverse_scalar_dispersion_2026_10_30/results.json').read_text())
   seeds=[CD(p['q']) for row in raw['pole_rows'] if row['L']==self.L for p in row['poles'] if p['kind']!='massless']
   slope=self.model.co[0]
  ps=[dict(q=CD(0),residue=CD(1/slope))]
  for q in seeds:
   best=q;error=abs(self.A(q))
   for _ in range(5):
    q-=self.A(q)/self.derivative(q);e=abs(self.A(q))
    if e<error:best=q;error=e
   assert error<1e-12
   ps.append(dict(q=best,residue=1/self.derivative(best)))
  cert=json.loads((ROOT/'experiments/scalar_cut_boundary_2026_10_30/results.json').read_text())
  expected=next(x['total_A_zeros_counting_multiplicity'] for x in cert['global_counts'] if x['L']==self.L and x['matched_to_L10']==self.matched)
  assert len(ps)==expected
  return ps

def integrate(f,a,b,label,complex_value=False):
 with warnings.catch_warnings(record=True) as caught:
  warnings.simplefilter('always')
  re,e1=quad(lambda t:float(np.real(f(t))),a,b,epsabs=1e-10,epsrel=1e-11,limit=500)
  if complex_value:im,e2=quad(lambda t:float(np.imag(f(t))),a,b,epsabs=1e-10,epsrel=1e-11,limit=500)
  else:im=e2=0.
 for w in caught:WARNINGS.append(dict(label=label,message=str(w.message)))
 return CD(re+1j*im),e1+e2

def cut_integral(response,R,kernel,label):
 edges=sorted(set([float(x) for x in 4*response.model.z]+[response.model.switch,float(R)]))
 edges=[s for s in edges if s<=R];total=CD(0);error=0.
 for lo,hi in zip(edges[:-1],edges[1:]):
  def fun(t):
   if t==0 or t==1:return CD(0)
   ss=lo+(hi-lo)*np.sin(PI*LD(t)/2)**2;jac=(hi-lo)*PI*np.sin(PI*LD(t))/2
   return response.rho(ss)*kernel(ss)*jac
  ans,err=integrate(fun,0,1,label+f' [{lo},{hi}]',True);total+=ans;error+=err
 return total,error

def tail_integral(response,R,kernel,label):
 # x=R/s; x=0 is not sampled by quad. The transformed kernels have finite limits.
 def fun(x):
  if x==0:return CD(0) # single endpoint value does not alter the integral
  s=LD(R)/LD(x);return response.rho(s)*kernel(s)*R/LD(x)**2
 return integrate(fun,0,1,label,True)

def pair(z):return [float(z.real),float(z.imag)]
def normerr(a,b,scale=None):return float(abs(a-b)/(max(abs(a),LD(1)) if scale is None else scale))
def reconstruct(response,ps,certificate,q):
 q=CD(q);q0=CD(1);R=certificate['R'];K=F(certificate['rho_K_upper']);p=certificate['rho_decay_power']
 kernel=lambda s:1/(s+q)-1/(s+q0)
 cut,ce=cut_integral(response,R,kernel,f'{response.L} {response.matched} subtracted {q}')
 tail,te=tail_integral(response,R,kernel,f'{response.L} {response.matched} subtracted tail {q}')
 pole=sum(x['residue']*(1/(q-x['q'])-1/(q0-x['q'])) for x in ps)
 exact=1/response.A(q)-1/response.A(q0)
 h=max(F(1),F(str(float(abs(q.real)+abs(q.imag)))))
 dq=F(str(float(abs((q-q0).real)+abs((q-q0).imag))))
 bound=K*dq/((p+1)*(1-h/R)**2*R**(p+1))
 assert abs(tail)<float(bound)
 reconstructed=pole+cut+tail
 return dict(L=response.L,matched_to_L10=response.matched,q=pair(q),q0=pair(q0),R=R,exact_difference=pair(exact),pole_sum=pair(pole),finite_cut=pair(cut),numerical_UV_tail=pair(tail),UV_tail_abs_bound=str(bound),UV_tail_abs_bound_decimal=float(bound),quadrature_error_estimate=ce+te,relative_closure_error=normerr(exact,reconstructed),finite_cut_only_relative_error=normerr(exact,pole+cut))

def unsubtracted(response,ps,certificate,q):
 q=CD(q);R=certificate['R'];K=F(certificate['rho_K_upper']);p=certificate['rho_decay_power']
 kernel=lambda ss:1/(ss+q)
 cut,ce=cut_integral(response,R,kernel,f'{response.L} {response.matched} absolute {q}')
 tail,te=tail_integral(response,R,kernel,f'{response.L} {response.matched} absolute tail {q}')
 pole=sum(x['residue']/(q-x['q']) for x in ps)
 constant=CD(0) if p==5 else 1/(C*response.model.constant)
 exact=1/response.A(q);value=constant+pole+cut+tail
 h=F(str(float(abs(q.real)+abs(q.imag))))
 bound=K/(p*(1-h/R)*R**p);assert abs(tail)<float(bound)
 return dict(L=response.L,matched_to_L10=response.matched,q=pair(q),R=R,constant_F_infinity=pair(constant),direct_F=pair(exact),pole_sum=pair(pole),finite_cut=pair(cut),numerical_UV_tail=pair(tail),UV_tail_abs_bound=str(bound),UV_tail_abs_bound_decimal=float(bound),quadrature_error_estimate=ce+te,relative_closure_error=normerr(exact,value))

def moments(response,ps,certificate):
 assert response.matched and response.L!=10
 R=certificate['R'];K=F(certificate['rho_K_upper']);a1=response.model.d1;a2=response.model.d2;ainf=C*response.model.constant
 targets=[LD(0),1/a2,-a1/a2**2,a1*a1/a2**3-ainf/a2**2];out=[]
 for n in range(4):
  cut,ce=cut_integral(response,R,lambda s:s**n,f'{response.L} moment {n}')
  tail,te=tail_integral(response,R,lambda s:s**n,f'{response.L} moment tail {n}')
  pole=sum(x['residue']*x['q']**n for x in ps)
  value=pole+(-1)**n*(cut+tail);scale=max(LD(1),abs(pole)+abs(cut)+abs(tail)+abs(targets[n]));bound=K/F((4-n)*R**(4-n))
  assert abs(tail)<float(bound)
  out.append(dict(L=response.L,n=n,pole_moment=pair(pole),finite_cut_moment=pair(cut),numerical_tail_moment=pair(tail),exact_UV_tail_bound=str(bound),UV_tail_bound_decimal=float(bound),expected_coefficient=float(targets[n]),reconstructed=pair(value),quadrature_error_estimate=ce+te,absolute_closure_error=float(abs(value-targets[n])),relative_to_component_scale=normerr(value,targets[n],scale),partial_sum_without_tail=pair(pole+(-1)**n*cut)))
 return out

def main():
 certificates=[];recon=[];absolute=[];moment_rows=[];asym=[];poles=[]
 for L in [3,5,10,30,100]:
  model=m.Model(L);R=256*(5*L+1);v=uv.exterior(L,'A',F(R))
  for matched in [False,True]:
   response=Response(model,matched);b=bounds(L,matched,v);certificates.append(b);ps=response.poles()
   poles.append(dict(L=L,matched_to_L10=matched,poles=[dict(q=pair(x['q']),residue=pair(x['residue'])) for x in ps]))
   for q in [.5,2.,1+2j]:
    recon.append(reconstruct(response,ps,b,q))
    if not (matched and L==10):absolute.append(unsubtracted(response,ps,b,q))
   if matched and L!=10:moment_rows+=moments(response,ps,b)
   p=b['rho_decay_power'];z,c=uv.m.spec(L);S3=sum(ci*zi**3 for zi,ci in zip(z,c));den=model.d2 if p==5 else C*model.constant
   alpha=-5*C*m.asld(S3)/(4*den**2)
   for multiple in [1,100,10000]:
    s=LD(R)*multiple;ratio=response.rho(s)*s**p/alpha
    asym.append(dict(L=L,matched_to_L10=matched,s=float(s),decay_power=p,leading_positive_coefficient=float(alpha),ratio_to_leading=float(ratio)))
 result=dict(scope='Specified raw and Lref10 matched auxiliary scalar coefficient inverse; no source amplitude or physical positivity claim.',baseline='ab6f85ff8e1d150e39c8f10782615df8b140213f',versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,longdouble_bits=np.finfo(LD).nmant),exact_tail_certificates=certificates,poles=poles,once_subtracted_reconstructions=recon,absolute_reconstructions=absolute,matched_inverse_moment_sum_rules=moment_rows,UV_density_asymptotic_checks=asym,quadrature_warnings=WARNINGS)
 (HERE/'results.json').write_text(json.dumps(result,indent=2))
 summary=dict(reconstructions=len(recon),absolute_checks=len(absolute),max_absolute_reconstruction_error=max(x['relative_closure_error'] for x in absolute),moment_checks=len(moment_rows),max_relative_reconstruction_error=max(x['relative_closure_error'] for x in recon),max_moment_component_relative_error=max(x['relative_to_component_scale'] for x in moment_rows),warnings=len(WARNINGS))
 print(json.dumps(summary,indent=2))
 assert summary['max_absolute_reconstruction_error']<1e-8 and summary['max_relative_reconstruction_error']<1e-8 and summary['max_moment_component_relative_error']<1e-8
if __name__=='__main__':main()
