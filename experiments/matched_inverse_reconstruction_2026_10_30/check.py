"""Matched inverse response: complete bounded conformal Cauchy and slit diagnostics.
Run python experiments/matched_inverse_reconstruction_2026_10_30/check.py.
No global pole completeness or physical-health inference from closure residuals.
"""
from pathlib import Path
import importlib.util,json,platform
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import roots_legendre
ROOT=Path(__file__).resolve().parents[2]
def load(path,name):
 s=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
base=load('experiments/inverse_scalar_dispersion_2026_10_30/check.py','scalar')
th=load('experiments/threshold_independent_2026_10_30/check.py','threshold')
uv=load('experiments/uv_asymptotic_audit_2026_10_30/check.py','uv')
LD=np.longdouble;CD=np.clongdouble;PI=np.arccos(LD(-1));C=base.C
PROBES=[.5,2.,1+2j]
def asld(f):return LD(f.numerator)/LD(f.denominator)
def wmap(q):
 b=np.sqrt(1+CD(q)/4);return (b-1)/(b+1)
def qmap(w):return 16*w/(1-w)**2
def qprime(w):return 16*(1+w)/(1-w)**3

class Model:
 def __init__(self,L):
  self.L=L;self.z,self.c=base.spectrum(L);self.ref=th.coeff(10);self.co=th.coeff(L)
  self.d1=self.ref[0]-self.co[0];self.d2=self.ref[1]-self.co[1];self.switch=16*float(self.z[-1])
  zz,cc,norm,const,rows=uv.terms(L,'A')
  self.constant=asld((const[0]+const[1])/2);self.norm=asld(norm)
  self.k=np.array([r[0] for r in rows],dtype=int)
  self.b=np.array([asld(r[1]) for r in rows]);self.a=np.array([asld(r[2]) for r in rows]);self.s=np.array([asld(r[3]) for r in rows]);self.t=np.array([asld((r[4][0]+r[4][1])/2) for r in rows])
 def raw_uv(self,q,logq=None):
  q=np.asarray(q,dtype=CD);logq=np.log(q) if logq is None else logq
  vals=q[...,None]**(2-self.k)*(self.b*(self.s*logq[...,None]-self.t)+self.a*self.s)
  return C*(self.constant+self.norm*np.sum(vals,axis=-1))
 def A(self,q):
  q=np.asarray(q,dtype=CD)
  if np.all(abs(q)>=self.switch):raw=self.raw_uv(q)
  elif np.all(abs(q)<self.switch):raw=base.A(q,self.L)
  else:
   raw=np.empty(q.shape,dtype=CD);sel=abs(q)>=self.switch
   raw[sel]=self.raw_uv(q[sel]);raw[~sel]=base.A(q[~sel],self.L)
  return raw+self.d1*q+self.d2*q*q
 def derivative(self,q):
  q=CD(q);z,c=self.z,self.c;r=np.sqrt(q/(q+4*z));J=np.log(z)-2+2*np.arctanh(r)/r
  Jp=(1-4*z*r*np.arctanh(r)/q)/q;U=z-q/2
  return C/4*np.dot(c,-U*J+U*U*Jp+z*np.log(z)/2)+self.d1+2*self.d2*q
 def derivative_quad(self,q,n=512):
  z,c=self.z,self.c;U=z-q/2
  if q.real<-4 and abs(q.imag)<.1*abs(q.real):
   J=[];Jp=[]
   for mass in z:
    pts=[]
    if q.real<-4*mass:
     beta=np.sqrt(1+4*float(mass)/q.real);roots=[(1-beta)/2,(1+beta)/2];width=abs(q.imag)*float(mass)/(q.real*q.real*beta)
     pts=sorted(set(p for root in roots for k in [-1000,-100,-10,-1,0,1,10,100,1000] if 0<(p:=root+k*width)<1))
    def integrate(f):
     re=quad(lambda x:float(f(x).real),0,1,points=list(map(float,pts)),epsabs=1e-13,epsrel=1e-12,limit=500)[0]
     im=quad(lambda x:float(f(x).imag),0,1,points=list(map(float,pts)),epsabs=1e-13,epsrel=1e-12,limit=500)[0]
     return complex(re,im)
    J.append(integrate(lambda x:np.log(complex(mass)+q*x*(1-x))))
    Jp.append(integrate(lambda x:x*(1-x)/(complex(mass)+q*x*(1-x))))
   J=np.array(J,dtype=CD);Jp=np.array(Jp,dtype=CD)
  else:
   x,w=roots_legendre(n);x=(x+1)/2;w=w/2;t=x*(1-x)
   J=np.log(z[:,None]+q*t)@w;Jp=(t/(z[:,None]+q*t))@w
  return C/4*np.dot(c,-U*J+U*U*Jp+z*np.log(z)/2)+self.d1+2*self.d2*q
 def lower(self,s):
  if s>=self.switch:raw=self.raw_uv(CD(-s),logq=np.asarray(np.log(LD(s))-1j*PI,dtype=CD))
  else:raw=base.boundary(s,self.L)
  return raw-self.d1*s+self.d2*s*s
 def rho(self,s):return float((1/self.lower(s)).imag/PI)

def pole_list(model):
 L=model.L;old=json.loads((ROOT/'experiments/threshold_independent_2026_10_30/results.json').read_text())
 new=json.loads((ROOT/'experiments/scalar_conformal_count_2026_10_30/verification_results.json').read_text())
 seeds=[(CD(-next(r['s'] for r in old['matched_comparison'] if r['L']==L)),'subthreshold')]
 for r in new['extra_matched_scalar_roots']:
  if r['L']==L:
   q=CD(complex(*r['q']));seeds.append((q,'extra'))
   if q.imag:seeds.append((q.conjugate(),'extra'))
 if L==30:
  certificate=json.loads(Path(__file__).with_name('local_pole_certificate.json').read_text());q=CD(complex(*certificate['center_float']))
  seeds.extend([(q,'collar'),(q.conjugate(),'collar')])
 poles=[dict(q=CD(0),residue=CD(1/model.ref[0]),kind='massless')];diagnostics=[]
 for seed,kind in seeds:
  q=seed;best=q;res=abs(model.A(q))
  for _ in range(6):
   q=q-model.A(q)/model.derivative(q);now=abs(model.A(q))
   if now<res:best=q;res=now
  q=best;der=model.derivative(q);direct=model.derivative_quad(q);error=float(abs(direct/der-1))
  assert error<1e-8 and res<1e-13
  poles.append(dict(q=q,residue=1/der,kind=kind))
  diagnostics.append(dict(kind=kind,q=[float(q.real),float(q.imag)],residual=float(res),residue=[float((1/der).real),float((1/der).imag)],independent_derivative_relative_error=error))
 return poles,diagnostics

def complex_pair(z):return [float(np.real(z)),float(np.imag(z))]
def relative_error(a,b):return float(abs(a-b)/max(abs(a),LD(1e-12)))
def conformal_reconstruct(model,poles,q):
 poles=[p for p in poles if abs(wmap(p['q']))<.8]
 q0=CD(1);w=wmap(q);w0=wmap(q0);exact=1/model.A(q)-1/model.A(q0)
 pole=sum(p['residue']/qprime(wmap(p['q']))*(1/(w-wmap(p['q']))-1/(w0-wmap(p['q']))) for p in poles)
 outer=[]
 for n in [128,256,512]:
  x,weights=roots_legendre(n);theta=x.astype(LD)*PI;zz=LD(.8)*np.exp(1j*theta)
  val=np.sum(weights*(1/model.A(qmap(zz)))*(zz/(zz-w)-zz/(zz-w0)))/2;outer.append(val)
 return dict(L=model.L,q=complex_pair(q),exact_difference=complex_pair(exact),pole_difference=complex_pair(pole),outer_boundaries=[dict(nodes=n,value=complex_pair(v)) for n,v in zip([128,256,512],outer)],relative_closure_error=relative_error(exact,pole+outer[-1]),omitting_extra_poles_error=relative_error(exact,outer[-1]+sum(p['residue']/qprime(wmap(p['q']))*(1/(w-wmap(p['q']))-1/(w0-wmap(p['q']))) for p in poles[:2])))

def slit_reconstruct(model,poles,q,R):
 q=CD(q);q0=CD(1);exact=1/model.A(q)-1/model.A(q0)
 pole=sum(p['residue']*(1/(q-p['q'])-1/(q0-p['q'])) for p in poles)
 edges=sorted(set([float(x) for x in 4*model.z]+[model.switch,R]));edges=[x for x in edges if x<=R]
 cut=0j;error=0.
 for lo,hi in zip(edges[:-1],edges[1:]):
  def integrand(t):
   if t==0 or t==1:return 0j
   ss=lo+(hi-lo)*np.sin(np.pi*t/2)**2;jac=(hi-lo)*np.pi*np.sin(np.pi*t)/2
   return complex(model.rho(ss)*(1/(ss+q)-1/(ss+q0))*jac)
  re,e1=quad(lambda t:integrand(t).real,0,1,epsabs=1e-10,epsrel=2e-10,limit=400)
  im,e2=quad(lambda t:integrand(t).imag,0,1,epsabs=1e-10,epsrel=2e-10,limit=400)
  cut+=complex(re,im);error+=e1+e2
 outer=[]
 for n in [64,128,256]:
  x,weights=roots_legendre(n);theta=x.astype(LD)*PI;zz=R*np.exp(1j*theta)
  outer.append(np.sum(weights*(1/model.A(zz))*(zz/(zz-q)-zz/(zz-q0)))/2)
 return dict(L=model.L,q=complex_pair(q),R=R,exact_difference=complex_pair(exact),pole_difference=complex_pair(pole),cut_difference=complex_pair(cut),cut_quadrature_error_estimate=error,outer_boundaries=[dict(nodes=n,value=complex_pair(v)) for n,v in zip([64,128,256],outer)],relative_closure_error=relative_error(exact,pole+cut+outer[-1]),relative_error_without_outer=relative_error(exact,pole+cut),relative_error_omitting_collar_pair=relative_error(exact,pole+cut+outer[-1]-sum(p['residue']*(1/(q-p['q'])-1/(q0-p['q'])) for p in poles if p['kind']=='collar')),interpretation='Diagnostic closure with known poles; no proof that the uncounted collar contains no further poles.')

def main():
 mapped=[];slit=[];pole_data=[];wrapper=[]
 saved=json.loads((ROOT/'experiments/uv_asymptotic_audit_2026_10_30/results.json').read_text())
 for L in [3,5,10,30,100]:
  model=Model(L);poles,diagnostics=pole_list(model);pole_data.append(dict(L=L,common_massless_residue=float(poles[0]['residue'].real),isolated=diagnostics))
  for r in saved['decimal_160_digit_comparisons']:
   if r['L']==L and r['channel']=='A':
    R=r['R'];val=model.raw_uv(-CD(R),logq=np.asarray(np.log(LD(R))-1j*PI,dtype=CD)) if r['lower_cut'] else model.raw_uv(CD(R))
    ref=CD(LD(r['direct'][0])+1j*LD(r['direct'][1]))*C;err=relative_error(ref,val);assert err<1e-15
    wrapper.append(dict(L=L,R=R,lower_cut=r['lower_cut'],relative_error=err))
  for q in PROBES:
   mapped.append(conformal_reconstruct(model,poles,q))
   for f in [16,64]:slit.append(slit_reconstruct(model,poles,q,f*float(model.z[-1])))
 result=dict(scope='Explicitly Lref=10 matched inverse scalar response. Conformal bounded inventory is complete from prior validated count; slit reconstruction outside that region remains diagnostic.',versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,longdouble_bits=np.finfo(LD).nmant),pole_residues=pole_data,stable_uv_wrapper_checks=wrapper,conformal_cauchy=mapped,slit_cut_diagnostics=slit)
 out=Path(__file__).resolve().parent/'results.json';out.write_text(json.dumps(result,indent=2))
 print(json.dumps(dict(mapped_checks=len(mapped),slit_checks=len(slit),max_mapped_closure=max(x['relative_closure_error'] for x in mapped),max_slit_closure=max(x['relative_closure_error'] for x in slit),max_error_omitting_new_poles=max(x['omitting_extra_poles_error'] for x in mapped)),indent=2))
 assert max(x['relative_closure_error'] for x in mapped)<1e-8
 assert max(x['relative_closure_error'] for x in slit)<1e-8
if __name__=='__main__':main()
