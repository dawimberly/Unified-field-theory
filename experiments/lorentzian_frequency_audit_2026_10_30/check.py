"""Conditional frequency continuation of saved Euclidean poles.
Not a microscopic Lorentzian contour choice or an EFT stability verdict.
"""
from pathlib import Path
import json,cmath,math,platform
from fractions import Fraction as F
import numpy as np
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent

def pair(z):return [float(z.real),float(z.imag)]
def read(path):return json.loads((ROOT/path).read_text())

def inventories():
 scalar=read('experiments/infinite_scalar_dispersion_2026_10_30/results.json')
 spin=read('experiments/spin2_global_count_2026_10_30/results.json')
 rows=[]
 for r in scalar['poles']:
  for p in r['poles']:
   q=complex(*p['q'])
   if q==0:continue
   rows.append(dict(L=r['L'],matched=r['matched_to_L10'],channel='scalar',q=q,residue=complex(*p['residue']),source_factor=(q-2)**2/16))
 # Exact certificates supply spin-2 complex centers and signed real brackets.
 # Residues are unnecessary for mapping existence, so none is invented here.
 for r in spin['positive_real_root_brackets']:
  a,b=map(F,r['bracket_exact']);q=complex(float((a+b)/2))
  rows.append(dict(L=r['L'],matched=True,channel='spin2',q=q,localization_radius=float((b-a)/2)))
 for r in spin['complex_local_certificates']:
  q=complex(*r['center_decimal']);rad=float(F(r['radius']))
  for qq in [q,q.conjugate()]:rows.append(dict(L=r['L'],matched=True,channel='spin2',q=qq,localization_radius=rad))
 return rows

def frequency_roots(q,k):return [cmath.sqrt(k*k-q),-cmath.sqrt(k*k-q)]
def lower_root(q,k):
 roots=frequency_roots(q,k);w=next(x for x in roots if x.imag<0)
 return w

def classify(q,k):
 if q.imag:return 'complex_pair_member_one_UHP'
 if q.real>k*k:return 'imaginary_pair_one_UHP'
 if q.real==k*k:return 'frequency_zero_double_root'
 return 'real_boundary_pair'

def causal(pole,t,k=0):
 w=cmath.sqrt(k*k-pole['q'])
 return pole['residue']*cmath.sin(w*t)/w if t>0 else 0j

def real_axis(pole,t,k=0):
 w=lower_root(pole['q'],k)
 return 1j*pole['residue']*cmath.exp(-1j*w*abs(t))/(2*w)

def main():
 ps=inventories();mapping=[];counts=[];residue_checks=[]
 for p in ps:
  for k in [0.,1.,10.,100.]:
   roots=frequency_roots(p['q'],k)
   errors=[abs(k*k-w*w-p['q'])/max(1,abs(p['q']),k*k) for w in roots]
   assert max(errors)<1e-12
   upper=[w for w in roots if w.imag>0]
   mapping.append(dict(L=p['L'],matched=p['matched'],channel=p['channel'],q=pair(p['q']),k=k,frequency_roots=[pair(w) for w in roots],upper_half_plane_count=len(upper),upper_growth_rates=[w.imag for w in upper],classification=classify(p['q'],k),q_mapping_error=max(errors),source_factor=pair(p['source_factor']) if 'source_factor' in p else None))
   if 'residue' in p:
    # Independent symmetric near-pole frequency limit checks the chain-rule sign.
    for w in roots:
     if abs(w.imag)<1e-12:continue
     expected=-p['residue']/(2*w);h=1e-5
     value=sum(hh*p['residue']/(k*k-(w+hh)**2-p['q']) for hh in [h,-h])/2
     err=abs(value-expected)/max(1,abs(expected));assert err<1e-7
     residue_checks.append(dict(L=p['L'],matched=p['matched'],k=k,q=pair(p['q']),omega=pair(w),expected_frequency_residue=pair(expected),symmetric_limit=pair(value),normalized_error=err))
 for matched in [False,True]:
  for L in [3,5,10,30,100]:
   for k in [0.,1.,10.,100.]:
    r=[x for x in mapping if x['L']==L and x['matched']==matched and x['k']==k]
    counts.append(dict(L=L,matched=matched,k=k,scalar_UHP=sum(x['upper_half_plane_count'] for x in r if x['channel']=='scalar'),spin2_UHP=sum(x['upper_half_plane_count'] for x in r if x['channel']=='spin2')))
 assert [x['scalar_UHP']+x['spin2_UHP'] for x in counts if x['matched'] and x['k']==0]==[3,3,0,5,3]
 assert all(x['scalar_UHP']==x['spin2_UHP']==0 for x in counts if not x['matched'])
 # New contour diagnostics for saved off-real scalar poles only. A pole and its
 # conjugate can be combined afterwards; the inversion is linear.
 chosen=[]
 for L in [3,30]:
  p=next(x for x in ps if x['L']==L and x['matched'] and x['q'].imag>0)
  chosen.append(p)
 p=next(x for x in ps if x['L']==100 and x['matched'] and x['q'].imag==0 and x['q'].real>0);chosen.append(p)
 inversion=[];time_rows=[]
 for p in chosen:
  for t in [-1.,-.25,.25,1.]:
   # On the real frequency axis the rational integrand has no pole. Its evenness
   # permits an independent cosine quadrature, including t<0.
   def f(w):return p['residue']/(0-w*w-p['q'])/math.pi
   re,e1=quad(lambda w:f(w).real,0,np.inf,weight='cos',wvar=abs(t),epsabs=2e-10,limit=400,limlst=400)
   im,e2=quad(lambda w:f(w).imag,0,np.inf,weight='cos',wvar=abs(t),epsabs=2e-10,limit=400,limlst=400)
   numerical=complex(re,im);closed=real_axis(p,t);err=abs(numerical-closed)/max(1,abs(closed));assert err<1e-8
   inversion.append(dict(L=p['L'],q=pair(p['q']),t=t,real_axis_closed_form=pair(closed),real_axis_numerical_inverse=pair(numerical),quadrature_error_estimate=e1+e2,normalized_error=err))
  for t in [-2.,-1.,0.,.25,.5,1.,2.,4.]:
   time_rows.append(dict(L=p['L'],q=pair(p['q']),t=t,causal_growing=pair(causal(p,t)),real_axis_two_sided=pair(real_axis(p,t))))
  # Fundamental-solution jump: g(0+)=0 and g'(0+)=r_q.
  h=1e-6;der=causal(p,h)/h;err=abs(der-p['residue'])/max(1,abs(p['residue']));assert err<1e-8
 result=dict(scope='Formal continuation q=k^2-omega^2 of fixed Euclidean coefficient inverses. Upper-pole locations and toy pole-sector inversions do not choose a microscopic contour, retarded quantum theory or physical completion.',baseline='945106beb3436dca4eb039926eb9686b635b937e',versions=dict(python=platform.python_version(),numpy=np.__version__),frequency_maps=mapping,UHP_counts=counts,frequency_residue_checks=residue_checks,off_real_pole_sector_inverse_transforms=inversion,time_domain_examples=time_rows,interpretation='Causal inverse-Laplace response can grow exponentially. Real-axis inversion of off-real pole sectors decays both ways but is generally nonzero before the source. No full-kernel time evolution or nonlinear stability was computed.')
 (HERE/'results.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(dict(mapped_entries=len(mapping),residue_checks=len(residue_checks),inverse_transforms=len(inversion),max_inverse_error=max(x['normalized_error'] for x in inversion),max_residue_error=max(x['normalized_error'] for x in residue_checks),matched_k0_UHP_counts=[x['scalar_UHP']+x['spin2_UHP'] for x in counts if x['matched'] and x['k']==0]),indent=2))
if __name__=='__main__':main()
