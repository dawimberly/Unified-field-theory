"""Independent root quadrature, interval primitive audits and focused refinements.
Run after check.py: python experiments/scalar_conformal_count_2026_10_30/verify.py.
"""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,cmath
import numpy as np
from scipy.optimize import root,brentq
from scipy.special import roots_legendre
s=importlib.util.spec_from_file_location('conformal_count',Path(__file__).with_name('check.py'))
core=importlib.util.module_from_spec(s);s.loader.exec_module(core)
t=core.m.base

def exact_inside(v,iv):return F(iv[0])<=v<=F(iv[1])
def direct(q,L,n):
 z,c=t.base.spectrum(L);x,w=roots_legendre(n);x=(x+1)/2;w=w/2
 logs=np.log(z[:,None]+q*x*(1-x))@w;U=z-q/2
 return complex(t.C/4*np.dot(c,U*U*logs-U*z*np.log(z)))
def main():
 primitives=[]
 for a in [F(-123456789,12345),F(-1,7),F(0),F(1,3),F(987654321,23456)]:
  for b in [F(-9,11),F(1,19),F(1234567,76543)]:
   aa=core.rational(a);bb=core.rational(b)
   assert exact_inside(a+b,core.ia(aa,bb)) and exact_inside(a*b,core.im(aa,bb)) and exact_inside(-a,core.neg(aa))
   primitives.append(dict(a=str(a),b=str(b),sum_product_negation_enclosed=True))
 roots=[];isolations=[];ref=t.coeff(10)
 for L in [3,5,30,100]:
  co=t.coeff(L);d1=ref[0]-co[0];d2=ref[1]-co[1]
  fun=lambda q:t.base.A(q,L)+d1*q+d2*q*q
  if L in [3,5]:
   init=[-17,13] if L==3 else [-25,30]
   rr=root(lambda v:[float(fun(complex(*v)).real),float(fun(complex(*v)).imag)],init,tol=1e-11)
   assert rr.success;q=complex(*rr.x)
  else:q=complex(brentq(lambda q:float(fun(q).real),1,100,xtol=1e-12))
  w=(cmath.sqrt(1+q/4)-1)/(cmath.sqrt(1+q/4)+1)
  center=(F(str(w.real)),F(str(w.imag)))
  isolated=core.count(L,True,center=center,path_radius=F(1,200));assert isolated['zeros_excluding_massless']==1
  isolations.append(isolated)
  integrals=[direct(q,L,n)+complex(d1)*q+complex(d2)*q*q for n in [128,256,512]]
  deriv=[complex((fun(q+h)-fun(q-h))/(2*h)) for h in [1e-3,1e-4,1e-5]]
  assert max(abs(x) for x in integrals)<1e-10
  roots.append(dict(L=L,q=[q.real,q.imag],w=[w.real,w.imag],compact_residual=float(abs(fun(q))),independent_parameter_integral_nodes=[128,256,512],independent_parameter_integral_residuals=[abs(x) for x in integrals],derivatives=[[x.real,x.imag] for x in deriv],inverse_A_residue=[(1/deriv[-1]).real,(1/deriv[-1]).imag],abs_q_over_4=abs(q)/4))
 # New accuracy settings justify this focused repetition; no old broad scan is repeated.
 archived=json.loads((core.ROOT/'experiments/inverse_scalar_dispersion_2026_10_30/results.json').read_text())
 matched_archive=json.loads((core.ROOT/'experiments/threshold_independent_2026_10_30/results.json').read_text())
 known=[]
 for flag in [False,True]:
  for L in [3,5,10,30,100]:
   if flag:q=-next(v['s'] for v in matched_archive['matched_comparison'] if v['L']==L)
   else:
    ps=next(v['poles'] for v in archived['pole_rows'] if v['L']==L)
    if len(ps)==1:continue
    q=ps[1]['q']
   w=(cmath.sqrt(1+q/4)-1)/(cmath.sqrt(1+q/4)+1)
   iv=core.count(L,flag,center=(F(str(w.real)),F(0)),path_radius=F(1,200))
   assert iv['zeros_excluding_massless']==1
   known.append(iv)
 core.PREC=80;core.coefficients.cache_clear()
 refinement=[core.count(L,True,N=160) for L in [3,100]]
 assert [x['zeros_excluding_massless'] for x in refinement]==[3,2]
 # Empirical overlap of the independently derived conformal polynomial and compact kernel.
 overlaps=[]
 for L in [3,100]:
  for matched in [False,True]:
   co=core.coefficients(L,128,matched)
   for w in [.2+.3j,-.3+.4j,.5-.2j]:
    iv=core.evaluate(co,(core.rational(F(str(w.real))),core.rational(F(str(w.imag)))))
    poly=complex(*map(float,core.midpoint(iv)));q=16*w/(1-w)**2
    v=t.base.A(q,L)
    if matched:
     c=t.coeff(L);v+=(ref[0]-c[0])*q+(ref[1]-c[1])*q*q
    exact=complex(4*(1-w)**4*v/(t.C*w))
    rel=abs(poly-exact)/max(abs(exact),1)
    assert rel<1e-10
    overlaps.append(dict(L=L,matched_to_L10=matched,w=[w.real,w.imag],relative_overlap_error=rel))
 result=dict(interval_primitive_checks=primitives,extra_matched_scalar_roots=roots,validated_isolating_counts=isolations,known_subthreshold_isolating_counts=known,higher_order_precision_refinement=refinement,empirical_conformal_overlap=overlaps,limitations='No claim of physical UV completion, EFT prediction, constrained-CEV equivalence or global count as |w| tends to one.')
 Path(__file__).with_name('verification_results.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(dict(roots=roots,isolating_counts=[x['zeros_excluding_massless'] for x in isolations],refined_counts=[x['zeros_excluding_massless'] for x in refinement],max_overlap=max(x['relative_overlap_error'] for x in overlaps)),indent=2))
if __name__=='__main__':main()
