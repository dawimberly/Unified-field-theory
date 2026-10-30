"""Analytic absolute tail bounds inside |q|<4 for explicitly matched kernels.
Run python experiments/low_energy_remainders_2026_10_30/check.py.
Bounds are exact-arithmetic inequalities; displayed evaluations use floating point.
"""
from pathlib import Path
import importlib.util,json,platform
import numpy as np
import scipy
p=Path(__file__).resolve().parents[1]/'scalar_identity_matching_2026_10_30/check.py'
s=importlib.util.spec_from_file_location('series_reference',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def bound(r,L,ch,N):
 z,c=m.scalar.spectrum(L);f=abs(m.rational_coeff(N+1,ch))
 return m.C*np.longdouble(f.numerator)/f.denominator*r**(N+1)*np.sum(abs(c)*z**(1-N)/(1-r/(4*z)))
def main():
 rows=[];inverse=[];radii=[]
 for ch in ['A','K2']:
  ref=m.coefficients(10,ch);kernel=m.scalar.A if ch=='A' else m.spin.kernel
  for L in [3,5,10,30,100]:
   co=m.coefficients(L,ch);matched={**co,1:ref[1],2:ref[2]}
   for r in [.25,.5,1.,2.]:
    for N in [2,4]:
     B=bound(r,L,ch,N)
     for phase in [0.,np.pi/3,2*np.pi/3,np.pi]:
      q=r*np.exp(1j*phase)
      exact=kernel(q,L)+(ref[1]-co[1])*q+(ref[2]-co[2])*q*q
      trunc=sum(matched[n]*q**n for n in range(1,N+1))
      error=abs(exact-trunc)
      rows.append(dict(channel=ch,L=L,N=N,q=[q.real,q.imag],absolute_error=float(error),absolute_bound=float(B),error_bound_ratio=float(error/B)))
      if abs(trunc)>B:
       rel=float(B/abs(trunc));observed=float(abs(1-exact/trunc))
       inverse.append(dict(channel=ch,L=L,N=N,q=[q.real,q.imag],relative_inverse_error=observed,relative_inverse_bound=rel,absolute_inverse_bound=float(B/(abs(trunc)*(abs(trunc)-B)))))
   # Uniform disk bound on factored kernel K/q, excluding only massless q=0.
   for r in [.25,.5,1.,2.]:
    # Quadratic factor ref1+ref2 q has minimum modulus >= |ref1|-|ref2|r.
    lower=abs(ref[1])-abs(ref[2])*r
    eps=bound(r,L,ch,2)/r
    radii.append(dict(channel=ch,L=L,r=r,factored_quadratic_lower_bound=float(lower),factored_tail_upper_bound=float(eps),no_extra_zero_certified_by_formula=bool(lower>eps),uniform_relative_inverse_bound=float(eps/lower) if lower>eps else None))
 result=dict(scope='Explicit Lref=10 matching; absolute geometric tail inequality for |q|<4, not a physical UV error estimate. Floating-point evaluation is not interval certification.',versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),kernel_checks=rows,inverse_checks=inverse,disk_bounds=radii)
 (Path(__file__).resolve().parent/'results.json').write_text(json.dumps(result,indent=2))
 assert all(x['absolute_error']<=x['absolute_bound'] for x in rows)
 assert all(x['relative_inverse_error']<=x['relative_inverse_bound'] for x in inverse)
 print(json.dumps(dict(kernel_checks=len(rows),inverse_checks=len(inverse),max_error_bound_ratio=max(x['error_bound_ratio'] for x in rows),disk_bounds=len(radii),certified_by_formula=sum(x['no_extra_zero_certified_by_formula'] for x in radii)),indent=2))
if __name__=='__main__':main()
