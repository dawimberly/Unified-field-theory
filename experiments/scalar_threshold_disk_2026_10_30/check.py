"""Positive scalar coefficient audit and subthreshold-disk zero classification.
Run python experiments/scalar_threshold_disk_2026_10_30/check.py.
Exact rational inequalities are separated from floating-point endpoint/bracket checks.
"""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,platform,math
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
def load(path,name):
 s=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
base=load('experiments/threshold_independent_2026_10_30/check.py','threshold_reference')
LD=np.longdouble

def spec(L):
 z=list(map(F,[1,L,2*L,5*L+1]));w=[1/math.prod(z[i]-z[j] for j in range(4) if i!=j) for i in range(4)]
 return z,[x/w[0] for x in w]
def homogeneous(x,N):
 h=[F(1)]+[F(0)]*N
 for a in x:
  for k in range(1,N+1):h[k]+=a*h[k-1]
 return h
def normalized_resolvent(u,z):return math.prod(x/(x+u) for x in z)
def rational_ratio(n):return F((n-2)*(3*n*n+3*n-1),2*(2*n+3)*(3*n*n-3*n-1))
def series(L,matched,N):
 z,c=base.base.spectrum(L);a1,a2=base.coeff(10 if matched else L)
 # F(s)=A(-s)/s = constant - sum d_k s^k.
 d=[LD(0),-a2];mu=LD(1)/6 # mu_(n-2) for n=3
 for n in range(3,N+3):
  r=3*mu*(3*n*n-3*n-1)/(16*(n-2)*(4*n*n-1))
  d.append(base.C*r*np.dot(c,z**(2-n)))
  mu*=LD(n-1)/(2*(2*n-1))
 return -a1,d

def bracket_values(L,matched,s,N):
 c,d=series(L,matched,N)
 powers=LD(s)**np.arange(1,N+1,dtype=LD)
 upper=c-np.dot(np.array(d[1:N+1]),powers)
 tail=d[N+1]*LD(s)**(N+1)/(1-LD(s)/4)
 return dict(lower=float(upper-tail),upper=float(upper),tail_bound=float(tail))

def plus(a,b):return a[0]+b[0],a[1]+b[1]
def times(a,c):return (a[0]*c,a[1]*c) if c>=0 else (a[1]*c,a[0]*c)
def log_series(t,N=64):
 value=2*sum(t**(2*n+1)/F(2*n+1) for n in range(N))
 tail=2*t**(2*N+1)/(F(2*N+1)*(1-t*t))
 return value,value+tail
def log_bounds(x):
 k=0;y=x
 while y>=2:y/=2;k+=1
 assert 1<=y<2
 return plus(times(log_series(F(1,3)),k),log_series((y-1)/(y+1)))
def threshold_interval(L,matched=False):
 z,c=spec(L);logs=[log_bounds(x) for x in z];ans=(F(0),F(0))
 for x,ci,ln in zip(z,c,logs):
  if x==1:S=(F(0),F(0))
  else:
   N=128;partial=sum(F((-1)**n,2*n+1)/(x-1)**n for n in range(N))
   nxt=F((-1)**N,2*N+1)/(x-1)**N
   S=(min(partial,partial+nxt),max(partial,partial+nxt))
  term=plus(times(ln,(x+2)/2),times(plus(S,(F(-1),F(-1))),(x+2)**2/2))
  ans=plus(ans,times(term,ci))
 if matched:
  for LL,sign in [(10,1),(L,-1)]:
   zz,cc=spec(LL)
   for x,ci in zip(zz,cc):ans=plus(ans,times(log_bounds(x),sign*ci*(x/2+1)))
 return ans

def main():
 # Exact right-endpoint Riemann sums on [0,16], 256 intervals of width 1/16.
 cert=[]
 for L,normalized in [(3,True),(10,False)]:
  z,c=spec(L);f0=sum(ci/zi for ci,zi in zip(c,z))
  lower=sum(normalized_resolvent(F(i,16),z)*F(1,16) for i in range(1,257))
  if not normalized:lower*=f0
  margin=lower-F(34,35);assert margin>0
  cert.append(dict(L=L,quantity='integral f/f(0)' if normalized else 'integral f',exact_lower_bound_decimal=float(lower),comparison='>34/35',exact_positive_margin=True,margin_numerator_bits=margin.numerator.bit_length(),margin_denominator_bits=margin.denominator.bit_length()))
 # Exact inverse moments from the complete homogeneous identity.
 moments=[]
 for L in [3,5,10,30,100]:
  z,c=spec(L);f0=sum(ci/zi for ci,zi in zip(c,z));h=homogeneous([1/x for x in z],16)
  for k in range(1,17):
   direct=sum(ci/zi**k for ci,zi in zip(c,z));positive=f0*h[k-1]
   assert direct==positive and positive>0
   moments.append(dict(L=L,m=k,equality_exact=True,positive=True))
 # n=3,4,5 exceptional finite prefix; remaining n>=6 follows the analytic induction.
 upper_vars=[F(1,3),F(1,6),F(1,16)];h=homogeneous(upper_vars,3);prefix=[]
 for n in [3,4,5]:
  b=rational_ratio(n)*(1+h[n-2]);assert b<F(1,4)
  prefix.append(dict(n=n,ratio_upper_bound=str(b),strictly_below_one_quarter=True))
 assert F(4*5*6,2)<=3**4 # k=3 induction base
 # Classify using finite threshold values, not a new broad pole scan.
 old=json.loads((ROOT/'experiments/inverse_scalar_dispersion_2026_10_30/results.json').read_text())
 oldmatch=json.loads((ROOT/'experiments/threshold_independent_2026_10_30/results.json').read_text())
 rows=[];brackets=[];ratios=[]
 for matched in [False,True]:
  for L in [3,5,10,30,100]:
   a4=float(base.uniform(0,L,matched).real);interval=threshold_interval(L,matched)
   assert interval[1]<0 or interval[0]>0
   count=1 if interval[1]<0 else 0
   assert (a4<0)==bool(count)
   rows.append(dict(L=L,matched_to_L10=matched,A_at_threshold=a4,extra_zeros_in_open_disk= count,threshold_over_C_rational_enclosure=[float(x) for x in interval],threshold_over_C_enclosure_width=float(interval[1]-interval[0]),classification_basis='analytic coefficient theorem plus exact rational threshold-sign enclosure'))
   const,d=series(L,matched,512)
   ratios.append(dict(L=L,matched_to_L10=matched,largest_d_ratio_first_513=float(max(d[k+1]/d[k] for k in range(1,513))),all_positive=all(x>0 for x in d[1:])))
   assert all(x>0 for x in d[1:]) and max(d[k+1]/d[k] for k in range(1,513))<.25
   if not count:continue
   if matched:r=next(x['s'] for x in oldmatch['matched_comparison'] if x['L']==L)
   else:r=-next(x for x in old['pole_rows'] if x['L']==L)['poles'][1]['q']
   delta=1e-8
   below=bracket_values(L,matched,r-delta,512);above=bracket_values(L,matched,r+delta,512)
   assert below['lower']>0 and above['upper']<0
   brackets.append(dict(L=L,matched_to_L10=matched,archived_s=r,bracket_half_width=delta,lower_endpoint=below,upper_endpoint=above))
 result=dict(scope='Scalar zero classification only within |q|<4 for the specified raw or explicitly Lref=10 matched PV determinant. No global first-sheet theorem or physical-health claim.',versions=dict(python=platform.python_version(),numpy=np.__version__,longdouble_bits=np.finfo(LD).nmant),exact_first_ratio_certificates=cert,exact_inverse_moment_checks=moments,exact_higher_ratio_prefix=prefix,coefficient_ratio_samples=ratios,disk_classifications=rows,new_positive_series_root_brackets=brackets)
 (Path(__file__).resolve().parent/'results.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(dict(exact_moments=len(moments),riemann_certificates=len(cert),coefficient_prefix=len(prefix),classifications=rows,independent_series_brackets=len(brackets)),indent=2))
if __name__=='__main__':main()
