"""Exact cut-boundary/Sturm audit and global argument-principle scalar count.
Run python experiments/scalar_cut_boundary_2026_10_30/check.py.
"""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,math,platform
ROOT=Path(__file__).resolve().parents[2]
def load(path,name):
 s=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
iv=load('experiments/matched_inverse_reconstruction_2026_10_30/local_pole_certificate.py','intervals')
def trim(p):
 while len(p)>1 and p[-1]==0:p.pop()
 return p
def add(p,q):return trim([(p[i] if i<len(p) else 0)+(q[i] if i<len(q) else 0) for i in range(max(len(p),len(q)))])
def scale(p,a):return trim([a*x for x in p])
def mul(p,q):
 r=[F(0)]*(len(p)+len(q)-1)
 for i,a in enumerate(p):
  for j,b in enumerate(q):r[i+j]+=a*b
 return trim(r)
def sub(p,q):return add(p,scale(q,-1))
def sq(p):return mul(p,p)
def val(p,x):
 r=F(0)
 for a in reversed(p):r=r*x+a
 return r
def primitive(p):
 den=math.lcm(*(a.denominator for a in p));r=[int(a*den) for a in p];g=math.gcd(*r)
 return [F(a//g) for a in r] # positive rescaling preserves Sturm signs

def remainder(p,q):
 r=p[:]
 while r!=[0] and len(r)>=len(q):
  k=len(r)-len(q);a=r[-1]/q[-1]
  for j,b in enumerate(q):r[k+j]-=a*b
  trim(r)
 return r
def sturm(p):
 seq=[primitive(p),primitive([i*p[i] for i in range(1,len(p))])]
 while True:
  r=remainder(seq[-2],seq[-1])
  if r==[0]:break
  seq.append(primitive(scale(r,-1)))
 return seq

def variations(seq,x):
 signs=[1 if y>0 else -1 for p in seq if (y:=val(p,x))!=0]
 return sum(a!=b for a,b in zip(signs,signs[1:]))
def count(seq,a,b):
 assert val(seq[0],a)!=0 and val(seq[0],b)!=0
 return variations(seq,a)-variations(seq,b)
def isolate(seq,a,b):
 n=count(seq,a,b)
 if not n:return []
 if n==1 and b-a<F(1,10**14):return [(a,b)]
 mid=(a+b)/2
 assert val(seq[0],mid)!=0,'Exact dyadic root: implement explicit endpoint split before proceeding.'
 return isolate(seq,a,mid)+isolate(seq,mid,b)
def polynomial(L,n):
 z,c=iv.m.spec(L)
 X=[scale(mul(sq(sq([2*z[i],F(1)])),[-4*z[i],F(1)]),c[i]**2) for i in range(n)]
 if n==2:return sub(X[0],X[1])
 assert n==3
 return sub(sq(sub(sub(X[2],X[0]),X[1])),scale(mul(X[0],X[1]),4))
def interval_sum(items):
 out=iv.Z
 for v in items:out=iv.add(out,v)
 return out
def logmass(z):return iv.log_positive(z)
def real_cut(L,s,matched):
 """Enclose Re[A(-s-i0)/C] for a rational interval away from thresholds."""
 z,c=iv.m.spec(L);out=iv.Z
 for mass,ci in zip(z,c):
  beta2=iv.Z if s[0]==s[1]==4*mass else iv.sub((F(1),F(1)),iv.scale(iv.reciprocal(s),4*mass))
  if beta2[0]>=0:
   b=iv.sqrt_iv(beta2)
   if b[1]==0:piece=iv.Z
   else:
    ratio=iv.mul(iv.add((F(1),F(1)),b),iv.reciprocal(iv.sub((F(1),F(1)),b)))
    piece=iv.mul(b,iv.log_iv(ratio))
  else:
   assert beta2[1]<=0
   u=iv.sqrt_iv(iv.neg(beta2))
   # 2u atan(1/u); monotonic interval product, atan increasing.
   if u[0]==0:
    atan=(F(0),iv.PI[1]/2)
   else:atan=(iv.atan_positive(1/u[1])[0],iv.atan_positive(1/u[0])[1])
   piece=iv.scale(iv.mul(u,atan),2)
  I=iv.add(iv.add(logmass(mass),(F(-2),F(-2))),piece)
  U=iv.add((mass,mass),iv.scale(s,F(1,2)))
  term=iv.sub(iv.mul(iv.square(U),I),iv.scale(iv.mul(U,logmass(mass)),mass))
  out=iv.add(out,iv.scale(term,ci/4))
 if matched:
  for LL,sign in [(10,1),(L,-1)]:
   zz,cc=iv.m.spec(LL)
   T1=interval_sum(iv.scale(logmass(a),b*a) for a,b in zip(zz,cc))
   T0=interval_sum(iv.scale(logmass(a),b) for a,b in zip(zz,cc))
   out=iv.add(out,iv.scale(iv.mul(T1,s),F(sign,8)))
   out=iv.add(out,iv.scale(iv.mul(T0,iv.square(s)),F(sign,16)))
 return out

def cut_H(L,s):
 z,c=iv.m.spec(L);out=iv.Z
 for a,b in zip(z,c):
  if s[0]>=4*a:
   beta=iv.sqrt_iv(iv.sub((F(1),F(1)),iv.scale(iv.reciprocal(s),4*a)))
   U=iv.add((a,a),iv.scale(s,F(1,2)))
   out=iv.add(out,iv.scale(iv.mul(iv.square(U),beta),b))
  else:assert s[1]<=4*a
 return out

def ser(ivalue):return dict(exact=[str(x) for x in ivalue],decimal=[float(x) for x in ivalue])
def main():
 rows=[];thresholds=[]
 for L in [3,5,10,30,100]:
  z,c=iv.m.spec(L)
  for n in [2,3]:
   p=polynomial(L,n);seq=sturm(p);a=4*z[n-1];b=4*z[n]
   roots=isolate(seq,a,b);candidates=[]
   for box in roots:
    H=cut_H(L,box);values={str(m):real_cut(L,box,m) for m in [False,True]}
    assert all(v[0]>0 or v[1]<0 for v in values.values())
    candidates.append(dict(s=ser(box),H=ser(H),real_A_over_C={k:ser(v) for k,v in values.items()},imag_zero_not_assumed='Polynomial squaring may introduce candidates; every candidate has separated real part.'))
   rows.append(dict(L=L,active_species=n,open_interval=[str(a),str(b)],polynomial=[str(x) for x in primitive(p)],degree=len(p)-1,sturm_degrees=[len(x)-1 for x in seq],exact_candidate_count=count(seq,a,b),candidates=candidates))
  for s in [4*x for x in z]:
   H=cut_H(L,(s,s));values={str(m):real_cut(L,(s,s),m) for m in [False,True]}
   assert (H[0]>0 or H[1]<0) or all(v[0]>0 or v[1]<0 for v in values.values())
   thresholds.append(dict(L=L,s=str(s),H=ser(H),real_A_over_C={k:ser(v) for k,v in values.items()}))
 global_counts=[]
 for matched in [False,True]:
  for L in [3,5,10,30,100]:
   z,c=iv.m.spec(L)
   threshold=next(x for x in thresholds if x['L']==L and x['s']=='4')
   real4=[F(x) for x in threshold['real_A_over_C'][str(matched)]['exact']]
   initial_phase_pi=-1 if real4[1]<0 else 0
   assert real4[1]<0 or real4[0]>0
   # H at regulator thresholds: +,-,+. Thus lower Im A: -,+,-.
   ts=[x for x in thresholds if x['L']==L and x['s']!='4']
   for t,sgn in zip(ts,[1,-1,1]):
    h=[F(v) for v in t['H']['exact']];assert h[0]>0 if sgn==1 else h[1]<0
   branch_turns=0
   for n,turn in [(2,-1),(3,1)]:
    seg=next(x for x in rows if x['L']==L and x['active_species']==n)
    assert seg['exact_candidate_count']==1
    r=[F(x) for x in seg['candidates'][0]['real_A_over_C'][str(matched)]['exact']]
    if r[1]<0:branch_turns+=turn
    else:assert r[0]>0
   T2=interval_sum(iv.scale(logmass(a),b*a*a) for a,b in zip(z,c))
   assert T2[1]<0
   T1=interval_sum(iv.scale(logmass(a),b*a) for a,b in zip(z,c))
   if matched:
    zr,cr=iv.m.spec(10);T1=interval_sum(iv.scale(logmass(a),b*a) for a,b in zip(zr,cr))
   massless_slope=iv.scale(T1,F(-1,8));assert massless_slope[1]<0
   degree=0;leading=iv.scale(T2,F(3,8))
   if matched and L!=10:
    zz,cc=iv.m.spec(10)
    Tref=interval_sum(iv.scale(logmass(a),b) for a,b in zip(zz,cc))
    TL=interval_sum(iv.scale(logmass(a),b) for a,b in zip(z,c))
    leading=iv.scale(iv.sub(Tref,TL),F(1,16));degree=2
   assert leading[0]>0 or leading[1]<0
   final_principal_pi=-1 if leading[1]<0 else 0
   final_phase_pi=final_principal_pi+2*branch_turns
   total=degree+final_phase_pi-initial_phase_pi
   assert total>=1
   global_counts.append(dict(L=L,matched_to_L10=matched,initial_lower_phase_over_pi=initial_phase_pi,negative_real_crossing_turns=branch_turns,final_lower_phase_over_pi=final_phase_pi,outer_degree=degree,leading_scaled_coefficient=ser(leading),massless_scaled_slope=ser(massless_slope),total_A_zeros_counting_multiplicity=total,massless_zero_multiplicity=1,additional_first_sheet_zeros=total-1))
 result=dict(scope='Global first-sheet scalar zero counts for five raw and explicitly Lref10 matched kernels, conditional on the stated determinant formulas and proved UV asymptotics.',global_counts=global_counts,baseline='b1d6e710fe71ef0e549fc42f31bba20530a4194f',python=platform.python_version(),grid=str(iv.GRID),candidate_segments=rows,thresholds=thresholds,analytic_segments='H>0 for 4<s<4L (one species); H>0 for s>=4zmax by third divided difference, including its continuous endpoint. See REPORT.md.')
 Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(dict(segments=len(rows),candidate_roots=sum(x['exact_candidate_count'] for x in rows),thresholds=len(thresholds),cut_boundary_zero_exclusion=True,additional_counts=[x['additional_first_sheet_zeros'] for x in global_counts])))
if __name__=='__main__':main()
