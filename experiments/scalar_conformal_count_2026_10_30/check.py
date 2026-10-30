"""Exact conformal coefficients, directed-Decimal contour boxes and log tails.
Run python experiments/scalar_conformal_count_2026_10_30/check.py.
Counts are for |w|<4/5, not the entire first sheet.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext,ROUND_FLOOR,ROUND_CEILING
import importlib.util,json,platform,math,time
from functools import lru_cache
ROOT=Path(__file__).resolve().parents[2]
s=importlib.util.spec_from_file_location('disk_helpers',ROOT/'experiments/scalar_threshold_disk_2026_10_30/check.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
PREC=60

def down(f):
 with localcontext() as ctx:ctx.prec=PREC;ctx.rounding=ROUND_FLOOR;return f()
def up(f):
 with localcontext() as ctx:ctx.prec=PREC;ctx.rounding=ROUND_CEILING;return f()
def rational(q):return down(lambda:D(q.numerator)/D(q.denominator)),up(lambda:D(q.numerator)/D(q.denominator))
def ia(a,b):return down(lambda:a[0]+b[0]),up(lambda:a[1]+b[1])
def neg(a):return a[1].copy_negate(),a[0].copy_negate()
def im(a,b):
 return min(down(lambda:x*y) for x in a for y in b),max(up(lambda:x*y) for x in a for y in b)
def ca(a,b):return ia(a[0],b[0]),ia(a[1],b[1])
def cm(a,b):return ia(im(a[0],b[0]),neg(im(a[1],b[1]))),ia(im(a[0],b[1]),im(a[1],b[0]))
ZERO=(D(0),D(0))
def evaluate(co,w):
 out=(ZERO,ZERO)
 for c in reversed(co):out=ca(cm(out,w),(c,ZERO))
 return out
def gap(box):
 gaps=[]
 for lo,hi in box:
  gaps.append(lo if lo>0 else hi.copy_negate() if hi<0 else D(0))
 return max(gaps)
def pt(t,rot,r):
 x=r*(1-t*t)/(1+t*t);y=2*r*t/(1+t*t)
 for _ in range(rot):x,y=-y,x
 return x,y
def box_between(a,b):
 xr=[rational(a[0]),rational(b[0])];yi=[rational(a[1]),rational(b[1])]
 return (min(x[0] for x in xr),max(x[1] for x in xr)),(min(y[0] for y in yi),max(y[1] for y in yi))
def node_box(p):return rational(p[0]),rational(p[1])
def midpoint(box):return tuple((F(a)+F(b))/2 for a,b in box)
def cheb_averages(z,N):
 if z==1:return [F(1)]+[F(1,1-4*n*n) for n in range(1,N+1)]
 prev=[F(1)];cur=[F(1),-8/z];out=[F(1)]
 mu=[F(math.factorial(k)**2,math.factorial(2*k+1)) for k in range(N+1)]
 for n in range(1,N+1):
  out.append(sum(a*b for a,b in zip(cur,mu)))
  nxt=[F(0)]*(len(cur)+1)
  for k,a in enumerate(cur):nxt[k]+=2*a;nxt[k+1]-=16*a/z
  for k,a in enumerate(prev):nxt[k]-=a
  prev,cur=cur,nxt
 return out
@lru_cache(maxsize=None)
def coefficients(L,N,matched):
 z,c=m.spec(L);co=[F(0)]*(N+4)
 for mass,ci in zip(z,c):
  av=cheb_averages(mass,N);p=[mass,-2*mass-8,mass];p2=[F(0)]*5
  for i,a in enumerate(p):
   for j,b in enumerate(p):p2[i+j]+=a*b
  # Common -2 log(1-w) cancels exactly by the three PV moments.
  for n in range(1,N+1):
   for k,a in enumerate(p2):co[n-1+k]+=ci*a*(-2*av[n]/n)
 intervals=[rational(a) for a in co]
 zz,cc=m.spec(10 if matched else L)
 G=(F(0),F(0));S=(F(0),F(0))
 for mass,ci in zip(zz,cc):
  ln=m.log_bounds(mass);G=m.plus(G,m.times(ln,ci*mass));S=m.plus(S,m.times(ln,ci))
 logs=[m.times(G,-8),m.plus(m.times(G,16),m.times(S,64)),m.times(G,-8)]
 for k,v in enumerate(logs):intervals[k]=ia(intervals[k],(rational(v[0])[0],rational(v[1])[1]))
 return intervals

def tail(L,N,r):
 z,c=m.spec(L)
 return 2*sum(abs(ci)*(mass*(1+r)**2+8*r)**2 for mass,ci in zip(z,c))*r**N/(F(N+1)*(1-r))
def count(L,matched,N=128,r=F(4,5),center=(F(0),F(0)),path_radius=None):
 start=time.monotonic();co=coefficients(L,N,matched);rem=rational(tail(L,N,r))[1]
 # Exact-rational coefficient norm bound for sup |P_prime| on |w|<=r.
 M=sum(F(n)*max(abs(F(c[0])),abs(F(c[1])))*r**(n-1) for n,c in enumerate(co) if n)
 M2=sum(F(n*(n-1))*max(abs(F(c[0])),abs(F(c[1])))*r**(n-2) for n,c in enumerate(co) if n>1)
 derivative=[(down(lambda:c[0]*D(n)),up(lambda:c[1]*D(n))) for n,c in enumerate(co) if n]
 rad=r if path_radius is None else path_radius
 assert center[0]**2+center[1]**2 <= (r-rad)**2
 def point(t,rot):
  x,y=pt(t,rot,rad);return x+center[0],y+center[1]
 nodes=[];gaps=[];maxdepth=0;boxes=0
 def arc(t0,t1,rot,depth):
  nonlocal maxdepth,boxes
  a=point(t0,rot);b=point(t1,rot)
  image=evaluate(co,node_box(a));g=gap(image);boxes+=1
  dp=evaluate(derivative,node_box(a))
  dp_bound=sum(max(abs(F(v[0])),abs(F(v[1]))) for v in dp)
  distance=abs(b[0]-a[0])+abs(b[1]-a[1])
  excursion=rational(dp_bound*distance+M2*distance*distance/2)[1]
  g=down(lambda:g-excursion)
  if g>rem:
   endpoint=evaluate(co,node_box(b))
   hull=tuple((min(down(lambda:v[0]-excursion),e[0]),max(up(lambda:v[1]+excursion),e[1])) for v,e in zip(image,endpoint))
   g=gap(hull)
  if g>rem:
   nodes.append(midpoint(image));gaps.append(down(lambda:g-rem));maxdepth=max(maxdepth,depth);return
  if depth>=20:raise RuntimeError('contour interval did not resolve')
  tm=(t0+t1)/2;arc(t0,tm,rot,depth+1);arc(tm,t1,rot,depth+1)
 for rot in range(4):arc(F(0),F(1),rot,0)
 # Exact rational winding of the nodal image polygon; rectangle homotopies avoid zero.
 winding=0
 for a,b in zip(nodes,nodes[1:]+nodes[:1]):
  cross=a[0]*b[1]-b[0]*a[1]
  if a[1]<=0<b[1] and cross>0:winding+=1
  elif b[1]<=0<a[1] and cross<0:winding-=1
 return dict(L=L,matched_to_L10=matched,N=N,w_radius=str(r),path_center=[str(x) for x in center],path_radius=str(rad),zeros_excluding_massless=winding,accepted_arc_boxes=len(nodes),total_interval_boxes=boxes,maximum_bisection_depth=maxdepth,log_tail_bound=str(rem),minimum_component_gap_after_tail=str(min(gaps)),decimal_precision=PREC,polynomial_derivative_bound=float(M),polynomial_second_derivative_bound=float(M2),elapsed_seconds=time.monotonic()-start)
def arithmetic_audit():
 cases=0
 for a in [F(-123456789,12345),F(-1,7),F(0),F(1,3),F(987654321,23456)]:
  for b in [F(-9,11),F(1,19),F(1234567,76543)]:
   aa=(rational(a),rational(b));bb=(rational(b),rational(a/2))
   prod=cm(aa,bb);expected=(a*b/2,a*a/2+b*b)
   assert all(F(iv[0])<=value<=F(iv[1]) for value,iv in zip(expected,prod))
   mid=midpoint(prod);assert all(F(iv[0])<=value<=F(iv[1]) for value,iv in zip(mid,prod))
   cases+=1
 return dict(exact_complex_product_and_midpoint_cases=cases,directed_decimal_precision=PREC)

def main():
 arithmetic_audit()
 rows=[]
 for L in [3,5,10,30,100]:
  for matched in [False,True]:
   row=count(L,matched);rows.append(row);print(json.dumps(row),flush=True)
 result=dict(scope='Free scalar PV determinant or explicitly Lref=10 matched version, whole first-sheet conformal disk |w|<4/5. Not global first-sheet control or physical viability.',versions=dict(python=platform.python_version()),map=dict(q_of_w='16w/(1-w)^2',w_of_q='(sqrt(1+q/4)-1)/(sqrt(1+q/4)+1)',principal_sqrt='positive at q=0, Re sqrt>0 off cut',function='G(w)=4(1-w)^4 A(q(w))/(C w), removable at w=0'),counts=rows)
 (Path(__file__).resolve().parent/'results.json').write_text(json.dumps(result,indent=2))
if __name__=='__main__':main()
