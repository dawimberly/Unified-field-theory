"""Independent exact-rational local Rouche certificate for the L=30 collar pole.
No NumPy complex evaluation is used in this certificate. Standard library plus
previous exact scalar helper imports. Run this file before the revised check.py.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import importlib.util,math,json
ROOT=Path(__file__).resolve().parents[2]
s=importlib.util.spec_from_file_location('moments',ROOT/'experiments/scalar_threshold_disk_2026_10_30/check.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
GRID=10**30

def iv(a,b=None):
 b=a if b is None else b
 return F((a*GRID).__floor__(),GRID),F((b*GRID).__ceil__(),GRID)
def add(a,b):return iv(a[0]+b[0],a[1]+b[1])
def neg(a):return -a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 v=[x*y for x in a for y in b];return iv(min(v),max(v))
def scale(a,c):return mul(a,(c,c))
def square(a):return iv(F(0) if a[0]<=0<=a[1] else min(a[0]**2,a[1]**2),max(a[0]**2,a[1]**2))
def reciprocal(a):
 assert a[0]>0 or a[1]<0
 return iv(1/a[1],1/a[0])
def sqrt_point(x):
 assert x>=0
 k=math.isqrt((x*GRID*GRID).__floor__())
 return F(k,GRID),F(k if F(k*k,GRID*GRID)==x else k+1,GRID)
def sqrt_iv(a):return iv(sqrt_point(max(F(0),a[0]))[0],sqrt_point(a[1])[1])
def abs_upper(a):return max(abs(x) for x in a)
def distance(a):return a[0] if a[0]>0 else -a[1] if a[1]<0 else F(0)
Z=(F(0),F(0))
def complex_iv(a,b=F(0)):return iv(a),iv(b)
def ca(a,b):return add(a[0],b[0]),add(a[1],b[1])
def cn(a):return neg(a[0]),neg(a[1])
def cs(a,b):return ca(a,cn(b))
def cm(a,b):return sub(mul(a[0],b[0]),mul(a[1],b[1])),add(mul(a[0],b[1]),mul(a[1],b[0]))
def cscale(a,c):return scale(a[0],c),scale(a[1],c)
def ci(a):
 inv=reciprocal(add(square(a[0]),square(a[1])));return mul(a[0],inv),neg(mul(a[1],inv))
def cd(a,b):return cm(a,ci(b))
def log_positive(x):
 assert x>0
 return iv(*m.log_bounds(x)) if x>=1 else iv(*neg(m.log_bounds(1/x)))
def log_iv(a):
 assert a[0]>0
 return iv(log_positive(a[0])[0],log_positive(a[1])[1])
@lru_cache(None)
def atan_small(x):
 if x<0:return neg(atan_small(-x))
 assert x<=F(1,2)
 N=64;p=sum(F((-1)**n,2*n+1)*x**(2*n+1) for n in range(N));nxt=x**(2*N+1)/(2*N+1)
 return iv(p,p+nxt)
PI=sub(scale(atan_small(F(1,5)),16),scale(atan_small(F(1,239)),4))
def atan_positive(x):
 assert x>=0
 if x>1:return sub(scale(PI,F(1,2)),atan_positive(1/x))
 return add(atan_small(F(1,2)),atan_small((x-F(1,2))/(1+x/2)))
def argument(x,y):
 assert y>0
 if x>0:return atan_positive(y/x)
 if x<0:return sub(PI,atan_positive(y/(-x)))
 return scale(PI,F(1,2))
def clog(a):
 assert a[1][0]>0 # This certificate stays strictly above the logarithm cut.
 real=scale(log_iv(add(square(a[0]),square(a[1]))),F(1,2))
 corners=[argument(x,y) for x in a[0] for y in a[1]]
 return real,iv(min(v[0] for v in corners),max(v[1] for v in corners))
def csqrt(a):
 assert a[1][1]<0 # beta^2 lies strictly in the lower half-plane.
 modulus=sqrt_iv(add(square(a[0]),square(a[1])))
 real=sqrt_iv(scale(add(modulus,a[0]),F(1,2)))
 imag=neg(sqrt_iv(scale(sub(modulus,a[0]),F(1,2))))
 return real,imag

def matching():
 vals=[]
 for L in [10,30]:
  z,c=m.spec(L);g=Z;s=Z
  for mass,ci0 in zip(z,c):
   ln=log_positive(mass);g=add(g,scale(ln,ci0*mass));s=add(s,scale(ln,ci0))
  vals.append((g,s))
 return scale(sub(vals[0][0],vals[1][0]),F(-1,8)),scale(sub(vals[0][1],vals[1][1]),F(1,16))
def point(center):
 q=complex_iv(*center);one=complex_iv(F(1));out=(Z,Z);der=(Z,Z);d1,d2=matching();z,c=m.spec(30)
 for mass,ci0 in zip(z,c):
  beta=csqrt(ca(one,cscale(ci(q),4*mass)))
  ratio=cd(ca(beta,one),cs(beta,one));lnratio=clog(ratio)
  I=ca(complex_iv(F(0)),cm(beta,lnratio));I=ca(I,(add(log_positive(mass),(F(-2),F(-2))),Z))
  Ip=cs(ci(q),cscale(cd(lnratio,cm(cm(q,q),beta)),2*mass))
  U=cs(complex_iv(mass),cscale(q,F(1,2)));lnmass=(log_positive(mass),Z)
  out=ca(out,cscale(cs(cm(cm(U,U),I),cscale(cm(U,lnmass),mass)),ci0/4))
  dd=ca(cs(cm(cm(U,U),Ip),cm(U,I)),cscale(lnmass,mass/2))
  der=ca(der,cscale(dd,ci0/4))
 out=ca(out,ca(cm((d1,Z),q),cm((d2,Z),cm(q,q))))
 der=ca(der,ca((d1,Z),cscale(cm((d2,Z),q),2)))
 return out,der,d2

def main():
 root=json.loads((Path(__file__).with_name('collar_search_diagnostic.json')).read_text())['roots'][0]
 center=tuple(F(str(x)) for x in root);trace=[]
 for _ in range(4):
  value,der,d2=point(center);a,b=[(v[0]+v[1])/2 for v in value];c,d=[(v[0]+v[1])/2 for v in der];den=c*c+d*d
  trace.append(dict(center=[float(x) for x in center],scaled_residual_upper=float(abs_upper(value[0])+abs_upper(value[1]))))
  re=center[0]-(a*c+b*d)/den;im=center[1]-(b*c-a*d)/den
  center=(F((re*GRID).__floor__(),GRID),F((im*GRID).__floor__(),GRID))
 r=F(1,10**10);value,der,d2=point(center)
 Q=abs(center[0])+abs(center[1])+r;y=center[1]-r;assert y>0
 M2=F(0);z,c=m.spec(30)
 for mass,ci0 in zip(z,c):
  lower=mass*y/Q;upper=mass+Q/4
  logbound=max(abs_upper(log_positive(lower)),abs_upper(log_positive(upper)))+4
  U=mass+Q/2
  M2+=abs(ci0)*(logbound/2+2*U/y+U*U/(y*y))/4
 M2+=2*abs_upper(d2)
 error=abs_upper(value[0])+abs_upper(value[1])+M2*r*r/2
 linear=max(distance(der[0]),distance(der[1]))*r
 assert error<linear
 result=dict(L=30,scope='Exact rational local Rouche certificate for A/C; conjugation gives the lower partner. Not a global count.',center=[str(x) for x in center],center_float=[float(x) for x in center],midpoint_newton_trace=trace,radius=str(r),rounding_grid_denominator=str(GRID),center_value_enclosure=[[float(x) for x in v] for v in value],center_derivative_enclosure=[[float(x) for x in v] for v in der],scaled_second_derivative_bound=float(M2),rouche_error_upper=float(error),linear_term_lower=float(linear),strict_margin=float(linear-error),exact_inequality=True,zeros_inside=1)
 Path(__file__).with_name('local_pole_certificate.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
