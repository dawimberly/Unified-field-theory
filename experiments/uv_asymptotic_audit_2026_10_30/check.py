"""Large-q scalar/spin-2 expansion and exact rational exterior bounds.
Run python experiments/uv_asymptotic_audit_2026_10_30/check.py.
The scaled kernels below are A/C, K2/C, C=(4pi)^-2.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
import importlib.util,json,math,platform
ROOT=Path(__file__).resolve().parents[2]
s=importlib.util.spec_from_file_location('exact_helpers',ROOT/'experiments/scalar_threshold_disk_2026_10_30/check.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
N=64

def log_bounds(x):
 k=0;y=x
 while y>=2:y/=2;k+=1
 return m.plus(m.times(m.log_series(F(1,3),128),k),m.log_series((y-1)/(y+1),128))

def add(a,b):return a[0]+b[0],a[1]+b[1]
def scale(a,c):return (a[0]*c,a[1]*c) if c>=0 else (a[1]*c,a[0]*c)
def absolute(a):return max(abs(x) for x in a)
def minimum(a):return a[0] if a[0]>0 else -a[1] if a[1]<0 else F(0)
def expansion(N):
 b=[F(1)]
 for j in range(1,N+1):b.append(b[-1]*F(2*(3-2*j),j))
 h=[F(1)]+[x/2 for x in b[1:]];inv=[F(1)]
 for j in range(1,N+1):inv.append(-sum(h[k]*inv[j-k] for k in range(1,j+1)))
 log=[F(0)]+[sum(k*h[k]*inv[n-k] for k in range(1,n+1))/n for n in range(1,N+1)]
 a=[2*sum(b[k]*log[n-k] for k in range(n+1)) for n in range(N+1)]
 return b,a
B,A=expansion(N)
def channel(ch):return ([F(1,4),F(-1),F(1)],F(1,4),F(3,8)) if ch=='A' else ([F(1,16),F(1,2),F(1)],F(1,15),F(-1,8))
def terms(L,ch):
 z,c=m.spec(L);logs=[log_bounds(x) for x in z];p,norm,constant=channel(ch)
 T=[];S=[]
 for k in range(N+3):
  S.append(sum(ci*x**k for ci,x in zip(c,z)))
  v=(F(0),F(0))
  for ci,x,ln in zip(c,z,logs):v=add(v,scale(ln,ci*x**k))
  T.append(v)
 rows=[]
 for k in range(3,N+3):
  bk=sum(p[j]*B[k-j] for j in range(3) if 1<=k-j<=N)
  ak=sum(p[j]*A[k-j] for j in range(3) if 1<=k-j<=N)
  rows.append((k,bk,ak,S[k],T[k]))
 return z,c,norm,scale(T[2],constant),rows

def remainder(L,ch,R):
 z,c=m.spec(L);p,norm,_=channel(ch);value=F(0)
 for mass,ci in zip(z,c):
  ratio=8*mass/R;assert ratio<=F(1,2)
  log_upper=log_bounds(F(R,mass))[1]+4 # |arg q|<pi<4
  tail=(F(5,4)*log_upper+F(5,6))*ratio**(N+1)/(1-ratio)
  pref=(R/2+mass)**2 if ch=='A' else (R/4+mass)**2
  value+=abs(ci)*pref*tail
 return norm*value

def exterior(L,ch,R):
 z,c,norm,const,rows=terms(L,ch);lr=log_bounds(R)[1]+4
 delta=norm*sum((abs(b)*(abs(sk)*lr+absolute(tk))+abs(a*sk))*R**(2-k) for k,b,a,sk,tk in rows)+remainder(L,ch,R)
 # Lower cut imaginary part divided by pi: leading = norm*b3*S3/s <0.
 _,b3,a3,S3,_=rows[0]
 imag_error=norm*sum(abs(b*sk)*R**(2-k) for k,b,a,sk,tk in rows[1:])
 for mass,ci in zip(z,c):
  ratio=8*mass/R;pref=(R/2+mass)**2 if ch=='A' else (R/4+mass)**2
  imag_error+=norm*abs(ci)*pref*F(5,4)*ratio**(N+1)/(1-ratio)
 return dict(R=int(R),raw_no_zero_margin=minimum(const)-delta,raw_difference_bound=delta,constant=const,leading_cut_over_pi=norm*b3*S3/R,cut_sign_margin=norm*abs(b3*S3)/R-imag_error)

def coefficient_log_uncertainty(L,ch,R):
 z,c,norm,const,rows=terms(L,ch)
 return (const[1]-const[0])/2+norm*sum(abs(b)*(tk[1]-tk[0])*R**(2-k)/2 for k,b,a,sk,tk in rows)

def dec(x):return D(x.numerator)/D(x.denominator)
def atan(x):
 p=x;out=x;n=1
 while True:
  p*=-x*x;new=out+p/D(2*n+1)
  if new==out:return out
  out=new;n+=1

def reference(q,L,ch,cut=False):
 z,c=m.spec(L);pi=16*atan(D(1)/5)-4*atan(D(1)/239);real=D(0);imag=D(0)
 for mass,ci in zip(z,c):
  x=dec(mass);cc=dec(ci)
  if cut:
   beta=(1-4*x/q).sqrt();I=x.ln()-2+beta*((1+beta)/(1-beta)).ln();im=-pi*beta;qq=-q
  else:
   beta=(1+4*x/q).sqrt();I=x.ln()-2+beta*((beta+1)/(beta-1)).ln();im=D(0);qq=q
  if ch=='A':
   u=x-qq/2;real+=cc*(u*u*I-u*x*x.ln())/4;imag+=cc*u*u*im/4
  else:
   u=x+qq/4;real+=cc*(u*u*I-(x*x-x*qq/8)*x.ln())/15;imag+=cc*u*u*im/15
 return real,imag

def approximation(q,L,ch,cut=False):
 z,c,norm,const,rows=terms(L,ch);qq=-q if cut else q;pi=16*atan(D(1)/5)-4*atan(D(1)/239)
 logq=q.ln();real=dec((const[0]+const[1])/2);imag=D(0)
 for k,b,a,sk,tk in rows:
  factor=dec(norm)*qq**(2-k);tval=dec((tk[0]+tk[1])/2)
  real+=factor*(dec(b)*(dec(sk)*logq-tval)+dec(a*sk))
  if cut:imag-=factor*dec(b*sk)*pi
 return real,imag

def main():
 certificates=[];comparisons=[];refG={}
 for LL in [10]:
  z,c=m.spec(LL);gg=(F(0),F(0));ss=(F(0),F(0))
  for x,ci in zip(z,c):gg=add(gg,scale(log_bounds(x),ci*x));ss=add(ss,scale(log_bounds(x),ci))
  refG=gg;refS=ss
 for L in [3,5,10,30,100]:
  z,c=m.spec(L);w0=1/math.prod(z[0]-z[j] for j in range(1,4))
  assert sum(ci*x**3 for ci,x in zip(c,z))==1/w0<0
  gg=(F(0),F(0));ss=(F(0),F(0))
  for x,ci in zip(z,c):gg=add(gg,scale(log_bounds(x),ci*x));ss=add(ss,scale(log_bounds(x),ci))
  da1=scale(add(refG,scale(gg,-1)),F(-1,8));da2=scale(add(refS,scale(ss,-1)),F(1,16))
  for ch in ['A','K2']:
   for multiple in [16,64,256,1024]:
    R=multiple*z[-1];v=exterior(L,ch,R)
    matched_margin=None
    if ch=='A' and L!=10:
     matched_margin=minimum(da2)-absolute(da1)/R-(absolute(v['constant'])+v['raw_difference_bound'])/R**2
    certificates.append(dict(L=L,channel=ch,R=int(R),raw_exterior_zero_excluded=v['raw_no_zero_margin']>0,raw_no_zero_margin=float(v['raw_no_zero_margin']),raw_difference_bound=float(v['raw_difference_bound']),constant_interval=[float(x) for x in v['constant']],negative_lower_cut_proved=v['cut_sign_margin']>0,cut_sign_margin=float(v['cut_sign_margin']),leading_cut_over_pi=float(v['leading_cut_over_pi']),matched_exterior_zero_excluded=matched_margin>0 if matched_margin is not None else (v['raw_no_zero_margin']>0 if L==10 else None),matched_dominance_margin=float(matched_margin) if matched_margin is not None else None))
   with localcontext() as ctx:
    ctx.prec=160
    for multiple in [16,64,256]:
     R=multiple*z[-1]
     for cut in [False,True]:
      direct=reference(dec(R),L,ch,cut);pred=approximation(dec(R),L,ch,cut)
      err=sum((a-b)**2 for a,b in zip(direct,pred)).sqrt();bound=remainder(L,ch,R)
      unc=coefficient_log_uncertainty(L,ch,R)
      assert err<dec(bound+unc)
      comparisons.append(dict(L=L,channel=ch,R=int(R),lower_cut=cut,direct=[str(x) for x in direct],uv_series=[str(x) for x in pred],absolute_error=str(err),analytic_remainder_bound=float(bound),coefficient_log_midpoint_uncertainty=float(unc)))
 result=dict(scope='Raw free PV determinant or explicit scalar matching; exterior/cut bounds do not count zeros in the remaining finite collar or establish a physical UV completion.',versions=dict(python=platform.python_version()),N=N,first_scalar_uv_coefficients=dict(B3=str(terms(3,'A')[-1][0][1]),D3=str(terms(3,'A')[-1][0][2])),first_spin2_uv_coefficients=dict(B3=str(terms(3,'K2')[-1][0][1]),D3=str(terms(3,'K2')[-1][0][2])),exact_rational_certificates=certificates,decimal_160_digit_comparisons=comparisons)
 (Path(__file__).resolve().parent/'results.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(dict(comparisons=len(comparisons),certificates=len(certificates),raw_exterior_successes=sum(x['raw_exterior_zero_excluded'] for x in certificates),negative_cut_successes=sum(x['negative_lower_cut_proved'] for x in certificates),scalar_leading=result['first_scalar_uv_coefficients'],spin2_leading=result['first_spin2_uv_coefficients']),indent=2))
if __name__=='__main__':main()
