"""Exact flat metric/source and covariant counterterm audit.
No loop reruns, new physical action, Lorentzian contour or equivalence assumption.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations
import json,platform
HERE=Path(__file__).resolve().parent
D=4
I=[[F(i==j) for j in range(D)] for i in range(D)]
def tr(a):return sum(a[i][i] for i in range(D))
def dot(a,b):return sum(a[i][j]*b[i][j] for i in range(D) for j in range(D))
def scale(a,c):return [[c*x for x in r] for r in a]
def add(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def mm(a,b):return [[sum(a[i][k]*b[k][j] for k in range(D)) for j in range(D)] for i in range(D)]
def outer(a,b):return [[x*y for y in b] for x in a]
def P(h,j):return dot(h,j)-tr(h)*tr(j)/2

def gamma(Q,h):
 return [[[sum(F(Q[b])*h[a][c]+F(Q[c])*h[a][b]-F(Q[a])*h[b][c] for _ in [0])/2 for c in range(D)] for b in range(D)] for a in range(D)]
def EH_ordered(G,H):
 return sum(G[a][mu][b]*H[b][mu][a]-G[a][mu][mu]*H[b][a][b] for a in range(D) for b in range(D) for mu in range(D))
def Ricci(Q,h):
 q=sum(x*x for x in Q);tau=tr(h);qh=[sum(Q[a]*h[a][b] for a in range(D)) for b in range(D)]
 return [[(q*h[a][b]+Q[a]*Q[b]*tau-Q[a]*qh[b]-Q[b]*qh[a])/2 for b in range(D)] for a in range(D)]
def scalarR(Q,h):return sum(x*x for x in Q)*tr(h)-sum(Q[a]*Q[b]*h[a][b] for a in range(D) for b in range(D))

def curvature_audit():
 rows=[]
 seed=[[F((i+2*j)%5-2,7) for j in range(D)] for i in range(D)]
 seed=scale(add(seed,list(map(list,zip(*seed)))),F(1,2))
 couplings=[(F(1),F(0),F(0)),(F(0),F(1),F(0)),(F(0),F(0),F(1)),(F(2,7),F(-3,11),F(5,13)),(F(0),F(1),F(2))]
 for Q in [[F(0),F(0),F(0),F(1)],[F(1),F(1),F(0),F(1)],[F(1),F(2),F(2),F(0)]]:
  q=sum(x*x for x in Q);theta=add(I,scale(outer(Q,Q),-1/q))
  transverse=mm(mm(theta,seed),theta);TT=add(transverse,scale(theta,-tr(transverse)/3))
  v=[F(1,3),F(-1,2),F(2,5),F(3,7)];gauge=add(outer(Q,v),outer(v,Q))
  hs=[I,seed,theta,TT,gauge]
  for kappa,alpha,beta in couplings:
   errors=[]
   for h in hs:
    for j in hs:
     G=gamma(Q,h);H=gamma(Q,j)
     direct=kappa*(EH_ordered(G,H)+EH_ordered(H,G))+2*alpha*scalarR(Q,h)*scalarR(Q,j)+2*beta*dot(Ricci(Q,h),Ricci(Q,j))
     th=mm(mm(theta,h),theta);tj=mm(mm(theta,j),theta)
     scalar=dot(theta,h)*dot(theta,j)/3
     spin2=dot(th,tj)-scalar
     K2=-kappa*q/2+beta*q*q/2
     A=3*kappa*q/2+(9*alpha+3*beta)*q*q
     predicted=K2*spin2+(2*A/3)*scalar
     assert direct==predicted
     errors.append(direct-predicted)
   rows.append(dict(Q=[str(x) for x in Q],q=str(q),couplings=[str(x) for x in [kappa,alpha,beta]],bilinear_cases=len(errors),maximum_error='0'))
 # Independent invariant directions and the minimal scalar regulator difference.
 da2=F(1,16);dk2=F(1,240)
 beta=2*dk2;alpha=(da2-6*dk2)/9
 assert (alpha,beta)==(F(1,240),F(1,120))
 assert (3*F(-1,12)/2,-F(-1,12)/2)==(F(-1,8),F(1,24))
 return rows,dict(convention='Gamma_ct=integral sqrt(g)[kappa R+alpha R^2+beta Ricci^2], q=Q^2, D=4; topological/boundary terms omitted',delta_A='3*kappa*q/2+(9*alpha+3*beta)*q^2',delta_K2='-kappa*q/2+beta*q^2/2',inverse_map=dict(kappa='-2*delta_k1',beta='2*delta_k2',alpha='(delta_a2-6*delta_k2)/9'),minimal_scalar_difference_units='C times difference of sum c log z (curvature squared); C times difference of sum c z log z (Einstein)',minimal_scalar_coefficients=dict(alpha=str(alpha),beta=str(beta),kappa=str(F(-1,12))),flat_tadpole='All displayed terms have vanishing flat first variation under compact-support/periodic boundary conditions; curvature-squared terms begin at q^2 in Hessian.')

# Two nilpotent variables eps,eta: coefficients (1,eps,eta,eps*eta).
def ja(a,b):return tuple(x+y for x,y in zip(a,b))
def js(a,c):return tuple(c*x for x in a)
def jm(a,b):return (a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[0]*b[2]+a[2]*b[0],a[0]*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0])
def detjet(g):
 ans=(F(0),)*4
 for p in permutations(range(D)):
  sign=(-1)**sum(p[i]>p[j] for i in range(D) for j in range(i+1,D));term=(F(sign),F(0),F(0),F(0))
  for i in range(D):term=jm(term,g[i][p[i]])
  ans=ja(ans,term)
 return ans

def direct_contact(h,j,p,lam):
 g=[[(I[a][b],h[a][b],j[a][b],F(0)) for b in range(D)] for a in range(D)]
 detg=detjet(g);assert detg[0]==1
 delta=(F(0),)+detg[1:]
 sqrtg=ja(ja((F(1),F(0),F(0),F(0)),js(delta,F(1,2))),js(jm(delta,delta),F(-1,8)))
 hj=add(mm(h,j),mm(j,h));trace_inverse=(F(D),-tr(h),-tr(j),tr(hj))
 lag=ja(js(trace_inverse,p/2),(-2*p+lam,F(0),F(0),F(0)))
 return jm(sqrtg,lag)[3]

def contact_audit():
 seed=[[F((i+j)%3-1,5) for j in range(D)] for i in range(D)]
 hs=[I,[[F((1 if i==0 else -1 if i==1 else 0) if i==j else 0) for j in range(D)] for i in range(D)],[[F((i,j) in [(0,1),(1,0)]) for j in range(D)] for i in range(D)],seed,add(seed,I)]
 p=F(7,11);rows=[]
 for lam in [F(0),p]:
  failures_of_universal_zero=0
  for i,h in enumerate(hs):
   for k,j in enumerate(hs):
    # Inverse-metric derivative contact and coordinate-reparametrization contact.
    inverse_part=-p*tr(h)*tr(j)/2
    reparam_part=p*dot(h,j)
    matter=inverse_part+reparam_part
    total=matter-lam*P(h,j)/2
    direct=direct_contact(h,j,p,lam);assert total==direct
    assert matter==p*P(h,j)
    if total:failures_of_universal_zero+=1
    rows.append(dict(lambda_over_p=str(lam/p),h=i,j=k,inverse_metric_contact=str(inverse_part),inverse_to_covariant_contact=str(reparam_part),stationary_or_cosmological_contact=str(-lam*P(h,j)/2),direct_nilpotent_second_variation=str(direct),total_contact=str(total)))
  assert failures_of_universal_zero>0
 return rows

def main():
 curv,dictionary=curvature_audit();contacts=contact_audit()
 result=dict(scope='Exact Euclidean flat-source counterexamples and covariant quadratic dictionary, not a literal proof about an inaccessible final publication or a new physical model.',baseline='fe85718a13098c30b01be0cd0019dc67a8d0194d',python=platform.python_version(),curvature_bilinears=curv,covariant_matching_dictionary=dictionary,canonical_contact_audit=contacts,source_convention_warning='Tadpole stationarity alone does not null the fixed-field seagull. Metric-dependent operator subtractions must be differentiated separately. No published error or quantum equivalence is asserted.')
 (HERE/'results.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(dict(exact_curvature_bilinears=sum(x['bilinear_cases'] for x in curv),exact_contact_bilinears=len(contacts),maximum_error='0',curvature_dictionary=dictionary),indent=2))
if __name__=='__main__':main()
