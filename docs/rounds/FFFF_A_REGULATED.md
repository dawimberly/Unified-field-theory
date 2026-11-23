# FFff(A) regulated extrapolation

Continues FFFF_A_CORNER_CORRECTION.md and FFFF_A_LONG_PAIR.md. This is a fixed-1/6, nonzero-ε calculation with z=μ² in the extraction. Both global long-pair models are subtracted and added back analytically. The effective variables P,Q use cutoff Taylor subtractions through degree one on both axes. The code keeps the finite-dimensional Weyl factors, including the radial Γ(3ε−1).

Use the conditioned moment T=N M. With the normalized bridge covariances and reduced vectors a and b in the code, the connected contribution to T is multiplied by d/S³. This factor follows from the explicit 2P/S and 2Q/S covariances before cancelling the spin denominators, and is needed to recover the earlier exact corner. The two tails are integrated as a sum of Dirichlet integrals. Their leading coefficient is 2/ε². The full regulated integral supports that leading coefficient after the finite remainder is included. The Abelian Taylor integrals of the short constant and linear terms are cut^(ε−1)/(ε−1) and cut^ε/ε respectively.

At ε=.001875 the three cuts 0.5,1,2 give ε²I≈1.99380. Quadratic extrapolations of ε(I−2/ε²) give intercepts −3.33400, −3.33324, −3.33343, supporting −10/3 for the next projective coefficient. Slopes are near 14.93 for the finite coefficient. This is numerical evidence for the full graph, not an analytic evaluation of its remainder. The cutoff spread and truncated ε fit do not supply rigorous error bars.

For the bare graph multiply by −Γ(3ε−1) exp(3γ_E ε−3ε L)/[32d(d−1)]. The resulting numerical inference in units h²z/K is

- Triple pole consistent with 1/576.
- Double pole consistent with 5/1152−L/192.
- Simple pole approximately 0.03724−5L/384+L²/128. The finite coefficient is more sensitive to extrapolation than the higher poles.

Using the independently derived combined FFff higher-pole targets, the implied benchmarks for FFff(B) are 5/192 at triple order and 11/1152−5L/64 at double order. These are predictions to check, not an evaluation of FFff(B). Its simple pole is still needed for the Newton test.

## Runnable check

Python 3, NumPy and SciPy. The three cutoffs test the numerical subtraction assembly.

```python
import numpy as np
from scipy.special import beta,gamma
from scipy.stats import qmc
def H(P,Q,f,g,v,e):
 d=4-2*e;k=3*e
 q,r=P/f,Q/g
 S,W=1+P+Q,1+q+r
 t=[0,P,P+v,P+v+Q]
 def C(i,j):return 2*t[min(i,j)]*(S-t[max(i,j)])/S
 N=d*(d*(1+Q)*(1+P)+2*P*Q)/(4*S*S)
 a=[0,S-P,S-t[2],S-t[3]]
 b=[0,-P,-t[2],S-t[3]]
 def T(i,j):
  return d*C(i,j)*N+d/S**3*(d*(P*(1+P)*a[i]*a[j]+Q*(1+Q)*b[i]*b[j])-2*P*Q*(a[i]*b[j]+b[i]*a[j]))
 A=-(4-d)*(T(1,1)+T(2,2)+T(3,3))/2-d*(4-d)*N/3
 B=2/3*(q*T(1,1)+r*(T(2,2)+T(3,3)+T(2,3))+v*(T(1,1)+T(2,2)+T(1,2))+(1-v)*T(3,3)+N*((d-1)*(q*q*(1-2*f)+r*r*(1-2*g))+d*(v*v+(1-v)**2)))
 return S**(-d/2)*(W**(1-k)*A+(k-1)*(W**(-k)*B-2*(d-1)*N*(q**(2-k)+r**(2-k))/3))
def tail(e):
 d=4-2*e;k=3*e
 # Long Q tail: P^(e-2) Q^(-2e) and rational spin factor.
 def D(i,j):
  a=e-1+i;b=1-2*e+j;c=4-e
  return gamma(a)*gamma(b)*gamma(c-a-b)/gamma(c)
 return (k-1)*2*(d-1)/3*d/4*beta(2-e,2-e)*beta(2*e,2*e)*(d*(D(0,0)+D(1,0)+D(0,1)+D(1,1))+2*D(1,1))
def run(e,seed,power,cut=1.):
 z=qmc.Sobol(5,scramble=True,seed=seed)
 out=[]
 # Real complex-step derivative; mixed derivative via centered difference
 # of the exact complex-step first derivative.
 h=1e-5
 def dp(P,Q,f,g,v):return H(P+1e-27j,Q,f,g,v,e).imag/1e-27
 def dq(P,Q,f,g,v):return H(P,Q+1e-27j,f,g,v,e).imag/1e-27
 for _ in range(2**(power-14)):
  x,y,a,b,v=z.random(2**14).T
  P,Q=x/(1-x),y/(1-y);f,g=a*(1-a),b*(1-b)
  X,Y=P<cut,Q<cut
  wp=P**(e-2)/(1-x)**2;wq=Q**(e-2)/(1-y)**2
  c=H(0,0,f,g,v,e);cp=dp(0,0,f,g,v);cq=dq(0,0,f,g,v)
  # Mixed derivative of powers P^(2-3e) vanishes at zero.
  # They contaminate a finite difference, so remove all four one-variable tails.
  def regular(p,q):
   d=4-2*e;k=3*e;S=1+p+q
   N=d*(d*(1+p)*(1+q)+2*p*q)/(4*S*S)
   return H(p,q,f,g,v,e)+(k-1)*2*(d-1)/3*N*S**(-d/2)*((p/f)**(2-k)+(q/g)**(2-k))
  d=4-2*e
  # Before averaging v and split fractions the mixed coefficient is tested elsewhere.
  # Compute using positive complex-step derivative and difference, adaptive to g.
  step=h*g
  dp0=lambda q:regular(1e-25j,q).imag/1e-25
  mixed=( -3*dp0(0)+4*dp0(step)-dp0(2*step))/(2*step)
  U=H(P,0,f,g,v,e);V=H(0,Q,f,g,v,e)
  D=dp(0,Q,f,g,v);E=dq(P,0,f,g,v)
  M=H(P,Q,f,g,v,e)-X*(V+P*D)-Y*(U+Q*E)+(X&Y)*(c+P*cp+Q*cq+P*Q*mixed)
  a0=cut**(e-1)/(e-1); a1=cut**e/e
  out.append(np.mean((f*g)**(1-e)*(wp*wq*M+wq*(V-Y*(c+Q*cq))*a0+wq*(D-Y*(cp+Q*mixed))*a1+wp*(U-X*(c+P*cp))*a0+wp*(E-X*(cq+P*mixed))*a1+c*a0*a0+(cp+cq)*a0*a1+mixed*a1*a1)))
 return np.mean(out)+2*tail(e)
for e in (.0075,.00375,.001875):
 for cut in (.5,1.,2.):
  I=run(e,12,18,cut)
  print(e,cut,e*e*I,e*(I-2/(e*e)),flush=True)

```

- [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177).
