# CcFf regulated pole cross-check

Continues [the analytical long-pair model](ffz_ccff_long_pair.md). This checks the complete CcFf graph integrand at nonzero ε in the equal-mass fixed-1/6 benchmark. It does not reuse the individual short-face integrals. The global long scalar-pair model is subtracted from the integrand and its analytic integral is added back.

## Implementation

The anchored cycle positions along its perimeter are 0,P,P+1−v,P+1−v+Q. The Brownian-bridge covariance is 2 min(t_i,t_j)[1−max(t_i,t_j)/S]. For the spin-conditioned covariance use

M_ij = d C_ij + 4Q l_i l_j/[S(1+P)],
l=(0,P,P+1−v,−v).

This formula avoids dividing by Q at the fermion face. The regulated integrand has effective-length weights P^(ε−1)Q^(ε−2) f_r^(−ε)g_r^(1−ε), where f_r=f(1−f), g_r=g(1−g). For Q<1 subtract H(P,0)+Q∂_QH(P,0), obtaining the derivative by complex step. Subtract the P=0 overlap for P<1 from each term. The added-back effective-length integrals are 1/ε, 1/(ε−1), or 1/ε for the scalar face, fermion constant or fermion linear term, respectively. The independent cutoff in effective lengths differs from the earlier local subtraction in physical pair totals. Only the global tail model is shared.

Using three scrambled quasi-Monte Carlo streams, ε times the full projective integral tends to −2, supporting the predicted double pole +1/(144ε²) and the corresponding scale logarithm −L/(48ε). At ε=.0025 the estimates range from −1.97664 to −1.97648, consistent with a smooth finite correction. After subtracting the double pole, positive-ε extrapolation gives a simple-pole constant about −0.0247 at L=0. This constant is numerical. The stream spread is a sampling diagnostic rather than a rigorous error bound.

In units h²z/K the checked CcFf contribution is therefore

1/(144ε²) + [approximately −0.0247 − L/48]/ε.

Adding it to the scalar and QFf subtotal gives a double pole 5/(144ε²) and a simple pole approximately [−0.1506−5L/48]/ε. Both FFff graphs and the matter counterterms still enter the projected Newton sum. The absence of the full projected higher poles and the comparison with 1/24 remain open.

## Runnable check

Python 3 with NumPy and SciPy. A power of 22 increases the point count eightfold compared with 19. The run below prints the two coefficients for three independent streams. Use one BLAS thread to reproduce the local timing.

```python
import numpy as np
from scipy.stats import qmc
from scipy.special import beta,gamma
from math import exp
def H(P,Q,f,g,v,e):
 d=4-2*e; k=3*e
 p,q=P/f,Q/g
 S,W=1+P+Q,1+p+q
 ts=[0,P,P+1-v,P+1-v+Q]
 ls=[0,P,P+1-v,-v]
 def M(i,j):
  return 2*d*ts[min(i,j)]*(S-ts[max(i,j)])/S+4*Q*ls[i]*ls[j]/(S*(1+P))
 diag=[M(i,i) for i in range(4)]
 A=-(6-d)/2*diag[1]-(4-d)/2*(diag[2]+diag[3])-d*(4-d)*(1+p)/3
 B=2/3*(p*diag[1]+d*p*p*(1-2*f)+v*(diag[3]+d*v)
     +(1-v)*(diag[1]+diag[2]+M(1,2)+d*(1-v))
     +q*(diag[2]+diag[3]+M(2,3))+(d-1)*q*q*(1-2*g))
 tail=4*d*e/3*p**(1-3*e)
 return d*(1+P)/2*S**(-1-d/2)*(W**(-k)*(A+k*B/W)-tail)
def cterms(P,f,g,v,e):
 return H(P,0,f,g,v,e),H(P,1e-25j,f,g,v,e).imag/1e-25
def run(e,seed,power):
 sums=np.zeros(6)
 sob=qmc.Sobol(5,scramble=True,seed=seed)
 batch=2**15
 for _ in range(2**power//batch):
  x,y,a,b,v=sob.random(batch).T
  P,Q=x/(1-x),y/(1-y)
  f,g=a*(1-a),b*(1-b)
  X,Y=P<1,Q<1
  wp=P**(e-1)/(1-x)**2*f**(-e)
  wq=Q**(e-2)/(1-y)**2*g**(1-e)
  C,D=cterms(P,f,g,v,e)
  C0,D0=cterms(0,f,g,v,e)
  F=H(P,Q,f,g,v,e)-Y*(C+Q*D)
  F0=H(0,Q,f,g,v,e)-Y*(C0+Q*D0)
  sums+=np.array([np.mean(val) for val in (
   wp*wq*(F-X*F0),wq*F0,wp*g**(1-e)*(C-X*C0),
   wp*g**(1-e)*(D-X*D0),g**(1-e)*C0,g**(1-e)*D0)])
 m,p,c,d,c0,d0=sums/(2**power//batch)
 fac=beta(1-e,1-e)/e
 tail=2*(4-2*e)**2*e/3*beta(2*e,2*e)*beta(2-e,2-e)*beta(e-1,4-2*e)*beta(1-2*e,2)
 I=m+fac*p+(c+fac*c0)/(e-1)+(d+fac*d0)/e+tail
 graph=-I*gamma(3*e)*exp(3*.5772156649015329*e)/(8*(4-2*e)*(3-2*e))
 return e*I,e*(graph-1/(144*e**2))
for seed in (11,22,33):
    rows=[run(e,seed,19) for e in (.01,.005,.0025)]
    c0,c1=(8*np.array(rows[2])-6*np.array(rows[1])+np.array(rows[0]))/3
    assert abs(c0+2)<.0001
    assert abs(c1+.0247)<.0004
    print(seed,"projective residue",c0,"simple pole",c1)
```

Next, use the two independent parallel fermion pairs in FFff(A) to extend the Gaussian spin weighting to two endpoint differences. Isolate the small-pair overlaps before integrating the four-vertex spin trace in FFff(B).

- [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177).
