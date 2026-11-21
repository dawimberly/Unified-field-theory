# FFZ QCc endpoint correction

This revises the QCc doubled graph in [the scalar endpoint checkpoint](ffz_scalar_endpoint_residues.md). The same gauge fixes the single edge at 1. Define p,q as pair totals, r=g(1−g) for the q pair, W=1+p+q and S=1+pf(1−f)+qr. Work in d=4−2ε and the original fixed-1/6 convention.

## Audit and correction

1. The homogeneous evanescent term in A is −d(4−d)W/3 in this gauge. It scales with the sum of the unnormalized lengths.
2. At p=0 the nonroot covariance is C=2qr/(1+qr). The vertex term is −d(6−d)C/2. The mass-insertion term is B=(2d/3)[WC+1+q²(1−2r)]. This follows by summing the insertion on the edge of length 1 and the pair edges with separate lengths qg,q(1−g).
3. The face numerator after radial integration is −d(6−d)C/2+dε[−2W/3+2(WC+1+q²(1−2r))/W]. The measure is q^(ε−1)W^(−3ε)S^(−2+ε) dq dg.
4. When q→∞ and one pair fraction tends to zero with u=qr fixed, the evanescent term grows as 4dεq/3. This region makes the ε expansion of the face nonuniform. It was omitted in the previous double-pole calculation.
5. Subtract the model (4dε/3)q^(1−3ε)(1+qr)^(−2+ε) from the face before expanding, using the measure q^(ε−1)dq dg. The model is integrable at q=0 near ε=0.
6. Integrating the model gives (4dε/3) Beta(1−2ε,1+ε) Beta(2ε,2ε). The second beta function includes both fraction endpoints. The finite limit is 16/3. The once-subtracted regular face has limit −4, so the full face has limit +4/3.
7. The small-p, small-q corner remains 2dε/3, with projective pole 8/(3ε). Together with the two corrected faces this gives a projective residue +16/3.
8. The mixed remainder cannot supply an additional projective pole from the long-pair region after also subtracting its p=0 limit. At fixed p the tail residue depends smoothly on p, with the difference proportional to p near zero. With both pair totals large and both effective lengths fixed, the measure scales as p^(ε−2)q^(ε−2) times a numerator at most linear in p+q and W^(−3ε). Any single-long tail contributes at most a finite term after multiplying its evanescent coefficient. Regions with large effective length are suppressed by S.
9. With Γ(3ε) and the Wick weight +1/4 the corrected bare QCc doubled pole is **+h²z/(108Kε²)**. It has no triple pole. Its log-dependent simple pole is −h²z L/(36Kε), L=ln(z/μ²), from z^(−3ε). This leaves its constant simple-pole coefficient open. CCcc has one more inverse W from its radial factor. The new endpoint mechanism therefore does not furnish a CCcc face residue.
10. Combining this correction with QQ, QCc tadpole and CCcc gives the known part of the bare scalar sector, in units h²z/K:
    5/(216ε²) + [c_QCc − 5/144 − (5/72)L]/ε.
    Here c_QCc is the outstanding constant for the doubled QCc graph. The scalar sum alone does not change the two Newton moment constraints and cannot settle the projected baseline pole.

## Runnable checks

The following Python check needs NumPy and SciPy. The nonzero-ε face calculation subtracts both the small-q corner and the new tail. Convergence to 4/3 checks the analytical endpoint calculation. Finite ε samples are not used to infer the remaining QCc simple-pole constant.

```python
import numpy as np
from scipy.special import roots_legendre
def face(e,n=120):
    d=4-2*e
    x,w=roots_legendre(n); x=(1+x)/2; w=w/2
    u=(x/(1-x))[:,None]
    r=(x*(1-x))[None,:]
    q=u/r
    S=1+u
    W=1+q
    C=2*u/S
    V=(-d*(6-d)/2*C + d*e*(-2*W/3+2*(W*C+1+q*q*(1-2*r))/W))*S**(-2+e)*W**(-3*e)
    # Subtract both the small-q corner and the large-q/fraction corner.
    F=V-(q<1)*2*d*e/3
    tail=4*d*e/3*q**(1-3*e)*S**(-2+e)
    # Integrable at r=0 after subtraction. The model is also integrable at u=0.
    num=(q**e*(F-tail)/u/(1-x[:,None])**2)
    from scipy.special import beta
    analytic=4*d*e/3*beta(1-2*e,1+e)*beta(2*e,2*e)
    return w@num@w+analytic
for e in (.015,.0075,.00375,.001875):
    print(e,face(e,220),face(e,400))
from fractions import Fraction as F
assert abs(face(.001875,400)-4/3)<.04
corner, face0 = F(8,3), -4+F(16,3)
corrected = (corner+2*face0)/3/4/48
assert corrected == F(1,108)
assert corrected+F(1,72)==F(5,216)
assert -3*(corrected+F(1,72))==-F(5,72)
assert F(1,16)-F(1,36)-F(5,72)==-F(5,144)
```

Next, subtract the long-pair tail at arbitrary spectator p and its overlap with p=0, then evaluate the finite projective term including the ε derivative of the corrected faces. Use the result in the projected sum with the Yukawa graphs and counterterms. The baseline simple pole has not yet been compared with 1/24.

- [Heinrich et al., Expansion by regions with pySecDec](https://arxiv.org/abs/2108.10807).
- [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177).

### Integrated tail with the spectator pair retained

The q-tail model at arbitrary p is obtained by replacing 1+qr with 1+pf(1−f)+qr. After integrating q and its fraction, the residual dependence is [1+pf(1−f)]^(−1−ε). Its integral against p^(ε−1)dp df is Beta(ε,1) Beta(1−ε,1−ε). The sum of the two symmetric tail models therefore has closed form

(8dε/3) Beta(1−2ε,1+ε) Beta(2ε,2ε) Beta(ε,1) Beta(1−ε,1−ε)
= 32/(3ε)+80/3+O(ε).

This gives the new tail contribution to both projective coefficients and provides a reusable global subtraction. The finite projective remainder after removing these two models still requires integration.
