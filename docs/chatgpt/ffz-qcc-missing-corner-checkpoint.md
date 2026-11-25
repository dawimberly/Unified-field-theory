# FFZ QCc endpoint correction

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**ASSUMED** This checkpoint records the FFZ QCc endpoint correction. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** The $1/24$ target is a normalization chosen to match Einstein gravity; it is not an independently predicted gravitational constant. A computed $c_Y=1/24$ in the stated matter prescription is a separate conditional calculation. Matching a baseline simple pole to that target remains a test, not a completed derivation.


**DERIVED** This revises the QCc doubled graph in [the scalar endpoint checkpoint](ffz-scalar-endpoint-residues.md). The same gauge fixes the single edge at 1. Define p,q as pair totals, $r=g(1-g)$ for the q pair, $W=1+p+q$ and $S=1+pf(1-f)+qr$. Work in $d=4-2\epsilon$ and the original fixed $1/6$ convention.

## Audit and correction

1. **DERIVED** The homogeneous evanescent term in A is $-d(4-d)W/3$ in this gauge. It scales with the sum of the unnormalized lengths.
2. **DERIVED** At $p=0$ the nonroot covariance is $C=2qr/(1+qr)$. The vertex term is $-d(6-d)C/2$. The mass-insertion term is $B=(2d/3)[WC+1+q^{2}(1-2r)]$. This follows by summing the insertion on the edge of length 1 and the pair edges with separate lengths qg,q(1−g).
3. **DERIVED** The face numerator after radial integration is $-d(6-d)C/2+d\epsilon [-2W/3+2(WC+1+q^{2}(1-2r))/W]$. The measure is $q^{\epsilon -1}W^{-3\epsilon }S^{-2+\epsilon } dq dg$.
4. **DERIVED** When $q\to \infty$ and one pair fraction tends to zero with $u=qr$ fixed, the evanescent term grows as $4d\epsilon q/3$. This region makes the $\epsilon$ expansion of the face nonuniform. It was omitted in the previous double-pole calculation.
5. **DERIVED** Subtract the model $(4d\epsilon /3)q^{1-3\epsilon }(1+qr)^{-2+\epsilon }$ from the face before expanding, using the measure $q^{\epsilon -1}dq dg$. The model is integrable at $q=0$ near $\epsilon =0$.
6. **DERIVED** Integrating the model gives $(4d\epsilon /3) \operatorname{B} (1-2\epsilon ,1+\epsilon ) \operatorname{B} (2\epsilon ,2\epsilon )$. The second beta function includes both fraction endpoints. The finite limit is $16/3$. The once-subtracted regular face has limit −4, so the full face has limit $+4/3$.
7. **DERIVED** The small-p, small-q corner remains $2d\epsilon /3$, with projective pole $8/(3\epsilon )$. Together with the two corrected faces this gives a projective residue $+16/3$.
8. **DERIVED** The mixed remainder cannot supply an additional projective pole from the long-pair region after also subtracting its $p=0$ limit. At fixed p the tail residue depends smoothly on p, with the difference proportional to p near zero. With both pair totals large and both effective lengths fixed, the measure scales as $p^{\epsilon -2}q^{\epsilon -2}$ times a numerator at most linear in $p+q$ and $W^{-3\epsilon }$. Any single-long tail contributes at most a finite term after multiplying its evanescent coefficient. Regions with large effective length are suppressed by S.
9. **DERIVED** With $\Gamma (3\epsilon )$ and the Wick weight $+1/4$ the corrected bare QCc doubled pole is **$+h^{2}z/(108K\epsilon ^{2})$**. It has no triple pole. Its log-dependent simple pole is $-h^{2}z L/(36K\epsilon ), L=\ln\left(z/\mu ^{2}\right)$, from $z^{-3\epsilon }$. This leaves its constant simple-pole coefficient open. CCcc has one more inverse W from its radial factor. The new endpoint mechanism therefore does not furnish a CCcc face residue.
10. **DERIVED** Combining this correction with QQ, QCc tadpole and CCcc gives the known part of the bare scalar sector, in units $h^{2}z/K$:
$$
\begin{gathered}
5/(216\epsilon ^{2}) + [c_QCc - 5/144 - (5/72)L]/\epsilon . \\
Here c_QCc is the outstanding constant for the doubled QCc graph. The scalar \sum alone does not change the two Newton moment constraints and cannot settle the projected baseline pole.
\end{gathered}
$$

## Runnable checks

**DERIVED** The following Python check needs NumPy and SciPy. The nonzero $\epsilon$ face calculation subtracts both the small-q corner and the new tail. Convergence to $4/3$ checks the analytical endpoint calculation. Finite $\epsilon$ samples are not used to infer the remaining QCc simple-pole constant.

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

**ASSUMED** Next, subtract the long-pair tail at arbitrary spectator p and its overlap with $p=0$, then evaluate the finite projective term including the $\epsilon$ derivative of the corrected faces. Use the result in the projected sum with the Yukawa graphs and counterterms. The baseline simple pole has not yet been compared with $1/24$.

- **CITED** [Heinrich et al., Expansion by regions with pySecDec](https://arxiv.org/abs/2108.10807). — eq./page UNVERIFIED.
- **CITED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177). — eq./page UNVERIFIED.

### Integrated tail with the spectator pair retained

**DERIVED** The q-tail model at arbitrary p is obtained by replacing $1+qr$ with $1+pf(1-f)+qr$. After integrating q and its fraction, the residual dependence is $[1+pf(1-f)]^{-1-\epsilon }$. Its integral against $p^{\epsilon -1}dp df$ is $\operatorname{B} (\epsilon ,1) \operatorname{B} (1-\epsilon ,1-\epsilon )$. The sum of the two symmetric tail models therefore has closed form

**DERIVED** $(8d\epsilon /3) \operatorname{B} (1-2\epsilon ,1+\epsilon ) \operatorname{B} (2\epsilon ,2\epsilon ) \operatorname{B} (\epsilon ,1) \operatorname{B} (1-\epsilon ,1-\epsilon )$
**DERIVED** $= 32/(3\epsilon )+80/3+O(\epsilon )$.

**DERIVED** This gives the new tail contribution to both projective coefficients and provides a reusable global subtraction. The finite projective remainder after removing these two models still requires integration.
