# FFZ Newton pole projection and fourth-order momentum coefficient

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**ASSUMED** This checkpoint records the FFZ Newton pole projection and fourth-order momentum coefficient. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Continue [the two-loop momentum slope](ffz-two-loop-momentum-slope.md) and [the fixed-1/6 conversion](ffz-fixed-one-sixth-checkpoint.md). The earlier analytic improvement conversion gave zero additive drift through cubic order in the stated source subtraction. The present analysis treats the pure-metric Newton pole and the next term in the even response separately.

## A single three-loop Newton residue controls the original family

**DERIVED** After matter subdivergence subtraction the order-$h^{2}$ R simple pole, in units $h^{2}/(K \epsilon )$, has group form $z(c_{0}+c_X a_X+c_Y a_Y)$. One curvature insertion makes the dependence affine. Independent groups and mass-independent local subtraction give common coefficients for common h. There is one dimension-two mass invariant z in each group, and the simple pole cannot contain ln z after all subdivergences. The moment constraints imply $\sum z a_X=0, \sum z a_Y=-\sum z$. The residual pole is therefore $7M^{2}(c_{0}-c_Y)$, independent of t. In particular the $t=0$ member, for which $a_X=0$ for every group, suffices. For a single-group computation choose $a_X=0, a_Y=-1 (w=0)$. Its pole directly gives c₀−c_Y. This test bypasses the physical source-split scan, while its nonzero result would apply to all t in the original equal-Yukawa family.

**DERIVED** One-loop F_R has pole $+zw/(12K\epsilon )$ in the parent's Euclidean convention. Thus if $D=c_{0}-c_Y$ is nonzero, a boundary change $\delta w_g=h^{2}r_g$ can cancel the aggregate pole at a chosen scale while maintaining $\sum \delta w=0$. An example is $r=(28D,0,-28D)$, giving $\sum z r=-84M^{2}D$. Putting the correction entirely in Y leaves both X moments null. The old mass-weighted moment itself then equals $-84h^{2}M^{2}D$. This is a change of source assignments compared with the original family. Cancellation through this order at the reference scale does not establish all-order invariance of the corrected constraints. A distinct finite change of the renormalized pure-metric coefficient by $\kappa h\sum z/K$ shifts its beta by $4\kappa h^{2}\sum z/K$ at this order, since $\beta _h=3h^{2}$ and $\beta _z=hz$. A gravitational matching condition must transform with that change. Neither manipulation computes D in the fixed prescription. The additive matter beta from the previous run is likewise insufficient to fix D.

## Finite two-loop $q^{4}$ coefficient of the even response

**DERIVED** Let q be Euclidean and $v=tw_g$, $z=m_g^{2}, L=\ln\left(z/\mu ^{2}\right), K=16\pi ^{2}$. The one-loop coefficient of $(q^{2})^{2}$ is $-v^{2}/(2160Kz^{2})$, by expanding the Codello–Zanusso potential form factor. The new order-h matter coefficient is

**DERIVED** $-h v^{2}/(2160K z^{2}) [L/2 - 1.49932545598]$.

**DERIVED** At $\mu =M$ the three-group contribution is $+0.00400063873331 h t^{2}/(KM^{4})$. The local finite $R\Box ^{2}R$ coefficient, if independently specified, contributes at this derivative order.

**DERIVED** The prior finite-momentum sunset construction uses effective resistances R_same,R_pair,R_end. Square each resistance and use the same location measures to obtain S,Q,E. The factor $1/2$ from the exponential Taylor term is included in the overall coefficient. After integrating the radial Schwinger scale the normalized sunset integrand is proportional to $(2+2\epsilon)(2S-Q)/\mathrm{scale}-E$. At the corner with only line a long, $S=1/60, Q=0, E=1/30$, so its subtracted endpoint is $(1+2\epsilon )(1+t)^{-2+\epsilon }/30$ in the x-sector parameterization below. The other endpoints vanish. The corresponding pole in units h times the one-loop $q^{4}$ coefficient is $+1/(2\epsilon )$. The figure-eight term supplies $-1/\epsilon$, and the matter counterterm $+1/(2\epsilon )$.

**DERIVED** For a separate finite subtraction check, the exact endpoint supplies finite $(5-\ln 2)/2$ at $L=0$. The two-bubble square contributes $-11/6$, including the product of the two $q^{2}$ pieces. The matter counterterm has zero finite part because $\Gamma (2+\epsilon )e^{\gamma _E \epsilon }=1+\epsilon +$… cancels the mass-insertion term $-1/2$. Adding the finite projective integral agrees with dimensional extrapolation. The RG gives an independent log check. At the reference family point the two source variations have equal squared weights and only X runs, so the numerator runs at 3h relative to itself. The denominator $z^{2}$ runs at 2h, leaving a relative h. The explicit logarithm therefore has coefficient $+1/2$. A t-independent improvement translation leaves both variations unchanged.

## Runnable check

**DERIVED** Python 3 with SciPy and NumPy. The three-node Gauss rule integrates the source-location polynomials exactly. The outer quadrature and the dimensional extrapolation are numerical.

```python
from itertools import permutations
from math import exp, log
import numpy as np
from scipy.integrate import dblquad, quad
from scipy.special import gamma
eg=.5772156649015329
nodes,weights=np.polynomial.legendre.leggauss(3)
nodes,weights=(nodes+1)/2,weights/2
def kernels(a,b,c,d):
    l=a*nodes
    ra=l-l*l*(b+c)/d
    r=b*nodes
    pair=l[:,None]+r-((b+c)*l[:,None]**2+(a+c)*r*r+2*c*l[:,None]*r)/d
    return a*np.dot(weights,(a-l)*ra**2), a*b*(weights@ (pair**2) @weights), a*np.dot(weights,ra**2)
def sunset(e,L,projective_only=False):
    def fun(t,x):
        lengths=(1,x,x*t)
        total=0.
        for a,b,c in permutations(lengths):
            d=a*b+a*c+b*c
            v=a+b+c
            S,Q,E=kernels(a,b,c,d)
            term=(d/x)**(-2+e)*((2+2*e)*(2*S-Q)/v-E)/v**(2+2*e)
            if a==1:
                term-=(1+2*e)/30*(1+t)**(-2+e)
            total+=term
        return x**(e-1)*total
    regular=dblquad(fun,0,1,0,1,epsabs=3e-9)[0]
    if projective_only:
        return regular
    endpoint=2*(1+2*e)/30*quad(lambda t:(1+t)**(-2+e),0,1)[0]/e
    return 15*exp(2*e*(eg-L))*gamma(2+2*e)*(endpoint+regular)
def figure(e,L):
    return -30*exp(2*e*(eg-L))*(gamma(e)*gamma(2+e)/30+gamma(1+e)**2/36)
def ct(e,L):
    return exp(e*(eg-L))*gamma(2+e)*(.5/e-.5)
def extrapolate(f):
    vals=[f(.016/2**j) for j in range(5)]
    for order in range(1,5):
        vals=[(2**order*b-a)/(2**order-1) for a,b in zip(vals,vals[1:])]
    return vals[0]
c=extrapolate(lambda e:sunset(e,0)+figure(e,0)+ct(e,0))
print(c)
for L in (0,log(2),log(4)):
    result=extrapolate(lambda e:sunset(e,L)+figure(e,L)+ct(e,L))
    print(L,result)
    assert abs(result-c-L/2)<1e-7
print('aggregate q4 coefficient',-sum(w*w/z**2*(c+log(z)/2) for w,z in zip((2,-3,1),(1,2,4)))/2160)
assert abs(c-(15*sunset(0,0,True)+(5-log(2))/2-11/6)) < 1e-7
for split in (0,1,-3):
    weights_ = (2,-3,1)
    z_ = (1,2,4)
    ax = [split*w for w in weights_]
    ay = [(w-2)/2-a for w,a in zip(weights_,ax)]
    assert sum(m*a for m,a in zip(z_,ax)) == 0
    assert sum(m*a for m,a in zip(z_,ay)) == -7
assert sum((28,0,-28)) == 0
assert sum(m*r for m,r in zip((1,2,4),(28,0,-28))) == -84
```

## Next calculation

**DERIVED** Compute a single order-$h^{2}$ pure-metric R residue at $a_X=0,a_Y=-1$ using the original d-dimensional Weyl metric continuation. Track the curvature-free conformal part and one Y mass-source insertion, including the two-loop matter counterterms. This projected point suffices to determine whether the original family needs an interacting shift or pure-metric input. For the momentum branch, compare the two-loop threshold discontinuity with the two-particle cut or compute the whole momentum-subtracted kernel with the known $q^{2}$ and $q^{4}$ limits.

## Primary papers

- **CITED** [Frolov, Fursaev and Zelnikov, Statistical Origin of Black Hole Entropy in Induced Gravity](https://arxiv.org/abs/hep-th/9607104). — eq./page UNVERIFIED.
- **CITED** [Codello and Zanusso, On the non-local heat kernel expansion](https://arxiv.org/abs/1203.2034), Eqs. (4.10)–(4.12).
- **CITED** [Barvinsky, Gusev, Vilkovisky and Zhytnikov, Asymptotic behaviours of the heat kernel in covariant perturbation theory](https://doi.org/10.1063/1.530428). — eq./page UNVERIFIED.
