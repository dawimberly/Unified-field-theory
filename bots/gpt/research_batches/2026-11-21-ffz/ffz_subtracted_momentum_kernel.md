# FFZ momentum-subtracted two-loop even response

Continues [the Y-improvement Newton residue](ffz_three_loop_Y_improvement_checkpoint.md), [the q⁴ check](ffz_newton_projection_q4_checkpoint.md), and [the finite momentum slope](ffz_two_loop_momentum_slope.md). Earlier fixed-1/6 matter beta and Newton moment results carry forward. The baseline curvature residue still needs a curved vacuum graph sum.

## Subtracted form factor at finite momentum

For each group let u=q²/z with Euclidean q, L=ln(z/μ²), v=tw_g and K=16π². The order-h contribution to the even curvature-squared response, subtracting its value at q=0, is h v²/(72K) times H(u,L). Let r=u x(1−x), g=ln(1+r) and all integrals below be over 0≤x≤1:
A=∫g, B=∫r/(1+r), C=∫g², D=∫g r/(1+r).

H = J + (1−ln 2)(B−A) + A−D+C/2+A² + L(3A−B).

Here J is a finite two-dimensional ordered-sector integral with internal location averages, implemented below. Use the six assignments (a,b,c) of (1,x,xt) and V=a+b+c. Denote the three effective resistance kernels from the momentum-slope checkpoint by R_same,R_pair,R_end. For each source-location integral put s=u R/V. The two sunset-resolvent terms have insertions −s/[V(1+s)] with weights −2 and +1. The explicit endpoint-source term has −ln(1+s). Use the same-line measure dl(a−l), the distinct-line measure dl dr, and the endpoint measure dl. Multiply by (D_graph/x)^(−2), subtract (B−A)/(1+t)² in each sector with a=1, divide by x, then integrate over x,t.

## Derivation and checks

Take the q=0 difference before integrating the common Schwinger scale. The resolvent terms have radial Γ(1+2ε) and the endpoint has Γ(2ε). The source-location integral controlling the corner a=1 is the integral of 1−(1+r)^(−1−2ε)+[(1+r)^(−2ε)−1]/(2ε). Its expansion is B−A+ε[2(A−D)+C]. The angular corner factor after summing the two sectors is 1+ε(1−ln2). Hence the sunset has pole (B−A)/ε and finite J+(1−ln2−2L)(B−A)+2(A−D)+C.

The bubble-square contribution has pole −2A/ε and finite 4LA+C+A². The matter counterterm is minus (3+z∂_z) times the subtracted dimensionally regulated bubble divided by ε, with q fixed. It contributes (3A−B)/ε and finite −L(3A−B)−3C/2−A+D. All poles cancel for every u in the regular region. This gives the displayed finite H, with no fit to individual momentum points. An independent check of its log follows from applying the one-loop matter flow to the one-loop response v² A/(36K). The numerator variation contributes 3h A, the mass flow contributes −h B. The explicit derivative ∂lnμ H = −2(3A−B) cancels that flow in the stated normalization.

At L=0 the numerical response yields the prior q² and q⁴ coefficients by symmetric finite differences. It also passes quadrature convergence checks from u=−3.8 to u=10. This representation gives the whole spacelike response and subthreshold timelike response. It assumes a common pure-metric subtraction and no extra assignment-dependent higher-derivative matching. The subtraction removes the common local R² constant.

## Threshold domain

Every pair of marked points in a theta graph belongs to a cycle of total length at most V. The effective resistance with all edges present is bounded by that of the cycle alone, which is the parallel combination of the two arcs and is at most V/4. The same bound holds for the bubble. The radial denominators and logarithms have no singularity for real u>−4 or in the complex Taylor disk |u|<4. Subtractions preserve analyticity on compact subsets away from the cut. The bound is approached when one line takes nearly the whole Schwinger length and the two insertions are nearly opposite on it. At free order the cut begins exactly at u=−4, since the bubble logarithm has a discontinuity proportional to the interval length sqrt(1−4/|u|). The order-h formula isolates the integrable sectors below that cut. Continuing the sum above the cut with the proper prescription, and resolving mass-shell shifts or singular perturbation terms at threshold, remains open.

## Curved Newton implementation step

At ξ=1/6 in the fixed-dimensional continuation, the quadratic Weyl probe may be implemented on a scalar line with two Schwinger lengths s₁,s₂. The vertex insertion samples an internal coordinate whose conditional mean is (s₂x+s₁y)/(s₁+s₂) and whose variance per component is 2s₁s₂/(s₁+s₂). The average of the coordinate squared is therefore the mean squared plus 2d s₁s₂/(s₁+s₂). The external interaction vertices can carry the same quadratic coordinate insertions from their Weyl weights. For the remaining graph integral the origin fixes one vertex and the Schwinger edges define a reduced conductance Laplacian. The covariance of each Cartesian coordinate vector is twice its inverse. This gives a graph prescription for scalar propagator and vertex-position insertions at any loop order. The fermion propagator still needs endpoint derivatives and spinor numerators, and the full sum requires matter counterterms. The published flat vacuum result fixes the Y part only. At the projected point the two higher poles must cancel, and the baseline simple pole must equal 1/24 for the original Newton family to survive.

## Runnable check

Python 3 with SciPy and NumPy. Set OPENBLAS_NUM_THREADS=1 if needed to keep the nodal sums inexpensive.

```python
from itertools import permutations
from math import log
import numpy as np
from scipy.special import roots_legendre

def rule(n):
    nodes, weights = roots_legendre(n)
    return (nodes+1)/2, weights/2

def response(u,L=0,n=32,m=32):
    loc, w = rule(m)
    edge, wx = rule(n)
    x, t = edge[:,None], edge[None,:]
    one = np.ones((n,n))
    lengths = (one,x*one,x*t)
    r = u*loc*(1-loc)
    g = np.log1p(r)
    A,B,C,D = [w@a for a in (g,r/(1+r),g*g,g*r/(1+r))]
    value = np.zeros((n,n))
    for idx in permutations(range(3)):
        a,b,c = (lengths[i] for i in idx)
        v = a+b+c
        det = a*b+a*c+b*c
        l = a[...,None]*loc
        rr = l-l*l*(b+c)[...,None]/det[...,None]
        s = u*rr/v[...,None]
        same = a*np.sum(w*(a[...,None]-l)*(-s/(1+s)),axis=-1)/v
        end = -a*np.sum(w*np.log1p(s),axis=-1)
        left,right = l[..., :,None],(b[...,None]*loc)[...,None,:]
        rp=left+right-((b+c)[...,None,None]*left**2+(a+c)[...,None,None]*right**2
                       +2*c[...,None,None]*left*right)/det[...,None,None]
        s = u*rp/v[...,None,None]
        pair=a*b*np.sum(w[:,None]*w*(-s/(1+s)),axis=(-1,-2))/v
        kernel = (-2*same+pair+end)/(det/x)**2
        if idx[0]==0:
            kernel -= (B-A)/(1+t)**2
        value += kernel/x
    J = wx@value@wx
    return J+(1-log(2))*(B-A)+A-D+C/2+A*A+L*(3*A-B)

p1=-.3073178731149105/3
p2=1.4993254559801423/30
d=.001
plus,minus=response(d),response(-d)
assert abs((plus-minus)/(2*d)-p1)<1e-7
assert abs((plus+minus)/(2*d*d)-p2)<1e-7
for u in (-3.8,-2,1,4,10):
    coarse=response(u,n=42,m=48)
    fine=response(u,n=56,m=64)
    assert abs(coarse-fine)<2e-7,(u,coarse,fine)
    print(u,fine)
```

## Sources

- [Schmidt and Schubert, Worldline Green functions for multiloop diagrams](https://doi.org/10.1016/0370-2693(94)90944-X).
- [Bunch and Parker, Feynman propagator in curved spacetime](https://doi.org/10.1103/PhysRevD.20.2499).
- [Martin, Three-loop effective potential for softly broken supersymmetry](https://doi.org/10.1103/PhysRevD.109.015019).
