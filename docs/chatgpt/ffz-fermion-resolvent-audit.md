# FFZ fermion curvature insertion audit

**Conversation date range** 2026-11-22

**DERIVED** This checkpoint records the FFZ fermion curvature insertion audit. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Continuation of [the locality audit](ffz-double-pole-locality-audit.md). This checks the spinor insertion used in [the curved projective kernels](ffz-ccff-projective-kernel.md) against the [two-loop radial reduction](ffz-two-loop-flat-space-checkpoint.md).

## Off-diagonal resolvent

**DERIVED** In the flat Weyl frame δD=−κm |u|². The first derivative of D^−1 is −S δD S. Write each flat spinor resolvent as (m−γ·∂)G, then transfer the derivative of the second factor to its other endpoint. The odd Clifford component is κz γ·(∂y−∂x)Q, where Q=∫du G(x−u)|u|²G(u−y). The even Clifford part vanishes between the two opposite chiral projections on every kinetic line.

**DERIVED** At fixed total Schwinger length s, the convolution in Q has extra factor s and its conditional position moment is μ²+2ds a(1−a), with μ=ax+(1−a)y. Differentiating before anchoring gives two terms, from the heat kernel and from μ². With Δ=x−y their combined numerator at fixed a is
**DERIVED** Δ[μ²+2ds a(1−a)] + 2s(1−2a)μ.
**DERIVED** The a integral is Δ[x²+y²+x·y+(d−1)s]/3. Thus the insertion multiplies the flat kinetic component by 2κzs[x²+y²+x·y+(d−1)s]/3. In particular the split-point derivative lowers d to d−1 in the variance term. This agrees with the formula used in both four-Yukawa kernels.

## Independent anchored radial check

**DERIVED** Set y=0 only after the differentiation. The flat radial kinetic factor at fixed s is b_s=r K_s/(2s). Using ∫s b_s=rG/2 and ∫s² b_s=−r ∂z G/2 yields δb/κ=(zr/3)[r²G−(d−1)∂zG]. The scaling relation 2z∂zG=(d−2)G+rG′ reduces it to
**DERIVED** zr³G/3+(d−1)r²b/6−(d−1)(d−2)rG/6.
**DERIVED** The older radial expression reduces to this same polynomial using b′=−zG−(d−1)b/r and ∂z b=−rG/2. This tests the relative sign and the dimension-dependent variance against the differently organized lower-loop derivation.

## Separate test for the scalar evanescent insertion

**DERIVED** At linear order in R, the difference between the two prescriptions ξ=1/6 and ξ=ξ_d is exactly
**DERIVED** R(1/6−ξ_d)(∂_x+∂_u)V_flat
**DERIVED** with independent scalar squared masses x,u. The derivatives keep the fermion mass and the cubic and Yukawa vertices fixed. The coefficient 1/6−ξ_d=(4−d)/[12(d−1)] multiplies the flat mass derivative, so only its leading triple pole is required for a three-loop double-pole comparison. This provides a flat-space diagnostic of the direct scalar curvature-mass terms, separate from the spinor line test. That comparison has yet to be evaluated.

## Exact checks

```python
from fractions import Fraction as F
from itertools import product

def integrate(poly):
    return sum(F(c,k+1) for k,c in enumerate(poly))

for x,y,s,d in product((-3,1,7),(-1,0,4),(1,3),(2,3,4,7)):
    delta=x-y
    h=(y*y,2*y*delta+2*d*s,delta*delta-2*d*s)
    term=(delta*h[0]+2*s*y,
          delta*h[1]+2*s*(delta-2*y),
          delta*h[2]-4*s*delta)
    assert integrate(term)==F(delta,3)*(x*x+y*y+x*y+(d-1)*s)

for d,r,z,g,b in product((2,3,4,7),(F(1,2),F(3)),(2,7),(F(3,5),),(F(4,9),)):
    dzg=((d-2)*g-r*b)/(2*z)
    dzb=-r*g/2
    bp=-z*g-(d-1)*b/r
    old=-(d-1)*r*r*b/6-r**3*bp/3+d*(d-1)*dzb/3+(d-1)*r*g/3
    new=z*r*(r*r*g-(d-1)*dzg)/3
    assert old==new
    assert 4*d*(d-1)*(F(1,6)-F(d-2,4*(d-1)))==F(d*(4-d),3)
print("Off-diagonal spinor insertion and anchored radial check passed")
```

- **DERIVED** [Bunch and Parker, Feynman propagator in curved spacetime](https://doi.org/10.1103/PhysRevD.20.2499).
