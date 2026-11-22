# FFff(A) analytic double pole and remaining discrepancy

The regulated FFff(A) extrapolation is now checked analytically through the subleading projective pole. Use physical pair totals q,r with unit cutoffs temporarily and the same two global long-pair subtraction models.

- The sum of the two tail models is 2/ε²+2/(3ε)+O(1). The script derives the two coefficients from the Dirichlet integrals with Gamma recurrences.
- The angular-average mixed short-short Taylor term is −14ε/3+O(ε²) from the finite-d corner check. It gives −14/(3ε).
- For the once-soft linear terms, take d=4 only after the two tail models are subtracted. In the power-divergent spectator variable, integrate (F−F(0))/r² on the whole positive axis. This is the cutoff-independent finite distribution when the spectator linear coefficient vanishes at ε=0. The unit-cutoff add-back is −F(0), equivalent to subtracting that constant also over r>1. With g and f the respective split-fraction products, the two once-soft integrals are −2(2g+1) and 2(24f+1)/3. Their averaged sum is 2/3. The once-soft constants carry no endpoint logarithmic pole from their own axes. The remaining mixed-subtracted terms have no projective pole under these subtractions.
- Together this gives 2/ε²−10/(3ε) in A projective units. Restoring the corrected full trace prefactor −1/2 and radial factor gives bare triple 1/144 and double 5/288−L/48.

Independently, symmetric nonzero-ε runs of the earlier A code with 2²¹ points gave centered projective slopes −3.35281, −3.33820, −3.33455 for positive ε=.01,.005,.0025 and their negatives. A second scramble and cutoff give −3.35275, −3.33819, −3.33455. This agrees with the analytic −10/3 limit.

The full B extraction was rerun with the narrower symmetric grid ε=0,±.00125,±.0025,±.005, 4096 points per sector, seed 73 and cutoff .9. The fit for ε²I gives (−47.99990,224.06719,190.26085) at powers (0,1,2), again supporting the B double near −1/96. It is stable against the change of grid and disagrees with the budget target −1/288.

## Current pole consistency limit

A's numerical extrapolation is no longer a credible source of the gap. With the stated kernels and all counterterms listed in COUNTERTERM_BUDGET.md, the projected triple cancels and the double residual is approximately −1/144. The same bookkeeping implies an uncanceled +L/48 in the simple pole, since each three-loop bare graph has the common exp(−3εL). This prevents interpreting the simple-pole constant as physical matching. A finite change of the Newton coefficient cannot repair this check.

An independent line identity supplies a route to an audit of the lower-loop/finite-d matching and an alternate derivation of the B curvature insertion. In midpoint and relative coordinates write the scalar Weyl correction as δG=κ(2zQ−c G₂), c=d(4−d)/3, where G₂ is the convolution of two scalar kernels. The kinetic spinor vector obeys δb=−2κz∂_r Q, hence b_curved=−∂_r G_curved−κc∂_r G₂ to linear order. A parallel mixed scalar–fermion pair can therefore be replaced by a relative derivative of G_curved²/2 plus an explicitly evanescent term. Relative derivatives are at fixed midpoint. This suggests reducing the mixed pairs by integration by parts with the vertex profiles left in place, as a check independent of the Gaussian Wick implementation. The reduction has not yet been carried out for the full graph.

## Runnable symbolic check

```python
import sympy as s
p,q,f,g,x=s.symbols('p q f g x',positive=True)
S=1+p+q
K=(4*(1+p)*(1+q)+2*p*q)/S**4
# q and r are physical lengths below, p,q above are effective for K.
r,u=s.symbols('r u',positive=True)
Fq=12*f*r*r*(2*g**3*r*r-g*g*r*r-8*g*g*r-8*g*g+g*r-2*g+2)/(g*r+1)**5
Fr=4*(18*f**3*g*u**4-9*f*f*g*u**4-90*f*f*g*u**3-96*f*f*g*u*u-32*f*f*u*u+9*f*g*u**3+36*f*g*u*u+60*f*g*u-40*f*u+18*g*u*u+12*g-8)/(3*(f*u+1)**5)
Fq=s.cancel(Fq+2*r*r*f*s.diff(K,p).subs({p:0,q:g*r}))
Fr=s.cancel(Fr+2*u*u*g*s.diff(K,q).subs({p:f*u,q:0}))
# Average the short-pair fraction.
Fq=Fq.subs(f,s.Rational(1,6))
Fr=Fr.subs(g,s.Rational(1,6))
print('qface',s.factor(Fq))
print('rface',s.factor(Fr))
# Power-endpoint continuation gives ∫(F-F0) dr/r², since a
# constant integrated over all scales is zero in this prescription.
lhs=s.factor(s.integrate(s.cancel((Fq-Fq.subs(r,0))/r**2),(r,0,s.oo)))
rhs=s.factor(s.integrate(s.cancel((Fr-Fr.subs(u,0))/u**2),(u,0,s.oo)))
print(lhs,rhs,s.simplify(lhs+rhs.subs(f,g)))
avg=s.simplify(s.integrate((lhs+rhs.subs(f,g)).subs(g,x*(1-x)),(x,0,1)))
print('face sum',avg)
ep=s.symbols('ep',positive=True); d=4-2*ep;k=3*ep
G=s.gamma
def D(i,j):
 a=ep-1+i;b=1-2*ep+j;c=4-ep
 ga=G(1+ep)/(ep-1) if i==0 else G(1+ep)
 return ga*G(b)*G(c-a-b)/G(c)
tail2=(k-1)*2*(d-1)/3*d/4*G(2-ep)**2/G(4-2*ep)*G(1+2*ep)**2/G(1+4*ep)*(d*(D(0,0)+D(1,0)+D(0,1)+D(1,1))+2*D(1,1))
print('tail leading and linear',s.simplify(tail2.subs(ep,0)),s.simplify(s.diff(tail2,ep).subs(ep,0)))

```

- [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177).
