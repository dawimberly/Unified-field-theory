# FFZ CcFf long scalar-pair subtraction

Continues [the short-pair checkpoint](ffz_ccff_short_pair.md). Use the same p,q,P,Q,S,W notation and d=4−2ε. The overall bare graph weight is −z|y|⁴/2.

At p→∞ with P=pf(1−f) fixed, the evanescent mass-profile term supplies a large-p numerator proportional to p. The fermion spin weight after removing 1/q is d(1+P)/(2S). The resulting global model in the same projective measure is

H_tail = (2d²ε/3) (1+P) S^(−1−d/2) p^(1−3ε).

The linear coefficient in the bracket follows from −d(4−d)W/3 + 3εB/W, whose leading large-p value is (4dε/3)p. For the competing long fermion pair, the measure is q^(ε−2) and an endpoint split fraction contributes a further 1/q at fixed Q. The corresponding linear growth is therefore integrable at q=∞. If both effective lengths grow, the determinant S also suppresses the integrand.

Change p,q to P,Q in the tail model. The fraction integrations give Beta(2ε,2ε) for the scalar pair and Beta(2−ε,2−ε) for the fermion pair. Integrate Q first, using analytic continuation of the power-divergent short-q constant. The Q integral gives Beta(ε−1,4−2ε) (1+P)^(−4+2ε). The remaining P integral gives Beta(1−2ε,2). Thus the full model is

(2d²ε/3) Beta(2ε,2ε) Beta(2−ε,2−ε) Beta(ε−1,4−2ε) Beta(1−2ε,2)
= −8/(3ε)−20/3+O(ε).

In particular the Q endpoint must be continued before interpreting the large-p subtraction. Treating it as a positive ordinary integral would give the wrong residue.

The earlier q-linear geometric face gives +2/ε and the small-p, small-q evanescent overlap gives −4/(3ε). The p=0 face is evanescent and its regular part supplies no pole apart from the already included short-q linear overlap. Subtracting the large-p model removes the linear growth at fixed P. In the short-q linear remainder the model's homogeneous large-p part is removed as well. The resulting candidate projective residue is −2. It implies bare CcFf double-pole coefficient +1/144 in units h²z/K and the scale-dependent part −L/48 of the simple coefficient. This assembly is analytical and still needs a full regulated-integrand check before being used as a final coefficient.

## Geometric part of the finite coefficient

For the cubic-vertex term, keep finite d in the conditioned variance. Omitting W^(−3ε), the numerator including qN is proportional to

−λ_C d [dP(1+P)(1+Q)+2P²Q] / S^(2+d/2).

The P,Q integrations with the original powers are Dirichlet integrals. The script verifies the residue +2 and gives the finite coefficient +2 for this model. Restoring W^(−3ε) adds a finite short-q face correction from the q derivative of both the geometric term and ln W. That correction and the remaining evanescent finite terms are needed for the graph's simple pole.

## Runnable algebra

Python 3 and SymPy.

```python
import sympy as s
e=s.symbols('e'); d=4-2*e; G=s.gamma
# e times the long-scalar-pair model, with poles shifted by recurrences
tail=-2*d*d/3*G(1+2*e)**2/G(1+4*e)*G(1+e)*G(4-2*e)/((1-e)*G(3-e))*G(2-e)**2/G(4-2*e)/((1-2*e)*(2-2*e))
assert s.simplify(tail.subs(e,0))==-s.Rational(8,3)
print('tail finite',s.simplify(s.diff(tail,e).subs(e,0)))
# geometric model (drop W^(-3e)); each monomial uses a Dirichlet integral.
P,Q=s.symbols('P Q'); lam=(6-d)/2
poly=s.Poly(-lam*d*(d*P*(1+P)*(1+Q)+2*P**2*Q),P,Q)
def gshift(a): # e*Gamma(e+a) for a=-1,0; Gamma otherwise
 if a==-1: return G(1+e)/(e-1)
 if a==0: return G(1+e)
 return G(e+a)
result=0
for (i,j),coef in poly.terms():
 result+=coef*G(e+i)*gshift(j-1)*G(5-3*e-i-j)/G(4-e)
result*=G(1-e)**2/G(2-2*e)*G(2-e)**2/G(4-2*e)
assert s.simplify(result.subs(e,0))==2
print('geom model finite',s.simplify(s.diff(result,e).subs(e,0)))
assert (2-s.Rational(4,3)+tail.subs(e,0))*(-s.Rational(1,2))/3/48==s.Rational(1,144)
assert s.simplify(s.diff(tail,e).subs(e,0))==-s.Rational(20,3)
assert s.simplify(s.diff(result,e).subs(e,0))==2
```

Next apply the global tail and the two short-q Taylor subtractions at nonzero ε, verify the full double pole, and extract the finite projective part. The three-loop Newton cancellation remains open.

- [Heinrich et al., Expansion by regions with pySecDec](https://arxiv.org/abs/2108.10807).
