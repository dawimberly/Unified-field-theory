# FFff(A) long-pair endpoint

Use the gauge and conditioned two-spin kernel in FFFF_A_CORNER_CORRECTION.md. At r→∞ with Q=rg fixed and q finite, B grows as 2(d−1)N r²/3. Including the radial projective factor, the leading full numerator is (k−1)2(d−1)N S^(−d/2) r^(2−k)/3, k=3ε. A is subleading in this limit. This gives a global tail model with the same short-q continuation as the original graph.

Change q,r to P,Q. In this model the angular integrals are Beta(2−ε,2−ε) and Beta(2ε,2ε). The remaining integrand is P^(ε−2)Q^(−2ε) times N S^(−d/2). For the double projective pole the first beta tends to 1/6 and the second has residue 1. In the short-P derivative at d=4,
 
N S^(−2) = [4(1+P)(1+Q)+2PQ]/(1+P+Q)^4,

so ∫dQ ∂_P(N S^(−2)) evaluated at P=0 equals −3. The combined long-r, short-q model residue is therefore (−2)(1/6)(−3)=+1 at projective order 1/ε². The long-q counterpart has the same residue by exchange of the pairs after integrating the scalar lengths. Together these two add-back models contribute +1/(576ε³) in units h²z/K, including Γ(3ε−1), the signed Wick weight and the conversion from κ to R. The earlier zero from the short-short corner does not extend to the long-pair models.

This is a contribution to the Laurent series, with the full subtracted remainder still to be checked. In particular the full FFff(A) pole and the FFff(B) pole remain open. The next check is to show that the sum of the globally subtracted remainder and add-backs is cutoff independent and to extract the remaining double and simple poles.

```python
import sympy as s
P,Q=s.symbols("P Q", nonnegative=True)
N=(4*(1+P)*(1+Q)+2*P*Q)/(1+P+Q)**4
derivative=s.factor(s.diff(N,P).subs(P,0))
assert s.integrate(derivative,(Q,0,s.oo))==-3
one=-2*s.beta(2,2)*(-3)
assert s.expand_func(one)==1
assert s.Rational(-1,384)*s.Rational(-1,3)*2==s.Rational(1,576)
```
