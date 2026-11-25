# FFZ CcFf long scalar-pair subtraction

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**ASSUMED** This checkpoint records the FFZ CcFf long scalar-pair subtraction. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues [the short-pair checkpoint](ffz-ccff-short-pair.md). Use the same p,q,P,Q,S,W notation and $d=4-2\epsilon$. The overall bare graph weight is $-z|y|^{4}/2$.

**DERIVED** At $p\to \infty$ with $P=pf(1-f)$ fixed, the evanescent mass-profile term supplies a large-p numerator proportional to p. The fermion spin weight after removing $1/q$ is $d(1+P)/(2S)$. The resulting global model in the same projective measure is

**DERIVED** $H_tail = (2d^2\frac{\epsilon }{3}) (1+P) S^(-1-\frac{d}{2}) p^(1-3\epsilon ).$

**DERIVED** The linear coefficient in the bracket follows from $-d(4-d)W/3 + 3\epsilon B/W$, whose leading large-p value is $(4d\epsilon /3)p$. For the competing long fermion pair, the measure is $q^{\epsilon -2}$ and an endpoint split fraction contributes a further $1/q$ at fixed Q. The corresponding linear growth is therefore integrable at $q=\infty$. If both effective lengths grow, the determinant S also suppresses the integrand.

**DERIVED** Change p,q to P,Q in the tail model. The fraction integrations give $\operatorname{B} (2\epsilon ,2\epsilon )$ for the scalar pair and $\operatorname{B} (2-\epsilon ,2-\epsilon )$ for the fermion pair. Integrate Q first, using analytic continuation of the power-divergent short-q constant. The Q integral gives $\operatorname{B} (\epsilon -1,4-2\epsilon ) (1+P)^{-4+2\epsilon }$. The remaining P integral gives $\operatorname{B} (1-2\epsilon ,2)$. Thus the full model is

**DERIVED** $(2d^{2}\epsilon /3) \operatorname{B} (2\epsilon ,2\epsilon ) \operatorname{B} (2-\epsilon ,2-\epsilon ) \operatorname{B} (\epsilon -1,4-2\epsilon ) \operatorname{B} (1-2\epsilon ,2)$
**DERIVED** $= -8/(3\epsilon )-20/3+O(\epsilon )$.

**DERIVED** In particular the Q endpoint must be continued before interpreting the large-p subtraction. Treating it as a positive ordinary integral would give the wrong residue.

**DERIVED** The earlier q-linear geometric face gives $+2/\epsilon$ and the small-p, small-q evanescent overlap gives $-4/(3\epsilon )$. The $p=0$ face is evanescent and its regular part supplies no pole apart from the already included short-q linear overlap. Subtracting the large-p model removes the linear growth at fixed P. In the short-q linear remainder the model's homogeneous large-p part is removed as well. The resulting candidate projective residue is −2. It implies bare CcFf double-pole coefficient $+1/144$ in units $h^{2}z/K$ and the scale-dependent part $-L/48$ of the simple coefficient. This assembly is analytical and still needs a full regulated-integrand check before being used as a final coefficient.

## Geometric part of the finite coefficient

**DERIVED** For the cubic-vertex term, keep finite d in the conditioned variance. Omitting $W^{-3\epsilon }$, the numerator including qN is proportional to

**DERIVED** $-\lambda _C d [dP(1+P)(1+Q)+2P^{2}Q] / S^{2+d/2}$.

**DERIVED** The P,Q integrations with the original powers are Dirichlet integrals. The script verifies the residue +2 and gives the finite coefficient +2 for this model. Restoring $W^{-3\epsilon }$ adds a finite short-q face correction from the q derivative of both the geometric term and ln W. That correction and the remaining evanescent finite terms are needed for the graph's simple pole.

## Runnable algebra

**DERIVED** Python 3 and SymPy.

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

**ASSUMED** Next apply the global tail and the two short-q Taylor subtractions at nonzero $\epsilon$, verify the full double pole, and extract the finite projective part. The three-loop Newton cancellation remains open.

- **CITED** [Heinrich et al., Expansion by regions with pySecDec](https://arxiv.org/abs/2108.10807). — eq./page UNVERIFIED.
