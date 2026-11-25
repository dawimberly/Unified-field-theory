# FFZ finite QCc projective remainder and bare scalar sum

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**ASSUMED** This checkpoint records the FFZ finite QCc projective remainder and bare scalar sum. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** The $1/24$ target is a normalization chosen to match Einstein gravity; it is not an independently predicted gravitational constant. A computed $c_Y=1/24$ in the stated matter prescription is a separate conditional calculation. Matching a baseline simple pole to that target remains a test, not a completed derivation.


**ASSUMED** Continues [the missing-endpoint correction](ffz-qcc-missing-corner-checkpoint.md). Use the same fixed $1/6$ continuation and bare graph normalization. Separate the projective QCc integral into the vertex-position term G and the terms proportional to $\epsilon$, denoted E. Write $W=1+p+q, P=p f(1-f), Q=q g(1-g), S=1+P+Q$. The anchor remains the quartic vertex throughout.

**DERIVED** The sum of the two nonroot coordinate variances per Cartesian dimension is $2(P+Q+2PQ)/S$. Hence G has numerator $-d(6-d)(P+Q+2PQ)/S$. First omit $W^{-3\epsilon }$ from this term. The P,Q integrals are Dirichlet integrals and give the following, including the two split-fraction integrals:

**DERIVED** $G_{0}=-d(6-d) \operatorname{B} (1-\epsilon ,1-\epsilon )^{2} [2\Gamma (1+\epsilon )\Gamma (\epsilon )\Gamma (2-3\epsilon )+2\Gamma (1+\epsilon )^{2}\Gamma (1-3\epsilon )]/\Gamma (3-\epsilon )$.

**DERIVED** $G_{0}=-8/\epsilon -32+O(\epsilon )$. Restoring $W^{-3\epsilon }$ contributes a finite term from each once-subtracted face. On $p=0$ the geometric measure is $-8r dq/(1+qr)^{3}$, with $r=g(1-g)$. Its correction comes from $-3\epsilon \ln\left(1+q\right)$, times the p pole. The two faces give a total

**DERIVED** $G_{1}=48\int dg \int _{0}\infty dq r \ln\left(1+q\right)/(1+qr)^{3}$
  $=24\int dg [-\ln r/(1-r)^{2}-1/(1-r)]$.

**DERIVED** The overlap has zero geometric numerator. The difference after face subtraction has no $1/\epsilon$ pole and vanishes for this correction at $\epsilon =0$.

**DERIVED** For E use the two globally integrated long-pair models from the preceding checkpoint. They contribute $32/(3\epsilon )+80/3$. After subtracting the models, the evanescent integrand contributes only through the small-p and small-q faces and their overlap. A model evaluated on its own short-pair face vanishes. On the opposite face keep the model subtraction for the long pair. With the $q<1$ cutoff for the overlap, the remaining q integration, divided by $d\epsilon$, is

**DERIVED** $\int _{0}\infty dq/q { [-2W/3 + 2(W C+1+q^{2}(1-2r))/W - 4q/3]/(1+qr)^{2} - (4/3)1_(q<1) }$,

**DERIVED** where $W=1+q$ and $C=2qr/(1+qr)$. Its exact partial fraction integral is $2/3 - [4(r+2)/(3(r-1))] \ln r$. The two faces together contribute eight times its g average. The overlap, $2d\epsilon /(3\epsilon ^{2})$, contributes $8/(3\epsilon )-4/3$. The mixed remainder after the two short-face and two long-tail subtractions is multiplied by $\epsilon$ and contributes zero to this Laurent order. Here the short-face difference supplies a power of the short pair parameter and the long-tail remainder loses the linear growth in the long pair. In the joint long limit the spectator pair measure suppresses the overlap. The geometric Dirichlet integral is treated separately, so it does not rely on these evanescent remainder estimates.

**DERIVED** Summing gives the projective coefficients $16/(3\epsilon )$ and $I_{0}=-17.5025158190$. Multiplication by $\Gamma (3\epsilon )$, the $MS-bar$ scale factors, and 1/[16d(d−1)] gives the bare doubled QCc Laurent coefficients, in units $h^{2}z/K$,

**DERIVED** $1/(108\epsilon ^{2}) + [-0.01958384305 - L/36]/\epsilon$.

**DERIVED** The exact constant is defined by the convergent one-dimensional integrals above. It includes the $7\epsilon /6$ correction from 1/[d(d−1)]. The bare sum over all scalar-only graphs is now

**DERIVED** $5/(216\epsilon ^{2}) + [-0.05430606527 - 5L/72]/\epsilon$.

**DERIVED** This sum omits Yukawa graphs and counterterms. The finite QCc result has an independent integration check of its individual face integrals below. A direct Laurent extrapolation of the entire graph would be a stronger cross-check.

## Counterterm check for assembly

**DERIVED** The flat supersymmetric mass renormalization factor for z follows from $\beta _z/z=h-h^{2}+O(h^{3})$. For $z_B=z Z_z$, the relation $(-2\epsilon h+3h^{2})\partial ₕ Z_z+(h-h^{2})Z_z=0$ gives $Z_z=1+h/(2\epsilon )+h^{2}[1/(2\epsilon ^{2})-1/(4\epsilon )]$. This supplies the mass-only substitutions in the one-loop baseline term. Curvature source matching and the remaining two-loop insertions must be tracked separately before asserting the projected higher-pole cancellation.

## Runnable check

**DERIVED** Python 3, SciPy and SymPy.

```python
from scipy.integrate import quad
from math import log
import sympy as s

def geom(r):
    return -log(r)/(1-r)**2-1/(1-r)

def evan(r):
    return 2/3-4*(r+2)*log(r)/(3*(r-1))

def evan_direct(r):
    def integrand(q):
        S,W=1+q*r,1+q
        C=2*q*r/S
        Y=(-2*W/3+2*(W*C+1+q*q*(1-2*r))/W-4*q/3)/S**2
        return (Y-(4/3 if q<1 else 0))/q
    return quad(integrand,0,1)[0]+quad(integrand,1,float("inf"))[0]

for r in (.001,.03,.25):
    val=2*quad(lambda q:r*log(1+q)/(1+q*r)**3,0,float("inf"))[0]
    assert abs(geom(r)-val)<1e-9
    assert abs(evan(r)-evan_direct(r))<1e-9
g=24*quad(lambda x:geom(x*(1-x)),0,1,epsabs=1e-11)[0]
v=8*quad(lambda x:evan(x*(1-x)),0,1,epsabs=1e-11)[0]
I0=-32+g+v+80/3-4/3
c=(I0/3+(16/9)*(7/6))/192
assert abs(c+0.0195838430501)<1e-10
assert abs(c-5/144+0.0543060652723)<1e-10

h,e=s.symbols("h e")
Z=1+h/(2*e)+h*h*(1/(2*e*e)-1/(4*e))
eq=s.expand((-2*e*h+3*h*h)*s.diff(Z,h)+(h-h*h)*Z)
assert all(eq.coeff(h,i)==0 for i in (0,1,2))
print("QCc",c,"scalar sum",c-5/144,"and mass pole recursion passed")
```

**ASSUMED** Next, apply the signed spin traces to the unreduced Yukawa graphs using the spinor line insertion, then combine those bare poles with the scalar sum and counterterms. Compare the resulting baseline simple pole with $1/24$ only after the projected higher poles cancel.

- **CITED** [Heinrich et al., Expansion by regions with pySecDec](https://arxiv.org/abs/2108.10807). — eq./page UNVERIFIED.

### Mass-only one-loop counterterm insertion

**DERIVED** In the baseline the one-loop curvature term is $z \exp [(\gamma _E-L)\epsilon ] \Gamma (\epsilon )/[6K(1-\epsilon )]$. Inserting $z_B=z Z_z$ raises Z_z to $1-\epsilon$, so at h² the binomial expansion, after dividing by $1-\epsilon$, gives $1/(2\epsilon ^{2})-3/(8\epsilon )$. Its Laurent coefficients in the same units $h^{2}z/K$ are

**DERIVED** $1/(12\epsilon ^{3}) - [1/16+L/12]/\epsilon ^{2} + [\zeta (2)/24+L/16+L^{2}/24]/\epsilon$.

**DERIVED** The linear-h term from the same substitution is $\Gamma (\epsilon ) \exp [(\gamma _E-L)\epsilon ] h z/(12K\epsilon )$, which recovers the two-loop checkpoint's counterterm normalization for $w=2$. This calculation accounts only for the mass substitution and uses no pure-metric finite redefinition.

```python
L=s.symbols("L")
binomial=s.simplify(s.diff((1-e)*(Z-1)-e*(1-e)*(Z-1)**2/2,h,2).subs(h,0)/2/(1-e))
assert s.cancel(binomial-1/(2*e**2)+3/(8*e))==0
P=s.expand((1/e-L+e*(s.zeta(2)+L**2)/2)*binomial/6)
assert P.coeff(e,-3)==s.Rational(1,12)
assert s.expand(P.coeff(e,-2)+s.Rational(1,16)+L/12)==0
assert s.expand(P.coeff(e,-1)-s.zeta(2)/24-L/16-L**2/24)==0
```
