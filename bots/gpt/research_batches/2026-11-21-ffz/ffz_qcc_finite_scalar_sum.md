# FFZ finite QCc projective remainder and bare scalar sum

Continues [the missing-endpoint correction](ffz_qcc_missing_corner_checkpoint.md). Use the same fixed-1/6 continuation and bare graph normalization. Separate the projective QCc integral into the vertex-position term G and the terms proportional to ε, denoted E. Write W=1+p+q, P=p f(1−f), Q=q g(1−g), S=1+P+Q. The anchor remains the quartic vertex throughout.

The sum of the two nonroot coordinate variances per Cartesian dimension is 2(P+Q+2PQ)/S. Hence G has numerator −d(6−d)(P+Q+2PQ)/S. First omit W^(−3ε) from this term. The P,Q integrals are Dirichlet integrals and give the following, including the two split-fraction integrals:

G₀=−d(6−d) Beta(1−ε,1−ε)² [2Γ(1+ε)Γ(ε)Γ(2−3ε)+2Γ(1+ε)²Γ(1−3ε)]/Γ(3−ε).

G₀=−8/ε−32+O(ε). Restoring W^(−3ε) contributes a finite term from each once-subtracted face. On p=0 the geometric measure is −8r dq/(1+qr)³, with r=g(1−g). Its correction comes from −3ε ln(1+q), times the p pole. The two faces give a total

G₁=48∫dg ∫₀∞ dq r ln(1+q)/(1+qr)³
  =24∫dg [−ln r/(1−r)²−1/(1−r)].

The overlap has zero geometric numerator. The difference after face subtraction has no 1/ε pole and vanishes for this correction at ε=0.

For E use the two globally integrated long-pair models from the preceding checkpoint. They contribute 32/(3ε)+80/3. After subtracting the models, the evanescent integrand contributes only through the small-p and small-q faces and their overlap. A model evaluated on its own short-pair face vanishes. On the opposite face keep the model subtraction for the long pair. With the q<1 cutoff for the overlap, the remaining q integration, divided by dε, is

∫₀∞ dq/q { [−2W/3 + 2(W C+1+q²(1−2r))/W − 4q/3]/(1+qr)² − (4/3)1_(q<1) },

where W=1+q and C=2qr/(1+qr). Its exact partial fraction integral is 2/3 − [4(r+2)/(3(r−1))] ln r. The two faces together contribute eight times its g average. The overlap, 2dε/(3ε²), contributes 8/(3ε)−4/3. The mixed remainder after the two short-face and two long-tail subtractions is multiplied by ε and contributes zero to this Laurent order. Here the short-face difference supplies a power of the short pair parameter and the long-tail remainder loses the linear growth in the long pair. In the joint long limit the spectator pair measure suppresses the overlap. The geometric Dirichlet integral is treated separately, so it does not rely on these evanescent remainder estimates.

Summing gives the projective coefficients 16/(3ε) and I₀=−17.5025158190. Multiplication by Γ(3ε), the MS-bar scale factors, and 1/[16d(d−1)] gives the bare doubled QCc Laurent coefficients, in units h²z/K,

1/(108ε²) + [−0.01958384305 − L/36]/ε.

The exact constant is defined by the convergent one-dimensional integrals above. It includes the 7ε/6 correction from 1/[d(d−1)]. The bare sum over all scalar-only graphs is now

5/(216ε²) + [−0.05430606527 − 5L/72]/ε.

This sum omits Yukawa graphs and counterterms. The finite QCc result has an independent integration check of its individual face integrals below. A direct Laurent extrapolation of the entire graph would be a stronger cross-check.

## Counterterm check for assembly

The flat supersymmetric mass renormalization factor for z follows from β_z/z=h−h²+O(h³). For z_B=z Z_z, the relation (−2εh+3h²)∂ₕ Z_z+(h−h²)Z_z=0 gives Z_z=1+h/(2ε)+h²[1/(2ε²)−1/(4ε)]. This supplies the mass-only substitutions in the one-loop baseline term. Curvature source matching and the remaining two-loop insertions must be tracked separately before asserting the projected higher-pole cancellation.

## Runnable check

Python 3, SciPy and SymPy.

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

Next, apply the signed spin traces to the unreduced Yukawa graphs using the spinor line insertion, then combine those bare poles with the scalar sum and counterterms. Compare the resulting baseline simple pole with 1/24 only after the projected higher poles cancel.

- [Heinrich et al., Expansion by regions with pySecDec](https://arxiv.org/abs/2108.10807).

### Mass-only one-loop counterterm insertion

In the baseline the one-loop curvature term is z exp[(γ_E−L)ε] Γ(ε)/[6K(1−ε)]. Inserting z_B=z Z_z raises Z_z to 1−ε, so at h² the binomial expansion, after dividing by 1−ε, gives 1/(2ε²)−3/(8ε). Its Laurent coefficients in the same units h²z/K are

1/(12ε³) − [1/16+L/12]/ε² + [ζ(2)/24+L/16+L²/24]/ε.

The linear-h term from the same substitution is Γ(ε) exp[(γ_E−L)ε] h z/(12Kε), which recovers the two-loop checkpoint's counterterm normalization for w=2. This calculation accounts only for the mass substitution and uses no pure-metric finite redefinition.

```python
L=s.symbols("L")
binomial=s.simplify(s.diff((1-e)*(Z-1)-e*(1-e)*(Z-1)**2/2,h,2).subs(h,0)/2/(1-e))
assert s.cancel(binomial-1/(2*e**2)+3/(8*e))==0
P=s.expand((1/e-L+e*(s.zeta(2)+L**2)/2)*binomial/6)
assert P.coeff(e,-3)==s.Rational(1,12)
assert s.expand(P.coeff(e,-2)+s.Rational(1,16)+L/12)==0
assert s.expand(P.coeff(e,-1)-s.zeta(2)/24-L/16-L**2/24)==0
```
