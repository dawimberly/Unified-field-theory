# FFZ three-loop quartic Newton sector

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**ASSUMED** This checkpoint records the FFZ three-loop quartic Newton sector. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues [the graph inventory](ffz-newton-graphs-and-pole-mass-checkpoint.md) and the verified [Yukawa trace script](ffz-yukawa-spin-traces-check.md). This checkpoint computes the bare QQ sector only, at fixed $\xi =1/6$ in $d=4-2\epsilon$. No change to the matter beta or the Newton moment constraints follows from an individual vacuum graph.

**ASSUMED** Use the flat resolvent from the original benchmark. For one endpoint at the origin its variation is $\delta G/\kappa=(4-d)r^2G/6-r^3G'/3$. The quartic interaction at the other endpoint contributes −(4−d)r². For n cross propagators, radial integration by parts gives $[(n(4-d))/6+(d+2)/3]$ times $\int r^{2}Gⁿ$ from their variation. The doubly self-contracted QQ graph has $n=2$ plus a local tadpole at each endpoint. The origin tadpole has zero variation at $\xi =1/6$. The other tadpole scales as $(z \exp\left(2\sigma \right))^{d/2-1}$, adding −(d−2)r². The total integrated curvature insertion therefore vanishes identically. The four-cross-line basketball has $n=4$ and no tadpoles. Its coefficient is $(2d-2)/3$ in units $\kappa$, giving $F_{R,QQ}=-|y|^4J_4/(48d)$ for $J_{4}=\int dᵈr r^{2}G(r)^{4}$. This uses the independently counted flat Wick factor −|y|⁴/8.

## Bare pole check

**DERIVED** The same J₄ equals minus the Laplacian in external momentum of the four-line banana at zero momentum. Equivalently direct Gaussian integration yields the Schwinger factor $2d R_{graph} U^{-d/2}$, where U is the sum of the four triple products and R_graph is the parallel resistance $1/\sum 1/s_i$. Use 24 ordered sectors with lengths (1,x,xy,xyt) times the largest length. The radial factor is $\Gamma (-1+3\epsilon )$ and the mass power $z^{1-3\epsilon }$, with the usual loop normalization $K^{-3}$. The ordered projective integrand is $x^{-1+2\epsilon }y^\epsilon t (1+t+yt+xyt)^{-3+\epsilon }[1+x+xy+xyt]^{1-3\epsilon }$. The x endpoint is the only pole in that projective integral. At $\epsilon =0$ its integrated residue before the factor 24 is $1/12$, and so the leading double pole in F_R,QQ is $+h^{2} z/(72K\epsilon ^{2})$. This is the bare contribution. Other graphs and counterterms must be combined before comparing with the projected cancellation.

## Runnable check

**DERIVED** Python 3, SymPy and SciPy.

```python
import sympy as s
from scipy.integrate import dblquad
d,t,y=s.symbols("d t y",positive=True)
def cross(n):
    return n*(4-d)/6+(d+2)/3
assert s.expand(cross(2)-(4-d)-(d-2))==0
assert s.simplify((cross(4)-(4-d))/(4*d*(d-1))) == 1/(6*d)
endpoint=s.integrate(s.integrate(t/(1+t+y*t)**3,(y,0,1)),(t,0,1))
assert endpoint==s.Rational(1,12)
quad=dblquad(lambda y,t:t/(1+t+y*t)**3,0,1,0,1)[0]
assert abs(quad-float(endpoint))<1e-13
assert s.Rational(-1,48)*2*s.Rational(-1,3)*24*endpoint/2==s.Rational(1,72)
print("Quartic curvature and double pole passed")
```

**ASSUMED** Next, retain the same ordered sectors to subtract the endpoint and extract the basketball simple pole. Then apply the position-space variation and covariance to the cubic and Yukawa graphs. A diagnostic on the eventual complete result is the known Y-improvement residue $1/24$ and the vanishing of the projected higher poles.

- **CITED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005). — eq./page UNVERIFIED.
- **CITED** [Bunch and Parker, Feynman propagator in curved spacetime](https://doi.org/10.1103/PhysRevD.20.2499). — eq./page UNVERIFIED.

## Simple-pole extraction

**DERIVED** With $L=\ln\left(z/\mu ^{2}\right)$, the same bare contribution has simple-pole coefficient $0.0625000000000-L/24$ in units $h^{2}z/K$. The constant is numerically consistent with $1/16$. To reproduce it subtract the $x=0$ numerator before setting $\epsilon =0$ and include the derivative of the endpoint numerator. The radial factor after the $MS-bar$ measure is $-1/(3\epsilon )+L-1$. Run the following after the first check.

```python
from math import log
from scipy.integrate import tplquad
angular = dblquad(lambda y,t:t/(1+t+y*t)**3*log(y*(1+t+y*t)),0,1,0,1)[0]
regular = tplquad(lambda y,t,x:t*((1+x+x*y+x*y*t)/(1+t+y*t+x*y*t)**3
                                  -1/(1+t+y*t)**3)/x,0,1,0,1,0,1)[0]
assert abs((angular/2+regular)/3+1/24 - 1/16)<1e-11
```

**DERIVED** The remaining graph classes, including the cubic and Yukawa sectors, and the counterterms are required to decide the projected pole. Audit priority is the evanescent numerator and the normalization of the subtracted bare sector before assembling those classes.
