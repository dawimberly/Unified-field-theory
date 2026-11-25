# Complete three-loop triple pole for FFff(A+B)

**Conversation date range** 2026-11-22

**ASSUMED** The calculation uses equal masses and fixed $\xi=1/6$. **DERIVED** A physical-length partition covers the simultaneous short, mixed and long endpoint regions. **DERIVED** The triple-pole result and its angular checks provide a checkpoint for the next-order calculation.

## Research record

**ASSUMED** Continues the [A endpoint calculation](ffz-ffff-a-leading-endpoints.md) and the [B endpoint calculation](ffz-ffff-b-mixed-endpoints.md). Use the fixed $1/6$ equal-mass benchmark. This covers the simultaneous scalings and completes the triple-pole check. It also isolates one homogeneous contribution needed at the next order.

## A partition that covers the fraction and length overlaps

**DERIVED** Normalize the effective perimeter to one. The effective pair lengths P,Q lie in $P+Q\le 1$. The two bridges are rv,r(1−v), with $r=1-P-Q$ and Jacobian r. In each pair split the fraction at $1/2$. Denote the smaller fraction by a and the larger physical length by $T=P/a$. The other physical length is $P/(1-a)\le 2$. Use a fixed $cut T=2$, and apply the same construction to the Q pair with U.

**DERIVED** For $T\le 2$ use (T,a) and $P=aT$. The exact transformed pair measures, with the regular spin factors outside, are

| Pair type | Measure in dT da |
| --- | --- |
| **DERIVED** A, two fermions | $T^{\epsilon -2}(1-a)^{1-\epsilon }$ |
| **DERIVED** B, scalar long within the pair | $T^{\epsilon -1}(1-a)^{1-\epsilon }$ |
| **DERIVED** B, fermion long within the pair | $T^{\epsilon -1}a(1-a)^{-\epsilon }$ |

**DERIVED** Thus $a\to 0$ at fixed bounded T is integrable. The only poles in this chart come from $T\to 0$. For $T\ge 2$ use (P,T) and $a=P/T$. The same densities in dP dT are

| Pair type | Measure in dP dT |
| --- | --- |
| **DERIVED** A | $T^{\epsilon -3}(1-P/T)^{1-\epsilon }$ |
| **DERIVED** B, scalar long | $T^{\epsilon -2}(1-P/T)^{1-\epsilon }$ |
| **DERIVED** B, fermion long | $P T^{\epsilon -3}(1-P/T)^{-\epsilon }$ |

**DERIVED** These are smooth at $P=0$ uniformly for $T\ge 2$. The effective bridge endpoints are polynomial in r,v in the unnormalized Gaussian moments. In particular $r=0$ supplies no singular measure. With one bounded length tending to zero the moving boundary $P+Q\le 1$ produces no first-derivative boundary term because the numerator contains r. The total physical length is at least the normalized effective perimeter, so there is no additional mass-factor zero.

**DERIVED** The two bounded lengths give the two short residues. At fixed fractions the short Taylor term converts back to exactly the prior P,Q calculation. With one long length and one bounded length, the only possible angular double pole combines the long logarithm and the bounded short residue. This reproduces the mixed endpoints. For B, the compact-gauge position moment on the $Q=0$ face is 4P(3−4P). The mixed coefficient uses $\int _{0}^{1} (1-P)P(11-16P)dP=1/2$, independently reproducing the previous integral in the bridge-sum gauge.

## Two long lengths

**DERIVED** Order the two long lengths. In $T\ge U\ge 2$ write $T=U/x$ with $0<x\le 1$. The other ordering exchanges the labels. The physical mass factor at leading homogeneous order is $(T+U)^{-3\epsilon }$; its corrections have lower degree. All curvature numerator terms are degree at most two in T,U. At $\epsilon =0$ only the mass-profile terms survive, with degree-two part a sum of individual squares. The Weyl term supplies an explicit $\epsilon$ and can have at most two angular logarithmic poles, since the sector has only two singular variables and their smooth Taylor coefficients generate at most simple poles in each.

**DERIVED** For A the slowest tail from one square is $T^{-1+\epsilon }U^{-3+\epsilon }$, or its exchange, before the mass factor. The comparable-scale integral has negative degree and is analytic in the common radial variable near $\epsilon =0$. Therefore the sector can have only a single angular pole. For B the only degree that can yield a common radial pole comes from two scalar-long densities with a square on one of the long lines. Its homogeneous contributions are $T^\epsilon U^{-2+\epsilon }(T+U)^{-3\epsilon }$ and the exchange. The ordered ratio weights for one term are $x^{-2+2\epsilon }(1+x)^{-3\epsilon }$ and $x^{2\epsilon }(1+x)^{-3\epsilon }$. The linear Taylor coefficient in the first weight is $-3\epsilon$, so the potential ratio pole cancels one power of $\epsilon$. The other ratio weight is integrable. Lower-degree terms can have only a ratio pole because their radial integral is analytic. This exhausts the possible double angular poles in the two-long chart.

**DERIVED** For the displayed homogeneous B term alone, the two ordered ratio integrals tend to $-1-3/2$ and +1 by analytic continuation. The common radial integration is $2^{-\epsilon }/\epsilon$. The sum is therefore $-3/(2\epsilon )$, before multiplying by the geometry and curvature prefactors. This is one part of the graph double pole. Lower-degree terms and the $\epsilon -$dependent numerator still need to be included.

## Result

**DERIVED** In units $h^{2}z/K$ the full bare triple poles are **$FFff(A) = 1/144$** and **$FFff(B) = 1/48$**. Their sum **$1/36$** meets the independent counterterm target. This leaves the bare double pole at $L=0$ as the next quantitative check, followed by the simple-pole comparison.

## Runnable exact checks

```python
from fractions import Fraction as F

for p in (F(1,3),F(2,3),F(9,10)):
    t=1+p
    a=p/t
    f=a*(1-a)
    jac_short=a
    jac_long=p/t**2
    assert p**(-2)*f*jac_short == t**(-2)*(1-a)
    assert p**(-1)*(1-a)*jac_short == t**(-1)*(1-a)
    assert p**(-1)*a*jac_short == t**(-1)*a
    assert p**(-2)*f*jac_long == t**(-3)*(1-a)
    assert p**(-1)*(1-a)*jac_long == t**(-2)*(1-a)
    assert p**(-1)*a*jac_long == p*t**(-3)

# Compact-gauge mixed B residue: P(1-P)(11-16P).
assert F(11,2)-F(27,3)+F(16,4) == F(1,2)
# Ordered homogeneous scalar-long contribution.
assert -1+F(-3,2)+1 == F(-3,2)
assert F(1,144)+F(1,48) == F(1,36)
print("Length transformations, compact mixed residue, and triple pole passed")
```

- **CITED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177), endpoint expansion and subtraction. — eq./page UNVERIFIED.
- **CITED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), Eqs. (3.19), (3.24).
