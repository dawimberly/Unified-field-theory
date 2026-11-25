# FFff(B) homogeneous two-long tail audit

**Conversation date range** 2026-11-22

**ASSUMED** This checkpoint records FFff(B) homogeneous two-long tail audit. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues the [two-long chart](ffz-ffff-double-long-double-pole.md) and the [double-pole discrepancy](ffz-ffff-double-pole-mismatch.md). This independently integrates the homogeneous squared-length and evanescent cross terms through the simple pole. Work at $L=0$ and use long cutoff c for both lengths.

**DERIVED** Combine the two orderings of T,U using $x=min(T,U)/max(T,U)$. The radial integral over the smaller length is $c^{-\epsilon }/\epsilon$. The combined ratio weight for one square is
**DERIVED** $J_s = ∫₀¹ [x^(-2+2\epsilon )+x^(2\epsilon )](1+x)^(-3\epsilon ) dx.$
**DERIVED** Integrating the derivative of $x^{-1+2\epsilon }(1+x)^{-3\epsilon }$, first in the convergent domain and then by continuation, gives
**DERIVED** $J_s = 2^(-3\epsilon )/(2\epsilon -1) + [3\epsilon /(2\epsilon -1)] H + J₂,$
**DERIVED** where $H=\int _{0}^{1} x^{-1+2\epsilon }(1+x)^{-1-3\epsilon } dx$ and J₂ is the second, regular term of J_s. Subtracting H's value at $x=0$ gives $H=1/(2\epsilon )-\ln 2+O(\epsilon )$. Direct differentiation gives $J_{2}=1+(1-6 \ln 2)\epsilon +O(\epsilon ^{2})$. Together this yields $J_s=-3/2-4\epsilon +O(\epsilon ^{2})$. The cross term has ratio weight $J_c=2\int _{0}^{1} x^{-1+2\epsilon }(1+x)^{-3\epsilon } dx=1/\epsilon +O(\epsilon )$. Its constant vanishes.

**DERIVED** The angular factor before the spin normalization is
**DERIVED** $(c^-\epsilon /\epsilon ) [(2dk/9) J_s - d(4-d) J_c/9]$,
**DERIVED** where $d=4-2\epsilon$ and $k=-1+3\epsilon$. This includes both square numerators and the simplex weight $1/6$. Its pole and constant are $4/(9\epsilon )-2/3-(4/9)\ln c$. Restoring $d(d+2)/4$ and the radial and curvature prefactor yields
**$-1/(54\epsilon ^{2}) + [-11/324+(\ln c)/54]/\epsilon$**
**UNVERIFIED** in units $h^{2}z/K$. The double pole matches the earlier homogeneous contribution exactly. The simple pole is new for this subtraction term. The nonhomogeneous remainders and other regions are separate, so this check leaves the total mismatch unresolved.

## Reproduction

**DERIVED** Python 3. The algebra uses exact rational coefficients with ln 2 and ln c collected separately.

```python
from fractions import Fraction as F

# J_s slopes. Each pair is (rational, coefficient of ln 2).
boundary = (-2,3)
h_term = (-3,3)
regular = (1,-6)
assert tuple(sum(x) for x in zip(boundary,h_term,regular)) == (-4,0)
js0,js1 = -F(3,2),-4
square0 = -F(8,9)*js0
square1 = -F(8,9)*js1+F(28,9)*js0
cross0,cross1 = -F(8,9),F(4,9)
angular = square0+cross0
constant = square1+cross1
assert (angular,constant)==(F(4,9),-F(2,3))
assert -6*angular/144==-F(1,54)
assert -(6*constant+(-5+6*F(25,6))*angular)/144==-F(11,324)
assert 6*angular/144==F(1,54)
print("Homogeneous two-long pole coefficients verified")
```

- **CITED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177). — eq./page UNVERIFIED.
