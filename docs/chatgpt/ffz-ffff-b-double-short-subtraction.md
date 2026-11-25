# FFZ(B) double-short subtraction through the double pole

**Conversation date range** 2026-11-22

**DERIVED** This checkpoint records FFZ(B) double-short subtraction through the double pole. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Continues [the physical-length partition](ffz-ffff-complete-triple-pole.md) and [the cutoff identity](ffz-ffff-double-pole-cutoff-identity.md). The result is the explicit double-short subtraction term of the FFff(B) bounded chart. The single-face integrals of that chart are separate.

**DERIVED** In effective perimeter gauge, use T,U∈[0,c], c≥2, with the smaller fraction within each pair denoted by a∈[0,1/2]. Subtract the corner value of the smooth numerator in the T^(ε−1) U^(ε−1) weight. This value is independent of T,U and its analytic add-back has the factor c^(2ε)/ε². At the corner the compact simplex constraint is automatically satisfied. The subtraction is defined over the full square, with the original numerator extended by zero outside its allowed simplex. This specifies the cutoff convention for the quoted term.

**DERIVED** At T=U=0 the sum of the scalar-long and fermion-long measures is (1−a)^(−ε). Its integral is J(ε)=[1−2^(ε−1)]/(1−ε). The spin normalization is d(d+2)/4 = 6−5ε+O(ε²). Both physical bridges are fermion lines. The curvature-mass contribution vanishes at the corner since the scalar lengths vanish. The Weyl vertex weight is −2ε M_XX, and the bridge average of M_XX vanishes at d=4. Its first contribution is of order ε².

**DERIVED** The profile average over bridge division v is
**DERIVED** $B_avg = (\frac{2}{3})[d+\frac{2}{3}-4d/(d+2)] = \frac{4}{3}-(\frac{28}{27})\epsilon +O(\epsilon ^2).$
**DERIVED** Multiplication by k=−1+3ε gives −4/3+(136/27)ε. Combining with the spin normalization, J(ε)² and c^(2ε) gives the angular subtraction
**DERIVED** −2/ε² + [47/9−4 ln(c/2)]/ε + O(1).

**DERIVED** The radial prefactor at L=0 is Γ(−1+3ε)e^(3γε)/[4d(d−1)]. Relative to its leading term −1/(144ε), its slope is 3+7/6=25/6. Consequently the graph contribution from this add-back is

**DERIVED** 1/(72ε³) + [7/324 + ln(c/2)/36]/ε² + O(1/ε)

**DERIVED** in units h²z/K. At c=2 its double-pole coefficient is **7/324**. The cutoff logarithm must be carried along with the bounded single-face and mixed terms. This is a component of the bounded chart, so the full double-pole target 1/72 remains untested.

## Exact check

```python
from fractions import Fraction as F

mean_m = F(16,6)-F(16,6)
assert mean_m == 0
profile0 = F(2,3)*(4+F(2,3)-F(16,6))
profile1 = -2*F(2,3)*(1-F(8,36))
assert (profile0,profile1) == (F(4,3),F(-28,27))
k_profile0 = -profile0
k_profile1 = 3*profile0-profile1
angular0 = 6*k_profile0/4
# At c=2, the log(2) terms from J(e)^2 and c^(2e) cancel.
angular1 = (6*k_profile1-5*k_profile0)/4+2*angular0
assert (angular0,angular1) == (-2,F(47,9))
assert -angular0/144 == F(1,72)
assert -(angular1+F(25,6)*angular0)/144 == F(7,324)
print("FFff(B) double-short add-back verified through the double pole")
```

- **DERIVED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177), endpoint subtraction and add-back.
