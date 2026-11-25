# FFZ FFff(B): one long constituent and one short mixed pair

**Conversation date range** 2026-11-22

**ASSUMED** This checkpoint records the FFZ FFff(B): one long constituent and one short mixed pair. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Continuation of [the four-Yukawa curvature checkpoint](ffz-ffff-curvature-corner.md). Equal masses, $d=4-2\epsilon$, fixed $\xi =1/6$. This evaluates the mixed long-constituent/short-pair regions of FFff(B) in the chart where the two bridge lengths sum to one.

## Exact numerator in this chart

**DERIVED** Let the effective cyclic lengths be P,v,Q,1−v. In each doubled mixed pair let a,b be the fermion fractions. The two physical lengths in the P pair are $P/a$ (scalar) and $P/(1-a)$ (fermion). Write $S=1+P+Q$. The spin numerator divided by the four kinetic denominators is exactly
**DERIVED** $N = d(d+2)(1-a)(1-b)/(4 S^2).$
**DERIVED** This follows because all three products of off-diagonal increment covariance entries in the four-cycle trace equal $4PQv(1-v)/S^{2}$. Combining N with the Gaussian determinant gives $S^{-d/2-2}$. The pair measures are $P^{\epsilon -1}[a(1-a)]^{-\epsilon }$ and the corresponding expression with Q,b.

## Mixed endpoint residue

**DERIVED** First keep P and v in a compact interior interval and take Q→0 at fixed b. For the long-constituent endpoint of the P pair, the quadratic position moment sums to M_11 because vertex 0 is the anchor. For the insertion |x_1|², the increment-position covariance divided by the effective lengths has entries $2/(1+P)$ on the P edge and $-2P/(1+P)$ on the other three. All three products in the trace consequently have the same relative tilt derivative, $4P(1-P)/(1+P)$. Adding d C_11 yields $M_11=4P(3-P)/(1+P)$ at $d=4$.

**DERIVED** At a→0 it is the scalar that becomes long. In the mass-profile term multiply by the trace weight (1−a) before expanding. The residue of the resulting logarithmic endpoint, with the common factor $2/3$ stripped, is P M_11−dP². The −dP² comes from multiplying the scalar length squared P²/a² by (1−a). At a→1 the fermion length is $P/(1-a)$; the trace weight suppresses the quadratic divergence and leaves a logarithm of coefficient (d−1)P². The sum is P M_11−P². The mass factor $W_{total}^{-3\epsilon }$ combines with $[a(1-a)]^{-\epsilon }$ so the long logarithm gives $1/(2\epsilon )$. The Q pole gives $1/\epsilon$. This calculation uses endpoint Taylor coefficients at finite $\epsilon$ followed by $\epsilon \to 0$. Terms proportional to $\epsilon$ at these two endpoints contribute to lower poles.

**DERIVED** After removing the interior cutoffs in this coefficient, the P integral is
**DERIVED** $\int _{0}\infty dP P(11-5P)/(1+P)^{5} = 1/2$.
**DERIVED** The v integral is unity. The short-pair b integral gives $1/2$. Include the $2/3$ from the mass profile, spin factor 6, long pole factor $1/2$, two choices for the long pair, and the radial curvature normalization ($+1/144$ after including $k=-1$). The mixed regions give **$+1/144$** multiplying $h^{2}z/(K \epsilon ^{3})$. Adding the previous double-short region gives a **localized subtotal +1/48**. This subtotal does not certify completeness of the bare triple pole. The proof still needs the complement, particularly simultaneous scaling of the fractions and bridge sum.

## Reproduction

**DERIVED** Python 3 standard library. This checks the signed long residues and integrates the projective coefficient exactly $via P=t/(1-t)$.

```python
from fractions import Fraction as F

for p in (F(1,9), F(1,2), F(1), F(3), F(20)):
    perimeter = 1+p
    u0, u = 2/perimeter, -2*p/perimeter
    tilt = -perimeter*(u0*u + u*u)
    moment = 8*p/perimeter + tilt
    assert moment == 4*p*(3-p)/perimeter
    scalar_log = p*moment - 4*p*p
    fermion_log = 3*p*p
    assert scalar_log + fermion_log == p*p*(11-5*p)/perimeter

# P=t/(1-t) gives 11 t(1-t)^2 - 5 t^2(1-t).
integral = 11*(F(1,2)-F(2,3)+F(1,4)) - 5*(F(1,3)-F(1,4))
assert integral == F(1,2)
mixed = F(1,144)*6*F(2,3)*F(1,2)*F(1,2)*2*integral
assert mixed == F(1,144)
assert F(1,72)+mixed == F(1,48)
print("FFff(B) mixed-endpoint pole and subtotal verified")
```

- **CITED** [Martin, Effective potential at three loops, Eq. (3.19)](https://doi.org/10.1103/PhysRevD.96.096005).
- **CITED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177). — eq./page UNVERIFIED.
