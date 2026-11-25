# FFz(A) mixed long-face and overlap double pole

**Conversation date range** 2026-11-22

**DERIVED** This checkpoint records FFz(A) mixed long-face and overlap double pole. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Continues [the FFff(A) leading endpoints](ffz-ffff-a-leading-endpoints.md). Use effective perimeter one, physical cut 2 and L=0. These two add-backs are computed without imposing the counterterm target.

## Long face with both short Taylor terms removed

**DERIVED** The long logarithm comes only from the long fermion variance. At ε=0 the integrated geometry for the other effective pair length q is
**DERIVED** $G(q)=∫₀^(1-q) (1-P-q)[4(1-P)(1-q)+2Pq]dP=(1-q)^3(\frac{4}{3}+q).$
**DERIVED** The bounded partner has physical length U and smaller fraction a, so q=aU. Its smooth weight for one assignment is 1−a and its length weight is U^−2. Subtract G(0)+qG′(0) because the short expansion has a constant and a linear term. The remainder is q²+(5/3)q³−q⁴. The two orientations per pair, the two pair choices, (2/3)×3 from the long profile, and the long pole factor 1/2 give the coefficient 8 multiplying the finite integral. Multiplying by k=−1 and the radial-curvature-graph prefactor +1/288 yields the double-pole term.

## Mixed overlap

**DERIVED** Write A=d²/4 and B=d/2 for the coefficients in n=A(1−P)(1−q)+BPq. The integrated linear geometry coefficient is (B−5A)/6=−3+(19/6)ε. The remaining factors relative to their leading values are:

- **DERIVED** Long variance d−1 gives slope −2/3.
- **DERIVED** k=−1+3ε gives slope −3.
- **DERIVED** The short fraction integral 2∫₀^½ a(1−a)^(1−ε) da gives slope 13/12−ln 2.
- **DERIVED** The two analytic length integrals give cutoff factor 2^−ε, adding −ln 2.

**DERIVED** The leading angular overlap is +2/ε². Its relative slope is −131/36−2 ln 2. The radial prefactor adds slope 25/6. The resulting overlap graph double pole is 19/5184−ln 2/72 in units h²z/K. The Weyl term is bounded in the long length before multiplication by W and therefore has no long logarithm in graph A.

**DERIVED** The complementary mixed short face and the bounded chart remain to be integrated. The full FFff(A) double pole, and its comparison with −209/2592, are open.

## Exact checks

```python
from fractions import Fraction as F

def moment(k):
    a = F(1,2)**(k+1)/(k+1)-F(1,2)**(k+2)/(k+2)
    return a * F(2)**(k-1)/(k-1)

angular = 8*(moment(2)+F(5,3)*moment(3)-moment(4))
long_face = -angular/288
assert (F(2)-5*4)/6 == -3
assert (-1+5*4)*F(1,6) == F(19,6)
slope = -F(19,18)-F(2,3)-3+F(13,12)
assert slope == -F(131,36)
assert 2*(slope+F(25,6))/288 == F(19,5184)
print("A mixed long-face double pole:", long_face)
```

- **DERIVED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005).
- **DERIVED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177).
