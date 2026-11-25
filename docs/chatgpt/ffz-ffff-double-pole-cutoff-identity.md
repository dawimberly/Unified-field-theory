# FFZ cutoff identity for the four-Yukawa double pole

**Conversation date range** 2026-11-22

**DERIVED** This checkpoint records the FFZ cutoff identity for the four-Yukawa double pole. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Continues [the two-long chart calculation](ffz-ffff-double-long-double-pole.md). Keep effective perimeter one but replace the physical cut 2 by any c≥2. At L=0, the complete two-long chart contributes the following to the bare double pole in units h²z/K.

- **DERIVED** A: −23/(4320c²)+1/(1296c³).
- **DERIVED** B: 31/(4320c)−1/54.

**DERIVED** For A, the finite partner density integrates to 1/(2c²)−Q/(3c³). The geometry integrals are ∫r n=23/60 and ∫r n(P+Q)=1/6. Insert these into 4∫r n[c^−2−(P+Q)/(3c³)] and multiply by −1/288. For B, the two partner densities add to U^−2, which integrates to 1/c. The fixed-other-length angular residue is 31/(180c). The homogeneous scalar-long and evanescent cross terms carry only c^−ε and c^−ε/ε² times an explicit ε respectively, so their angular simple-pole residues are independent of c. The graph normalization then gives the formula above.

## Check for assembling the other charts

**DERIVED** Moving c₁→c₂ transfers the band c₁<T<c₂ with its other length unbounded from the two-long chart to the mixed chart, plus its exchange. Within the band the other endpoint may still be long. The coefficient of its logarithm is integrated against the finite partner densities. The bounded square with both lengths in the band has no projective pole. The bounded/mixed boundary terms with a short other length cancel between those charts because the same short endpoint Taylor coefficients occur on each side. These facts give the cutoff change of the sum of bounded and mixed terms as the negative of the displayed change of the two-long chart. This identity is independent of the counterterm prediction.

**DERIVED** At c=2 the two-long coefficient is −419/25920. At c=3 it is −73/4374. If the predicted full double pole 1/72 holds, the remaining charts at c=2 must give 779/25920. These are checks for the eventual direct assembly. The remaining charts' cutoff-independent constant and the simple pole have not been obtained here.

## Exact checks

```python
from fractions import Fraction as F
from math import factorial

def moment(i,j):
    return F(factorial(i)*factorial(j),factorial(i+j+3))

def weighted(i,j):
    return 4*(moment(i,j)-moment(i+1,j)-moment(i,j+1))+6*moment(i+1,j+1)

assert weighted(0,0)==F(23,60)
assert 2*weighted(1,0)==F(1,6)

def chart(c):
    c=F(c)
    return -F(23,4320)/c**2+F(1,1296)/c**3+F(31,4320)/c-F(1,54)

def band(c1,c2):
    lo,hi=F(c1),F(c2)
    a=(4*weighted(0,0)*(lo**-2-hi**-2)
       -8*weighted(1,0)*(lo**-3-hi**-3)/3)
    b=F(31,180)*(1/lo-1/hi)
    return -a/288+b/24

assert chart(2)==F(-419,25920)
assert chart(3)==F(-73,4374)
assert chart(2)-chart(3)==band(2,3)
assert chart(2)-chart(5)==band(2,5)
assert F(1,72)-chart(2)==F(779,25920)
print("Cutoff and finite-band identities passed")
```

- **DERIVED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177).
