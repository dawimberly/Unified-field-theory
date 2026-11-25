# FFZ scalar curvature-mass audit

**Conversation date range** 2026-11-22

**ASSUMED** This checkpoint records the FFZ scalar curvature-mass audit. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues [the fermion-resolvent audit](ffz-fermion-resolvent-audit.md). This concerns only the explicit scalar mass insertion from the difference between $\xi =1/6$ and $\xi =\xi _d$. It uses bare three-loop graphs, with matter insertions separately subtracted from the RG result.

## Independent flat-space target

**DERIVED** Give the X and Y scalars independent infinitesimal mass shifts s_X and s_Y. Hold the fermion and the cubic and Yukawa vertices fixed in taking these derivatives. For one source s, write $t=3$ for X and $t=0$ for Y, so $\beta _s/s=t h+O(h^{2})$. Normalize the pole by $zs/K$. The one-loop bare linear pole is $-1/\epsilon$, since the scalar tadpole is $-z/(K\epsilon )$. The vacuum counterterm with opposite sign has leading coefficient +1.

**DERIVED** In this normalization the local two-loop double pole before the overall vacuum counterterm is $-(1+t)/4$. The one-loop matter insertion contributes $-(1+t)/2$, leaving $+(1+t)/4$ as the bare two-loop double pole. At three loops the local triple pole is $-(4+t)(1+t)/24$. The one-loop matter insertion's triple pole is $-[1/2+t(t+3)/8+t/4]$: the second-order bare mass factors for z and s are $1/2$ and $t(t+3)/8$, and their first-order cross gives $t/4$. The two-loop matter insertion contributes $+(1+t)(4+t)/8$. For $t=0$ and 3 these insertions sum to zero separately. Hence the sum of bare mass derivatives has triple-pole coefficient $-(1/6+7/6)=-4/3$ in units $h^{2}z/K$.

**DERIVED** The curvature mass shift multiplies this by $(1/6-\xi _d)R = [\epsilon /18+O(\epsilon ^{2})]R$ at $d=4$. **The explicit scalar curvature-mass insertion in the full bare graph inventory must therefore give −2/27 at the double pole**, conditional on the benchmark's flat source continuation. This target follows from flat-space counterterm recursion.

## Extracting the same term from the four-Yukawa kernels

**DERIVED** For FFff(A) only the two scalar bridges contribute. At the double-short corner the relevant expression is $-(8/3)\epsilon W r^{2}n$. In effective lengths P,Q, write $W=1+P(1/f-1)+Q(1/g-1)$. The PQ coefficient of W r²n is $54-12/f-12/g$. Its fraction-weighted integral is $-5/2$. Converting the $\epsilon$ times angular double pole by the graph and radial prefactor $+1/288$ gives **$5/216$**. The Weyl term on this scalar mass insertion has no long logarithm in A, since both long lines are fermions.

**DERIVED** For B the double-short corner vanishes. The two mixed overlaps integrate the scalar-long logarithm proportional to $-(8/3)\epsilon (1-P)$ against the Jacobian 1−P, giving $+1/54$ after radial conversion. The two-long overlap gives $+1/27$ from the scalar-long cross term. Thus the scalar mass piece is **$1/18$** for B, and the two graphs sum to **$17/216$**.

**DERIVED** The other bare graph classes must supply **$-11/72$** in this isolated sector for the direct calculation to agree with the flat target. The remaining graphs' isolated scalar term has not been compared yet. The full double-pole discrepancy remains open.

## Exact check

```python
from fractions import Fraction as F

def bare(t):
    two_local = -F(1+t,4)
    two_bare = two_local+F(1+t,2)
    three_local = -F((4+t)*(1+t),24)
    ct_one = -F(1,2)-F(t*(t+3),8)-F(t,4)
    ct_two = two_bare*F(4+t,2)
    assert ct_one+ct_two==0
    return three_local-ct_one-ct_two

assert (bare(0),bare(3)) == (F(-1,6),F(-7,6))
flat_curvature = (bare(0)+bare(3))/18
corner_a=F(54,36)-F(24,6)
assert corner_a==F(-5,2)
a=F(-8,3)*corner_a/288
b=F(1,54)+F(1,27)
assert (flat_curvature,a,b,flat_curvature-a-b)==(
    F(-2,27),F(5,216),F(1,18),F(-11,72))
print("Flat-source target and four-Yukawa scalar insertion checked")
```

- **CITED** [Jack, Jones and Pickering, The soft scalar mass beta-function](https://arxiv.org/abs/hep-ph/9803405). — eq./page UNVERIFIED.
- **CITED** [Martin, Three-loop effective potential for softly broken supersymmetry](https://doi.org/10.1103/PhysRevD.109.015019). — eq./page UNVERIFIED.
