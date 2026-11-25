# FFff graph-by-graph scalar insertion audit

**Conversation date range** 2026-11-22

**DERIVED** This checkpoint records FFff graph-by-graph scalar insertion audit. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** The [previous scalar comparison](ffz-scalar-evanescent-comparison.md) used the combined scalar insertion from the four-Yukawa kernels. Here the two graph coefficients are determined separately by differentiating their flat propagator reductions. This checks the scalar leg assignments, relative Wick weights and signs in this isolated sector.

**DERIVED** Use the bare integrals in the momentum-algebra reductions of Martin's (3.24) and (3.19), rather than substituting the UV-subtracted basis. Keep only the leading Laurent term of each master: tadpole −x, sunset −(x+y+z)/2, basketball E equal to one third of the sum of the six pair mass products, five-line G equal to minus one third of its central squared mass minus one sixth of the other four squared masses. The corresponding divided differences are 1, 1/2 and 1/3 for tadpole, sunset and G. Set the common fermion squared mass to one. Denote the scalar squared mass by s and set it to one only after differentiating. The linear curvature shift is εR/18 at the order needed.

**DERIVED** For A, the scalar legs are the doubled central propagators x,w in (3.24). The derivative of the reduced leading triple pole is the sum of 2/3 from K, −4/3 from G, +2 from the doubled sunset terms and +2 from the undoubled sunset terms. The E and double-tadpole terms are independent of s at this order. After the Wick weight 1/8, the result is 5/12 before multiplying by εR/18, giving **5/216** at the double pole.

**DERIVED** For B, the central x,w are fermions, and the scalar legs are y,z in (3.19). The derivatives of the K, G and E terms sum to zero. The two tadpole difference terms multiplying doubled sunsets give 1/2, and those multiplying undoubled sunsets give 3/2. The product of the two tadpole differences has zero first derivative at equal masses. After the Wick weight 1/2, the result is 1 before multiplying by εR/18, giving **1/18** at the double pole.

**DERIVED** Both values match the direct curved graph values individually. This excludes a compensating error between these graphs in the explicit scalar curvature-mass term. The graph normalization under a space-dependent mass profile and the full pole mismatch remain open.

## Reproduction

**DERIVED** The expressions below are the leading triple-pole polynomials evaluated with only the scalar squared mass varied. The symmetric difference differentiates these polynomials exactly because they have degree at most two.

```python
from fractions import Fraction as F

def E(a,b,c,d):
    return (a*b+a*c+a*d+b*c+b*d+c*d)/3

def G(w,a,b,c,d):
    return -w/3-F(1,6)*(a+b+c+d)

def a(s):
    return ((8*s-8-2*s*s)/6-E(F(1),F(1),F(1),F(1))
            +(2*s-4)*G(s,1,1,1,1)-2*(2-s)+2*(s+2)-4)/8

def b(s):
    return ((8-8*s+2*s*s)/12+(s-2)*G(F(1),1,s,1,s)
            +(2-s)*(s-1)/2+E(F(1),s,F(1),s)/2
            +(s-1)*(2+s)/2+(s-1)**2)/2

def derivative(fn):
    h=F(1,7)
    return (fn(1+h)-fn(1-h))/(2*h)

assert (derivative(a),derivative(b))==(F(5,12),F(1))
assert (derivative(a)/18,derivative(b)/18)==(F(5,216),F(1,18))
print("Both graphwise scalar insertion coefficients passed")
```

- **DERIVED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), Eqs. (3.19), (3.24).
- **DERIVED** [Martin and Robertson, Evaluation of the general 3-loop vacuum Feynman integral](https://doi.org/10.1103/PhysRevD.95.016008).
