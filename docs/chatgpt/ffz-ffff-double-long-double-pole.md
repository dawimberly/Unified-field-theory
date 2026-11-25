# FFZ double pole from the two-long chart

**Conversation date range** 2026-11-22

**ASSUMED** This checkpoint records the FFZ double pole from the two-long chart. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues [the complete triple-pole calculation](ffz-ffff-complete-triple-pole.md). Use effective perimeter one, the same physical length $cut T,U\ge 2$, bridge Jacobian $r=1-P-Q$, and $L=0$. The intermediate coefficients below are for the angular simple pole. Multiplication by the overall radial pole converts them to the graph double pole. The two-long chart has no angular double pole after the explicit evanescent factors are included.

## A

**DERIVED** At $\epsilon =0$ the density on either long orientation is $T^{-3}(1-P/T)$. Only the fermion variance 3T² can make a long logarithm. Its pole factor is $1/(2\epsilon )$. The other long length has an integrable density with $\int 1/8-Q/24$ for one orientation. Include the multiplicity of two orientations in each pair and the factor $2/3$ in the mass profile. With $n=4(1-P)(1-Q)+2PQ$, the angular B-profile residue for the chart is
**DERIVED** $4\int dP dQ r n [1/4-(P+Q)/24] = 16/45$.
**DERIVED** At $\epsilon =0$ the Weyl numerator has no degree-two term in this chart, because both long lines are fermions. The graph double-pole coefficient is therefore $-(16/45)/288 = -1/810$.

## B

**DERIVED** There are three terms. Strip the spin constant 6 and keep the phase-space Jacobian r in every integral.

1. **DERIVED** The fixed-other-length residue of the mass profile has coefficient $(2/3)(M_{ii}+M_{jj}+M_{ij}-P)$ for the P pair, summing scalar-long and fermion-long orientations. Integrating the density of the other pair gives $1/2$, since its two densities sum to $U^{-2}$. Including the long pole factor $1/2$ and the Q counterpart gives $31/360$.
2. **DERIVED** The homogeneous scalar-long squared-length contribution adds $-4/3$ to the mass-profile angular residue. For each of the two squares the preceding checkpoint obtained the analytic ratio and radial coefficient $-3/2$, which multiplies 4 times $2/3$. Integration of r over the simplex gives $1/6$.
3. **DERIVED** The evanescent scalar curvature-mass term cannot be dropped. In the scalar-long/scalar-long sector W A contains the cross term $-16\epsilon TU/3$. With the two leading densities, the ordered integrals give a combined $1/\epsilon ^{2}$, from a common-scale $1/\epsilon$ and two ordered ratio factors $1/(2\epsilon )$. The angular residue is therefore $-8/9$ after the simplex integral. Terms in W A with only one long power or one squared long length have at most one angular pole after extracting the explicit $\epsilon$ and give no contribution here.

**DERIVED** The kernel contains W A + kB with $k=-1$ at this order. Thus the stripped angular residue is $-31/360+4/3-8/9=43/120$. Multiplying by 6 and by the radial and curvature prefactor $-1/144$ yields $-43/2880$.

## Result and use

**DERIVED** The full two-long chart contributes **$-419/25920$** to the double pole of FFff$(A+B)$ in units $h^{2}z/K$. This chart value uses the chosen cut at 2. The bounded/bounded and the two mixed charts must use the same partition when assembling the full double pole. Those charts and the full simple pole remain open.

## Reproduction

**DERIVED** Python 3 standard library. The integrations are exact because the substituted geometry integrands are polynomials of degree at most five in each quadrature variable. The tilt weights $(2,2,-1)/3$ follow from the three equal products of off-diagonal increment covariances in the four-cycle trace at $d=4$.

```python
from fractions import Fraction as F

def quad(fn):
    return sum(w*fn(F(j,4)) for j,w in enumerate((7,32,12,32,7)))/90

def simplex(fn):
    return quad(lambda p: (1-p)*quad(lambda y: fn(p,(1-p)*y)))

def b_moment(p,q,v,i,j):
    r=1-p-q
    t=(0,p,p+r*v,p+r*v+q)
    ui=[2*((edge<i)-t[i]) for edge in range(4)]
    uj=[2*((edge<j)-t[j]) for edge in range(4)]
    def tilt(a,b):
        return -(ui[a]*uj[b]+ui[b]*uj[a])/2
    change=(2*(tilt(0,1)+tilt(2,3)+tilt(0,3)+tilt(1,2))
            -(tilt(0,2)+tilt(1,3)))/3
    return 8*t[min(i,j)]*(1-t[max(i,j)])+change

def b_edges(p,q,v):
    m=lambda i,j: b_moment(p,q,v,i,j)
    return m(1,1)+m(2,2)+m(3,3)+m(2,3)

a_profile = 4*simplex(lambda p,q: (1-p-q)*(4*(1-p)*(1-q)+2*p*q)
                                 *(F(1,4)-(p+q)/24))
b_profile = simplex(lambda p,q: (1-p-q)*quad(lambda v:b_edges(p,q,v)-p-q))/6
volume = simplex(lambda p,q: 1-p-q)
b_homogeneous = -8*volume
b_evanescent = -F(16,3)*volume
a_pole = -a_profile/288
b_pole = -6*(-b_profile-b_homogeneous+b_evanescent)/144
assert (a_profile,b_profile,volume)==(F(16,45),F(31,360),F(1,6))
assert (a_pole,b_pole,a_pole+b_pole)==(F(-1,810),F(-43,2880),F(-419,25920))
print('Two-long chart rational checks passed')

```

- **CITED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177). — eq./page UNVERIFIED.
- **CITED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005). — eq./page UNVERIFIED.
