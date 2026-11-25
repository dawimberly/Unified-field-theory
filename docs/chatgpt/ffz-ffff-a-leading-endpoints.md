# FFZ FFff(A) leading endpoint calculation

**Conversation date range** 2026-11-22

**ASSUMED** This checkpoint records the FFZ FFff(A) leading endpoint calculation. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues [the four-Yukawa kernels](ffz-ffff-curvature-corner.md) and checks the independent [matter-counterterm targets](ffz-matter-counterterms-yukawa-targets.md). Equal masses, fixed $\xi =1/6$. Normalize the effective perimeter to one, with lengths (P,rv,Q,r(1−v)) and $r=1-P-Q$. The bridge change of variables gives Jacobian r. Define $f=a(1-a), g=b(1-b)$, the products of physical length fractions in the two parallel pairs.

**DERIVED** The spin numerator with the kinetic denominators restored is $fg n/(PQ)$, where $n=[d^{2}(1-P)(1-Q)+2dPQ]/4$. The determinant is one after factoring the physical pair totals. Hence the pair measure including spin denominators is $P^{\epsilon -2} Q^{\epsilon -2} f^{1-\epsilon } g^{1-\epsilon }$, with r n in the numerator. The two short poles select the PQ Taylor coefficient. The mass-profile numerator is multiplied by $k=-1$ at $d=4$. The Weyl weight is proportional to $\epsilon$ and contributes to lower orders in the localized double-endpoint terms.

## Two short pairs

**DERIVED** The position moments in the check are unnormalized n times M_ij. For cumulative positions t_i, the Gaussian moment contributes 8 min(t_i,t_j)[1−max(t_i,t_j)]n. The additional terms are the directional derivative of the parallel-pair spin polynomial under the covariance tilt. The script expands these moments in the quotient where $P^{2}=Q^{2}=0$, retaining the actual PQ coefficient. Physical pair totals $P/f$ and $Q/g$ multiply their position moments. Squared pair lengths cannot enter this Taylor coefficient. The fraction integrals give $1/36$ for the bridge terms and $1/6$ for each term linear in a pair total. The resulting function of v has degree at most three. Exact quadrature gives **zero** for the integrated localized double-short triple pole.

## One long constituent and one short pair

**DERIVED** At a→0 the constituent length $P/a$ makes the fermion variance term 3P²/a² at $d=4$. Multiplication by f turns this into a logarithm. The opposite endpoint contributes the same amount. The long-constituent residue before multiplying by $2/3$ is therefore 6P² summed over the two constituent choices, with pole factor $1/(2\epsilon )$ from the physical mass-scale factor. The short Q pole differentiates r n with respect to Q, giving $(1-P)(-8+6P)$. Its P integral after cancelling P² is −3. The b fraction contributes $1/6$ and there are two choices for which pair contains the long line. The complete factor at this localized order is −2. Multiplication by $-1/288$ (radial residue $-1/3$, mass-profile factor −1, curvature conversion $1/48$ and graph coefficient $-1/2$) gives **+1/(144 ε³)** in units $h^{2}z/K$.

**DERIVED** Adding the verified contributions of A and B gives **$1/36$** at the triple pole, exactly the counterterm target. This agreement applies to the computed regions. The complementary scaling regions still need uniform control, in particular when large physical lengths arise simultaneously from both fractions. Consequently the agreement is a check rather than a completed graph pole extraction. The double- and simple-pole coefficients are still open.

## Reproduction

**DERIVED** Python 3 standard library. The quadrature integrates the degree-three polynomial exactly. All arithmetic uses rational inputs.

```python
from fractions import Fraction as F

class Jet:
    def __init__(self,x): self.c = x if isinstance(x,tuple) else (x,0,0,0)
    def __add__(self,y):
        y = y if isinstance(y,Jet) else Jet(y)
        return Jet(tuple(a+b for a,b in zip(self.c,y.c)))
    __radd__ = __add__
    def __neg__(self): return Jet(tuple(-v for v in self.c))
    def __sub__(self,y): return self+-Jet(y) if not isinstance(y,Jet) else self+-y
    def __rsub__(self,y): return -self+y
    def __mul__(self,y):
        y = y if isinstance(y,Jet) else Jet(y)
        out=[0]*4
        for i,x in enumerate(self.c):
            for j,z in enumerate(y.c):
                if not i&j: out[i|j]+=x*z
        return Jet(tuple(out))
    __rmul__ = __mul__

p,q=Jet((0,1,0,0)),Jet((0,0,1,0))
r=1-p-q
n=4*(1-p)*(1-q)+2*p*q

def integrand(v):
    lengths=(p,r*v,q,r*(1-v))
    t=[sum(lengths[:i]) for i in range(4)]
    # d=4. Set b_i/P and c_i/Q before taking the short limit.
    a=[2*((i>0)-t[i]) for i in range(4)]
    c=[2*((i>2)-t[i]) for i in range(4)]
    def m(i,j):
        x,y=t[min(i,j)],t[max(i,j)]
        return (8*x*(1-y)*n
                +4*p*(1-q)*a[i]*a[j]+4*q*(1-p)*c[i]*c[j]
                -2*p*q*(a[i]*c[j]+c[i]*a[j]))
    def edge(i,j): return m(i,i)+m(j,j)+m(i,j)
    bridge=(lengths[1]*(edge(1,2)+4*lengths[1]*n)
           +lengths[3]*(edge(3,0)+4*lengths[3]*n))
    return (F(2,3)*r*(bridge*F(1,36)
                         +(p*edge(0,1)+q*edge(2,3))*F(1,6))).c[3]

vals=[integrand(F(j,4)) for j in range(5)]
short=sum(x*y for x,y in zip(vals,(7,32,12,32,7)))/90
long=F(2,3)*6*F(1,6)*F(1,2)*2*(-8+F(14,2)-F(6,3))
assert short == 0 and long == -2
assert -long/288 == F(1,144)
assert F(1,48)-long/288 == F(1,36)
print('FFff(A) endpoint checks passed')

```

- **CITED** [Martin, Effective potential at three loops, Eq. (3.24)](https://doi.org/10.1103/PhysRevD.96.096005).
- **CITED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177). — eq./page UNVERIFIED.
