# FFZ(B) mixed-region long face after short subtraction

**Conversation date range** 2026-11-22

**DERIVED** This checkpoint records FFZ(B) mixed-region long face after short subtraction. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Continues [the bounded-chart integration](ffz-ffff-b-bounded-double-pole.md) using effective perimeter one, cutoff 2, and L=0. This computes the long-face add-back with its short limit subtracted.

**DERIVED** At ε=0 the sum of the scalar-long and fermion-long logarithm coefficients in the mass profile is (2/3)(M_ii+M_jj+M_ij−P) for long effective length P. The analytic long-length integral contributes 1/(2ε). For the partner's bounded physical length U and smaller fraction a, the sum of the two assignments gives the unit fraction density and short effective length s=aU. Integrate the logarithm coefficient over P from 0 to 1−s and the bridge-division fraction from 0 to 1, retaining bridge Jacobian r=1−P−s. Strip 2/3 in defining the resulting polynomials.

**DERIVED** For long edges 0 and 2 respectively the polynomials are
**DERIVED** G₀(s)=1/2−s/6−5s²/2+7s³/2−4s⁴/3,
**DERIVED** G₂(s)=1/2+23s/6−33s²/2+39s³/2−22s⁴/3.

**DERIVED** After subtracting G_i(0), integrate with dU/U da over U∈[0,2] and a∈[0,1/2]. A term s^k integrates to 1/[2k(k+1)]. The two answers are −11/80 and +17/80. Multiplying their sum by 1/72, from (2/3), long-pole factor 1/2, spin 6, mass-profile factor −1 and radial conversion −1/144, gives **+1/960** to the double pole in units h²z/K.

**DERIVED** This term is insensitive to the ε slope of the integrand because the short-subtracted face is finite. The complementary short face with its long subtraction and the simultaneous endpoint's ε slope are still needed for the full mixed coefficient.

## Exact check

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

def polynomial(fn):
    ts = [F(j,4) for j in range(5)]
    rows = [[t**k for k in range(5)]+[fn(t)] for t in ts]
    for k,row in enumerate(rows):
        pivot=row[k]
        rows[k]=[x/pivot for x in row]
        for j,r in enumerate(rows):
            if j!=k:
                rows[j]=[a-r[k]*b for a,b in zip(r,rows[k])]
    return [row[-1] for row in rows]

def face(short,edge):
    def at(x):
        long=(1-short)*x
        p,q=(long,short) if edge==0 else (short,long)
        def mv(v):
            if edge==0: return b_moment(p,q,v,1,1)
            return sum(b_moment(p,q,v,i,j) for i,j in ((2,2),(3,3),(2,3)))
        return (1-short)*(1-p-q)*(quad(mv)-long)
    return quad(at)

polys=[polynomial(lambda s:face(s,i)) for i in (0,2)]
for i,coeff in zip((0,2),polys):
    assert coeff[0] == F(1,2)
    assert sum(c*F(2,7)**k for k,c in enumerate(coeff)) == face(F(2,7),i)
# Integrate (aU)^k over a in [0,1/2] and dU/U in [0,2].
totals=[sum(c/F(2*k*(k+1)) for k,c in enumerate(row) if k)
        for row in polys]
coef=sum(totals)/72
assert totals == [F(-11,80),F(17,80)]
assert coef == F(1,960)
print('Mixed long-face coefficient verified')

```

- **DERIVED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177).
