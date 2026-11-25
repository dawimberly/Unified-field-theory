# FFZ double-pole locality audit

**Conversation date range** 2026-11-22

**ASSUMED** This checkpoint records the FFZ double-pole locality audit. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues [the double-pole mismatch](ffz-ffff-double-pole-mismatch.md) and independently checks the [counterterm targets](ffz-matter-counterterms-yukawa-targets.md). The mismatch survives two checks on the target. This audit uses the fixed $1/6$ external-metric continuation and mass-independent subtraction.

## Check 1: baseline higher poles from baseline data alone

**DERIVED** For the one-loop coefficient R, the scalar heat-kernel term vanishes at $\xi =1/6$. The four-component Dirac coefficient is $tr(R/6-R/4)=-R/3$. With the fermionic determinant sign and $\Gamma (-1+\epsilon )$, the positive Newton pole is $z/(6K\epsilon )$. This provides a direct normalization check on the one-loop input. The two-loop baseline pole after matter insertions is $hz/K$ times $[1/(24\epsilon ^{2})-1/(24\epsilon )]$.

**DERIVED** For the normalized local pole polynomial $C(h,\epsilon )$, the bare-scale derivative operator is $-2\epsilon +\beta _h\partial _h+\beta _z/z$, where $\beta _h=-2\epsilon h+3h^{2}+O(h^{3})$ and $\beta _z/z=h-h^{2}+O(h^{3})$. Requiring its Laurent-pole part to vanish yields $6 C_23=4 C_12$ and $6 C_22=4 C_11-C_01$. Thus the baseline higher poles at order $h^{2}$ are **1/36 and −1/18**, respectively, without importing the improvement-source computation. Subtracting the known counterterm double pole $-5/48$ and the earlier graph subtotal $5/144$ gives the remaining bare target **$1/72$** again.

## Check 2: the mismatch produces a nonlocal simple pole

**DERIVED** The three-loop bare terms carry the common factor $\exp\left(-3\epsilon L\right)$ with $L=\ln\left(z/\mu ^{2}\right)$. Hence the residual double-pole discrepancy $\Delta =59/648$ necessarily produces an uncanceled **$-59 L/(216\epsilon )$** in the simple pole of the one-mass baseline, given the existing matter insertions. A mass-independent local pure-metric subtraction cannot compensate this scale logarithm.

**DERIVED** The three-group projection does not hide it. In the benchmark masses $z_g=M^{2}(1,2,4)$, the improvement weights satisfy $\sum z_g a_Y,g=-7M^{2}$ and the baseline sum is $\sum z_g=7M^{2}$. The discrepancy is in the improvement-independent coefficient. It leaves the projected simple-pole logarithm proportional to $7 \ln\left(M^{2}/\mu ^{2}\right)+10 \ln 2$, in units $h^{2}M^{2}/K$. Thus the mismatch has to be repaired in the complete calculation before the projected finite matching or simple coefficient can be interpreted.

## Independent numerical-free check of the spin-conditioned covariance

**DERIVED** The script computes quadratic insertions by brute-force Wick pairings of all six Gaussian factors at integer dimensions. It compares the result to the covariance-tilt formula for the two squared-vector spin weights and the signed four-cycle trace. This checks the moment contraction independently of the integrations. It does not verify the Weyl variation or the completeness of the graph sum.

```python
from fractions import Fraction as F
from itertools import product

c01,c12,c11=F(1,6),F(1,24),-F(1,24)
t3,t2=4*c12/6,(4*c11-c01)/6
assert (t3,t2)==(F(1,36),-F(1,18))
target=t2+F(5,48)-F(5,144)
delta=F(17,162)-target
assert (target,delta,-3*delta)==(F(1,72),F(59,648),-F(59,216))
assert sum((1,2,4))==7 and sum(z*k for z,k in ((1,0),(2,1),(4,2)))==10

def wick(items,C):
    if not items:
        return F(1)
    i,a=items[0]
    return sum(C[i][j]*wick(items[1:s]+items[s+1:],C)
               for s,(j,b) in enumerate(items[1:],1) if a==b)

def brute(pairings,i,j,d,C):
    return sum(sign*wick(((a,u),(b,u),(c,v),(f,v),(i,w),(j,w)),C)
               for sign,a,b,c,f in pairings
               for u,v,w in product(range(d),repeat=3))

def tilt(pairings,i,j,d,C):
    def dm(a,b):
        return C[a][i]*C[b][j]+C[a][j]*C[b][i]
    total=F(0)
    for sign,a,b,c,f in pairings:
        for weight,x,y,z,t in ((d*d,a,b,c,f),(d,a,c,b,f),(d,a,f,b,c)):
            total+=sign*weight*(d*C[i][j]*C[x][y]*C[z][t]
                               +dm(x,y)*C[z][t]+C[x][y]*dm(z,t))
    return total

A=((1,0,0,2,2),)
B=((1,0,1,2,3),(-1,0,2,1,3),(1,0,3,1,2))
for lengths in ((F(1,4),)*4,(F(1,10),F(1,5),F(3,10),F(2,5))):
    K=[[2*((lengths[i] if i==j else 0)-lengths[i]*lengths[j])
        for j in range(4)] for i in range(4)]
    vectors=[[int(i==j) for i in range(4)] for j in range(4)]
    vectors += [[int(i<j) for i in range(4)] for j in (1,2,3)]
    C=[[sum(x[i]*y[j]*K[i][j] for i,j in product(range(4),repeat=2))
        for y in vectors] for x in vectors]
    for xy in (A,B):
        for i,j,d in product(range(4,7),range(4,7),(1,2,4)):
            assert brute(xy,i,j,d,C)==tilt(xy,i,j,d,C)
print("Baseline recurrence, nonlocal discrepancy and Gaussian moments verified")
```

- **CITED** [Martin, Three-loop effective potential for softly broken supersymmetry](https://doi.org/10.1103/PhysRevD.109.015019). — eq./page UNVERIFIED.
