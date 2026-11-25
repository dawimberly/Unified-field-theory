# FFZ(B) complete bounded-region double pole

**Conversation date range** 2026-11-22

**ASSUMED** This checkpoint records FFZ(B) complete bounded-region double pole. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues [the double-short subtraction](ffz-ffff-b-double-short-subtraction.md) with physical length cut 2 and $L=0$. At $\epsilon =0$ a single-face pole selects the finite integral of the other face after subtracting the common corner. Calculate both faces with vertex 0 anchored. Their integrands are different in this representation.

**DERIVED** On either face let T be the longer length of the surviving pair and a its smaller length fraction. The opposite constituent has length $s=aT/(1-a)$. At $\epsilon =0 \sum$ the two assignments of which field occupies T, with respective weights 1−a and a. The summed pair contribution to the profile before $2/3$ is
**DERIVED** $(T+s)(M_{ii}+M_{jj}+M_{ij})+(4-a)T^{2}+(3+a)s^{2}$.
**DERIVED** Add both fermion bridges and multiply by $r=1-aT$. The angular face residue comes from integrating the bridge-averaged profile minus $4/3$ with measure $dT/T da$ over $0<T<2$ and $0<a<1/2$. No simplex boundary clipping is needed on a single face at cut 2. The pole of the collapsed pair supplies the fraction factor $1/2$. With spin $6, k=-1$ and the radial conversion $-1/144$ the overall face factor is $1/48$.

**DERIVED** Call the two integrals I₀ and I₂ according to the surviving cyclic edge. Exact polynomial integration gives
**DERIVED** $I_{0}=-671/60+20 \ln 2$,
**DERIVED** $I_{2}=-386/45+16 \ln 2$.
**DERIVED** Adding $(I_{0}+I_{2})/48$ to the corner coefficient $7/324$ gives the **full bounded-region** double pole for graph B:
**−10111/25920+(3/4) ln 2 ≈ 0.12977550887675**
**DERIVED** in units $h^{2}z/K$. This is a chart coefficient at cut 2. The mixed regions and the bounded region of graph A are still required for the full two-graph result.

## Reproducible exact integration

**DERIVED** The first interpolation integrates the polynomial in T after corner subtraction and division by T. The second multiplies by b² with $b=1-a$ to remove the physical-length denominators, integrates each monomial of the recovered polynomial divided by b², and retains ln 2 symbolically. The Gaussian tilt uses rational arithmetic throughout.

```python
from fractions import Fraction as F
from math import log

def quad(fn,hi=1):
    return hi*sum(w*fn(hi*F(j,4)) for j,w in enumerate((7,32,12,32,7)))/90

def numerator(T,a,v,edge):
    p,q=(a*T,0) if edge==0 else (0,a*T)
    r=1-p-q
    lengths=(p,r*v,q,r*(1-v))
    t=[sum(lengths[:i]) for i in range(4)]
    def m(i,j):
        x=[F(2)*((k<i)-t[i]) for k in range(4)]
        y=[F(2)*((k<j)-t[j]) for k in range(4)]
        def tilt(k,l):
            return -(x[k]*y[l]+x[l]*y[k])/2
        return (8*t[min(i,j)]*(1-t[max(i,j)])
                +(2*(tilt(0,1)+tilt(2,3)+tilt(0,3)+tilt(1,2))
                  -tilt(0,2)-tilt(1,3))/3)
    def edge_m(i):
        j=(i+1)%4
        return m(i,i)+m(j,j)+m(i,j)
    bridge=sum(lengths[i]*(edge_m(i)+3*lengths[i]) for i in (1,3))
    s=a*T/(1-a)
    return F(2,3)*r*(bridge+(T+s)*edge_m(edge)
                       +(4-a)*T*T+(3+a)*s*s)

def interpolate(xs,ys):
    out=[F(0)]*len(xs)
    for i,(x,y) in enumerate(zip(xs,ys)):
        row=[y]
        for j,t in enumerate(xs):
            if i==j: continue
            nxt=[F(0)]*(len(row)+1)
            for k,v in enumerate(row):
                nxt[k]-=t*v/(x-t)
                nxt[k+1]+=v/(x-t)
            row=nxt
        out=[a+b for a,b in zip(out,row)]
    return out

def integrate_face(edge):
    def integral_t(b):
        a=1-b
        # (v averaged numerator - corner average) / T is degree ≤3.
        xs=[F(j,2) for j in range(1,5)]
        ys=[(quad(lambda v:numerator(t,a,v,edge))-F(4,3))/t
            for t in xs]
        coeff=interpolate(xs,ys)
        assert sum(x*F(1,3)**j for j,x in enumerate(coeff)) == (
            quad(lambda v:numerator(F(1,3),a,v,edge))-F(4,3))*3
        return b*b*sum(value*F(2)**(j+1)/(j+1) for j,value in enumerate(coeff))
    xs=[F(j+8,16) for j in range(8)]
    coeff=interpolate(xs,[integral_t(b) for b in xs])
    assert sum(x*F(5,7)**j for j,x in enumerate(coeff)) == integral_t(F(5,7))
    rational=coeff[0]+sum(coeff[j]*(1-F(1,2)**(j-1))/(j-1)
                          for j in range(2,len(coeff)))
    return rational,coeff[1]

rows=[integrate_face(edge) for edge in (0,2)]
r=sum(row[0] for row in rows)
l=sum(row[1] for row in rows)
assert rows == [(F(-671,60),F(20)),(F(-386,45),F(16))]
assert F(7,324)+r/48 == F(-10111,25920) and l/48 == F(3,4)
print("bounded double pole:",F(7,324)+r/48,l/48,
      float(F(7,324)+r/48)+float(l/48)*log(2))

```

- **CITED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177). — eq./page UNVERIFIED.
