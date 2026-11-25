# FFZ(B) complete double pole at L=0

**Conversation date range** 2026-11-22

**ASSUMED** This checkpoint records FFZ(B) complete double pole at $L=0$. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues [the mixed long face](ffz-ffff-b-mixed-long-face.md). Effective perimeter one, physical cut 2. This completes the complementary short face and the overlap. The graph coefficients below are in units $h^{2}z/K$.

## Mixed short face

**DERIVED** At $\epsilon =0 \sum$ the two long assignments. After subtracting its constant and logarithmic long tails the short-face mass profile is $(2/3)r$ times
**DERIVED** $T^{-2} [B_{bridge} + s M + (3+P/T)s^{2}]$,
**DERIVED** where $M=M_{ii}+M_{jj}+M_{ij}$ for the long pair, B_bridge includes the two bridge lengths and their variance terms, and $s=P/(1-P/T)$. Substitute $h=1/T$ and $x=1-P/2$. Integration over $T\ge 2$ gives $(2/3)r[B_{bridge}/2+4P(1/x-1)+(P-M)\ln x]$. The resulting integrals over P and bridge division are $2543/180-(59/3)\ln 2$ and $113/10-(47/3)\ln 2$ for edges 0 and 2. Multiply their sum by $1/48$, including the short fraction weight, spin, radial residue and $k=-1$.

## Mixed overlap

**DERIVED** On the short face $W=T+1+O(1/T)$. The long logarithm coefficient of the profile before multiplying by k is $(2/3)[M-P+\epsilon (4P-12)]$, where M retains its d dependence. The 4P comes from the density expansion and −12 from $W^{-3\epsilon }$ multiplying the scalar quadratic tail. The W A term adds $-\epsilon [\sum _i M_{ii}+(8/3)(1-P)]$ to the logarithm. The integrated $\epsilon$ slopes of the full bracket $W A+kB$ are $23/9$ and $10/3$ with r included. Multiply by the spin normalization, the short fraction $\int J(\epsilon )$, and $2^{-\epsilon }/2$ from the product of the analytic short and long cut factors. The angular overlap is $-1/\epsilon ^{2}+[26/3+2 \ln 2]/\epsilon$. Radial conversion with slope $25/6$ gives the overlap double-pole coefficient $-1/32-\ln 2/72$.

## Assembly

**DERIVED** Add the [long face](ffz-ffff-b-mixed-long-face.md), the short face and the overlap. The mixed regions give $1079/2160-(3/4)\ln 2$. Adding the [bounded region](ffz-ffff-b-bounded-double-pole.md) and the [two-long region](ffz-ffff-double-long-double-pole.md) cancels ln 2 exactly:
**DERIVED** FFff(B) double pole: $c_{-2}(L=0)=245/2592$.
**DERIVED** The combined target $1/72$ therefore requires FFff(A) to give $-209/2592$. That independent graph evaluation and the two simple poles remain open.

## Reproduction

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

def geometry(p,v,edge):
    r=1-p
    t=(0,p,p+r*v,p+r*v) if edge==0 else (0,0,r*v,r*v+p)
    def moments(i,j):
        x=[F(2)*((k<i)-t[i]) for k in range(4)]
        y=[F(2)*((k<j)-t[j]) for k in range(4)]
        def tilt(k,l): return -(x[k]*y[l]+x[l]*y[k])/2
        a=tilt(0,1)+tilt(2,3)+tilt(0,3)+tilt(1,2)
        b=tilt(0,2)+tilt(1,3)
        c=2*t[min(i,j)]*(1-t[max(i,j)])
        return 4*c+(2*a-b)/3,-2*c+(-a+2*b)/9
    def pair(i):
        j=(i+1)%4
        return tuple(sum(moments(*ij)[k] for ij in ((i,i),(j,j),(i,j)))
                     for k in (0,1))
    m,m1=pair(edge)
    vertices=sum(moments(i,i)[0] for i in range(4))
    bridge=sum(length*(pair(i)[0]+3*length) for i,length in
               ((1,r*v),(3,r*(1-v))))
    return m,m1,vertices,bridge

def finite_short(edge):
    # x=1-p/2. Split into rational and log(x) terms.
    def parts(x):
        p=2*(1-x)
        m,bridge=(quad(lambda v:geometry(p,v,edge)[i]) for i in (0,3))
        return (1-p)*(bridge/2+4*p*(1/x-1)), (1-p)*(p-m)
    xs=[F(j+5,10) for j in range(6)]
    rational=interpolate(xs,[x*parts(x)[0] for x in xs])
    logs=interpolate(xs,[parts(x)[1] for x in xs])
    at=F(7,11)
    assert sum(c*at**j for j,c in enumerate(rational)) == at*parts(at)[0]
    assert sum(c*at**j for j,c in enumerate(logs)) == parts(at)[1]
    # Integral dx from 1/2 to 1. Include dp=2 dx and 2/3 profile.
    plain=sum(rational[j]*(1-F(1,2)**j)/j for j in range(1,6))
    plain+=sum(c*(F(1,2)**(j+1)-1)/(j+1)**2 for j,c in enumerate(logs))
    ln=rational[0]+sum(c*F(1,2)**(j+1)/(j+1) for j,c in enumerate(logs))
    return F(4,3)*plain,F(4,3)*ln

def slope(edge):
    def f(p,v):
        m,m1,verts,_=geometry(p,v,edge)
        return (1-p)*(2*(m-p)-F(2,3)*(m1+4*p-12)
                      -verts-F(8,3)*(1-p))
    return quad(lambda p:quad(lambda v:f(p,v)))

short=[finite_short(e) for e in (0,2)]
slopes=[slope(e) for e in (0,2)]
# Corner angular coefficient -1/e^2, next 3/2*sum(slope)-1/6+2 ln2.
corner=-(F(3,2)*sum(slopes)-F(1,6)-F(25,6))/144
mixed=F(1,960)+sum(r for r,l in short)/48+corner
mixed_ln=sum(l for r,l in short)/48-F(1,72)
b=mixed-F(10111,25920)-F(43,2880)
assert slopes == [F(23,9),F(10,3)]
assert (mixed,mixed_ln,b)==(F(1079,2160),F(-3,4),F(245,2592))
assert F(1,72)-b == F(-209,2592)
print('FFff(B) double pole and log cancellation verified')

```

- **CITED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177). — eq./page UNVERIFIED.
