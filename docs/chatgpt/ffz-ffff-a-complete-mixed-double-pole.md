# FFff(A) mixed-region double pole

**Conversation date range** 2026-11-22

**ASSUMED** This checkpoint records FFff(A) mixed-region double pole. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues [the long-face and overlap calculation](ffz-ffff-a-mixed-long-terms.md). Effective perimeter one, physical cut 2 and $L=0$. This finishes the mixed regions for A by computing the short-face contribution after the long subtraction.

**DERIVED** Write $r=1-P-Q$. The long pair has effective length Lp, physical lengths T and $Lp/(1-Lp/T)$. Its density for one assignment at $\epsilon =0$ is $T^{-3}(1-Lp/T)$. Use the unnormalized moments n M_ij. Write E_long and E_short for their sums on the two pair edges, and H for their bridge sum with bridge variance terms. After removing the long logarithmic tail, integrating T from 2 to infinity gives the following profile, before the factor $2/3$:
**DERIVED** r [(1/8−Lp/24)(H + 6 q E_short)
   +(E_long−3Lp n)/2 + 3n(−ln(1−Lp/2)−Lp/2)].
**DERIVED** Here q is the short effective length and only its first Taylor coefficient is used. The factor 6 accounts for the short pair's physical total q/[a(1−a)] in the fraction integral. The two fermion assignments give $2\int a(1-a) da=1/6$ on $0\le a\le 1/2$ for the other terms and $2\int da=1$ for the q E_short term. The moving simplex boundary gives no derivative term because r vanishes there.

**DERIVED** The first Taylor coefficients, integrated over Lp and bridge split, are $-139/210+3 \ln 2$ and $-703/210+3 \ln 2$ with the two possible long edges and the same anchor. Multiplication by $-1/1296$ gives the short-face term $421/136080-\ln 2/216$. Including the long face and overlap, the total mixed-region coefficient is **1283/272160−ln 2/54** in units $h^{2}z/K$. The bounded region still needs evaluation before checking the full graph against $-209/2592$.

## Reproduction

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

def quad(fn):
    return sum(w*fn(F(j,4)) for j,w in enumerate((7,32,12,32,7)))/90

def parts(L,v,edge):
    short=Jet((0,1,0,0))
    p,q=(Jet(L),short) if edge==0 else (short,Jet(L))
    r=1-p-q
    n=4*(1-p)*(1-q)+2*p*q
    lens=(p,r*v,q,r*(1-v))
    t=[sum(lens[:i]) for i in range(4)]
    a=[2*((i>0)-t[i]) for i in range(4)]
    c=[2*((i>2)-t[i]) for i in range(4)]
    def m(i,j):
        return (8*t[min(i,j)]*(1-t[max(i,j)])*n
                +4*p*(1-q)*a[i]*a[j]+4*q*(1-p)*c[i]*c[j]
                -2*p*q*(a[i]*c[j]+c[i]*a[j]))
    def e(i):
        j=(i+1)%4
        return m(i,i)+m(j,j)+m(i,j)
    bridge=sum(lens[i]*(e(i)+4*lens[i]*n) for i in (1,3))
    weight=F(1,8)-L/24
    plain=r*(weight*(bridge+6*short*e(2-edge))
             +F(1,2)*(e(edge)-3*L*n)-F(3,2)*L*n)
    return plain.c[1],(-3*r*n).c[1]

def integrate(edge):
    def val(x,k): return quad(lambda v:parts(2*(1-x),v,edge)[k])
    xs=[F(j+6,12) for j in range(7)]
    plain,lp=[interpolate(xs,[val(x,k) for x in xs]) for k in (0,1)]
    for k,row in enumerate((plain,lp)):
        assert sum(z*F(7,11)**j for j,z in enumerate(row)) == val(F(7,11),k)
    r=sum(z*(1-F(1,2)**(j+1))/(j+1) for j,z in enumerate(plain))
    r+=sum(z*(F(1,2)**(j+1)-1)/(j+1)**2 for j,z in enumerate(lp))
    l=sum(z*F(1,2)**(j+1)/(j+1) for j,z in enumerate(lp))
    return 2*r,2*l
rows=[integrate(i) for i in (0,2)]
coef=tuple(-sum(row[i] for row in rows)/1296 for i in (0,1))
assert rows==[(F(-139,210),F(3)),(F(-703,210),F(3))]
assert coef==(F(421,136080),F(-1,216))
print('Mixed-region coefficient:',coef[0]-F(53,25920)+F(19,5184),coef[1]-F(1,72))

```

- **CITED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177). — eq./page UNVERIFIED.
