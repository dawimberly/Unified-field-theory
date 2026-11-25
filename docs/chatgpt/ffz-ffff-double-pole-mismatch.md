# FFZ four-Yukawa double-pole mismatch

**Conversation date range** 2026-11-22

**DERIVED** The bounded, mixed and two-long regions give $1/96$ for FFff(A). **DERIVED** The $A+B \sum$ is $17/162$ against the independent $1/72$ target, leaving $59/648$. **UNVERIFIED** The physical simple pole remains uninterpretable until that inconsistency is resolved.

## Research record

**DERIVED** This completes the bounded FFff(A) calculation at physical cut 2 and $L=0$, following [the A mixed regions](ffz-ffff-a-complete-mixed-double-pole.md).

**DERIVED** Expand the bounded smooth numerator in T,U through first order in each variable. The $T^{\epsilon -2} U^{\epsilon -2}$ factors select the first Taylor coefficient of either length for a simple endpoint pole. For the other length keep the analytically integrated constant ($-1/2$ at cutoff 2), the finite logarithmic factor ln 2 on the linear term, and ordinary integration on higher powers. The summed corner at $\epsilon =0$ is zero. Its first-order $\epsilon$ slope is $-182/27$, including the expansion of $W^{-3\epsilon }$ and the fraction factors. The logarithms of the two fractions vanish separately by the zero marginal integrals of the leading corner. The script uses exact Taylor arithmetic in T,U, polynomial integration, and an exact derivative in $\epsilon$.

**DERIVED** The two bounded single-face contributions to the angular simple pole, with vertex 0 anchored, are $7246/945+(8/3)\ln 2$ and $1012/945+(8/3)\ln 2$. Add the corner slope and convert by $+1/288$. The bounded graph double pole has logarithm $+\ln 2/54$, which cancels the mixed region's logarithm exactly. Adding the rational terms from the bounded, mixed and two-long regions gives **$FFff(A) = 1/96$** in units $h^{2}z/K$.

**DERIVED** Combining with [FFff(B) = 245/2592](ffz-ffff-b-complete-double-pole.md) gives **$17/162$**, versus the independent target **$1/72$**. The discrepancy is **$59/648$**. This is an inconsistency between the current calculation and the pole equations. The physical simple-pole test cannot be interpreted until the mismatch is resolved. An independent audit should rederive the curvature numerators and subgraph residues, and check the matter insertions and the earlier bare subtotal in one common convention.

## Reproduction

**DERIVED** Python 3 standard library. The polynomial integrations are exact. The full run may take a minute.

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

def kernel(T,U,a,b,v,e=0):
    p,q=a*T,b*U
    r=1-p-q
    d=4-2*e
    n=F(1,4)*d*d*(1-p)*(1-q)+F(1,2)*d*p*q
    lens=(p,r*v,q,r*(1-v))
    t=[sum(lens[:i]) for i in range(4)]
    u=[2*((i>0)-t[i]) for i in range(4)]
    w=[2*((i>2)-t[i]) for i in range(4)]
    def m(i,j):
        return (2*d*t[min(i,j)]*(1-t[max(i,j)])*n
               +F(1,4)*d*d*(p*(1-q)*u[i]*u[j]+q*(1-p)*w[i]*w[j])
               -F(1,2)*d*p*q*(u[i]*w[j]+w[i]*u[j]))
    def edge(i):
        j=(i+1)%4
        return m(i,i)+m(j,j)+m(i,j)
    s=a*T*(1/(1-a))
    z=b*U*(1/(1-b))
    profile=F(2,3)*(
        sum(lens[i]*(edge(i)+d*lens[i]*n) for i in (1,3))
        +(T+s)*edge(0)+(U+z)*edge(2)+(d-1)*n*(T*T+s*s+U*U+z*z))
    weyl=-e*sum(m(i,i) for i in range(4))-F(1,3)*d*(4-d)*r*n
    W=r+T+s+U+z
    return r*(W*weyl+(-1+3*e)*profile),W

D=Jet((0,1,0,0))
E=Jet((0,0,1,0))
def exact_face(flip):
    def f(x,k):
        return (x*x if k==0 else x*x)*2*x*face(1-x,flip)[k]
    xs=[F(j+8,16) for j in range(9)]
    out=[]
    for k in (0,1):
        coef=interpolate(xs,[f(x,k) for x in xs])
        at=F(7,11)
        assert f(at,k)==sum(y*at**j for j,y in enumerate(coef))
        plain=coef[0]+sum(y*(1-F(1,2)**(j-1))/(j-1)
                          for j,y in enumerate(coef) if j>=2)
        out.append((plain,coef[1]))
    return out

# For fixed surviving pair and v the short derivative is affine in
# the short effective fraction b and its total-length multiplier 1/(1-b).
# Recover the two terms by b=0 and b=1/2 with exact weights.
def face(a,flip):
    def fn(T):
        def fv(v):
            def fb(b):
                x,y,aa,bb=(D,T,b,a) if flip else (T,D,a,b)
                return kernel(x,y,aa,bb,v)[0].c[1]
            # fb = b*C+S/(1-b); integrated with 2(1-b).
            return (fb(F(1,2))+fb(F(0)))/3
        return quad(fv)
    xs=[F(j,3) for j in range(7)]
    coef=interpolate(xs,[fn(t) for t in xs])
    assert fn(F(1,3))==sum(x*F(1,3)**j for j,x in enumerate(coef))
    return -coef[0]/2+sum(x*F(2)**(j-1)/(j-1) for j,x in enumerate(coef) if j>=2),coef[1]

def raw(a,b,k):
    def f(v):
        h,W=kernel(D,E,a,b,v)
        if k==1: return (1-a)*(1-b)*h.c[3]
        zs=[F(j,10) for j in range(6)]
        slope=interpolate(zs,[kernel(D,E,a,b,v,e)[0].c[3] for e in zs])[1]
        lw=Jet((0,W.c[1],W.c[2],W.c[3]-W.c[1]*W.c[2]))
        return (1-a)*(1-b)*(slope-3*(h*lw).c[3])
    return quad(f)

def marginal(a):
    return quad(lambda b:raw(a,b/2,1))
# The two fraction log terms can be dropped if either marginal vanishes.
for a in (F(0),F(1,8),F(1,4),F(3,8),F(1,2)):
    assert marginal(a)==0
    assert quad(lambda b:raw(b/2,a,1))==0
ans=quad(lambda a:quad(lambda b:raw(a/2,b/2,0)))
print(ans)
assert ans==F(-182,27)

rows=[exact_face(flip) for flip in (False,True)]
assert rows==[[(F(7246,945),F(8,3)),(0,0)],[(F(1012,945),F(8,3)),(0,0)]]
bounded=(F(8258,945)+ans)/288
assert F(16,3)/288==F(1,54)
a=bounded+F(1283,272160)-F(1,810)
total=a+F(245,2592)
assert (a,total,total-F(1,72))==(F(1,96),F(17,162),F(59,648))
print("A, total, discrepancy:",a,total,total-F(1,72))

```

- **CITED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177). — eq./page UNVERIFIED.
- **CITED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005). — eq./page UNVERIFIED.
