# Triple-pole mismatch resolved by FFff(A) spin normalization

The FFff(A) numerical implementation integrates the product of two coordinate dot products. Its numerator N uses E[Δ₁²Δ₂²]/(16s₁s₂s₃s₄), with no factor of 2 from either Weyl trace. In the signed flat Wick check each two-vertex trace is twice the dot product (with its oriented sign). The two separate traces factorize. The Wick and expansion coefficient −8/64=−1/8 multiplies the **full** trace product. Therefore the coordinate integrand used in FFFF_A_REGULATED.md must instead carry −1/2. This correction multiplies every inferred FFff(A) Laurent coefficient by four. It leaves the projective integral and the estimated extrapolation errors in those stripped units intact. The FFff(B) implementation already accounts for its own trace factor of two and keeps its +1 coordinate prefactor.

## Revised result

The inferred A triple is 1/144 and the B triple is exactly 1/48 for the derived kernel. They sum to 1/36, the independent bare triple target in the projected comparison. This explains the entire 1/192 gap without a change of model. A still relies on the regulated numerical extraction of its projective pole. B now has exact integrated sector coefficients −8,−16,−16,−8 for the four two-log families. The four power/log families integrate to zero exactly.

The corrected inferred A double is 5/288−L/48. The budget target for B becomes **−1/288−L/16**. B's double and simple constants still need independent evaluation. For the bare B simple pole, homogeneity predicts L/96+3L²/32 after using the revised double target. At L=0, a B simple coefficient near −0.04183 would be needed to cancel the projected Newton term, given the present numerical A and CcFf coefficients. This is a target only. The Newton simple-pole test remains open.

## Twenty-check batch

1. Read the sector-decomposition Laurent extraction reference.
2. Read Martin's flat two-fermion-loop normalization.
3. Recomputed B's two-log subgroup structure in ordered coordinates.
4. Integrated the first five-line B log family exactly.
5. Integrated the first paired-cycle B log family exactly.
6. Integrated the reflected paired-cycle B log family exactly.
7. Integrated the reflected five-line B log family exactly.
8. Integrated the nonzero-looking B power/log family and showed its sum vanishes.
9. Checked the reflected power/log family vanishes.
10. Checked the box nested in the first power endpoint gives zero.
11. Checked the box nested in the reflected power endpoint gives zero.
12. Compared A's polynomial numerator and endpoint-subtracted numerical integrand at ε=.01.
13. Repeated that comparison at ε=.1.
14. Audited the signed flat one-loop trace convention against the flat Wick script.
15. Audited the two-loop trace product and its symmetry coefficient.
16. Rechecked the Y pole recursion through cubic loop order.
17. Verified cancellation of the two counterterm triple poles.
18. Reassembled the corrected bare triple budget.
19. Recomputed the independent B double target and induced scale logarithms.
20. Recomputed the B simple coefficient needed at the projected finish line.

For checks 12–13 use the A polynomial in the B ordered sector script with edges 02,02,13,13,03,12 and trace polynomial d² c_ab c_cd+2d c_ac c_bd (the two pairs have identical within-pair directions). Compare after restoring W^(−3ε), undoing the effective-length measure, and subtracting the global tails, as in FFFF_A_REGULATED.md. The integrated A pole is distinct from this pointwise check.

## Exact log/log integration

Python 3 with SymPy.

```python
import sympy as s
from itertools import permutations
from collections import Counter
d=s.symbols('d')
x=s.symbols('a b c e f g')
E=((0,2),(1,2),(1,3),(0,3),(0,2),(1,3))
V=s.eye(4)[:,1:]
W=s.Matrix.vstack(*(V[i,:]-V[j,:] for i,j in E))
lap=W.T*s.diag(*[1/t for t in x])*W
U=s.cancel(lap.det()*s.prod(x))
C=(2*V*lap.adjugate()*V.T*s.prod(x)).applyfunc(s.cancel)
D=s.Matrix.vstack(*(s.eye(4)[i,:]-s.eye(4)[j,:] for i,j in E[:4]))
c=D*C*D.T
def tr(M):
 a,b,c,d1,e,f=(M[i,j] for i,j in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)))
 return d*d*(a*f-b*e+c*d1)+2*d*b*e
N=s.expand(tr(c))
cv=C*D.T
def moment(i,j):
 a,b=cv[i,:].T,cv[j,:].T
 dc=a*b.T+b*a.T
 return s.expand(d*C[i,j]*N+tr(c+dc)-tr(c)-tr(dc))
T=[moment(i,i) for i in range(4)]
A=s.expand(-(4-d)*sum(T)/2-d*(4-d)*N*U*sum(x[4:])/3)
B=s.expand(sum(s.Rational(2,3)*x[k]*(T[i]+T[j]+moment(i,j)+(d-int(k<4))*x[k]*N*U) for k,(i,j) in enumerate(E)))
# Fix largest length at 1. W=sum of lengths never vanishes on any sector.
# N,U and A,B have respective degrees 8,3,12,13.
assert [s.Poly(p,*x).total_degree() for p in (U,N,A,B)]==[3,8,12,13]
# After radial integration numerator is W*A+(k-1)B.
J=s.Poly(sum(x)*A+(5-3*d/2)*B,*x)
J4=s.Poly(J.as_expr().subs(d,4),*x)
up=s.Poly(U,*x)
from itertools import combinations
def powers(poly,subset):
 return min(sum(m[i] for i in subset) for m in poly.monoms())
table={}
for n in range(1,6):
 for ids in combinations(range(6),n):
  u=powers(up,ids)
  spin=sum(i<4 for i in ids)
  E0=n-1+powers(J,ids)-5*u-spin
  E4=n-1+powers(J4,ids)-5*u-spin
  table[ids]=(E0,E4,u)
print('singular subsets',[(tuple(i+1 for i in ids),p) for ids,p in table.items() if p[0]<0])
counts=Counter()
for perm in permutations(range(6)):
 z=[table[tuple(sorted(perm[i:]))] for i in range(1,6)]
 counts[(sum(a<-1 for a,b,c in z),sum(a==-1 for a,b,c in z),sum(a<0<=b for a,b,c in z))]+=1
print('sector multiplicities (power, log, evanescent improved):',sorted(counts.items()))
print('monomials',len(J.terms()),len(J4.terms()))
from collections import defaultdict
group=defaultdict(Counter)
z=s.symbols('z0:3')
for order in permutations(range(6)):
 sets=[tuple(sorted(order[i:])) for i in range(1,6)]
 dat=[table[ids] for ids in sets]
 axes=[i for i,p in enumerate(dat) if p[0]<0]
 if len(axes)!=2 or any(dat[i][0]!=-1 for i in axes):continue
 reg=[i for i in range(5) if i not in axes]
 def coeff(poly):
  vals=[powers(poly,ids) for ids in sets]
  out=0
  for mon,c in poly.terms():
   powers0=[sum(mon[j] for j in ids)-v for ids,v in zip(sets,vals)]
   if all(powers0[i]==0 for i in axes):
    out+=c*s.prod(v**powers0[i] for v,i in zip(z,reg))
  return out
 f=s.factor(coeff(J4)*s.prod(v**dat[i][1] for v,i in zip(z,reg))/coeff(up)**5/s.prod(dat[i][2] for i in axes))
 group[tuple(tuple(i+1 for i in sets[j]) for j in axes)][f]+=1
for g,terms in group.items():
 print(g,len(terms),s.factor(sum(f*n for f,n in terms.items())))
def integrate(expr):
 for a in reversed(z):
  expr=s.factor(s.integrate(expr,(a,0,1)))
 return expr
values=[]
for key,terms in group.items():
 value=integrate(sum(f*n for f,n in terms.items()))
 values.append(value)
 print(key,value)
assert sorted(values)==[-16,-16,-8,-8]
assert sum(values)==-48

```

## Exact power/log integration

```python
import sympy as s
from itertools import permutations
from collections import Counter
d=s.symbols('d')
x=s.symbols('a b c e f g')
E=((0,2),(1,2),(1,3),(0,3),(0,2),(1,3))
V=s.eye(4)[:,1:]
W=s.Matrix.vstack(*(V[i,:]-V[j,:] for i,j in E))
lap=W.T*s.diag(*[1/t for t in x])*W
U=s.cancel(lap.det()*s.prod(x))
C=(2*V*lap.adjugate()*V.T*s.prod(x)).applyfunc(s.cancel)
D=s.Matrix.vstack(*(s.eye(4)[i,:]-s.eye(4)[j,:] for i,j in E[:4]))
c=D*C*D.T
def tr(M):
 a,b,c,d1,e,f=(M[i,j] for i,j in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)))
 return d*d*(a*f-b*e+c*d1)+2*d*b*e
N=s.expand(tr(c))
cv=C*D.T
def moment(i,j):
 a,b=cv[i,:].T,cv[j,:].T
 dc=a*b.T+b*a.T
 return s.expand(d*C[i,j]*N+tr(c+dc)-tr(c)-tr(dc))
T=[moment(i,i) for i in range(4)]
A=s.expand(-(4-d)*sum(T)/2-d*(4-d)*N*U*sum(x[4:])/3)
B=s.expand(sum(s.Rational(2,3)*x[k]*(T[i]+T[j]+moment(i,j)+(d-int(k<4))*x[k]*N*U) for k,(i,j) in enumerate(E)))
# Fix largest length at 1. W=sum of lengths never vanishes on any sector.
# N,U and A,B have respective degrees 8,3,12,13.
assert [s.Poly(p,*x).total_degree() for p in (U,N,A,B)]==[3,8,12,13]
# After radial integration numerator is W*A+(k-1)B.
J=s.Poly(sum(x)*A+(5-3*d/2)*B,*x)
J4=s.Poly(J.as_expr().subs(d,4),*x)
up=s.Poly(U,*x)
from itertools import combinations
def powers(poly,subset):
 return min(sum(m[i] for i in subset) for m in poly.monoms())
table={}
for n in range(1,6):
 for ids in combinations(range(6),n):
  u=powers(up,ids)
  spin=sum(i<4 for i in ids)
  E0=n-1+powers(J,ids)-5*u-spin
  E4=n-1+powers(J4,ids)-5*u-spin
  table[ids]=(E0,E4,u)
print('singular subsets',[(tuple(i+1 for i in ids),p) for ids,p in table.items() if p[0]<0])
counts=Counter()
for perm in permutations(range(6)):
 z=[table[tuple(sorted(perm[i:]))] for i in range(1,6)]
 counts[(sum(a<-1 for a,b,c in z),sum(a==-1 for a,b,c in z),sum(a<0<=b for a,b,c in z))]+=1
print('sector multiplicities (power, log, evanescent improved):',sorted(counts.items()))
print('monomials',len(J.terms()),len(J4.terms()))
from collections import defaultdict
group=defaultdict(Counter)
z=s.symbols('z0:3')
for order in permutations(range(6)):
 sets=[tuple(sorted(order[i:])) for i in range(1,6)]
 dat=[table[ids] for ids in sets]
 axes=[i for i,p in enumerate(dat) if p[0]<0]
 if len(axes)!=2 or all(dat[i][0]==-1 for i in axes):continue
 power=next(i for i in axes if dat[i][0]==-2)
 reg=[i for i in range(5) if i not in axes]
 def coeff(poly,derivative=None):
  vals=[powers(poly,ids) for ids in sets]
  out=0
  for mon,c in poly.terms():
   powers0=[sum(mon[j] for j in ids)-v for ids,v in zip(sets,vals)]
   if all(powers0[i]==int(i==derivative) for i in axes):
    out+=c*s.prod(v**powers0[i] for v,i in zip(z,reg))
  return out
 f=s.factor((coeff(J4,power)-5*coeff(J4)*coeff(up,power)/coeff(up))*s.prod(v**dat[i][1] for v,i in zip(z,reg))/coeff(up)**5/s.prod(dat[i][2] for i in axes))
 group[tuple(tuple(i+1 for i in sets[j]) for j in axes)][f]+=1
for g,terms in group.items():
 print(g,len(terms),s.factor(sum(f*n for f,n in terms.items())))
def integrate(expr):
 for a in reversed(z):
  expr=s.factor(s.integrate(expr,(a,0,1)))
 return expr
values=[]
for key,terms in group.items():
 value=integrate(sum(f*n for f,n in terms.items()))
 values.append(value)
 print(key,value)
assert values==[0]*4

```

## Pole budget

```python
from fractions import Fraction as F
import sympy as s
e,h=s.symbols('e h')
# Bare-parameter recursion for Y source in the flat normalization.
C=1/e+h*(1/(4*e**2)-1/(4*e))+h*h*(1/(6*e**3)-1/(3*e**2)+1/(4*e))
eq=s.expand(-2*e*C+(-2*e*h+3*h*h)*s.diff(C,h)+(h-h*h)*C+2-h+3*h*h/2)
for k in range(3):assert eq.coeff(h,k)==0
# Newton curvature conversion and budget.
assert F(1,6)/6==F(1,36)
assert -F(1,3)/6==-F(1,18)
assert F(1,4)/6==F(1,24)
assert F(1,12)+2*F(-1,24)==0
assert -F(1,16)-F(1,24)==-F(5,48)
assert F(1,72)*2+F(1,144)==F(5,144)
assert -F(1,18)+F(5,48)-F(5,144)==F(1,72)
assert F(48,16*48*3)==F(1,48)
assert F(1,48)+F(1,576)==F(13,576)
assert F(1,36)-F(13,576)==F(1,192)
assert F(1,36)-F(1,48)==F(1,144)
assert F(1,144)/F(1,576)==4
print("Y, counterterm and combined pole checks passed")
# Coordinate numerator for A used one dot per cycle; each full Weyl trace has two.
wick_A=-F(8,64)
assert wick_A==-F(1,8)
assert 2*2*wick_A==-F(1,2)
A3=4*F(1,576); B3=F(1,48)
assert A3+B3==F(1,36)
A2=4*F(5,1152)
B2=F(1,72)-A2
assert B2==-F(1,288)
assert -3*B3==-F(1,16)
assert -3*B2==F(1,96)
assert F(9,2)*B3==F(3,32)
import mpmath as mp
C0=-mp.mpf(15)/2+2*mp.sqrt(3)*mp.im(mp.polylog(2,mp.exp(1j*mp.pi/3)))
ct=(C0+6)/18+(1-mp.zeta(2))/24
subtotal=-.05121964552-.07466151326-.0247+4*.03724051498984436+float(ct)
print("A,B triples and B double target:",A3,B3,B2)
print("CT simple:",ct,"B simple needed at L=0:",1/24-subtotal)

```

- [Martin, Effective potential at three loops](https://arxiv.org/abs/1709.02397).
- [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177).
