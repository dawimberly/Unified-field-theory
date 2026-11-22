# FFff(B) ordered UV sectors and triple-pole check

This continues FFFF_B_PROJECTIVE.md. Ordered projective coordinates cover all 720 regions. Polynomial valuations of the full finite-d numerator give eight singular subsets, including the mixed scalar–fermion pairs, the fermion box, and the five-line regions. Of 720 ordered sectors, 40 contain two logarithmic endpoints and 72 contain one logarithmic plus one power endpoint. Each of the latter needs the linear Taylor coefficient at its power endpoint. The other sectors have at most one pole in the projective integral.

The code derives the trace numerator and the weighted moment numerators symbolically. Direct comparison with the earlier covariance kernel at random finite-d points passes. The double projective residue comes from the 112 sectors with two singular variables. Its integrals, before the factor 1/16 from the four fermion heat kernels, converge to −48 across three quadrature orders. In units h²z/K, applying Γ(3ε−1)≈−1/(3ε), 1/48 from curvature conversion, 1/16 from the four kinetic denominators and the stated signed Wick prefactor +1 gives **a candidate bare triple-pole coefficient +1/48**.

That is below the current budget target 5/192 by 1/192. With FFff(A) the full reported bare triple sum is 13/576 rather than the required 16/576. This discrepancy must be resolved before the simple pole can be interpreted. The independent arithmetic excludes the simple-pole extrapolation as an explanation. Audit the anchor choice in the pole extraction and the common spinor and matter-counterterm normalization before using the result for a projected cancellation claim. The finite-d kernel comparison tests two implementations of the same placement and does not by itself prove correct curvature matching.

## Runnable check

Save the Python block in FFFF_B_PROJECTIVE.md as fffb.py beside this script. Set OPENBLAS_NUM_THREADS=1 and use NumPy, SciPy and SymPy.

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
import numpy as np
from scipy.special import roots_legendre
def residue(n):
 nodes,weights=roots_legendre(n)
 nodes,weights=(nodes+1)/2,weights/2
 grid=np.array(np.meshgrid(nodes,nodes,nodes,indexing="ij")).reshape(3,-1)
 w=np.prod(np.array(np.meshgrid(weights,weights,weights,indexing="ij")),axis=0).ravel()
 total=0.
 for order in permutations(range(6)):
  sets=[tuple(sorted(order[i:])) for i in range(1,6)]
  dat=[table[k] for k in sets]
  singular=[i for i,(a,b,c) in enumerate(dat) if a<0]
  if len(singular)!=2:continue
  regular=[i for i in range(5) if i not in singular]
  def coefficient(poly,derivative=None):
   data=[(tuple(sum(m[i] for i in ids) for ids in sets),float(c)) for m,c in poly.terms()]
   mins=[min(m[i] for m,c in data) for i in range(5)]
   out=np.zeros(grid.shape[1])
   for m,c in data:
    if all(m[i]-mins[i]==(i==derivative) for i in singular):
     out+=c*np.prod(grid**np.array([m[i]-mins[i] for i in regular])[:,None],axis=0)
   return out
  u=coefficient(up); j=coefficient(J4)
  power=next((i for i in singular if dat[i][0]==-2),None)
  if power is not None:
   j=coefficient(J4,power)-5*j*coefficient(up,power)/u
  reg=np.prod(grid**np.array([dat[i][1] for i in regular])[:,None],axis=0)
  total+=w@(reg*j/u**5)/np.prod([dat[i][2] for i in singular])
 return total
for n in (12,24,40):
 r=residue(n)
 print(n,r,"bare triple", -r/2304)
from fffb import kernel
rng=np.random.default_rng(21)
lengths=rng.uniform(.1,2,(10,6))
for dim in (3.2,4.):
 _,a,b=kernel(lengths,dim,0)
 args=lengths.T
 j=s.lambdify((*x,d),J.as_expr(),'numpy')(*args,dim)
 u=s.lambdify(x,U,'numpy')(*args)
 expected=j/(16*np.prod(lengths[:,:4],axis=1)*u**(3+dim/2))
 np.testing.assert_allclose(expected,lengths.sum(axis=1)*a+(5-3*dim/2)*b,rtol=2e-10)
print('Polynomial versus direct covariance passed')

```

- [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177).
