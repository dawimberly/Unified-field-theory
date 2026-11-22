# FFff(B) full regulated extraction

This batch performs 28 regulated evaluations (seven values of ε including zero, four sampling/cutoff/anchor configurations) following the ordered-sector and normalization audits. The independent numerical extraction does not use the predicted double or simple pole.

Use the finite-dimensional polynomial J and U from FFFF_B_ORDERED_SECTORS.md. In each of 720 ordered sectors the normalized integrand is x^a U_reduced^(−5+ε) W^(−3ε) J_reduced. The endpoint exponents are a_i+ε u_i. Subtract a constant at each logarithmic axis and a constant and a linear term at each power axis. Add them back analytically using c^(a_i+ε u_i+n+1)/(a_i+ε u_i+n+1), n=0 or 1, where c is the subtraction cutoff. Include the overlap by applying the subtraction operator separately to each axis. The code evaluates ε² times the projective integral, treating every 1/ε analytically so the value at zero is well defined. Use positive and negative ε, with the regularized remainder convergent in the tested interval.

## Numerical results

The polynomial-fit coefficients of ε² I at powers 0,1,2 are:

| Anchor | Cutoff | Samples per sector | Intercept | Linear | Quadratic |
| --- | --- | --- | --- | --- | --- |
| 0 | 1 | 4096 | −48.00073 | 223.98233 | 191.77900 |
| 0 | 0.5 | 4096 | −47.99655 | 223.99751 | 191.47522 |
| 0 | 0.75 | 4096 | −47.99853 | 224.05950 | 190.46241 |
| 1 | 1 | 2048 | −48.00818 | 224.04255 | 193.18188 |

The exact intercept −48 was checked earlier. The linear coefficient clusters near 224, while the pole budget would require 208. In units h²z/K at L=0, multiply by exp(3γ_E ε)Γ(3ε−1)/[64d(d−1)]; this yields a double coefficient near −0.01042, consistent numerically with −1/96. The required −1/288 differs by −1/144. The numerical simple coefficient near −0.052 is provisional because the summed higher poles still fail. The observed discrepancy is much larger than the spread in these runs. These differences are diagnostic and do not provide rigorous error bounds.

With the current scalar, CcFf and corrected FFff(A) entries, the projected triple pole cancels but the projected double pole has a residual around −1/144. This prevents interpreting the simple pole as a failure of physical matching. The next audit should independently extract the finite-d first subleading pole for A, verify other bare double poles and check evanescent matter counterterms in this continuation. A successful top-pole check does not exclude errors at the next order.

## Reproduce

Save as b_full.py and run with OPENBLAS_NUM_THREADS=1 and arguments 12 23 1, 12 41 .5, or 12 61 .75. The anchor-1 run uses V=s.eye(4)[:,[0,2,3]] and arguments 11 23 1.

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
from scipy.special import gamma
from scipy.stats import qmc
from itertools import product
import sys
power=int(sys.argv[1]) if len(sys.argv)>1 else 11
root_seed=int(sys.argv[2]) if len(sys.argv)>2 else 23
cut=float(sys.argv[3]) if len(sys.argv)>3 else 1.
eps=np.array([.0,.005,.01,.02,-.005,-.01,-.02])
# We integrate the Laurent components after multiplying by eps^2.
# Include negative eps for symmetric coefficient estimates.
coef=s.lambdify(d,[c for m,c in J.terms()],"numpy")
coeff=np.array([coef(4-2*e) for e in eps]).T
orders=[]
for order in permutations(range(6)):
 sets=[tuple(sorted(order[i:])) for i in range(1,6)]
 dat=np.array([table[k] for k in sets])
 def exponents(poly):
  m=np.array(poly.monoms())
  out=np.stack([m[:,ids].sum(axis=1) for ids in sets],axis=1)
  return out-out.min(axis=0)
 orders.append((dat,exponents(up),exponents(J)))
pts=qmc.Sobol(5,scramble=True,seed=root_seed).random_base2(power)
def sector(dat,U,J):
 ax=[i for i in range(5) if dat[i,0]<0]
 if not ax:return None
 # weights, regular coordinates sampled uniformly. Integrate over a short
 # (0,cut) interval and choose ordinary Taylor add-backs analytically.
 D=dat[:,0,None]+dat[:,2,None]*eps
 def evaluate(zeros,deriv):
  arrays=[]
  vals=[]
  for ex,co in ((U,np.ones((len(U),1))),(J,coeff)):
   keep=np.ones(len(ex),dtype=bool)
   for i in zeros:keep&=ex[:,i]==int(i==deriv)
   p=np.ones((len(pts),keep.sum()))
   for i in range(5):
    if i not in zeros:p*=pts[:,i,None]**ex[keep,i][None,:]
   vals.append(p@co[keep])
  return vals
 def W(zeros,deriv=None):
  out=np.zeros(len(pts))
  for n in range(6):
   ids=set(range(n))
   if ids.intersection(zeros)!=({deriv} if deriv is not None else set()):continue
   m=np.ones(len(pts))
   for i in ids-set(zeros):m*=pts[:,i]
   out+=m
  return out[:,None]
 cache={}
 def F(zeros,deriv):
  key=(tuple(zeros),deriv)
  if key in cache:return cache[key]
  u,j=evaluate(zeros,None)
  w=W(zeros)
  b=u**(-5+eps)*w**(-3*eps)
  a=j
  if deriv is not None:
   du,dj=evaluate(zeros,deriv)
   a=dj+j*((-5+eps)*du/u-3*eps*W(zeros,deriv)/w)
  cache[key]=a*b
  return cache[key]
 result=np.zeros(len(eps))
 for symbol in product(*(range(1,2-dat[i,0]) for i in ax)):
  done=[i for i,j in zip(ax,symbol) if j>1]
  remaining=[i for i in ax if i not in done]
  der=next((i for i,j in zip(ax,symbol) if j==3),None)
  M=F(done,der).copy()
  for b in product(*(range(1-dat[i,0]) for i in remaining)):
   rem=[i for i,j in zip(remaining,b) if j]
   if not rem:continue
   newder=next((i for i,j in zip(remaining,b) if j==2),der)
   # only one power endpoint per sector
   v=np.prod([((pts[:,i]<cut)*pts[:,i]**(j-1)) if j else np.ones(len(pts)) for i,j in zip(remaining,b)],axis=0)
   M+=(-1)**len(rem)*F(sorted(done+rem),newder)*v[:,None]
  fac=np.ones(len(eps))
  singular_count=0
  for i,j in zip(ax,symbol):
   if j>1:
    expo=D[i]+j-1
    if dat[i,0]+j-1==0:
     fac*=cut**expo/dat[i,2]
     singular_count+=1
    else:fac*=cut**expo/expo
  reg=[i for i in range(5) if i not in done]
  wt=np.prod(pts[:,reg,None]**D[reg][None,:,:],axis=1)
  result+=np.mean(wt*M,axis=0)*fac*eps**(2-singular_count)
 return result
total=np.zeros(len(eps))
for pre in orders:
 out=sector(*pre)
 if out is not None:total+=out
# separately integrate regular sectors, small but required for simple pole
for dat,U,J in orders:
 if np.any(dat[:,0]<0):continue
 def poly(ex,c):
  q=np.ones((len(pts),len(ex)))
  for i in range(5):q*=pts[:,i,None]**ex[None,:,i]
  return q@c
 u=poly(U,np.ones((len(U),1)));j=poly(J,coeff)
 w=1+np.cumprod(pts,axis=1).sum(axis=1,keepdims=True)
 D=dat[:,0,None]+dat[:,2,None]*eps
 total+=eps**2*np.mean(j*u**(-5+eps)*w**(-3*eps)*np.prod(pts[:,:,None]**D[None,:,:],axis=1),axis=0)
print('eps^2 projective',list(zip(eps,total)),flush=True)
print('fit ascending',np.polynomial.polynomial.polyfit(eps,total,4),flush=True)

```

- [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177).
