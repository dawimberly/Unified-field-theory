# FFZ CcFf curved projective kernel

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**ASSUMED** This checkpoint records the FFZ CcFf curved projective kernel. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues [the QFf pole](ffz-qff-curvature-pole.md) using the fixed $1/6$ benchmark. This checkpoint develops and checks the CcFf integrand. Its Laurent poles are still open.

1. **DERIVED** Label the cubic vertices 0,1 and the Yukawa vertices 2,3. The scalar edges are 01 twice, 03 and 12. The fermion edges are 23 twice. A cubic pair insertion on the scalar line of the Yukawa sunset gives a positive correction to that propagator, with coefficient $z|y|^{2}/2$. The sunset coordinate factor is $-|y|^{2}$ times the product of the two kinetic vectors. Thus the graph factor is $-z|y|^{4}/2$ with a positive coordinate inner product of the vectors.
2. **DERIVED** A fermion segment of Schwinger length s has a flat kinetic vector $K_s \Delta /(2s), \Delta =x-y$. The Weyl-flattened insertion is $\kappa m(m-\gamma \cdot \partial _x)(m+\gamma \cdot \partial _y)Q$. In the parallel pair only the odd vector component survives the chiral projection.
3. **DERIVED** In Q, the split fraction a has mean $ax+(1-a)y$ and variance 2sa$(1-a)$ per direction. The vector insertion uses $z(\partial _y-\partial _x)Q$. Perform these derivatives before setting a vertex coordinate to zero.
4. **DERIVED** Integrate the split fraction. The kinetic vector is multiplied at order $\kappa$ by $2zs[x^{2}+y^{2}+x\cdot y+(d-1)s]/3$. Relative to the scalar insertions this loses one s in the variance term and has no direct curvature-mass insertion.
5. **DERIVED** Denote coordinate covariance by C. Let $c=Var(\Delta )$ per direction and $v_i=Cov(x_i,\Delta )$. The flat spin factor is $N=dc/(4s_{4}s_{5})$, where s₄,s₅ are the fermion lengths. Isserlis contraction gives $E[\Delta ^{2} x_i\cdot x_j]/E[\Delta ^{2}]=M_{ij}=d C_{ij}+2v_i v_j/c$. A sampled Gaussian independently checks this condition.
6. **DERIVED** For each vertex use $-\lambda _i M_{ii}$ for its Weyl weight, $\lambda _C=(6-d)/2$ and $\lambda _F=(4-d)/2$. For each scalar edge use $-d(4-d)s/3$ for the curvature mass, and retain it before taking $d\to 4$.
7. **DERIVED** The remaining growing-length term is $(2z/3)\sum _e s_e(M_{ii}+M_{jj}+M_{ij}+d_e s_e)$, with $d_e=d$ on a scalar edge and d−1 on a fermion edge. These formulae keep the spin weight in the position moments, including at vertices outside the fermion pair.
8. **DERIVED** Write $s=\rho u$. N scales as $\rho ^{-1}$, the vertex and curvature mass weights as $\rho A$, and the mass-profile weight as $z\rho ^{2}B$. Combining with $\rho ^{k-1} d\rho , k=6-3d/2$, gives the exact radial factor $\Gamma (k) z^{-k}$ and projective numerator $N(u)[A(u)+kB(u)]/U^{d/2}$. The overall coordinate to curvature conversion is $\kappa =R/[4d(d-1)]$. Restore the same three-loop measure as the scalar checkpoint.
9. **DERIVED** The determinant still factors into two parallel-pair totals times the effective perimeter. For pair totals p,q assign q to the fermions. The extra N is proportional to $1/q$ when $q\to 0$ at fixed fractions. Thus the projective q measure behaves as $q^{-2+\epsilon }$, with both its leading and linear terms needed in endpoint subtraction. Applying the scalar pair subtraction alone would be insufficient.
10. **UNVERIFIED** The runnable check tests the split polynomial, conditionally weighted covariance and the flat momentum numerator, then evaluates the finite-d kernel using two different anchors. At $d=1$ it tends toward zero and at $d=1.25$ the anchors agree within about .003 in stripped units. At $d=1.5$ the spread grows to .018, so these exploratory quadratures do not verify high accuracy near a pole. The exact algebra checks are the validation of the integrand at this stage.

**CITED** For the Laurent computation, use analytic subtraction of both the leading and linear fermion-pair face terms and track the overlap with the scalar pair. The flat normalization uses Martin's K_SSSSFF class. The graph's pole and the projected Newton comparison remain open. — eq./page UNVERIFIED.

## Runnable check

**DERIVED** Python 3 with NumPy, SciPy and SymPy.

```python
import sympy as s
import numpy as np
a,x,y,t,d=s.symbols("a x y t d")
mu=a*x+(1-a)*y
H=mu**2+2*d*t*a*(1-a)
term=s.integrate((x-y)*H+2*t*(1-2*a)*mu,(a,0,1))
assert s.factor(term-(x-y)*(x*x+y*y+x*y+(d-1)*t)/3)==0
# Parallel fermion pair: same endpoint difference on both edges.
# Verify the conditional covariance directly by sampling the Gaussian.
rng=np.random.default_rng(122)
R=np.array([[1,.2,.1],[.2,.8,.3],[.1,.3,1.1]])
v=np.array([0,1,-1])
C=R@R.T
c=v@C@v
b=C@v
draw=rng.multivariate_normal(np.zeros(3),C,size=(700000,4))
delta=draw@v
w=np.sum(delta*delta,axis=1)
measured=np.einsum("b,bni,bnj->ij",w,draw,draw)/w.sum()
assert np.max(np.abs(measured-(4*C+2*np.outer(b,b)/c)))<.025
# Independent flat numerator reduction uses a common bubble momentum p.
A1,A2,PP,z1,z2,B=s.symbols("A1 A2 PP z1 z2 B")
# k·(p−k) = [p² - k² - (p−k)²]/2.
bubble=(PP+z1+z2)*B/2-(A1+A2)/2
# In the CcFf convention the reduced bubble is minus twice this form.
assert s.expand(-2*bubble-(A1+A2-(PP+z1+z2)*B))==0
print("Spinor split polynomial, conditioned Gaussian and flat numerator verified")

# CcFf graph and noninteger-dimensional projective integrals.
from scipy.stats import qmc
from scipy.special import ndtri,gamma
ends=np.array([(0,1),(0,1),(0,3),(1,2),(2,3),(2,3)])
def stripped(u,d,root):
    V=np.delete(np.eye(4),root,axis=1)
    incidence=V[ends[:,0]]-V[ends[:,1]]
    Q=np.einsum("ei,be,ej->bij",incidence,1/u,incidence)
    C=2*V@np.linalg.inv(Q)@V.T
    U=np.prod(u,axis=1)*np.linalg.det(Q)
    v=C[:,:,2]-C[:,:,3]
    c=v[:,2]-v[:,3]
    M=d*C+2*v[:,:,None]*v[:,None,:]/c[:,None,None]
    lam=np.array([(6-d)/2]*2+[(4-d)/2]*2)
    a=-np.einsum("i,bii->b",lam,M)-d*(4-d)*u[:,:4].sum(axis=1)/3
    b=np.zeros(len(u))
    for l,(i,j) in enumerate(ends):
        dim=d-(l>=4)
        b+=2*u[:,l]*(M[:,i,i]+M[:,j,j]+M[:,i,j]+dim*u[:,l])/3
    k=6-3*d/2
    N=d*c/(4*u[:,4]*u[:,5])
    return gamma(k)*N*(a+k*b)/U**(d/2)

def estimates(d,power=18,seed=12):
    u=qmc.Sobol(6,scramble=True,seed=seed).random_base2(power)
    u=ndtri((1+u)/2)**2
    u/=u.sum(axis=1,keepdims=True)
    f=np.prod(u,axis=1)**.5*gamma(.5)**6/gamma(3)
    return [sum(np.sum(w*stripped(x,d,root)) for x,w in zip(
                 np.array_split(u,32),np.array_split(f,32)))/len(u)
            for root in (0,2)]
for d0 in (1.,1.25,1.5):
    res=estimates(d0)
    print(d0,res)
```

- **CITED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), Eq. (3.17).
- **CITED** [Bunch and Parker, Feynman propagator in curved spacetime](https://doi.org/10.1103/PhysRevD.20.2499). — eq./page UNVERIFIED.
