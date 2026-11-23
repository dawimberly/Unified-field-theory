# FFff(B) full projective curvature kernel

This implements the earlier FFff(B) Gaussian reduction for the actual graph. In the edge order 02,12,13,03,02,13 the first four lines are fermions; the last two are scalars parallel to the first and third fermion lines. The standard signed coefficient for the trace divided by two is +|y|⁴. This follows from the 32 Wick matchings, the four vertex factors, 2!², the minus sign of the vacuum cumulant, and the negative signed trace. The flat equal-mass momentum reduction is the K topology of Martin (3.19).

For each projective edge tuple the code evaluates the reduced graph Laplacian, then obtains both the unweighted trace and the Weyl-weighted moments by exact polynomial polarization. The linear curvature numerator is A+(k−1)B, k=6−3d/2, and the radial factor is Γ(k−1). The scalar mass insertions and the Yukawa Weyl weights remain at finite d. At z=1 restore the loop measure (4π)^(−3d/2), |y|⁴, and curvature conversion 1/[4d(d−1)] when extracting the bare Laurent coefficients. The checks use the stripped values and do not depend on that final overall sign and normalization.

In d=1, G=e^(−|r|)/2, each flat kinetic vector is sign(r)G, and the trace divided by two is the product of the four one-dimensional signs. Summing the 24 position orders by the exact exponential gap integral gives 1/432. The independent Schwinger integral gives approximately 0.103111867 before restoring (4π)^(−3/2), consistent with (4π)^(3/2)/432. Both anchors reproduce the same flat result and give a near-zero curvature insertion at d=1. At d=1.2 they give 0.01109276 and 0.01109361 for the stripped curvature term. At d=1.4 they give 0.02885908 and 0.02885793. The Weyl-probe checks use κ directly and never divide by d−1 at d=1.

The short-subgraph inventory finds two mixed scalar–fermion pair subgraphs and a four-fermion box. Several larger subsets have nonnegative superficial degree as well. The inventory deliberately includes disconnected or reducible candidates. For pole extraction the subtraction must track nested and overlapping UV regions and use the actual numerator scaling to remove superficial candidates with zero leading coefficient. The convergent-dimension checks do not evaluate the Laurent coefficients near d=4. The independent FFff(B) targets are still the higher-pole budget 5/192 and 11/1152−5L/64.

## Runnable check

Python 3 with NumPy/SciPy. Use OPENBLAS_NUM_THREADS=1.

```python
import numpy as np
from scipy.special import gamma, ndtri
from scipy.stats import qmc
# Fermion trace order: 02,12,13,03, all arrows F to f.
edges=np.array([(0,2),(1,2),(1,3),(0,3),(0,2),(1,3)])
def kernel(s,d,root):
    V=np.delete(np.eye(4),root,axis=1)
    W=V[edges[:,0]]-V[edges[:,1]]
    lap=np.einsum('ei,be,ej->bij',W,1/s,W)
    C=2*V@np.linalg.inv(lap)@V.T
    U=np.prod(s,axis=1)*np.linalg.det(lap)
    D=np.eye(4)[edges[:4,0]]-np.eye(4)[edges[:4,1]]
    c=D@C@D.T
    def trace(M):
        ab,ac,ad,bc,bd,cd=[M[:,i,j] for i,j in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))]
        return d*d*(ab*cd-ac*bd+ad*bc)+2*d*ac*bd
    n=trace(c)
    cv=C@D.T
    # Trace is quadratic in c. Polarize it to take its exact derivative.
    def moment(i,j):
        x,y=cv[:,i,:],cv[:,j,:]
        dc=x[:,:,None]*y[:,None,:]+y[:,:,None]*x[:,None,:]
        return d*C[:,i,j]*n+trace(c+dc)-n-trace(dc)
    diag=[moment(i,i) for i in range(4)]
    A=-(4-d)/2*sum(diag)-d*(4-d)/3*n*s[:,4:].sum(axis=1)
    B=np.zeros(len(s))
    for k,(i,j) in enumerate(edges):
        B+=2*s[:,k]*(diag[i]+diag[j]+moment(i,j)+(d-(k<4))*s[:,k]*n)/3
    pref=1/(16*np.prod(s[:,:4],axis=1)*U**(d/2))
    return pref*n,pref*A,pref*B
def integrate(d,seed,power):
    k=6-3*d/2
    sob=qmc.Sobol(6,scramble=True,seed=seed)
    totals=np.zeros((2,2))
    for _ in range(2**(power-13)):
        u=ndtri((1+sob.random(2**13))/2)**2
        u/=u.sum(axis=1,keepdims=True)
        w=np.prod(u,axis=1)**.5*gamma(.5)**6/gamma(3)
        for idx,root in enumerate((0,2)):
            n,A,B=kernel(u,d,root)
            totals[idx]+=np.array([gamma(k-2)*np.mean(w*n),gamma(k-1)*np.mean(w*(A+(k-1)*B))])
    return totals/2**(power-13)
if __name__=='__main__':
    rng=np.random.default_rng(5)
    u=rng.uniform(.05,1,(20,6))
    for d in (1.,2.5,3.8):
        a=kernel(u,d,0)
        b=kernel(u,d,2)
        assert np.allclose(a[0],b[0],rtol=1e-12)
        assert np.allclose(np.array(kernel(2*u,d,0)),np.array(a)*np.array([2**(-2-1.5*d),2**(-1-1.5*d),2**(-1.5*d)])[:,None])
    from fractions import Fraction as F
    # At d=z=1, direct position ordering gives flat trace/2 integral 1/432.
    for d in (1.,1.2,1.4):
        print(d,integrate(d,12,19),flush=True)
    from itertools import permutations,combinations
    exact=F(0)
    for order in permutations(range(4)):
        pos={v:i for i,v in enumerate(order)}
        term=F(1,64)
        for i,j in edges[:4]: term*=1 if pos[i]>pos[j] else -1
        for cut in (1,2,3):
            term/=sum((pos[i]<cut)!=(pos[j]<cut) for i,j in edges)
        exact+=term
    assert exact==F(1,432)
    z=integrate(1.,12,19)
    np.testing.assert_allclose(z[:,0]/(4*np.pi)**1.5,float(exact),rtol=.01)
    assert np.max(abs(z[:,1]))<.00015
    # Candidate proper UV sectors before spin cancellations.
    for n in range(2,6):
        for subset in combinations(range(6),n):
            es=edges[list(subset)]
            vertices=np.unique(es)
            rank=np.linalg.matrix_rank((np.eye(4)[es[:,0]]-np.eye(4)[es[:,1]]))
            loops=n-rank
            deg=4*loops-sum(1+(k>=4) for k in subset)
            if loops and deg>=0:
                print(subset,loops,deg)

```

- [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), Eq. (3.19).
