# FFZ scalar Newton projection with exact scale integration

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**ASSUMED** This checkpoint records the FFZ scalar Newton projection with exact scale integration. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** The $1/24$ target is a normalization chosen to match Einstein gravity; it is not an independently predicted gravitational constant. A computed $c_Y=1/24$ in the stated matter prescription is a separate conditional calculation. Matching a baseline simple pole to that target remains a test, not a completed derivation.


**ASSUMED** Continues [the cubic insertion checkpoint](ffz-cubic-newton-checkpoint.md). The previous QCc tadpole pole remains $-h^{2}z/(36K\epsilon )$. The open scalar graph calculation reduces to two distinct integrals.

## Identical CCcc integrands

**DERIVED** The two CCcc matchings have identical undirected edge multisets on the same four labeled vertices. Both have double edges 02 and 13, and single edges 03 and 12. Only the assignment of X and Y species differs. In the baseline all scalar propagators have the same mass and curvature coupling, and all four vertices have the cubic Weyl weight. Consequently the line and vertex insertion integrands agree pointwise up to ordering their Schwinger lengths. The combined signed coefficient is $-(5/8)$|y|⁴z². This reduction applies to the baseline. An explicit Y-improvement insertion would distinguish the assignments.

## Projective formula for QCc and CCcc

**ASSUMED** Set $s_e=\rho u_e$ with $\sum u_e=1$. Let $C=2 L_g(u)^{-1}$ padded at the root, and $U=\det L_g(u)$ ∏u_e. With N internal lines, $k=N-3d/2$. Define

- **DERIVED** $A=-d[\sum ᵢ \lambda ᵢ C_{ii} + (4-d)/3]$.
- **DERIVED** $B=(2d/3)\sum _(e=ij) u_e(C_{ii}+C_{jj}+C_{ij}+u_e)$.

**DERIVED** The full linear Weyl insertion at fixed Schwinger lengths is $\rho A+z\rho ^{2}B$. Integration over $\rho$ with measure $\rho ^{k-1} \exp\left(-z\rho \right)$ gives the exact projective numerator $A+(k+1)B$ and radial factor $\Gamma (k+1)z^{-k-1}$. In units $\kappa =R/[4d(d-1)]$, the remaining integral is $\int_{\mathrm{simplex}}$ $U^{-d/2}[A+(k+1)B]$. Restore the three-loop factor $(4\pi )^{-3d/2}$, the dimensional scales of the original checkpoint, and the signed vertex factors. In $d=4-2\epsilon$ the radial gamma is $\Gamma (3\epsilon )$ for QCc and $\Gamma (1+3\epsilon )$ for CCcc. Thus the remaining QCc simple pole needs the finite part of its projective integral after endpoint subtraction, whereas CCcc needs its projective pole residues. A numerical integral from a convergent dimension is insufficient to fix either.

**DERIVED** For subtracting the endpoints, the two parallel pairs of lengths a,b and c,e can be combined exactly at the level of the graph determinant. Their effective lengths are $ab/(a+b)$ and $ce/(c+e)$. The determinant polynomial is $(a+b)(c+e)$ times the sum of the effective lengths and the remaining single-edge lengths. Writing total pair lengths $p_{1}=a+b, p_{2}=c+e$ and fractions $x_{1}=a/p_{1}, x_{2}=c/p_{2}$, the Jacobian contributes p₁p₂. The potentially singular individual-pair factors are $pᵢ^{-1+\epsilon }$ near four dimensions. Their joint soft limit must be treated separately when overlaps or additional degenerations are present. Keep the covariance numerator at finite d for these subtractions.

## Checks

**ASSUMED** Run Python 3 with NumPy and SciPy. The exact $d=1$ cancellation from the preceding position-space calculation is recovered with two independent scrambled quadratures. At $d=2$ both CCcc orderings agree numerically and the relabeling is checked exactly. Printed numbers are stripped integrals after radial integration with $z=1$, including the radial gamma but excluding the loop measure, Wick factors and $\kappa$. The conservative tolerances reflect endpoint sampling.

```python
import numpy as np
from scipy.special import gamma, ndtri
from scipy.stats import qmc
graphs = {
    "QCc": ((0,1),(0,1),(0,2),(0,2),(1,2)),
    "CCcc-A": ((0,2),(0,2),(1,3),(1,3),(0,3),(1,2)),
    "CCcc-B": ((0,2),(0,3),(1,2),(1,3),(0,2),(1,3)),
}
def weight(edges, q, d, u):
    n=max(max(e) for e in edges)+1
    vertices=np.eye(n)[:,1:]
    incidence=np.array([vertices[i]-vertices[j] for i,j in edges])
    lap=np.einsum("ei,be,ej->bij",incidence,1/u,incidence)
    cov=np.pad(2*np.linalg.inv(lap),((0,0),(1,0),(1,0)))
    U=np.prod(u,axis=1)*np.linalg.det(lap)
    lam=np.full(n,(6-d)/2)
    if q: lam[0]=4-d
    A=-d*(np.einsum("i,bii->b",lam,cov)+(4-d)/3)
    B=np.zeros(len(u))
    for e,(i,j) in enumerate(edges):
        B+=2*d*u[:,e]*(cov[:,i,i]+cov[:,j,j]+cov[:,i,j]+u[:,e])/3
    k=len(edges)-3*d/2
    return (A+(k+1)*B)/U**(d/2)

for d in (1.,2.):
    for name in ("QCc","CCcc-A","CCcc-B"):
        edges=graphs[name]
        vals=[]
        for seed in (13,37):
            vals1=[]
            for pts in np.array_split(qmc.Sobol(len(edges),scramble=True,seed=seed)
                                       .random_base2(18),16):
                ss=ndtri((1+pts)/2)**2
                u=ss/ss.sum(axis=1,keepdims=True)
                vals1.append((weight(edges,name=="QCc",d,u)*np.prod(u,axis=1)**.5).mean())
            vals.append(np.mean(vals1)*gamma(len(edges)+1-3*d/2)*gamma(.5)**len(edges)/gamma(len(edges)/2))
        print(d,name,vals)
        if d==1: assert max(map(abs,vals))<.001
        else: assert abs(vals[1]-vals[0])<.15*abs(np.mean(vals))
# Verify the exact relabeling, including the edge multiset.
def shape(edges):
    return sorted(tuple(sorted(edge)) for edge in edges)
perm=(0,1,2,3)
assert shape(graphs["CCcc-A"])==shape([(perm[i],perm[j]) for i,j in graphs["CCcc-B"]])
```

**ASSUMED** Next, subtract the small-pair endpoints with the factorized determinant and finite-d numerator, then verify the Laurent extraction against the QCc tadpole check before using it in the remaining Yukawa and counterterm assembly. The baseline simple pole and its comparison with $1/24$ remain open.

- **CITED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), scalar G and K topologies. — eq./page UNVERIFIED.
