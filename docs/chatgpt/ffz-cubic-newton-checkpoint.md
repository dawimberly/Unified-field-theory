# FFZ cubic-sector Newton insertion

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**ASSUMED** This checkpoint records the FFZ cubic-sector Newton insertion. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues [the quartic checkpoint](ffz-quartic-newton-checkpoint.md). Claude independently reproduced the quartic Schwinger normalization and numerically confirmed the simple-pole constant. The work here stays in the original fixed $1/6$ dimensional continuation.

## QCc with a tadpole

**DERIVED** The local scalar tadpole has zero linear curvature response at $\xi =1/6$, since its first heat-kernel coefficient is $(1/6-\xi )R$ identically in d. Work in the covariant frame, where the flat tadpole value T is constant at this order. The tadpole graph is therefore a mass derivative of the cubic sunset, taken at fixed cubic coupling. The flat Wick coefficient of the sunset is −z|y|²/2 and that of this QCc graph is +z|y|⁴. A single-line mass derivative inserts −G². Thus the QCc graph is obtained by applying 2|y|²T times the *single-line* mass derivative to the cubic sunset. At equal masses the single-line derivative is one third of the common mass derivative, including in the linear curvature term, by symmetry of the three scalar kernels with common curvature coupling.

**DERIVED** The cubic sunset curvature term is −z|y|² $J/(24d)$ with $J=\int dᵈr r^{2}G^{3}$. In taking the derivative, the z in the prefactor comes from the cubic vertex and must be held fixed. The resulting QCc term is −z|y|⁴T J'/(36d). Here J scales as $z^{d-4}$, and its residue is $J=2/(K^{2}\epsilon )+O(1)$. This follows directly from the short-distance $G=1/(4\pi ^{2}r^{2})$ at $d=4$ and the radial measure $r^{-1+4\epsilon }$. Then $T=-z/(K\epsilon )+O(1)$ implies the bare QCc tadpole simple pole

**DERIVED** $F_R = -h^2 z/(36K\epsilon ) + O(1).$

**DERIVED** The factor (d−4) from J' removes a possible double pole. The derivation uses the covariant tadpole property, so no extra Weyl weight for the coincident loop is multiplied into this differentiated covariant graph.

## General scalar insertion and independent check

**DERIVED** Let L_g be the reduced conductance Laplacian with edge conductances $1/s_e$, with one interaction vertex fixed. The coordinate covariance per Cartesian direction is $C=2 L_g^{-1}$, padded with zeros for the fixed vertex. For a split scalar line ij of total Schwinger length s, the conditional coordinate average at split fraction a is $a x_i+(1-a)x_j$ and the variance per coordinate is 2sa(1−a). The curvature insertion relative to the flat Gaussian integrand, in units $\kappa =R/[4d(d-1)]$, is therefore

**DERIVED** $2zd s(C_{ii}+C_{jj}+C_{ij}+s)/3 - d(4-d)s/3$.

**DERIVED** Sum this over lines and add the vertex term $-d \sum \lambda _i C_{ii}$. Under the Weyl field rescaling, $\lambda _Q=4-d$ and $\lambda _C=(6-d)/2$ from the number of scalar fields in each monomial. This expression retains the evanescent curvature mass term. The signed flat Wick weights, after extracting z|y|⁴ for QCc and z²|y|⁴ for CCcc, are $+1,+1/4$ and $-1/8,-1/2$ respectively. They follow from the slot multiplicities 16,4 and 8,32, the monomial factors $1/4$ and $1/2$, the expansion sign of −ln Z, and the 2!² denominator for two repeated conjugate cubic pairs.

**DERIVED** The code independently integrates each scalar graph in position space at $d=1, z=1$, using $G(r)=\exp(-|r|)/2$ and ordered coordinate gaps. It checks the full Gaussian implementation including the insertion against exact rational position integrals. At $d=1 \kappa$ is used as the Weyl probe amplitude, without dividing by d−1. All four variation integrals vanish in that limit. This is a normalization and cancellation check in a convergent dimension, with limited power to test the d-dependent pole. The QCc tadpole Laurent residue above is obtained by reduction and short-distance analysis. The other three Laurent residues require remaining projective endpoint subtractions.

## Runnable check

**DERIVED** Python 3, NumPy and SciPy.

```python
from itertools import permutations
from fractions import Fraction as F
import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc

graphs = (
 ("QCc tadpole", ((0,0),(0,1),(0,2),(1,2),(1,2)), ("Q","C","C"), F(1)),
 ("QCc doubled", ((0,1),(0,1),(0,2),(0,2),(1,2)), ("Q","C","C"), F(1,4)),
 ("CCcc parallel", ((0,2),(0,2),(1,3),(1,3),(0,3),(1,2)), ("C",)*4, -F(1,8)),
 ("CCcc crossed", ((0,2),(0,3),(1,2),(1,3),(0,2),(1,3)), ("C",)*4, -F(1,2)),
)

def position(edges,n,mark=None):
    total=F(0)
    for order in permutations(range(n)):
        rates=[sum((i in order[:k])!=(j in order[:k]) for i,j in edges)
                for k in range(1,n)]
        fac=F(1,2**len(edges))
        for r in rates: fac/=r
        if mark is not None:
            a,b=sorted((order.index(0),order.index(mark)))
            gaps=[F(1,r) for r in rates[a:b]]
            fac*=sum(gaps)**2+sum(g*g for g in gaps)
        total+=fac
    return total

def exact(edges, kinds):
    n=len(kinds)
    geom=-sum(F(3 if k=="Q" else 5,2) * position(edges,n,i)
              for i,k in enumerate(kinds) if i)
    for k,(i,j) in enumerate(edges):
        split=edges[:k]+edges[k+1:]+((i,n),(n,j))
        geom += 2*position(split,n+1,n)-position(split,n+1)
    return position(edges,n),geom

def gaussian(edges,kinds,power=19):
    n=len(kinds)
    u=qmc.Sobol(len(edges),scramble=True,seed=19).random_base2(power)
    s=ndtri((1+u)/2)**2/2
    V=np.eye(n)[:,1:]
    incidence=np.array([V[i]-V[j] for i,j in edges])
    mat=np.einsum('ei,be,ej->bij',incidence,1/s,incidence)
    cov=np.pad(2*np.linalg.inv(mat),((0,0),(1,0),(1,0)))
    base=2.**(-len(edges))*(4*np.pi)**((n-1)/2)/np.sqrt(np.linalg.det(mat))
    geom=-sum((1.5 if k=="Q" else 2.5)*cov[:,i,i] for i,k in enumerate(kinds))
    for l,(i,j) in enumerate(edges):
        geom+=2*s[:,l]*(cov[:,i,i]+cov[:,j,j]+cov[:,i,j]+s[:,l])/3-s[:,l]
    return np.array([base.mean(),(base*geom).mean()])

for name,edges,kinds,wick in graphs:
    flat,w=exact(edges,kinds)
    num=gaussian(edges,kinds)
    assert np.allclose(num,[float(flat),float(w)],rtol=.07,atol=.00035),(name,num,flat,w)
    print(name, 'flat',flat,'Weyl variation times Wick factor',wick*w, 'QMC raw',num)
assert -F(1,36)*F(1,4)*(-1)*(-2)*2 == -F(1,36)
```

**ASSUMED** Next, fix the three remaining scalar Laurent series by integrating the common scale and subtracting ordered-sector endpoints, using the tadpole reduction as a nontrivial check against the direct QCc implementation before extracting the other residues. The Yukawa insertions and matter counterterms still enter the projected full simple pole.

- **CITED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), Eqs. (3.5)–(3.9).
- **CITED** [Schmidt and Schubert, Worldline Green functions for multiloop diagrams](https://doi.org/10.1016/0370-2693(94)$90944-X$). — eq./page UNVERIFIED.
