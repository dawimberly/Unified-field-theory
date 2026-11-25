# FFZ Newton vacuum graph inventory and pole-mass response

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**ASSUMED** This checkpoint records the FFZ Newton vacuum graph inventory and pole-mass response. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues [the threshold and spinor checkpoint](ffz-threshold-and-spinor-checkpoint.md) and [the original two-loop benchmark](ffz-two-loop-flat-space-checkpoint.md). The absence of additive matter-source drift through order h³ in the stated fixed $1/6$ prescription carries forward for both Newton constraints.

## Three-loop Newton graph set

**DERIVED** For the ungauged theory at zero scalar background, Q denotes the scalar quartic X²X*², C and c the conjugate cubic pair X²Y* and X*²Y, and F and f the conjugate Yukawa pair $X\psi _X\psi _X$ and $X\,\psi _X$†ψ_X†. Each component carries a coupling magnitude |y|²/4, |my|/2 or |y|/2 respectively. The kinetic fermion propagator is the only fermion line between these vertices, since the supersymmetric mass matrix is off-diagonal and its XX helicity-flip propagator vanishes. A metric probe varying the same off-diagonal mass preserves the X versus Y flavor selection at a single insertion. Scalar and fermion line variations are attached after generating the graph.

**DERIVED** The complete order-|y|⁴ list from Wick matching has ten oriented, field-labeled classes up to relabeling vertices of the same type. The unsigned multiplicities are:

| Multiset | Classes | Multiplicities |
| --- | --- | --- |
| **DERIVED** QQ | 2 | 16, 4 |
| **DERIVED** QCc | 2 | 16, 4 |
| **DERIVED** QFf | 1 | 8 |
| **DERIVED** CCcc | 2 | 8, 32 |
| **DERIVED** CcFf | 1 | 8 |
| **DERIVED** FFff | 2 | 8, 32 |

**CITED** All connected matchings are one-particle irreducible in this model. A bridge would give a one-point subdiagram carrying a nonzero discrete charge. The enumeration separately tests each edge. Ignoring species and directions yields the Martin E, G, J, K and L vacuum types. The H tetrahedron is absent: with four cubic vertices the graph is bipartite between unconjugated and conjugated interaction sets, whereas H contains odd cycles. For a check of the scalar factors, the QQ connected matchings divide into 16 having one self-contraction at each vertex and 4 with all lines across, giving −|y|⁴/2 and −|y|⁴/8 in the Euclidean effective action. The QCc coefficients are +z|y|⁴ and +z|y|⁴/4. The CCcc coefficients are −z²|y|⁴/8 and −z²|y|⁴/2. These include the expansion factorials for repeated monomials. Fermion-containing multiplicities in the table count slot matchings before the signed spin trace. They cannot be used as scalar integral weights. — eq./page UNVERIFIED.

**DERIVED** The projected Newton pole requires the Weyl variation of these graphs plus the one- and two-loop counterterm graphs. The graph list narrows the integration set and identifies two independent checks: the flat sum at zero curvature must vanish in the supersymmetric vacuum and the projected curvature triple and double poles must vanish. Only the projected simple pole tests cancellation of the two-loop-source result with the three-loop vacuum term. The baseline simple-pole coefficient remains unknown.

## Threshold subtraction in the pole-mass variable

**DERIVED** In the original finite kernel let $u=q^{2}/z, L=\ln\left(z/\mu ^{2}\right), A(u)=\int dx \ln [1+u x(1-x)]$ and $B(u)=uA'(u)$. Its $q=0$ subtracted result is $v^{2}[2A(u)+h H(u,L)]/(72K)$. Let $\sigma _{1}=L+A(-1)$. The pole-mass relation is $z_{pole}=z(1+h \sigma _{1}/2)$. Changing variables at fixed q gives

**DERIVED** $H_pole(u_p,L_p)= H(u_p,L_p) + [L_p+A(-1)] B(u_p),$

**DERIVED** where $u_p=q^{2}/z_{pole}$ and $L_p=\ln\left(z_{pole}/\mu ^{2}\right)$. Because B∼$-\pi /$√$(u+4)$, the entire leading inverse-square-root term cancels in H_pole. The explicit scale logarithm is now 3L_p A, since the mass is RG invariant through this order. Its scale derivative cancels the separate running of the X member of the two source variations. This check fixes the sign of the added term.

**DERIVED** The low-momentum coefficients give a second check. At the reference scale the pole-mass order-h slope has bracket $3L_p/2 + c_{slope} + A(-1)/2$, multiplying h v²/(216K z_pole). The coefficient of (q²)² has bracket $3L_p/2 + c_{q4} + A(-1)$, multiplying −h v²/(2160K z_pole²). Thus moving to the pole mass changes both finite Taylor coefficients by known constants and makes the relative log coefficient identical for the two orders. At $L_p=0$ the brackets are approximately −0.40041819 and −1.68552609. The subleading cut remains to be evaluated.

## Runnable verification

**DERIVED** Python 3 standard library. The directed inventory preserves fields for curvature insertion. Each connected graph is checked for the loop count, bridges, and the expected number of equivalence classes.

```python
from itertools import permutations,product
from collections import Counter
# Fields: x,y,f and the corresponding uppercase conjugate slots.
slots={"Q":"xxXX","C":"xxY","c":"XXy","F":"xff","f":"XFF"}
def connected(n,edges):
    seen={0}
    while True:
        more=seen|{v for _,u,v in edges if u in seen}|{u for _,u,v in edges if v in seen}
        if more==seen:
            return len(more)==n
        seen=more
def canonical(vertices,edges):
    orders=(order for order in permutations(range(len(vertices)))
            if tuple(vertices[i] for i in order)==vertices)
    return min(tuple(sorted((s,order[u],order[v]) for s,u,v in edges)) for order in orders)
def graphs(vertices):
    families=[]
    for species in ("x","y","f"):
        src=[i for i,v in enumerate(vertices) for _ in range(slots[v].count(species))]
        snk=[i for i,v in enumerate(vertices) for _ in range(slots[v].count(species.upper()))]
        if len(src)!=len(snk):
            return {}
        families.append([tuple((species,u,v) for u,v in zip(src,perm))
                         for perm in permutations(snk)])
    out=Counter()
    for e in product(*families):
        edges=sum(e,())
        if not connected(len(vertices),edges):
            continue
        assert all(connected(len(vertices),edges[:i]+edges[i+1:]) for i in range(len(edges)))
        assert len(edges)-len(vertices)+1==3
        out[canonical(vertices,edges)]+=1
    return out
sets=("QQ","QCc","QFf","CCcc","CcFf","FFff")
expected={"QQ":20}
for vertices in sets:
    results=graphs(tuple(vertices))
    print(vertices,len(results),sum(results.values()))
    for edges,multiplicity in sorted(results.items()):
        counts=Counter(tuple(sorted((u,v))) for _,u,v in edges)
        loops=sum(n for (u,v),n in counts.items() if u==v)
        shape=loops,tuple(sorted(n for (u,v),n in counts.items() if u!=v))
        print(" ",shape,multiplicity,edges)
    if vertices in expected:
        assert sum(results.values())==expected[vertices]
checks = {"QQ": [4,16], "QCc": [4,16], "QFf": [8],
          "CCcc": [8,32], "CcFf": [8], "FFff": [8,32]}
assert sum(len(groups) for groups in checks.values()) == 10
for multiset, expected in checks.items():
    assert sorted(graphs(tuple(multiset)).values()) == expected
from math import pi, sqrt
on_shell = pi/sqrt(3)-2
assert abs(-.3073178731149105+on_shell/2+.4004181909978016)<1e-12
assert abs(-1.4993254559801423+on_shell+1.6855260917459245)<1e-12
for n in (1,2):
    assert (3-n)/2+n/2==1.5
```

## Next calculation

**DERIVED** Assign the signed Weyl spin traces to the four Yukawa-containing classes, generate the curvature numerator polynomials for all classes with the Gaussian graph covariance, and reduce their integrals and counterterms. The fixed-scheme test remains equality of the baseline simple pole with the Y improvement coefficient $1/24$. For the response, evaluate the remaining two-particle cut in the pole-mass representation with its leading threshold displacement removed.

## Sources

- **CITED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), Fig. 3.1 and the scalar and fermion vacuum terms. — eq./page UNVERIFIED.
- **CITED** [Martin, Three-loop effective potential for softly broken supersymmetry](https://doi.org/10.1103/PhysRevD.109.015019), quadratic mass basis and counterterm prescription. — eq./page UNVERIFIED.
