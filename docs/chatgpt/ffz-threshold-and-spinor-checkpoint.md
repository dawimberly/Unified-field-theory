# FFZ pole-mass threshold and curved spinor insertion

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**DERIVED** This checkpoint records the FFZ pole-mass threshold and curved spinor insertion. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Continues [the subtracted momentum kernel](ffz-subtracted-momentum-kernel.md) and [the Newton Y-improvement projection](ffz-three-loop-y-improvement-checkpoint.md). The matter beta and the two Newton moments retain their prior scheme-qualified results.

## Spinor line for the remaining curved vacuum residue

**DERIVED** Use the original Weyl frame σ(x)=−κx², κ=R/[4d(d−1)]. In the rescaled Dirac operator only m e^σ varies. Let x,y be the endpoints and let Q(x,y)=∫du G(x−u)u²G(u−y), with G the flat scalar kernel. The first-order Dirac resolvent insertion is exactly

**DERIVED** δS/κ = m(m−γ·∂_x)(m+γ·∂_y)Q.

**DERIVED** The derivatives act before imposing endpoints. For Schwinger segment lengths s₁,s₂, write s=s₁+s₂, a=s₂/s, b=s₁/s. The integrand of Q after integrating the insertion position is the full free heat kernel K_s(x−y) multiplied by (a x+b y)²+2d s a b. This formula supplies all scalar, vector, and antisymmetric gamma components of the insertion on an edge with arbitrary endpoints. In the graph representation the remaining Wick contractions reduce to Gaussian vertex moments. Use the chiral projections and spin multiplicities of the flat vacuum diagrams.

**DERIVED** An independent normalization follows at x=y=0. The vector term vanishes, and the mixed derivative contributes −d(d+2)ab K_s. Integrate ab over the split fraction to get 1/6 and integrate by parts in s. The result is δS(x,x)|₀/R=−m∫dp G(p)²/12, agreeing with the coincident Dirac heat-kernel coefficient R/6−R/4. In momentum space at y=0 the convolution gives δS/κ=−m SΔ_p S, as in the original two-loop benchmark. The coincident and noncoincident checks use different projections of this same insertion.

**DERIVED** The remaining Newton computation is a three-loop vacuum graph sum with these line operators and the scalar/vertex operators from the preceding checkpoint. It must cancel the known triple and double poles at a_X=0,a_Y=−1. The improvement-independent simple pole is still unknown, so cancellation against the Y value 1/24 remains an open test.

## On-shell threshold from the flat matter sector

**DERIVED** Set z=m², h=|y|²/(16π²), L=ln(z/μ²), and A(u)=∫₀¹ dx ln[1+u x(1−x)]. The scalar Y couples via the cubic term m* y Y* X²/2 plus conjugate. Its one-loop inverse propagator correction is minus |my|²/2 times the subtracted equal-mass bubble. The factor 1/2 is from the two identical X contractions. Thus the common mass-shell location (the multiplets are degenerate in this flat supersymmetric prescription) is

**DERIVED** $z_pole = z [1+(\frac{h}{2})(L+A(-1))]+O(h^2), with A(-1)=\pi /√3-2.$

**DERIVED** Its RG derivative vanishes through order h. The supersymmetric fermion check gives the same formula: only X has a Yukawa wavefunction correction, the off-diagonal supersymmetric mass stays in the quadratic term, and the pole condition uses one corrected X kinetic factor and the uncorrected Y factor. The two-particle threshold in invariant mass squared is 4 z_pole. The perturbation in the location is nonuniform when |4z−s|/z is of order h.

## Independent extraction from the finite curvature response

**DERIVED** The preceding checkpoint writes the subtracted response as h v² H(u,L)/(72K), with u=q²/z and v=t w_g. Put δ=u+4>0. The leading inverse-square-root coefficient is

**DERIVED** lim_(δ→0+) √δ H(−4+δ,L)/π = L+A(−1).

**DERIVED** It follows directly from the corner of the finite kernel. In either ordered sector with a=1 set b=x,c=x t and mark the long line at separation l. Near l=1/2 and x=0, the denominator is P+kx, where P=δ/4+4(l−1/2)² and k=1+t−t/(1+t). The singular same-line combination in J is the integral of [1/P−1/(P+kx)]/x. Integrating x gives ln(k/P)/P to leading order. The endpoint logarithm and the pair sectors do not furnish a leading 1/√δ: they have an extra radial integration or a vanishing short-edge source measure. Away from this corner, any nearly maximal resistance requires both a small auxiliary-edge length and a source location tuned to the antipode. Integrating these directions with a smooth graph measure gives a term subleading relative to δ^(−1/2).

**DERIVED** The term −D of the finite response contributes +ln P/P and cancels the logarithmically enhanced leading singularity. Including the explicit (1−ln 2)(B−A) and L(3A−B) terms leaves the coefficient L−1+ln 2+2∫₀¹ dt ln(k)/(1+t)². Change variables y=t/(1+t), giving k=[1−y(1−y)]/(1−y). Reflection symmetry of the first logarithm and direct integration of the second give A(−1)+1−ln 2 for the angular average. Finally ∫dl/P has leading coefficient π/√δ. The result is the same as expanding the one-loop form factor about z_pole, since its mass derivative at fixed q is −B and B∼−π/√δ.

**DERIVED** The signs matter. At L=0 the limiting coefficient is negative, π/√3−2≈−0.18620. Moderate-distance quadrature can still give positive values for √δ H/π because the subleading terms are sizeable. That quadrature is not a verification of the asymptotic coefficient. The analytical corner extraction and the independent mass-shell calculation supply the check. The full above-threshold discontinuity and its nonsingular vertex corrections remain open.

## Runnable checks

**DERIVED** Python 3 with SciPy and SymPy.

```python
import math
import sympy as s
from scipy.integrate import quad
d = s.symbols("d")
spinor = s.simplify((2*(2-d/2)-(d+2))*d/(6*4*d*(d-1)))
assert spinor == -s.Rational(1,12)
A = quad(lambda y: math.log1p(-y*(1-y)),0,1)[0]
angular = 2*quad(lambda t: math.log(1+t-t/(1+t))/(1+t)**2,0,1)[0]
assert abs(A-(math.pi/math.sqrt(3)-2))<1e-13
assert abs(angular-1+math.log(2)-A)<1e-13
for delta in (1e-3,1e-5,1e-7):
    u=-4+delta
    minus_B=-quad(lambda y: u*y*(1-y)/(1+u*y*(1-y)),0,1,
                  epsabs=1e-8,points=[.5])[0]
    assert abs(math.sqrt(delta)*minus_B/math.pi-1)<math.sqrt(delta)
print("Spinor coincidence, angular threshold coefficient and bubble limit passed")
```

## Next calculation

**DERIVED** Generate the flat Wess–Zumino three-loop component vacuum graphs, attach a single quadratic Weyl source to each line or vertex, project to a_X=0,a_Y=−1, and test the higher-pole cancellation before interpreting the simple pole. For the momentum branch, express the form factor in the pole-mass variable and compute the finite two-particle cut after removing the shift term.

## Sources

- **DERIVED** [Bunch and Parker, Feynman propagator in curved spacetime](https://doi.org/10.1103/PhysRevD.20.2499).
- **DERIVED** [Martin, Two-loop scalar self-energies in a general renormalizable theory at leading order in gauge couplings](https://doi.org/10.1103/PhysRevD.70.016005).
