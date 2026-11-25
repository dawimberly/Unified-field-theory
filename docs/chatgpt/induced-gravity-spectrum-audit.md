# Induced gravity: spectrum balance and its vacuum-energy cost

**Conversation date range** 2026-10 (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records induced gravity: spectrum balance and its vacuum-energy cost. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

## Result

**ASSUMED** A mass-matched free-field spectrum can cancel its regulated flat-space vacuum determinant while leaving a positive Einstein-Hilbert coefficient. The simplest example in the scanned class contains two real scalars and one neutral two-component Majorana fermion. This is a known boson/fermion balance, not a new particle model.

**ASSUMED** The cancellation does not remove higher-curvature terms and is sensitive to mass splitting unless a suitable symmetry protects it. A separate calculation using a known O'Raifeartaigh matter model produces a small decrease of the induced heavy-sector gravitational coupling, but also generates vacuum energy. At the tested parameters the larger decreases fail a necessary low-curvature consistency check. No observed gravitational constant, matter spectrum, stable cosmic vacuum, or microscopic origin of spacetime has been predicted.

## Fixed conventions and assumptions

**ASSUMED** Use hbar=c=1 and a common proper-time ultraviolet cutoff Lambda_UV=1. On a supplied smooth Euclidean metric, the one-loop scalar convention is W=+1/2 Tr log D, with D=-laplacian+m^2+xi R and positive R for a round sphere. The Euclidean Einstein-Hilbert term is -B integral sqrt(g) R, where B=1/(16 pi G)=Mbar^2/2.

**ASSUMED** Matter loops are evaluated at one loop. Gravitational loops, ultraviolet completion, metric emergence and a supergravity completion are absent. A chosen curvature coupling is not determined by flat-space particle data. Bare gravitational/vacuum contributions are set to zero in the illustrative examples; they could change the total coupling. Outputs are regulated coefficients, not automatically unique renormalized observables.

**ASSUMED** The scalar and fermion induced-Newton coefficients follow [1]. General induced-gravity finiteness and vacuum-energy issues are discussed in [2]. Curvature-squared coefficients below are the free-field heat-kernel magnitudes after statistics and gauge ghosts are included. Euler-density and total-derivative terms are omitted from the local bulk comparison at fixed topology without a boundary.

## Coefficient scan

**ASSUMED** Let Ns count real minimal scalars, Nf neutral Majorana fermions, and Nv massless gauge vectors. The relevant ultraviolet coefficients are

**ASSUMED** $a0=Ns-2Nf+2Nv,$

**ASSUMED** $a1=Ns(\frac{1}{6}-xi)+\frac{Nf}{6}-\frac{2Nv}{3},$

**ASSUMED** $aC=\frac{Ns}{120}+\frac{Nf}{40}+\frac{Nv}{10},$

**ASSUMED** $aR=(\frac{Ns}{2})(\frac{1}{6}-xi)^2.$

**ASSUMED** Here a0 controls the leading flat vacuum term, a1 the Einstein term, and aC,aR the displayed curvature-squared terms. The vector expression includes its gauge/ghost contributions and must not be used for a massive vector without modification.

**ASSUMED** The ranges Ns=0..24 and Nf,Nv=0..12 contain 4,224 nonempty tuples. Ninety have a0=0, and 42 of those have a1>0 at xi=0. These are coefficient patterns, not 42 successful theories. Vacuum cancellation at all proper times additionally requires appropriate mass matching, not just a0=0.

**ASSUMED** Every nonempty tuple has aC>0. Thus these ordinary free fields cannot cancel the Weyl-squared coefficient through boson/fermion counting alone. Having a nonzero higher-curvature term is not itself failure: it can be suppressed in a low-energy effective expansion. It is not a fundamental high-energy completion or proof of stability of an untruncated theory.

## Explicit matched spectrum

**ASSUMED** Choose Ns=2, Nf=1, Nv=0, xi=0, and equal masses m=0.1. With

**ASSUMED** $I0(m)=integral_1^infinity ds exp(-m^2 s)/s^3,$

**ASSUMED** $I1(m)=integral_1^infinity ds exp(-m^2 s)/s^2,$

**ASSUMED** $I2(m)=integral_1^infinity ds exp(-m^2 s)/s,$

**ASSUMED** the flat determinant contribution is proportional to 2I0-2I0=0. In contrast a1=1/2, aC=1/24, aR=1/36, so

**ASSUMED** B=a1 I1/(32 pi^2)=0.001503464734,

**ASSUMED** 1/G=16 pi B=0.075572380214.

**ASSUMED** These numbers are in arbitrary cutoff units. Choosing the cutoff to match measured G would be calibration. The common mass was imposed; the field content resembles a chiral multiplet, but no interacting locally supersymmetric gravity theory was built.

**ASSUMED** Using |C|=aC I2/(32 pi^2) and |D|=aR I2/(32 pi^2), a simple derivative power-counting estimate is (|C|+|D|)q^2/B. At q/m=0.001,0.01,0.1 the ratios are 5.91e-9,5.91e-7,5.91e-5. They show that these local coefficients can be small in their intended low-momentum range. The estimate is not a complete tensor dispersion relation or a quantum-gravitational validity bound.

## Mass splitting exposes the cancellation requirement

**ASSUMED** Keep the two scalar masses fixed and increase the Majorana mass to m(1+delta), keeping the cutoff and subtraction prescription fixed. Direct proper-time integration of the difference, evaluated stably with expm1, agrees with the independent difference of the individual integrals.

| Fractional fermion mass split | Magnitude of regulated vacuum contribution / Mbar^4 |
|---:|---:|
| **DERIVED** 0 | 0 |
| **DERIVED** 1e-8 | 1.3303e-7 |
| **DERIVED** 1e-6 | 1.3303e-5 |
| **DERIVED** 1e-4 | 0.0013303 |
| **DERIVED** 0.01 | 0.133710 |

**ASSUMED** This is an illustrative sensitivity calculation, not a bound from observations. Separate readjustment of the vacuum subtraction at every split would hide the issue. Exact global supersymmetry can protect a flat-space vacuum cancellation, while broken supersymmetry and coupling to gravity require further analysis [3,4].

**ASSUMED** An additional nondegenerate example uses squared masses (0.01,0.02,0.03,0.04) and signed mode counts (-2,6,-6,2). It obeys sum w=sum w m^2=sum w m^4=0, removing the flat vacuum quartic, quadratic and logarithmic divergences. Nevertheless its finite term

**ASSUMED** rho_finite=sum w m^4 log(m^2)/(64 pi^2)=2.6468390874e-7

**ASSUMED** does not vanish. A convergent proper-time integral independently agrees. Canceling divergent moments is not a solution of the small-vacuum-energy problem. No curvature-divergence cancellation or protecting symmetry was supplied for this engineered spectrum.

## A specified symmetry-breaking model

**ASSUMED** To avoid arbitrary mass changes, use the known rigid-SUSY matter superpotential

**ASSUMED** $W=fX+(\frac{lambda}{2})X phi1^2+m phi1 phi2$

**ASSUMED** with canonical field kinetic terms. This is an O'Raifeartaigh model of the type discussed in [3,4]. On its X=phi1=phi2=0 branch, define y=lambda f/m^2<1. The heavy real scalar squared masses are (m^2+lambda f,m^2-lambda f,m^2,m^2); the two heavy Majorana fermions have mass m. The tree vacuum energy is f^2. These relations are obtained from the potential, and the script independently verifies the real-scalar Hessian and the vanishing supertrace of squared masses.

**ASSUMED** The light X degrees of freedom and goldstino remain in the low-energy theory; they are not integrated out. Their possible contributions and interactions must not be mistaken for absent physics. The calculation fixes the tree reference X=0, not a fully solved curved-space or loop-stabilized configuration. We vary a parameter f; we do not simulate a dynamical ordering transition at fixed Hamiltonian. SUSY breaking is not a quantitative measure of nuclear ordering.

**ASSUMED** Set m=0.1, lambda=0.5 and minimal scalar curvature coupling. Integrating out only the heavy sector gives a coefficient B(y) that increases slightly with y, hence its isolated induced G decreases. Analytically, I1 as a function of mass squared has positive second derivative; therefore splitting a scalar pair symmetrically increases I1(m^2+delta)+I1(m^2-delta). This explains the direction of the effect without a fitted G(order) formula.

**ASSUMED** The same heavy determinants generate a regulated vacuum shift. Adding it to the REQUIRED tree energy f^2 gives the following diagnostic:

| y | Heavy-sector G(y)/G(0) | Curvature magnitude / lightest heavy mass squared |
|---:|---:|---:|
| **DERIVED** 0.001 | 0.9999999983 | 0.00002655 |
| **DERIVED** 0.01 | 0.9999998262 | 0.002679 |
| **DERIVED** 0.03 | 0.9999984360 | 0.02461 |
| **DERIVED** 0.1 | 0.9999825956 | 0.2947 |
| **DERIVED** 0.3 | 0.9998411884 | 3.409 |
| **DERIVED** 0.9 | 0.9982809521 | 214.43 |

**ASSUMED** The curvature diagnostic uses R=4rho/Mbar^2 as an Einstein-level estimate. It is a consistency gate, not a self-consistent solution of the full higher-curvature equations. Values near or above unity invalidate the assumed local low-curvature heavy-field expansion. A criterion such as |R|/m_min^2<0.01 is an illustrative conservative gate, not an experimental limit or a universal upper bound across all parameters. The tiny changes at small y can pass that gate; larger changes at these chosen parameters cannot be trusted as flat-background gravitational predictions.

**ASSUMED** This is conditional evidence that a concrete change in a quantum spectrum can shift an Einstein-like coefficient. It does not establish the total measured G, ordinary matter screening, a realistic cosmic vacuum, or the emergence of spacetime. Full local supersymmetry introduces additional gravitational-sector structure and is not represented by merely adding this rigid model to an induced Einstein term.

## Decision

**ASSUMED** Retain the regulated induced-gravity calculation as a consistency benchmark. Reject field-count cancellation alone as a complete explanation. Retain vacuum energy, light degrees of freedom and higher-curvature terms when testing any proposed coupling change. A successful extension needs protected matching conditions and a covariant treatment of the resulting background; larger parameter scans without those ingredients would not resolve the obstacle.

**ASSUMED** All 45 numerical checks pass. The checks validate algebra, quadrature and known-model consistency; they do not validate a new theory of our universe. The complete source and raw results are supplied.

## Sources

**ASSUMED** [1] Solodukhin, Newton constant, contact terms and entropy (2015), https://arxiv.org/abs/1502.03758

**ASSUMED** [2] Visser, Sakharov's induced gravity: a modern perspective (2002), https://arxiv.org/abs/gr-qc/0204062

**ASSUMED** [3] Intriligator and Seiberg, Lectures on Supersymmetry Breaking (2007), https://arxiv.org/abs/hep-ph/0702069

**ASSUMED** [4] Goodsell, Lectures on Supersymmetry, Part 3 (2019-2020), https://www.lpthe.jussieu.fr/~goodsell/SUSY/part3.pdf

**ASSUMED** Run induced_gravity_spectrum_audit.py with NumPy, SciPy and Matplotlib. All numerical scales are illustrative and dimensionless relative to the cutoff.
