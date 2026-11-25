# First reduction of the independent-quartic interaction calculation

**Conversation date range** $2026-10$ (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records first reduction of the independent-quartic interaction calculation. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

## What is established

**DERIVED** A factorized scalar-channel bridge between two stress insertions has zero transverse-traceless (TT) projection in a flat, isotropic, zero-temperature vacuum. This removes that particular contribution to the spin-2 response. It does not remove self-energy insertions or local metric contact terms, and it is not a proof that the complete composite-graviton response is unchanged.

**DERIVED** The diagram inventory is an organized list of required contribution classes. It is NOT a proven-complete list of fully contracted composite-theory Feynman diagrams. The independent quartic potential, full determinant-action vertices, counterterm prescription, measure and selected background must all be treated at consistent orders before making that claim.

## Tensor calculation

**DERIVED** For nonzero Euclidean four-momentum $q, \theta =I-qq/q^{2}$. The spin-zero projector is $P_{0,abcd}=\theta_{ab}\theta_{cd}/3$. The TT projector is $P_{2,abcd}=(\theta_{ac}\theta_{bd}+\theta_{ad}\theta_{bc})/2-P_{0,abcd}$.

**DERIVED** P2 annihilates both I and qq. It is idempotent and orthogonal to P0. Its rank is five in four-dimensional Euclidean off-shell tensor space; this is NOT a statement that a physical massless graviton has five polarizations. On-shell Lorentzian constraints and gauge equivalence are separate.

**DERIVED** The mixed scalar/stress loop in a rotationally invariant vacuum has the tensor form $S_{ab}=A(q^{2})\delta _{ab}+B(q^{2})q_a q_b$. Thus $P2 S=0$. A factorized bridge S_ab S_cd also vanishes when projected on both ends, irrespective of its scalar prefactor.

**DERIVED** The code checks this with independent Gaussian-Hermite quadrature of the shifted tensor numerator, not just by constructing the assumed $A I+B qq$ form. Twelve momentum orientations were used. The largest scalar-bubble TT residual is $4.52\times10^{-16}$, and the largest doubly projected bridge residual is $1.07\times10^{-16}$. These are numerical zeros at floating-point precision. Isotropic integration and a compatible regulator are assumptions. The argument does not justify an inconsistent loop-momentum shift with a noninvariant cutoff.

**DERIVED** A medium introduces an additional vector u, and P2(uu) need not vanish. The code verifies that counterexample. Finite temperature, finite density, an anisotropic condensate or a time-dependent state requires additional tensor structures. The vacuum cancellation must not be transferred to such states without a new calculation.

## Contributions that survive or remain unresolved

1. **DERIVED** Ordinary quartic scalar tadpoles shift the propagator mass. At leading large $N, \delta \chi =\lambda I_1(\chi )/(16\pi ^{2})$ in the ordinary scalar control with proper-time cutoff one. Such a shift changes the TT loop. At finite N the ordinary one-loop combinatorics include $a (1+2/N)$ factor; this run uses the leading-N limit only. Whether this mass shift is retained or subtracted depends on the composite theory's renormalization conditions; it is not a scheme-independent new prediction.
2. **DERIVED** A single metric variation of a local potential is trace-only and drops under TT projection. Its SECOND variation does not. For a normalized TT shear h, the second variation of $\sqrt{\det (I+t h)}$ at zero is $-1/2$. The numerical calculation gives approximately -0.50000006. It is a local contact term, not independently an observable graviton mass.
3. **DERIVED** Vacuum tadpoles, counterterms and background stationarity must be enforced together. An apparent quadratic term around a background that is not a solution cannot simply be called a physical mass.
4. **DERIVED** The determinant action's additional six-field terms from the preceding quartic expansion still require contractions and counterterms. Their derivative pieces cannot be discarded using the trace-only argument for a local potential.
5. **DERIVED** Additional composite vertices, regulator-induced local terms and the scalar response remain to be calculated. The separate spin-zero calculation is necessary even if a spin-2 contribution vanishes.

**DERIVED** The order labels in diagram_inventory.csv refer to the leading ordinary-matter response when appropriate, with external source normalization stripped. They do not establish a complete power counting for every composite scattering graph. In the composite large-N scaling V0 is proportional to N, and loop contractions can compensate explicit $1/N$ factors.

## Explicit surviving-loop control

**DERIVED** To check the self-energy insertion method, define a reduced off-diagonal TT scalar master integral

**DERIVED** $F(Q,\chi )=\int_0^1 dx \int_1^\infty ds \exp [-s(\chi +x(1-x)Q^2)]/s^3.$

**DERIVED** Then $dF/d\chi=-\int _0^{1} dx I_1[\chi +x(1-x)Q^{2}]$. This proper-time master integral omits source normalization, species multiplicity and local metric contact terms. Its cutoff prescription defines a loop block; it is not by itself the second metric variation of the complete regulated quantum action.

**DERIVED** The analytic derivative agrees with independent finite differences at nine mass/momentum combinations. Solving an ordinary leading-N scalar Hartree gap and expanding F to first order in lambda leaves a remainder proportional to $\lambda ^{2}$, with measured powers 1.99934 through 2.00000. That validates the mass-insertion control. It does not calculate the full composite pole, enforce its background condition or show that all other interaction corrections are absent.

## What the next calculation actually requires

**UNVERIFIED** The immediate unresolved task is the renormalized local/contact sector and the contractions of the added determinant-action vertices, at the same lambda and large-N orders as the surviving loop. Then the sum must satisfy the appropriate interacting Ward identities on the chosen stationary background. Only after that sum can a resummed spin-2 denominator, its residue and the scalar spectrum be interpreted physically.

**DERIVED** The published base-model simplifications use specific counterterms to remove tadpole attachments. They cannot be assumed unchanged after adding an independent quartic constituent potential. If one adopts a condition fixing the physical constituent mass, its tadpole subtraction must be carried consistently into all source/contact terms. Otherwise an apparent gravitational correction can be a mismatch of schemes.

## Relation to a future simulation

**DERIVED** See SIMULATION_SPEC.md. A numerical ordinary scalar theory would be a useful control for this selection rule but would not simulate the constrained composite-gravity theory. A TT component of a matter stress correlator is not automatically a gravitational propagator. The simulation requires a definition of the full action, measure, regulator and gauge treatment before it can be launched as a test of this candidate model.

## Reproducibility and novelty status

**ASSUMED** Run `python spin2_interaction_inventory.py`. Outputs are results.json and diagram_inventory.csv. All 115 implementation and identity checks pass. They are not experimental tests. This is ordinary tensor/perturbation groundwork; no originality, massless interacting graviton, physical change in G, or explanation of nuclear binding is claimed.

**DERIVED** Primary background sources:

- **CITED** Composite scalar action and corrected base-model context: https://arxiv.org/abs/1710.09367 — eq./page UNVERIFIED.
- **CITED** Regulator and counterterm dependence of the base model: https://arxiv.org/abs/2512.06132 — eq./page UNVERIFIED.
- **CITED** Background selection and contact variations: https://arxiv.org/abs/2011.06541 — eq./page UNVERIFIED.
- **CITED** Ordinary large-N scalar methods: https://arxiv.org/abs/hep-th/0306133 — eq./page UNVERIFIED.
