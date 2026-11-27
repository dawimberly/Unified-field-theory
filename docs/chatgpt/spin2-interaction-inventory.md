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

## Follow-up: first-order contact cancellation in the external-metric control

**ASSUMED** Use the ordinary Euclidean scalar control on a prescribed metric, with potential $\lambda(\phi_a\phi_a)^2/(4N)$ and a common covariant regulator for the coincident propagator, mass counterterm, and vacuum term. The mass is fixed by subtracting the flat momentum-independent tadpole, and the flat vacuum is made stationary. This calculation does not integrate over the metric or identify an external source with the composite metric.

**DERIVED** Let $G_g(x,x)$ be the regulated single-species free propagator and $G_0$ its flat value. Wick contraction gives $\langle(\phi_a\phi_a)^2\rangle=N(N+2)G_g^2$. The order-$\lambda$ vacuum functional before counterterms therefore has integrand $C G_g^2$, with $C=\lambda(N+2)/4$. The tadpole is $\lambda(1+2/N)G_0$. The fixed-mass counterterm contributes $-2C G_0G_g$ to the integrand, and the flat vacuum subtraction contributes $+CG_0^2$. The complete correction in this control is
$$
\Delta\Gamma[g]=C\int d^Dx\sqrt g\,[G_g(x,x)-G_0]^2+O(\lambda^2).
$$

**DERIVED** At the flat background its second variation is $2C\int d^Dx\,\delta G_g[h_1]\delta G_g[h_2]$. Terms with $\delta^2G_g$ or variations of $\sqrt g$ have zero coefficient there. Thus the mass-insertion and local contact terms are included in the cancellation. At nonzero momentum, rotational invariance permits only $\delta G_g[h]=[a(q^2)\delta_{\mu\nu}+b(q^2)q_\mu q_\nu]h_{\mu\nu}(q)$, which vanishes for transverse-traceless sources. The first-order TT correction to the prescribed-metric control is zero at all nonzero momenta with a symmetry-preserving regulator and the stated matching. The scalar Hessian is not eliminated by this argument.

**CITED** Carone, Erlich and Vaman report the large-$N$ composite spin-2 pole in their baseline scattering model. [Emergent Gravity from Vanishing Energy-Momentum Tensor, abstract](https://arxiv.org/abs/1610.08521) **UNVERIFIED** No interacting composite pole or pair of Newton residues is extracted here. The prescribed-metric cancellation cannot settle the determinant-action six-field contractions or the reduced-measure and source matching in [the constraint identity](matched-constraint-identity.md). The next calculation is to test whether that remaining reduced-theory sector adds a term to the TT kernel at this order, then evaluate source-saturated residues.

## Follow-up: reduced six-field vertex at one insertion

**ASSUMED** Use the flat, translationally and reflection invariant free reference state and a common regulator. Let $S=\sum_a\phi_a^2$, $D=\sum_a\partial_\mu\phi_a\partial^\mu\phi_a$, and $I=\langle\phi_a^2\rangle_0$ for a fixed species. Colons below order each block separately with respect to that state. This is an explicit local contraction prescription for a control calculation. Its equivalence to the full renormalized reduced measure is **UNVERIFIED**.

**DERIVED** The previously expanded deformation is $\Delta L_6=\lambda S^2(m^2S-D/2)/(4NV_0)$ at tree level. For its potential block the Wick product identity gives
$$
:S^2:\,:S:
=
:S^3:+8I:S^2:+4(N+2)I^2:S:.
$$
One cross-contraction gives $I\sum_a\partial_{\phi_a}(S^2)\partial_{\phi_a}S=8IS^2$. Two give $I^2\sum_a\partial_{\phi_a}^2(S^2)=4(N+2)I^2S$. More crossings are impossible because the smaller block has two fields. In the derivative block every cross-contraction is proportional to $\langle\phi_a\partial_\mu\phi_b\rangle_0=0$ by reflection symmetry. Hence $:S^2:\,:D:=:S^2D:$ for this prescription. A single local insertion generates no lower-field derivative term.

**DERIVED** The induced bilinear is proportional to $S$, so its contribution at this insertion can be absorbed by a matched mass counterterm. The induced quartic is a local scalar potential. If used in the factorized scalar bridge it has zero TT projection by the selection rule above. This leaves the Wick-ordered six-field vertices themselves, contributions from the consistent counterterm expansion of the full reduced action, and contractions joining distinct vertices. The one-insertion result does not remove any of those classes.

**DERIVED** Independent slot enumeration for $N=1,2,3,8$ gave double-cross coefficients $12,16,20,40$, respectively. A reproduction uses four ordered slots $(i,i,j,j)$, two slots $(k,k)$, and loops over $i,j,k$ and all ordered distinct pairs of the first four slots. Keep a pair if both selected labels equal $k$ and both remaining labels equal the same fixed external species. The count is the coefficient of that species' square in the identity.

**CITED** Carone, Claringbold and Vaman describe their baseline tadpole subtraction and give the source-saturated amplitude in Sec. III, Eqs. (3.1)–(3.3), arXiv:1710.09367v2. Their discussion after Eq. (3.10) also identifies a specific higher-vertex graph removed by the baseline tadpole prescription. [Primary paper](https://arxiv.org/pdf/1710.09367) **UNVERIFIED** That baseline subtraction does not establish the contraction prescription above for the independent quartic deformation. The first distinct-vertex TT graphs and the source and counterterm matching still have to be evaluated before extracting a pole shift.

## Follow-up: first two-vertex transverse-traceless contraction census

**ASSUMED** Continue the preceding flat Gaussian normal-ordering control. Count the connected two-source term with bilinear stress insertions $A,B$, the fully normal-ordered independent six-field insertion $X$, and one fully normal-ordered baseline quartic insertion $Y$. Exclude same-vertex contractions. The off-shell TT calculation below uses a rotationally invariant Euclidean regulator and nonzero momentum $q$. **CITED** The baseline stress-source Feynman rule is given by Carone, Claringbold and Vaman, arXiv:1710.09367v2, Eqs. (3.1)–(3.3); see also the higher-vertex discussion after Eq. (3.10). [Primary paper](https://arxiv.org/pdf/1710.09367)

**DERIVED** The single $X$ insertion cannot connect to two bilinear sources without internal contractions: it has six fields and the sources together supply four. Adding $Y$ admits precisely three connection patterns. Set $a=n_{AB}$, $b=n_{AX}$ and $c=n_{BX}$. Degree conservation gives $b+c=3-a$, $n_{XY}=3+a$, $n_{AY}=2-a-b$, and $n_{BY}=2-a-c$. Nonnegative integer multiplicities give $(a,b,c)=(0,2,1),(0,1,2),(1,1,1)$.

**DERIVED** In each $a=0$ pattern a source sends both lines to $X$. Its loop integrates separately at fixed $q$ and communicates with the rest of the graph only through its free derivative indices at $X$. Zero transmitted derivative indices give a rank-two tensor made from $\delta_{\mu\nu}$ and $q_\mu q_\nu$, killed by TT. One transmitted index gives a rank-three tensor made from $\delta$ and $q$; each term has a trace or at least one $q$ on the source indices and is also killed. Two transmitted indices survive the selection rule. After the loop shift, a rotational Gaussian with normalized variance $1/(2s)$ has $\langle (h_{\mu\nu}k_\mu k_\nu)k_\rho k_\sigma\rangle=h_{\rho\sigma}/(2s^2)$ for traceless $h$. This is the numerator tensor before the vertex signs, propagator integrals, and downstream contraction. It shows why the derivative six-field term cannot be discarded using the scalar-bridge argument.

**DERIVED** In the $a=1$ pattern all four slots of $Y$ go to $X$. Their subgraph receives zero net external momentum. As a two-leg insertion at $X$ it is local, giving a mass term or, when both uncontracted $X$ legs carry derivatives, a kinetic term. A single uncontracted derivative has an odd loop vector and vanishes. Holding the pole mass fixed does not by itself fix field normalization or the matching of its source. No net kinetic coefficient is inferred from this classification.

**DERIVED** Independent Gaussian slot census for an individual pairwise flavor contraction at $Y$: label the equal-species pairs of $A,B,X,Y$ as $(1),(2),(3,4,5),(6,7)$, with pair $5$ carrying the derivative at $X$. Enumerate perfect matchings of the fourteen slots, rejecting any pair within the same vertex. Require equal species across each line, and weight by $N$ for each independent equivalence class among the seven pair labels. For either $a=0$ pattern with both derivative legs at the doubled source, the resulting polynomial is $128N^3+256N^2$. For the $a=1$ pattern, summing over placements of those derivative legs before any momentum projection gives $192N^3+1152N^2+1536N$. The first polynomial is nonzero before the momentum and vertex sum. It is a Wick count, not a pole residue.

**DERIVED** With $V_0=O(N)$ and quartic couplings held fixed in the repo normalization, $X$ carries $\lambda/(NV_0)$ and the baseline $Y$ carries $1/V_0$. The counted block therefore scales at most as $O(\lambda)$, one power of $N$ below the free $O(N)$ two-source loop. Changes to the tuned pole at that order are not excluded by this counting. **UNVERIFIED** The calculation omits additional vertices, reducible dressing, and the matched counterterm and source sectors. The next explicit test is to contract the surviving two-index tensor through the actual baseline quartic tensor vertex and compute the low-$q^2$ coefficient with the same field and source normalization.

## Follow-up: finish the TT projection of the first mixed block

**CITED** The baseline stress coupling and the higher-vertex expansion conventions are given in Carone, Claringbold and Vaman, arXiv:1710.09367v2, Eqs. (3.3), (3.6), and (3.10). [Primary paper](https://arxiv.org/pdf/1710.09367) **ASSUMED** Continue the preceding normal-ordered flat Gaussian control with a rotationally invariant Euclidean regulator. This is a calculation of the specified two-vertex block, with no same-vertex contractions.

1. **DERIVED** The $\phi^6$ potential at $X$ sends no derivative indices through the doubled attachment in either $a=0$ topology. Its isolated source loop is annihilated by the TT projection.
2. **DERIVED** In the derivative $\phi^6$ term at $X$, both derivative legs have their indices contracted with each other in $S^2D$. When both go to the doubled source, the rank-four loop computed in the previous checkpoint is traced at $X$. Its would-be surviving tensor is proportional to $h_{\rho\sigma}$, so its contraction with $\delta^{\rho\sigma}$ is zero.
3. **DERIVED** The zero- and one-derivative attachments were already killed by the source projection. The two-derivative trace now completes the argument for both $a=0$ topologies, before the momentum integral at $Y$. Equivalently the contracted momentum factor is $k\cdot(q-k)$, a rotational scalar made from $k,q$. The integrated tensor on the doubled source then only has a trace and a longitudinal term. The preceding positive Wick count remains a count of the open-index terms and gives no TT contribution after closing the indices.
4. **DERIVED** In the $a=1$ graph both sources attach to each other and to $X$, while $Y$ attaches only to $X$. Momentum conservation leaves the four-line $X$–$Y$ loop independent of the external momentum. The $X$ insertion into the remaining loop is therefore a local scalar self-energy, with a mass term, a kinetic term, or a vanishing odd-vector term according to the placement of the derivatives at $X$.
5. **DERIVED** Isolate the kinetic case by leaving both derivatives of $D_X$ uncontracted. The four slots of $S_X^2$ must contract across to the four slots of $Y$. For a single-species radial propagator $G(z)$ and $u_\mu=\partial_\mu G(z)$, all differentiated legs at $Y$ see the same separation. Each such leg gives $u_\mu$, and each undifferentiated leg gives $G$. This factorization is in position space and does not assume the individual momentum-space momenta are collinear.
6. **DERIVED** Parametrize the formal Euclidean baseline polynomial, with $1/V_0$ stripped, as $b_1 D^2+b_2 K_{\mu\nu}K_{\mu\nu}+b_3DS+b_4 S^2$. The kinetic contraction is proportional to
   $$
   I_Y=\int d^Dz\,[ (b_1+b_2)|u|^4+b_3G^2|u|^2+b_4G^4].
   $$
   The common flavor factor is $8N(N+2)$, from the all-cross contraction of the two $S^2$ blocks. The expression strips the prefactor of $X$, the $1/V_0$ at $Y$, and the perturbation-expansion sign. The tensor contractions at $Y$ coincide in the all-cross term because each differentiated line is parallel in position space.
7. **DERIVED** For the formal polynomial of the tree expansion with Euclidean contractions, $b_1=1/8$, $b_2=-1/4$, $b_3=-m^2/4$, and $b_4=m^4/4$. These follow from $a_2-a_1U_2+U_2^2$ in [the tree-level expansion](quartic-vertex-derivation.md). They yield $-|u|^4/8-m^2G^2|u|^2/4+m^4G^4/4$. This is a reduction of that algebraic polynomial, conditional on using those coefficients in the Euclidean control. No Lorentzian normalization or renormalized value has been assigned.
8. **DERIVED** Independent enumeration of the all-cross pairings at $N=5$ gives $280$, equal to $8N(N+2)$, in either ordering of the target slots. For a smooth test function $G=e^{-r^2}$ in four dimensionless Euclidean dimensions, the three radial integrals of $G^4$, $G^2|u|^2$, and $|u|^4$ in units of $\pi^2/16$ are $1,2,6$. Midpoint integration agrees within $10^{-12}$. For the formal polynomial at $m=1$ the same units give $-1$. The Gaussian is a check on the contraction and integration, not the physical scalar propagator.
9. **DERIVED** An isolated local kinetic term can be tracked by field normalization in a matched external-metric stress response. If the renormalized kinetic coefficient is $Z$ at fixed pole mass and the source is the metric derivative of that same kinetic action, two propagators contribute $Z^{-2}$ and the two corresponding kinetic source vertices contribute $Z^2$. Their product is unity. This identity shows exactly what must be matched, including the contact sector if a full metric Hessian is evaluated. It does not extend automatically to an independently normalized reduced composite source.
10. **DERIVED** The first mixed block is reduced to local self-energy and its matching. The candidate two-index nonlocal graph from the previous checkpoint vanishes. **UNVERIFIED** The finite composite kernel at this order, extra determinant vertices, the reduced measure, and matching of the nonlinear composite source still remain. The immediate comparison is to derive the counterterm and matched source variations from the same reduced action and test the cancellation of the remaining local insertion there. No interacting pole shift or universal Newton residue is inferred.

## Follow-up: matched normalization of the reduced action and composite source

**CITED** The composite metric expansion is given by Carone, Claringbold and Vaman, arXiv:1710.09367v2, Eqs. (3.7)–(3.9). The need to transform the source together with the action under a field redefinition is discussed by Criado and Pérez-Victoria, JHEP 03 (2019) 038, §2. [Composite metric](https://arxiv.org/pdf/1710.09367) · [Field redefinitions](https://arxiv.org/pdf/1811.09413)

1. **ASSUMED** Use the nonsingular four-dimensional branch and a finite flat field-coordinate regulator for this change-of-variables check. Keep the clock and vacuum normalization $V_0$, canonical mass $m$, and canonical quartic coupling $\lambda$ fixed. Introduce a positive, constant normalization $Z$ of the matter fields.
2. **DERIVED** With $S=\sum_a\phi_a^2$ and $K_{\mu\nu}=\sum_a\partial_\mu\phi_a\partial_\nu\phi_a$, set $V_Z=V_0+Zm^2 S/2+\lambda Z^2S^2/(4N)$. The constrained metric and reduced density in this direction are
   $$
   g_Z[\phi]=\frac{V_0\eta+ZK}{V_Z},\qquad
   L_Z[\phi]=\frac{V_0^2}{V_Z}\sqrt{\det(I+Z\eta^{-1}K/V_0)}.
   $$
   These formulas follow by using the kinetic normalization $Z$ before eliminating the metric in the same tree-level determinant prescription.
3. **DERIVED** The change $\psi_a=\sqrt Z\,\phi_a$ gives $g_Z[\phi]=g_1[\psi]$ and $L_Z[\phi]=L_1[\psi]$ pointwise. This is an exact redundant direction in the prescribed regulated action.
4. **DERIVED** With $n_{\rm coord}$ regulated scalar integration coordinates, the Jacobian is $Z^{-n_{\rm coord}/2}$ and is independent of the composite-metric source. It cancels from the normalized generator with source $J^{\mu\nu}(g_{Z,\mu\nu}-\eta_{\mu\nu})$. All connected derivatives with respect to that source are invariant along this direction.
5. **DERIVED** Expand the density and source in homogeneous field degree $2j$. A normalization change $Z=1+z$ inserts $jzL_{2j}$ in the action and $jz(g-\eta)_{2j}$ in the source. The kinetic and mass terms are accompanied by the four- and six-field descendants $2zL_4$ and $3zL_6$, plus their matched source changes and the analogous higher terms.
6. **DERIVED** Writing $U_2=m^2S/2$ and $U_4=\lambda S^2/(4N)$, the bilinear part of $g-\eta$ is $(K-\eta U_2)/V_0$. Its quartic part is $[-KU_2+\eta(U_2^2-V_0U_4)]/V_0^2$. TT removes the terms proportional to $\eta$. In particular the nonlinear $KU_2$ source term is rescaled by $Z^2$, while the bilinear $K$ is rescaled by $Z$.
7. **DERIVED** In the free bilinear TT bubble the two source factors of $Z$ cancel the two propagator factors of $Z^{-1}$. Keeping the old source with the rescaled quadratic action would instead give a spurious relative factor $1-2z+O(z^2)$. The exact normalized-generator identity is the corresponding check with the descendants and nonlinear source restored.
8. **DERIVED** A fixed-mass local self-energy $\delta\Sigma(p^2)=z(p^2+m^2)$ can be removed from the quadratic block by a counterterm in this direction. The matching subtraction also gives the opposite four-field, six-field and source terms listed above. Discarding those terms would be a different matching prescription.
9. **DERIVED** An arithmetic check at twenty-four positive normalizations and regular diagonal metrics gave equality between the two evaluations of the density and all diagonal composite-metric components to the working precision. The check repeats by evaluating the displayed formulas first at $(Z,S,K)$ and then at $(1,ZS,ZK)$. The symbolic change of variables establishes the result for generic allowed matrices.
10. **DERIVED** The contribution along the complete matched normalization direction cannot shift a composite-metric pole or its matched residue. **UNVERIFIED** The radiatively generated reduced action and operator corrections have not been shown to lie entirely along that direction. Subtract it from both and calculate the residual four-field action and nonlinear source coefficients in the chosen regulated measure. Only then can a non-redundant shift of the pole or differences between distinct physical probes be tested. The local self-energy alone supplies neither result.
