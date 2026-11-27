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

## Follow-up: pair the local kinetic graph with its nonlinear source graph

**CITED** The baseline metric expansion is given by Carone, Claringbold and Vaman, arXiv:1710.09367v2, Eqs. (3.7)–(3.9). [Primary paper](https://arxiv.org/pdf/1710.09367) **ASSUMED** Continue the preceding four-dimensional Euclidean polynomial and flat Wick-ordering control, at fixed $V_0$. Use the same regulator and perturbation sign for matched insertions. The following comparison is for specific contractions and does not identify the full reduced measure with an auxiliary-metric measure.

1. **DERIVED** The TT-projected source is the projection of $K/(V_0+U_2+U_4)$. Expanding in field degree gives $K/V_0$, $-KU_2/V_0^2$, and $K(U_2^2/V_0^3-U_4/V_0^2)$ at degrees two, four, and six, respectively, before projection.
2. **DERIVED** The explicit degree-four source term proportional to $\lambda$ comes from the metric numerator's $\eta_{\mu\nu}$ piece and is trace-only. The first direct deformation of the TT source is the degree-six term $-KU_4/V_0^2$. Here $U_4=\lambda S^2/(4N)$.
3. **DERIVED** A single action insertion $X_6$ can join a quartic source and a bilinear source without self-contractions only if all six source lines go to $X_6$. The bilinear-source loop has just $q$ as an external vector once the derivative indices of $D_X$ are closed. Its TT contraction vanishes, for the same trace and longitudinal reasons as in the preceding mixed-block calculation.
4. **DERIVED** Conversely, the direct degree-six TT source $-KU_4/V_0^2$ can join a bilinear source and a single baseline quartic action vertex $Y$. With no same-vertex contractions, all four legs of $Y$ and both legs of the bilinear source must end at the degree-six source.
5. **DERIVED** Select the placement in which the two derivative legs of $K$ reach the bilinear source. The four legs of $S^2$ reach $Y$ and their bundle has zero net external momentum. The resulting insertion is a local correction proportional to the original tensor $K$. Placements with neither or only one leg of $K$ reaching the bilinear source vanish: its separate loop can only supply the previously excluded rank-two or rank-three structure, respectively.
6. **DERIVED** Define $J_Y$ as the integrated all-cross contraction of $S^2$ with $Y$, including the Wick count, the $Y$ coefficient, and its action-insertion sign. It is the very same local subgraph for the selected kinetic graph and source graph. In the notation of the preceding checkpoint it is $8N(N+2)I_Y/V_0$, with the common insertion sign attached.
7. **DERIVED** The selected kinetic graph gives $\delta L_2=-\lambda J_Y D/(8NV_0)$, since the six-field action contains $-DU_4/(2V_0)$. Relative to the canonical coefficient $D/2$, this is $z_Y=-\lambda J_Y/(4NV_0)$. The selected source graph gives $\delta g_{2,\mathrm{TT}}=-\lambda J_Y K_{\mathrm{TT}}/(4NV_0^2)$, a relative source coefficient $s_Y=-\lambda J_Y/(4NV_0)$.
8. **DERIVED** Therefore $s_Y-z_Y=0$ for this paired subgraph before integrating the common scalar loop. When the selected kinetic correction is used at fixed pole mass, its two-propagator normalization change cancels the two matched source corrections in the bilinear TT bubble. This check isolates the normalization term. A mass counterterm and the other source and action terms still need their own matching.
9. **DERIVED** Independent all-cross slot enumeration for $N=1,2,3,5,8$ produced $24,64,120,280,640$, respectively, agreeing with $8N(N+2)$. Arithmetic substitution in the two independently factored coefficient formulas also reproduced $s_Y=z_Y$ to working precision. These checks require no fitted $J_Y$.
10. **DERIVED** This removes the selected $X_6$–baseline-$Y$ normalization pair from the bilinear response and the single-$X_6$ quartic-source term from its TT channel. **UNVERIFIED** The complementary placement with the independent quartic in the action and an undeformed higher vertex or source, the remaining counterterms and source contacts, and additional distinct-vertex graphs have not been summed. The next calculation is the complementary placement and its subtraction at fixed canonical parameters. No physical pole shift or universal Newton residue is inferred from the paired cancellation.

## Follow-up: complementary quartic placement with a retained derivative pair

**CITED** The baseline composite source expansion is in Carone, Claringbold and Vaman, arXiv:1710.09367v2, Eqs. (3.7)–(3.9). [Primary paper](https://arxiv.org/pdf/1710.09367) **ASSUMED** Continue the four-dimensional, rotationally invariant Euclidean polynomial control. Put the independent scalar quartic in the action at $Y$, and retain an entire derivative bilinear at the undeformed six-field vertex or source at $X$. This selects a subset of contractions.

**DERIVED** In the all-cross bundle the four lines to $Y$ have zero total momentum. At separation $z$ an internal undifferentiated pair contributes $G^2$, and a derivative pair contributes $B_{\mu\nu}=\partial_\mu G\,\partial_\nu G$. A retained external derivative pair contributes $C_{\mu\nu}$. Denote $b=\operatorname{tr}B$ and $c=\operatorname{tr}C$. For this subset the coefficient linear in $C$ can be read by differentiating the action polynomial at the rank-one matrix $B$. The common all-cross flavor factor is $8N(N+2)$ as above.

**DERIVED** The determinant polynomials obey
$$
\delta a_1=c/2,\quad
\delta a_2=bc/4-\operatorname{tr}(BC)/2,\quad
\delta a_3=-b^2c/16+b\,\operatorname{tr}(BC)/4.
$$
The last identity uses $B^2=bB$. Angular averaging gives $\langle\operatorname{tr}(BC)\rangle=bc/4$ in this control. Thus $\langle\delta a_2\rangle=bc/8$ and $\langle\delta a_3\rangle=0$. An independent check averages four equally weighted rank-one matrices along the coordinate axes, which gives the same results for a general $C$ because the expressions are linear in its contraction with $B$ at fixed $b$.

**DERIVED** The selected $a_1U_2^2$ term contributes a relative kinetic correction equal to the relative source correction from $KU_2^2$, after stripping their shared four-line integral and the independent quartic insertion. The selected $-a_2U_2$ term also contributes an action term proportional to $-m^2c\,G^2|\nabla G|^2/(16V_0^2)$ with the shared flavor factor and insertion sign stripped. The direct baseline TT source has no term with this internal derivative structure. The selected $a_3$ term vanishes in exactly four Euclidean dimensions. The angular factor before setting $D=4$ is $-1/16+1/(4D)$, so an evanescent term times a regulator pole needs a separate check in dimensional regularization.

**UNVERIFIED** This does not establish a physical remainder. Contractions with retained derivatives taken from different bilinears, regulated contact terms, and matched source and action counterterms can contribute at the same order. The previous equality of one selected source and kinetic pair cannot be applied to this whole placement. The comparison now has an explicit term that must be accounted for in the complete sum.

## Follow-up: sum the retained and split derivative contractions

**CITED** Carone, Claringbold and Vaman, arXiv:1710.09367v2, Eqs. (3.7)–(3.10), specify the baseline source and higher vertices. [Paper](https://arxiv.org/pdf/1710.09367) **ASSUMED** Retain the Gaussian, all-cross Euclidean polynomial control above, with a scalar quartic at $Y$. Compute the local kinetic coefficient at $X$ with two external differentiated legs. This is a bare block and its direct source counterpart.

**DERIVED** Introduce a flavor vector $v$ for internal legs and external derivatives $w_{a\mu}$. At fixed separation put $u=\nabla G$, $A_\mu=\sum_a v_a w_{a\mu}$, $B=v^2uu^\top$, $C_{\mu\nu}=\sum_aw_{a\mu}w_{a\nu}$, and $M=uA^\top+Au^\top$. Then the derivative bilinear expands as $K=B+tM+t^2C$. The $t^2C$ term reproduces the retained pair. The quadratic $tM$ terms contain both split placements. Undifferentiated external legs cannot contribute to this kinetic coefficient. A single external derivative integrates to zero by reflection symmetry.

**DERIVED** Direct substitution into the trace polynomials gives the split parts
$$
[a_2(B+tM)]_{t^2}=-|u|^2|A|^2/2,\qquad
[a_3(B+tM)]_{t^2}
=\operatorname{tr}B\,[|u|^2|A|^2/4+(u\cdot A)^2/2].
$$
After the all-cross contraction, the rotationally invariant flavor weight replaces $v^2 v_a v_b$ by $(v^2)^2\delta_{ab}/N$ inside the average. Spatial rotation replaces $(u\cdot A)^2$ by $|u|^2|A|^2/D$. The flavor normalization for $(v^2)^2$ is the same $8N(N+2)$ as in the retained-pair calculation. For a check, symmetric finite differences of the two trace polynomials for ten test matrices gave the stated quadratic coefficients within $8\times10^{-12}$ absolute error.

**DERIVED** Combining both placements, the coefficient of $\operatorname{tr}C/V_0^2$, stripping the common flavor factor and independent quartic insertion at $Y$, has integrand
$$
\frac{m^4}{8}G^4
-\frac{m^2}{2}\left(\frac14-\frac{1}{2D}-\frac{1}{2N}\right)G^2|u|^2
+\left[-\frac1{16}+\frac1{4D}
+\frac1N\left(\frac14+\frac1{2D}\right)\right]|u|^4 .
$$
The explicit $D$ dependence here varies only the angular average of the stated algebraic polynomial. It is not the continued $D$-dimensional determinant action. In exactly four dimensions the mixed integrand is $-m^2(1/16-1/(4N))G^2|u|^2$, and the last coefficient is $3/(8N)$.

**DERIVED** The direct undeformed six-field TT source with this scalar action insertion still receives only the retained derivative pair. With at most one derivative endpoint on the doubled source, its loop is killed by TT as before. The $G^4$ term matches after normalization. At leading flavor order the split terms carry an extra $1/N$ and therefore do not cancel the retained $-m^2G^2|u|^2/16$ mismatch within this block. No cross-species Newton residue follows. **UNVERIFIED** The other action and operator contractions, contact subtractions, and finite counterterm matching are required for a statement about a physical pole. In particular, setting $D=4$ before regulating the physical integrals cannot settle the evanescent contributions.

## Follow-up: dimensional continuation and contact reduction

**CITED** The determinant action and static gauge convention are given in Carone, Claringbold and Vaman, arXiv:1710.09367v2, Eqs. (2.6)–(2.7). [Paper](https://arxiv.org/pdf/1710.09367) **ASSUMED** Use that tree-level dimension dependence in the formal Euclidean polynomial control, with the clock counterterm omitted at this stage. Match the scalar mass and project on the same TT source as in the preceding block.

1. **DERIVED** Set $p=D/2-1$, $U=U_2+U_4$. The continued density is $(V_0/p)\sqrt{\det(I+pK/V_0)}(1+U/V_0)^{-p}$, with Euclidean indices. The metric is $(\eta+pK/V_0)/(1+U/V_0)$ before Euclidean projection.
2. **DERIVED** The undeformed six-field action, with $V_0^{-2}$ stripped, is $p^2a_3-p^2a_2U_2+[p(p+1)/2]a_1U_2^2-[(p+1)(p+2)/6]U_2^3$. The last term has no external derivative pair.
3. **DERIVED** The relative TT source correction from $KU_2^2$ has coefficient unity after dividing by the baseline bilinear coefficient $p/V_0$ and stripping the common $V_0^{-2}$. The relative kinetic correction from $a_1U_2^2$ has coefficient $p(p+1)/2$. Their apparent equality in four dimensions leaves an evanescent difference in the continued block.
4. **DERIVED** Let $I_{40}=\int G^4$, $I_{22}=\int G^2|\nabla G|^2$, and $I_{04}=\int|\nabla G|^4$, all with the same regulator. Using the split and retained angular factors above, the kinetic-minus-source difference, divided by twice the common graph factor (which includes $V_0^{-2}$), is
   $$
   \Delta=\frac{m^4}{8}\left[\frac{p(p+1)}2-1\right]I_{40}
   -\frac{p^2m^2}{2}\left(\frac14-\frac{1}{2D}-\frac{1}{2N}\right) I_{22}
   +p^2\left[\frac{4-D}{16D}+\frac1N\left(\frac14+\frac1{2D}\right)\right]I_{04}.
   $$
5. **DERIVED** For $D=4-\varepsilon$, the leading-flavor $I_{40}$ prefactor bracket is $-3\varepsilon/4+\varepsilon^2/8$, and the leading-flavor angular factor of $I_{04}$ is $\varepsilon/(16D)$. Both can multiply UV poles. Taylor and direct angular-factor evaluations agreed to floating-point precision in ten trial dimensions. This does not evaluate either pole.
6. **ASSUMED** For the integration-by-parts check let a smooth massive regulated propagator obey $(-\nabla^2+m^2)G=\rho_\Lambda$ with vanishing boundary flux. Its smeared source and any regulator counterterms must use the same prescription.
7. **DERIVED** Integrating $\nabla\cdot(G^3\nabla G)$ gives $3I_{22}=\int G^3\rho_\Lambda-m^2I_{40}$. Replacing the regulated contact by a point delta multiplied by the unregulated singular propagator would not define the required integral.
8. **DERIVED** Consequently even at exactly four dimensions and leading flavor, the leftover $-m^2I_{22}/16$ contains both a scalar integral $m^4I_{40}/48$ and a contact term $-m^2\int G^3\rho_\Lambda/48$. A calculation that drops the contact would alter the kinetic-minus-source coefficient in this block.
9. **DERIVED** The isolated block cannot give a regulator-independent finite mismatch by setting $D=4$ before integration. There are two distinct evanescent prefactors and an explicit contact. These must be combined with subtractions of the same regulator and the other graphs before comparing the pole.
10. **UNVERIFIED** The Laurent coefficients after subdivergence subtraction and the full matched composite operator contact terms are still open. No physical pole shift or cross-species Newton residue is extracted.

## Follow-up: dimensional contact identity and three-loop masters

**CITED** Groote, Körner and Pivovarov, Phys. Rev. D 60, 061701 (1999), use configuration space for water-melon vacuum topologies. [Paper](https://doi.org/10.1103/PhysRevD.60.061701) **ASSUMED** Use the massive continuum Euclidean propagator $G_D(r)$ with Fourier measure $d^Dk/(2\pi)^D$. Define bare dimensional integrals by analytic continuation from $1<D<4/3$, where the following integrals converge. Renormalization scale factors are suppressed uniformly.

1. **DERIVED** At the origin $G_D(0)=(4\pi)^{-D/2}(m^2)^{D/2-1}\Gamma(1-D/2)$ in that strip. This follows from its Schwinger representation at zero separation.
2. **DERIVED** The radial boundary flux is $-\Omega_{D-1}r^{D-1}G_D'(r)\to1$, since $(-\nabla^2+m^2)G_D=\delta^{(D)}$ and the mass integral over the ball vanishes as its radius shrinks.
3. **DERIVED** Applying integration by parts outside the ball and taking its radius to zero therefore gives $3I_{22}=G_D(0)^3-m^2I_{40}$. By meromorphic continuation the same identity holds for the bare dimensional masters. It does not identify an independently renormalized contact in an arbitrary cutoff or source scheme.
4. **DERIVED** Near $D=4-\varepsilon$, the product $G_D(0)^3$ has leading Laurent term $-8m^6/[(4\pi)^6\varepsilon^3]$ in these conventions. Thus the contact appearing in the preceding block has explicit pole content. Keeping only a possible simple pole is insufficient at the bare level.
5. **DERIVED** Write $G_D(r)=\int_0^\infty dt\,(4\pi t)^{-D/2}e^{-m^2t-r^2/(4t)}$. For four copies set $t_i=sx_i$, $\sum x_i=1$, $X=\prod x_i$, $H=\sum x_i^{-1}$, and $U=XH$. Gaussian integration gives
   $$
   I_{40}=(4\pi)^{-3D/2}(m^2)^{3D/2-4}
     \Gamma(4-3D/2)\int_{\sum x=1} U^{-D/2}\,d^3x.
   $$
6. **DERIVED** Each gradient contributes $-r_\mu/(2t_i)$ before contractions. The normalized fourth radial moment under the combined Gaussian is $D(D+2)/(\sum t_i^{-1})^2$. Combining these factors gives
   $$
   I_{04}=\frac{D(D+2)}{16}(4\pi)^{-3D/2}(m^2)^{3D/2-2}
     \Gamma(2-3D/2)\int_{\sum x=1} X\,U^{-D/2-2}\,d^3x.
   $$
7. **DERIVED** These masters both have explicit scale-integration gamma poles near four dimensions. Their simplex boundaries encode additional UV subgraphs. Coefficients obtained by expanding only the gamma factors would omit those boundaries.
8. **DERIVED** At the convergent endpoint $D=1$, $G=e^{-m|r|}/(2m)$ gives $(I_{40},I_{22},I_{04})=(1,m^2,m^4)/(32m^5)$. The contact identity becomes $3m^2/(32m^5)=1/(8m^3)-1/(32m^3)$. Midpoint integration at $m=1.7$ reproduced the three exact values within $2\times10^{-6}$ relative error.
9. **DERIVED** Substituting $I_{22}$ into the preceding kinetic-minus-source expression reduces the isolated bare block to two simplex masters and one explicit tadpole cube. Its exact four-dimensional formal coefficients at leading flavor are $m^4/48$ on $I_{40}$ and $-m^2/48$ on the tadpole cube, with the evanescent contributions retained in the continued calculation.
10. **UNVERIFIED** The subtracted simplex coefficients and the corresponding composite-operator renormalization remain to be matched. The exact bare identity reduces the calculation and fixes the contact sign, but provides no renormalized graviton pole or residue by itself.

## Follow-up: published basketball poles and the evanescent subtraction

**CITED** Martin and Robertson, *Evaluation of the general 3-loop vacuum Feynman integral*, arXiv:1610.07720, Secs. II–III, define the scalar basketball master $E$ and give its bare Laurent coefficients in Eqs. (3.1), (3.5)–(3.8). [Paper](https://arxiv.org/pdf/1610.07720) **ASSUMED** Use their loop normalization throughout this comparison, with $D=4-2e$, $x=m^2$, and $L=\ln(x/Q^2)$. This $e$ is half the previous checkpoint's $\varepsilon$.

1. **DERIVED** Fourier transformation of $I_{40}$ gives four equal-mass propagators and one momentum-conserving delta function. This is the bare basketball $E$ after applying the three loop normalization factors.
2. **DERIVED** With $A_{\rm bare}/x=-1/e+(L-1)+O(e)$, the normalized contact is $A_{\rm bare}^3$, in the same three-loop convention.
3. **DERIVED** Substituting equal masses in the published basketball coefficients gives $E/x^2=2/e^3+(23/3-6L)/e^2+O(e^{-1})$. The six pair products contribute $2$ to the cubic coefficient. The double coefficient is $-6(L-1)+2-1/3$.
4. **DERIVED** In the leading-flavor block write $\Delta=x^2 c_E E-xc_T A_{\rm bare}^3+c_F I_{04}$, suppressing graph factors and using the same normalization on $I_{04}$. The coefficients are $c_T=p^2(1/4-1/(2D))/6$, $c_E=c_T+[p(p+1)/2-1]/8$, and $c_F=p^2(4-D)/(16D)$.
5. **DERIVED** Their expansions are $c_T=1/48-5e/96+O(e^2)$, $c_E=1/48-23e/96+O(e^2)$, and $c_F=e/32+O(e^2)$. In particular using the previous $\varepsilon$ in the published pole series without converting it would give the wrong evanescent contribution.
6. **DERIVED** The two scalar terms together give
   $$
   \frac{x^2c_E E-xc_T A_{\rm bare}^3}{x^4}
   =\frac{1}{16e^3}+
   \frac{-89/288-3L/16}{e^2}+O(e^{-1}).
   $$
   The cubic coefficient follows independently by combining the bare $E$ coefficient $2$ and the tadpole cube coefficient $-1$ with the opposite signs in $\Delta$. The displayed double coefficient excludes $I_{04}$.
7. **DERIVED** The exact Schwinger expression for the gradient master gives an alternative reduction target: in the unscaled Fourier convention it equals $D(D+2)(4\pi)^6/16$ times the equal-mass basketball in $D+4$ with all four propagators squared. The squared powers supply one $t_i$ each, and the dimension shift adds two powers of the inverse graph determinant. Their ratio gives the gradient factor $X/U^2$ in projective parameters. The mass dimension of the shifted integral is $3D-4$, as required.
8. **DERIVED** That doubled basketball is also the mixed derivative with respect to the four independent squared masses of the unsquared basketball in $D+4$. Four minus signs cancel. Differentiating the equal-mass specialization four times would instead include derivatives that strike the same propagator more than once.
9. **DERIVED** Subtracting the scalar master first and multiplying by its dimension-dependent coefficient is different from multiplying its bare Laurent series first and then subtracting poles. Terms of order $e$ can multiply subgraph poles, including the explicit tadpole cube. The bare scalar coefficients therefore provide a check for the subtraction, not a finite kinetic-minus-source coefficient.
10. **UNVERIFIED** The remaining input for this block is the Laurent series of the gradient master or its dimension recurrence, followed by the same subgraph and composite-source subtraction. Published scalar basketball subtractions do not by themselves fix counterterms for the full reduced operator or establish a graviton pole shift.

## Follow-up: reduce the gradient master in the original dimension

**CITED** Martin and Robertson, arXiv:1610.07720, Eqs. (3.5)–(3.6), give the scalar basketball poles used below. [Paper](https://arxiv.org/pdf/1610.07720) **ASSUMED** Use massive Euclidean dimensional regularization with a common factor per loop and scaleless polynomial integrals set to zero. Let $x=m^2$, $s=P^2$, $A$ be the bare tadpole, and $B(s)$ the equal-mass bubble.

1. **DERIVED** Group the four lines into two pairs with total momenta $P,-P$. Within each pair the gradient numerator is $(s-k^2-(P-k)^2)/2$, up to a Fourier sign common to the pair. The signs cancel between the pairs.
2. **DERIVED** Cancel the inverse propagators in this numerator to obtain the pair integral $(s+2x)B(s)/2-A$. Hence the gradient master is the momentum integral of its square.
3. **DERIVED** Define $M_j=\int_P s^j B(s)^2$, with $M_0=E$. Translation and rotational symmetry give $\int_P B=A^2$ and $\int_P s B=-2x A^2$. In the latter identity the mixed loop moment vanishes and $\int_k k^2/(k^2+x)=-x A$ in this prescription.
4. **DERIVED** The cross term in the square therefore cancels, giving $I_{04}=M_2/4+xM_1+x^2E$. This equality uses dimensionally regulated polynomial integrals. It is not a pointwise equality of the momentum integrands.
5. **DERIVED** Feynman parameterization writes $B$ as a constant times $\int_0^1 [x+s t(1-t)]^{D/2-2}dt$. Integrating the derivative of $(t-1/2)[x+s t(1-t)]^{D/2-2}$ and using the tadpole gamma recurrence gives
   $$
   s(s+4x)B'(s)=[(D-4)s/2-2x]B(s)-(D-2)A.
   $$
6. **DERIVED** Radial integration by parts with $s^{D/2+j}(s+4x)B^2$ yields
   $$
   (3D/2+j-3) M_{j+1}+x(2D+4j-4)M_j
   =2(D-2)A\int_P s^j B.
   $$
   These are dimensional identities. Endpoint polynomial terms in individual summands must be treated in the common dimensional prescription rather than assigned ordinary convergent integrals at positive integer dimension.
7. **DERIVED** The cases $j=0,1$ reduce $I_{04}$ to $E$ and $A^3$ alone. In particular, with $a=3D/2-3$,
   $$
   M_1=[2(D-2)A^3-(2D-4)x E]/a,\quad
   M_2=[-4(D-2)xA^3-2Dx M_1]/(a+1).
   $$
   This gives the original-dimension reduction without a recurrence from $D+4$.
8. **DERIVED** In $D=1$ at $m=1$, the formula gives $I_{04}=1/32$ using $E=1/32$ and $A=1/2$, in agreement with the direct position-space result. At $D=4$ the rational reduction coefficients multiplying $x^2 E$ and $x A^3$ approach $1/3$ and $1/6$. These are limits of the coefficients, since the masters themselves have poles.
9. **DERIVED** In $D=4-2e$, the gradient master has $I_{04}/x^4=1/(2e^3)+O(e^{-2})$ in the previous loop normalization. Multiplication by the evanescent $c_F=e/32+O(e^2)$ adds $1/(64 e^2)$ to the bare leading-flavor block. Including the scalar terms, the checked bare result is
   $$
   \Delta/x^4=1/(16e^3)+(-169/576-3L/16)/e^2+O(e^{-1}).
   $$
10. **DERIVED** The isolated gradient-master reduction closes at two published scalar inputs. **UNVERIFIED** A bare pole is not a renormalized composite-source correction. The local and subgraph counterterms still need to be combined with the remaining action and source graphs before assigning any shift to the physical graviton pole.

## Follow-up: close the isolated bare Laurent series

**CITED** Martin and Robertson, arXiv:1610.07720, Eqs. (3.5)–(3.8), give the basketball Laurent coefficients in terms of the tadpole and the finite subtracted basketball. [Paper](https://arxiv.org/pdf/1610.07720) **ASSUMED** Use that common loop normalization, $D=4-2e$, $x=m^2$, and $L=\ln(x/Q^2)$, keeping the same leading-flavor graph factors stripped.

1. **DERIVED** The preceding bubble recurrence gives $I_{04}=g x^2 E+h x A^3$, where $a=3D/2-3$, $r=-(2D-4)/a$, $t=2(D-2)/a$, $g=1+r[1-D/(2(a+1))]$, and $h=t[1-D/(2(a+1))]-(D-2)/(a+1)$.
2. **DERIVED** Combining all three terms yields $\Delta/x^4=v(E/x^2)+w(A/x)^3$, with $v=c_E+c_Fg$ and $w=-c_T+c_Fh$, using the coefficients in the earlier basketball checkpoint.
3. **DERIVED** Through cubic order, $v=1/48-11e/48+17e^2/192-e^3/256$ and $w=-1/48+11e/192-35e^2/768+5e^3/1024$, up to $O(e^4)$.
4. **DERIVED** Write $A/x=-1/e+b+eu+e^2t_A+O(e^3)$. Expansion of the tadpole gamma function gives $b=L-1$, $u=-1+L-L^2/2-\zeta(2)/2$, and $t_A=-1+L-L^2/2+L^3/6+(L-1)\zeta(2)/2+\zeta(3)/3$ (Eq. (2.10)).
5. **DERIVED** The tadpole cube has coefficients $(-1,3b,3u-3b^2,b^3-6bu+3t_A)$ in powers $e^{-3},e^{-2},e^{-1},e^0$.
6. **DERIVED** Equal masses in Eqs. (3.5)–(3.8) give the normalized basketball coefficients $(2,23/3-6L,6b^2-6u-5b+1/2,\mathcal E+12bu-5u-6t_A)$ in the same powers, where $\mathcal E$ is the paper's finite subtracted basketball divided by $x^2$.
7. **DERIVED** Convolution reproduces the previous cubic and double poles and gives the simple pole
   $$
   [\Delta/x^4]_{e^{-1}}
   =\frac{9L^2}{32}+\frac{169L}{192}
    +\frac{3\zeta(2)}{32}-\frac{2803}{2304}.
   $$
8. **DERIVED** The finite bare coefficient is
   $$
   \frac{\mathcal E-b^3+(18b-5)u-9t_A}{48}
   -\frac{11}{48}(6b^2-6u-5b+1/2)
   +\frac{11}{192}(3u-3b^2)
   +\frac{17}{192}(23/3-6L)
   -\frac{105b}{768}-\frac{13}{1024}.
   $$
   Here $t_A$ includes $+\zeta(3)/3$, so the explicit $\zeta(3)$ term in this coefficient is $-\zeta(3)/16$. Compared with the previous sign this shifts the finite bare coefficient by $-\zeta(3)/8$. The pole terms use only $b$ and $u$ and are unchanged.
9. **DERIVED** Numerical polynomial convolution at $L=0,\pm1$, with the zeta terms formally suppressed, gives simple-pole coefficients $-2803/2304,-127/2304,-4183/2304$, respectively in the order $0,+1,-1$. Restoring a unit test value for $\zeta(2)$ changes the coefficient by $3/32$. These are reproducibility checks with artificial inputs.
10. **DERIVED** The published subtractions enter through $\mathcal E$. They do not subtract the entire composite block. Its bare finite part includes evanescent prefactors multiplying lower-order pole coefficients. **UNVERIFIED** The full counterterm insertion and composite-source matching are still required to determine a physical pole or residue.

## Subtraction checks (independently verified)

**CITED** Martin and Robertson (M–R), [arXiv:1610.07720](https://arxiv.org/pdf/1610.07720), Eqs. (2.3), (2.6)–(2.13), (2.35)–(2.38), and (3.5)–(3.8). Use $\epsilon=e$, $x=m^2$, $y=L-1$, $a=A_{\rm bare}/x$, and $A=A_{\rm bare}$. Denote the finite, subtracted equal-mass basketball divided by $x^2$ by $\mathcal E$. Write $R'$ for subtraction of proper subgraphs, $K$ for pole extraction, and $R=(1-K)R'$.

1. **DERIVED** In the M–R loop normalization, the six pair and four triple terms in Eqs. (2.35)–(2.38), plus the overall subtraction, give
   $$
   E_{\rm div}/x^2=6a^2/\epsilon+(6/\epsilon^2-5/\epsilon)a
       +2/\epsilon^3-10/(3\epsilon^2)+1/(2\epsilon).
   $$
   Subtracting $v E_{\rm div}/x^2$ from $\Delta/x^4$ leaves the poles of $wa^3$ because $v(E-E_{\rm div})/x^2$ is finite. Their coefficients, from cubic to simple, are $1/48$, $1/192-L/16$, and $3L^2/32-L/64-1/768+\zeta(2)/32$. These follow by multiplying the tadpole series (Eqs. (2.6)–(2.10)) by the stated $w$ and checking the basketball cancellation against Eqs. (3.5)–(3.8).

2. **DERIVED** Replacing $a^3$ by $(a+1/\epsilon)^3$ and using the subtracted basketball produces the scalar reference $(\mathcal E-y^3)/48$. This is the $D=4$ coefficients times the $\overline{\rm MS}$ masters by construction (Eqs. (2.7)–(2.10), (2.35)–(2.38)). It is a reference for matching the reduced masters.

3. **DERIVED** With the contact product subtracted as three independent tadpoles, keep its fully contracted graph as an overall term. At the common dimension and loop normalization of Eqs. (2.3), (2.6)–(2.13), the resulting identity is
   $$
   R'I_{22}=-\frac{x}{3}R'E+\frac{(A+x/\epsilon)^3}{3}-\frac{x^3}{3\epsilon^3}.
   $$
   This can also be obtained by inserting the two-line and three-line contact corrections $x A^2/\epsilon+x^2 A/\epsilon^2$ into $I_{22}=(A^3-x I_{40})/3$, with the contact terms counted once. The M–R overall basketball polynomial in Eq. (2.38) then gives $K R'I_{22}=x^3[-1/\epsilon^3+10/(9\epsilon^2)-1/(6\epsilon)]$. Its coefficients contain no $L$. The finite part is $x^3(y^3-\mathcal E)/3$, the same combination as the scalar reference after accounting for the prefactor $-1/16$ on $I_{22}/x^3$ at $D=4$.

4. **DERIVED** For $I_{04}$, the two-line forest contribution in units of $x^4$ is $C(D)a^2/\epsilon$, with $C(D)=2/3+10/(3D)$. The two same-product hard pairs each give $1/D$ after the spectator average. The four split hard pairs each give $1/6+1/(3D)$. To obtain the split term use the free-index hard pole $\{\delta_{\mu\nu}(x/2+P^2/12)+P_\mu P_\nu/6\}/\epsilon$ from the Feynman-parameter shift, then contract in $D$ dimensions. For the same-product pairs the traced hard pole is $(2x+P^2/2)/\epsilon$. Its constant part integrates to zero against the remaining $r\cdot s$. The spectator moment follows from $\int k^2/(k^2+x)=-x A$ in the conventions of Eqs. (2.3) and (2.6)–(2.10).

5. **DERIVED (quotient check)** Inserting the proposed local three-line pole
   $$
   K\bar R S(p)=\frac{p^2}{96\epsilon}
       [p^4+24p^2x+70x^2-36x^2/\epsilon]
   $$
   in the spectator tadpole uses $\int (p^2)^n/(p^2+x)=(-x)^n A$. Summing four spectators gives $(\alpha/\epsilon^2+\beta/\epsilon)a$ in units of $x^4$, with $\alpha=3/2$ and $\beta=-47/24$. These coefficients also follow independently by requiring the $L$ and $L^2$ pieces of the quotient poles to cancel given $C(D)$ and the bare reduction. That condition by itself checks only the integrated polynomial at $p^2=-x$. The Fourier and locality calculation below supplies the two independent off-shell coefficients missing from that quotient check.

6. **DERIVED** With that three-line quotient, the subtraction check is the Laurent convolution
   $$
   R'I_{04}/x^4=g E/x^2+h a^3-C(D)a^2/\epsilon
        -(3/(2\epsilon^2)-47/(24\epsilon))a .
   $$
   Expand $g=(1/3,1/6,1/8,3/32)$, $h=(1/6,-1/24,-1/32,-3/128)$, and $C=(3/2,5/12,5/24,5/48)$ through cubic order. Use the tadpole and basketball coefficients from Eqs. (2.6)–(2.13), (2.35)–(2.38), and (3.5)–(3.8). The $L$ and zeta terms cancel from the pole, leaving $1/(2\epsilon^3)-13/(9\epsilon^2)+149/(288\epsilon)$. The finite part is
   $$
   \lim R I_{04}/x^4=\mathcal E/3+y^3/6+(17/24)y^2-(121/96)y+51/128 .
   $$

7. **DERIVED** Subtracting the scalar-master reference $\mathcal E/3+y^3/6$ leaves $(17/24)y^2-(121/96)y+51/128$, with value $303/128$ at $Q^2=x$. Both expressions have leading log $y^3/6$. A finite insertion of $17/36$ times the derivative-pair residue tensor yields $(17/36)(3/2)y^2=(17/24)y^2$ with finite spectator tadpoles. The remaining linear and constant terms have the form of a local three-line insertion into a tadpole and an overall local counterterm. $I_{22}$ needs no corresponding shift. **UNVERIFIED** The value $239/96$ does not follow from retaining $g(\epsilon)$ and $h(\epsilon)$ while fully subtracting the two scalar masters before taking the limit. Their positive powers then multiply finite series and vanish. A prescription for what evanescent products to retain would be needed to check that number. The proposed $\sigma=-749/864$ and $\delta_0=-1037/3456$ also require an explicit definition of their insertion polynomials and finite subgraph prescription. The finite scalar-master prescription assigns different shifts to derivative and non-derivative vertices. Compatibility with shared renormalization conditions is open.

### Open

- Match finite shifts consistently across every graph sharing each vertex.
- Determine where $c_F(\epsilon)$ enters the nested and overall counterterms and retain evanescent terms in the same prescription throughout.
- Renormalize the composite source and its coincident products, and include the remaining graphs.
- Impose the Ward/BRST identities and on-shell conditions before constructing the resummed TT propagator.

No graviton pole or residue follows from these checks.


### Off-shell differentiated sunset from two Fourier residues and locality

**CITED** The configuration-space product representation of the sunset is described by Groote, Körner and Pivovarov, [hep-ph/9805224](https://arxiv.org/pdf/hep-ph/9805224). Use M–R Eq. (2.3) for loop normalization and the two-line pole convention above.

**DERIVED** Set $q=p^2$ and write $G$ for the Euclidean massive propagator in position space. The three hard derivatives and the last, external derivative give
$$
S(p)=-\int d^D r\,e^{ip\cdot r}\,
   \nabla\cdot[(\nabla G)^2\nabla G],
$$
with the two M–R loop factors understood. This follows by representing the spectator derivative by $ip$ and integrating by parts. For $r\ne0$, the equation $(\nabla^2-x)G=0$ and the Coulomb singularity fix the first two terms
$$
\partial_\mu G=-b_D r_\mu r^{-D}[1+a_D x r^2+\cdots],
\qquad b_D=\frac{\Gamma(D/2)}{2\pi^{D/2}},
\quad a_D=\frac{1}{2(2-D)}.
$$
The regular branch starts at relative order $x^{D/2}r^D$ in the gradient and contributes to the $q x^2$ pole sector near four dimensions. Two-line counterterms inside the sunset also contribute only to $qx^2$. In each case their numerator is a polynomial in $p$ and the lone remaining loop momentum $k$, linear in the explicit $p$ from the spectator. Terms with an odd number of $k$ components vanish. The even terms contain two $k$ components, so the surviving scalar is proportional to $q\int k^2/(k^2+x)=-qxA$. This leaves the $q^3$ and $q^2x$ residues directly accessible from the first two short-distance terms.

**DERIVED** Apply the continued Fourier formula
$$
\int d^D r\,e^{ip\cdot r}(r^2)^{-\nu}
 =\pi^{D/2}2^{D-2\nu}\frac{\Gamma(D/2-\nu)}{\Gamma(\nu)}
 (p^2)^{\nu-D/2}.
$$
For the leading term $\nu=3D/2-1$. The coefficient of $q^3/\epsilon$ is
$$
-(16\pi^2)^2\,b_4^3(6)\,\pi^2 2^{-6}
       \frac{-1}{12\,\Gamma(5)}=\frac1{96}.
$$
The second term shifts $\nu$ down by one. In the residue the ratio of Fourier factors is $4(4)(-3)$, the ratio of divergence factors is $4/6$, and the relative gradient-cube coefficient is $3a_4=-3/4$. The product is $24$. Thus the $q^2x/\epsilon$ coefficient is $1/4$ and neither coefficient has a double pole.

**DERIVED** The remaining local polynomial must be proportional to $qx^2$, because rotational invariance and the explicit external derivative supply a factor of $q$. Determine its coefficient from the independently required absence of $L$ in the three-loop vacuum counterterm. With the two-line term $C(D)a^2/\epsilon$ already fixed, subtract $(\alpha/\epsilon^2+\beta/\epsilon)a$ from the reduced vacuum master after its two-line subtraction. The coefficient of $y$ in the double pole is $3/2-\alpha$. At $\alpha=3/2$, the coefficient of $y$ in the simple pole is $-47/24-\beta$ and the $y^2$ coefficient vanishes. This fixes $\alpha=3/2$ and $\beta=-47/24$ without supplying the proposed off-shell polynomial as input. For the four identical spectators the quotient operation is exactly $4A S(q=-x)$ on the counterterm polynomial. If the unknown $qx^2$ coefficients are denoted $B_2/\epsilon^2+B_1/\epsilon$, the resulting normalized quotient gives
$$
-4 B_2=\alpha,\qquad -4(1/96-1/4+B_1)=\beta.
$$
It follows that $B_2=-3/8$, $B_1=35/48$, and hence
$$
K\bar R S(p)=\frac{q}{96\epsilon}
   [q^2+24qx+70x^2-36x^2/\epsilon].
$$
This is an off-shell reconstruction from two directly calculated residues, locality and the vacuum quotient. The $qx^2$ term is determined by matching and has not been separately evaluated by direct off-shell integration. All of this concerns the selected UV subtraction. The finite shared-vertex and composite-source matching conditions remain open.


### Evanescent prefactors in the overall subtraction

**CITED** Use the proper-subgraph versus overall subtraction split of Martin and Robertson, [arXiv:1610.07720](https://arxiv.org/pdf/1610.07720), Eqs. (2.35)–(2.38). The general dependence on the renormalization prescription for structures that vanish at four dimensions is discussed by Dugan and Grinstein, [Phys. Lett. B 256, 239 (1991)](https://doi.org/10.1016/0370-2693(91)90680-O).

**ASSUMED** For this comparison carry each graph's prefactor through its already established proper-subgraph subtraction, without changing the inner tensor counterterms. Compare overall pole subtraction on the prefactored expression with subtraction on each master before multiplying. This fixes the scheme of this isolated block. An endpoint coefficient arising from the final rotational projection is carried outside the inner subgraphs, where the spectator momentum is still unintegrated.

**DERIVED** With $p=D/2-1=1-\epsilon$ and $B=1/4-1/(2D)$, the coefficients of the unreduced, dimensionless masters in $\Delta/x^4$ are
$$
\begin{aligned}
d_{40}&=\tfrac18[p(p+1)/2-1]=-3\epsilon/16+\epsilon^2/16,\\
d_{22}&=-p^2B/2=-1/16+5\epsilon/32-7\epsilon^2/64+\epsilon^3/128+O(\epsilon^4),\\
c_F&=p^2(4-D)/(16D)=\epsilon/32-3\epsilon^2/64+\epsilon^3/128+O(\epsilon^4).
\end{aligned}
$$
Let $P_j$ denote the overall pole of the corresponding normalized $R'I_j$, in units of $x^2,x^3,x^4$ respectively. In decreasing pole order the rows for $P_{40},P_{22},P_{04}$ are $(2,-10/3,1/2)$, $(-1,10/9,-1/6)$, and $(1/2,-13/9,149/288)$. The first row follows from M–R Eq. (2.38), and the last two are the forest results above.

**DERIVED** Only the $P_j$ rows can contribute a finite difference between the two prescriptions. For a coefficient with Taylor entries $(d_0,d_1,d_2,d_3)$ and a pole row $(P_3,P_2,P_1)$, the finite difference is $d_1P_1+d_2P_2+d_3P_3$. Multiplying these rows gives
$$
\begin{array}{c|rrrr}
 &\epsilon^{-3}&\epsilon^{-2}&\epsilon^{-1}&\epsilon^0\\ \hline
 d_{40}P_{40}&0&-3/8&3/4&-29/96\\
 d_{22}P_{22}&1/16&-65/288&169/576&-179/1152\\
 c_FP_{04}&0&1/64&-79/1152&809/9216\\ \hline
 {\rm sum}&1/16&-337/576&1123/1152&-3407/9216
\end{array}
$$
For instance the gradient finite term is $149/(32\cdot288)+39/(64\cdot9)+1/(128\cdot2)$. The gradient coefficient vanishes in four dimensions, but multiplying before overall pole subtraction leaves that finite local term. All three differences are independent of $L$.

**DERIVED** In this comparison the master-by-master result is $(\mathcal E-y^3)/48$. Keeping the full coefficients and subtracting only the poles at the overall step gives $(\mathcal E-y^3)/48-3407/9216$, in the same stripped normalization. To match the latter to the former requires an additive local finite shift $+3407/9216$ in this block. It specifies only the net shift in the selected projected channel. It does not fix the finite operator coefficients before projection, whose values must agree across all occurrences of the vertices.

**OPEN** The next condition is operator-level matching of the action and the composite source, including the Ward/BRST and on-shell conditions. The two overall prescriptions above leave the inner subtractions fixed. A different extension of the evanescent tensors to inner vertices needs its own nested finite insertions. The result supplies no graviton pole or residue.


### Forward stress matching and the independent derivative condition

**CITED** Brown and Collins, [Ann. Phys. 130, 215 (1980)](https://doi.org/10.1016/0003-4916(80)90232-8), construct the external-metric counterterms and distinguish flat-space renormalization from curvature-dependent terms. The following is a conditional external-metric control for the selected block.

**DERIVED** In the stripped normalization above, write $s$ and $z$ for the relative local bilinear source and kinetic corrections, each divided by twice the same common graph factor and by $x^4$. Thus $z-s=\Delta/x^4$ in this paragraph. The direct source has coefficient $1/8$ on $I_{40}/x^2$, independent of $\epsilon$. Its overall subtraction is $P_{40}/8$ and its finite coefficient is $s=\mathcal E/8$ in either overall prescription. The kinetic pole is obtained by adding $P_{40}/8$ to the previous pole row for the difference. Its entries are
$$
(5/16,\ -577/576,\ 1195/1152)
$$
in decreasing pole order. The kinetic finite coefficient with full prefactors and overall pole subtraction is
$$
z=\frac{7\mathcal E}{48}-\frac{y^3}{48}-\frac{3407}{9216}.
$$
In particular the scheme conversion $+3407/9216$ goes entirely to the kinetic side if the scalar source subtraction is held fixed. This is an allocation of the conversion **within this block**, not a condition on the sum of the operator insertions.

**DERIVED** For a source defined by variation of a covariant external metric, let $\Gamma(k^2)$ be the scalar inverse two-point function and $\Gamma_{\mu\nu}^T(k',k)$ the corresponding amputated translation-current vertex. In the Euclidean convention $q=k'-k$ and with tree kinetic vertex $k'_\mu k_\nu+k'_\nu k_\mu-\delta_{\mu\nu}(k'\cdot k+x)$, the flat-space Ward identity is
$$
q^\mu\Gamma^T_{\mu\nu}(k',k)
=k_\nu\Gamma(k'^2)-k'_\nu\Gamma(k^2).
$$
Varying the scalar arguments under a background diffeomorphism gives the two inverse-propagator contact terms. Differentiating at $q=0$ for a regular vertex yields
$$
\Gamma^T_{\mu\nu}(k,k)=2 k_\mu k_\nu\Gamma'(k^2)-\delta_{\mu\nu}\Gamma(k^2).
$$
At a simple scalar pole, $\Gamma(-x_{\rm pole})=0$ and the two external legs multiply the vertex by the inverse of $\Gamma'(-x_{\rm pole})$. The forward matrix element therefore has its tree-level $2k_\mu k_\nu$ normalization. The local source and kinetic renormalizations must agree in this forward component **after summing the graphs and counterterms** if the composite source is identified with this metric variation. For the selected block, the amount to be balanced in the full forward condition is the computed $z-s$. This identity supplies a condition on the sum and cannot be applied to each graph independently. The scalar improvement from $R\phi^2$ is proportional to $q_\mu q_\nu-\delta_{\mu\nu}q^2$ in a single insertion, so it cannot supply a TT forward shift.

**DERIVED** The same Ward identity permits an independent slope in the TT scalar form factor. To exhibit it within a local external-metric derivative expansion, add
$$
\lambda_R\int d^D r\,\sqrt g\,R_{\mu\nu}(g)
  \sum_a\nabla^\mu\phi_a\nabla^\nu\phi_a
$$
with an appropriately dimensionful coefficient. Its flat-space action vanishes, so it leaves the constituent mass and wavefunction conditions unchanged. In the convention $g=\delta+h$, the linearized Ricci tensor for a transverse-traceless $h$ is $-\partial^2 h_{\mu\nu}/2$. Thus the first variation supplies a TT bilinear proportional to $\lambda_R q^2 K_{\mu\nu}$ and vanishes in the forward limit. The background Ward identity is preserved by covariance. This is a counterexample to fixing the first TT momentum-dependent coefficient from the forward identity and scalar on-shell data alone. It applies to the external-source correlator and does not, by itself, imply a change to a dynamical pole.

**OPEN** For the reduced composite source the identification with a translation-normalized metric variation still has to be derived, or its finite forward normalization imposed separately. Even with that identification, a finite-momentum form factor or a microscopic matching calculation is required to fix the derivative coefficient. The remaining graph sum and quadratic source terms are still needed before constructing a resummed TT propagator.


### Leading TT pole matching

**CITED** Carone, Claringbold and Vaman, [arXiv:1710.09367](https://arxiv.org/pdf/1710.09367), extract the composite pole and compare scalar species in the resummed scattering amplitude.

**DERIVED** For an isolated simple TT pole, strip the on-shell stress-tensor kinematics and write $t=q^2$ and $\mathcal M_{ij}=F_i(t)F_j(t)/\mathcal K(t)+H_{ij}(t)$ with $H$ regular. At a simple zero $t_\star$ of the inverse exchange kernel, the residue is $\mathcal R_{ij}=F_i(t_\star)F_j(t_\star)/\mathcal K'(t_\star)$. An analytic invertible exchange-field rescaling sends $F_i$ to $f F_i$ and $\mathcal K$ to $f^2\mathcal K$ and leaves this ratio invariant. The derivative of $f^2$ multiplies the vanishing kernel. The relation $\mathcal R_{12}^2=\mathcal R_{11}\mathcal R_{22}$ tests factorization. Universality additionally requires $\mathcal R_{12}/\mathcal R_{11}=F_2/F_1=1$ for nonzero couplings in common stress-tensor normalization.

**DERIVED** At a massless pole with fixed kernel and forward vertices, the independently changeable source slope $F_i\to F_i+\rho_i t$ changes the regular term by $(\rho_i F_j(0)+\rho_j F_i(0))/\mathcal K'(0)$ and leaves the residue invariant. A curvature operator added to the microscopic action could also affect the kernel. The result here is for the external insertion. A purely external contact term in $H$ also has zero residue.

**DERIVED** Allow the pole to move by expanding $\mathcal K=Z_0t+Y_0t^2/2+\lambda(\kappa_0+\kappa_1t)+\cdots$ and $F_i=f_i+b_it+\lambda a_i+\cdots$. Put $u=-\kappa_0/Z_0$. The root is $t_\star=\lambda u+O(\lambda^2)$ and the coefficient of $\lambda$ in the residue is
$$
\frac{a_i f_j+f_i a_j+u(b_i f_j+f_i b_j)}{Z_0}
-\frac{f_i f_j(\kappa_1+uY_0)}{Z_0^2}.
$$
The baseline slope enters when the pole moves. An order-$\lambda$ source slope first enters at order $\lambda^2$. After massless retuning, only forward corrections $a_i$ and dynamical kernel slope $\kappa_1$ enter the first-order residue. This assumes a simple pole with no coincident threshold.

**CHECK** For $\mathcal K=2t+3t^2$, $F_1=4+5t$ and $F_2=7-2t$, the mixed residue is $14$ and the regular constant is $-15/2$. Replacing the coefficient $5$ by $6$ raises the constant by $7/2$ and retains the residue.

**OPEN** The summed forward matching and the constant and linear terms of the dynamical TT kernel still need calculation. The endpoint block alone does not supply $a_i$, $\kappa_0$ or $\kappa_1$. No graviton pole or residue is established for the deformation.
