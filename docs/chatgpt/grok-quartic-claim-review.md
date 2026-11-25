# Independent review of Grok quartic claims C1-C8

**Conversation date range** $2026-10$ (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records independent review of Grok quartic claims $C1-C8$. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Scope: algebra and ordinary scalar controls were independently checked. The complete interacting constrained composite theory, its regulator extension and its renormalized response were not provided or computed. Agreement on algebra does not establish a physical prediction or originality. Source conventions are mostly-minus Minkowski unless explicitly Wick rotated.

## C1: AGREE

**DERIVED** Classification: reproduced classical algebra, conditional on the tree-background parametrization.

**DERIVED** Let $p=D/2-1, A=\eta _{inverse} K$ and $K_{\mu\nu} =\sum \partial _\mu \phi \partial _\nu \phi$. Before counterterms, the standard general-D background gives

**DERIVED** $L(t)=(V0/p) \sqrt{\det [I+p t^2 A/V0]} [1+t^2 U2/V0+t^4 U4/V0]^(-p).$

**DERIVED** The cross term at $t^{6}$ is $[(p+1)U2-(p/2)tr(A)]U4/V0$. At $D=4$ it is $[2U2-tr(A)/2]U4/V0$. The normalization of A and the background is essential; the $D=4$ expression cannot be continued by changing only the potential exponent while ignoring the p in the determinant/background factor. This does not supply counterterm-dressed quantum vertices.

## C2: AGREE

**DERIVED** Classification: reproduced order$-1/V0$ identity, not a complete action or loop proof.

**DERIVED** With $U=U2+U4, T_{\mu\nu} =K_{\mu\nu} -\eta _{\mu\nu} [tr(A)/2-U]$, and

**DERIVED** $Pi_abcd=(\frac{p}{2})(\eta _ac \eta _bd+\eta _ad \eta _bc)-(\frac{1}{2})\eta _ab \eta _cd,$

**DERIVED** the order$-1/V0$ term is $-[p T_{\mu\nu} T^{\mu\nu} -(tr T)^{2}/2]/(4V0)$. Its expanded form is

**DERIVED** $-p[\operatorname{tr}(A^2)-(\operatorname{tr}A)^2/2]/(4V_0)-pU\operatorname{tr}(A)/(2V_0)+DU^2/(4V_0)$.

**CITED** This reproduces the C1 mixed terms. U4 enters this T as +eta U4, but the leading matter term -U4 remains separately in the action. A trace-only piece inside T Pi T does not eliminate cross-vertices or second metric variations. The base conventions are from Carone, Claringbold and Vaman, arXiv:1710.09367, Eqs. 2.8 and 3.4-3.5; the extension here is an algebraic substitution, not their published interacting-potential result.

## C3: DISAGREE

**DERIVED** Classification: dimensional identity reproduced; hard-cutoff formula corrected; PV statement conditional on source prescription.

**DERIVED** Define Minkowski $I_1=int d^D l/(l^{2}-m^{2}+i0)$ and $I_2=int d^D l/(l^{2}-m^{2}+i0)^{2}$, with a common normalization and regulator. Then $m^{2} I_2=(D/2-1)I_1$ in dimensional regularization. Substitution into the stated bracket makes it zero. Signs differ if I1 is defined Euclidean; the Wick rotation must be performed consistently.

**DERIVED** For a sharp Euclidean radial cutoff, after stripping angular/2pi factors, the corresponding integral is

**DERIVED** $\int_0^\Lambda dk k^{3} (k^{2}/2+m^{2})/(k^{2}+m^{2})^{2}$
**DERIVED** $= \Lambda ^{4}/[4(\Lambda ^{2}+m^{2})]$
**DERIVED** $= \Lambda ^{2}/4-m^{2}/4 + m^{4}/[4(\Lambda ^{2}+m^{2})]$.

**DERIVED** The claimed expression is only the large-cutoff expansion through the constant term. With the standard $d^{4}k/(2\pi )^{4}$ measure, multiply the radial answer by $1/(8\pi ^{2})$, plus any separately specified source/Wick factors. At $\Lambda =10,m=2$ the exact radial result is 24.0384615385, not 24.

**DERIVED** PV cancellation at $q=0$ holds for a specified sum using each species' own mass in its stress vertex and compatible source insertions, with cancellation of the relevant moments. We checked one such implementation. A physical-only $\phi ^{2}$ insertion is not automatically the same regulated observable. The interacting regulator action and composite operator must be defined. A nonzero hard-cutoff bubble alone is not an observable mass; covariant source/contact completion matters.

## C4: UNVERIFIED

**UNVERIFIED** Classification: conditional low-momentum estimate; central interacting-graviton conclusion unverified.

**DERIVED** $B(0)=0$ implies $B=O(q^{2})$ only if the massive vacuum response is analytic near $q=0$. B chi B is then $O(q^{4})$ only if chi is regular there. A massless/critical auxiliary response can change this power counting.

**UNVERIFIED** More directly, in an isotropic vacuum $B_{\mu\nu} =a \eta _{\mu\nu} +b q_\mu q_\nu$ has identically zero TT projection for nonzero non-null momentum, regardless of whether a(0) vanishes. That eliminates the factorized scalar bridge from the projected channel. It does not compute self-energy insertions, interaction-dependent determinant vertices, metric contacts, source normalization or counterterm shifts. Therefore the full replacement $m \to m_{eff}$ in the pole and residue is not established.

## C5: AGREE

**DERIVED** Classification: reproduced, scheme-dependent leading-N Hartree statement; not an interacting-CEV theorem.

**DERIVED** For Gaussian contractions $<\phi _a \phi _b\ge \delta _{ab} G$ at the chosen normal-ordering reference,

**DERIVED** $:(\phi^2)^2:=(\phi^2)^2-2(N+2)G\phi^2+N(N+2)G^2$.

**DERIVED** The ordinary one-loop tadpole is $\delta m^{2}=\lambda (N+2)G/N$, reducing to lambda G at leading N. Normal-ordering the quartic adds the opposite mass counterterm and cancels that tadpole at this order, at the selected reference. Normal-ordering only the quadratic potential subtracts a vacuum constant and does not cancel this quartic tadpole. This is not zero physical mass renormalization at every order or for every changed state.

**DERIVED** If the same dimensionally regulated G and the free CEV coefficient are used, $M_P^{D-2}=N G/6$, so $\lambda G=6\lambda M_P^{D-2}/N$. Writing $M_P^{2}$ requires the four-dimensional convention/limit with its regulator scale handled consistently. This algebraic relation cannot independently prove that the interacting gravitational residue retains the free CEV formula; doing so would assume C4. In an interacting PV theory G also requires its specified regulator coupling prescription.

## C6: UNVERIFIED

**ASSUMED** Classification: proposed interacting tuning; partial vacuum combinatorics reproduced.

**ASSUMED** The displayed expression is the known free-kernel tuning evaluated at a proposed effective mass. Applying it to the full interacting model requires the calculation missing in C4. At Gaussian order $<(\phi ^{2})^{2}\ge N(N+2)G^{2}$, whereas $<\phi ^{2}>^{2}=N^{2}G^{2}$. The quartic potential expectation equals $\lambda (N+2)G^{2}/4$, approaching $\lambda N G^{2}/4$ at leading N. Thus the proposed term has the correct leading-N size, but not the finite-N expression.

**DERIVED** It is not by itself the complete c2 condition. Mass counterterm insertions, vacuum diagrams, derivative-sector/c1 effects, local metric variations and the background condition must be included consistently. Subtracting the quartic expectation without specifying the mass/normal-ordering scheme can double count or miss contributions. The graviton tadpole is a metric variation, not merely a bare potential expectation.

## C7: DISAGREE

**DERIVED** Classification: resummation technique established; necessity and spin-2 consistency asserted too strongly.

**ASSUMED** A reduced auxiliary response may be written lambda/(1+lambda Pi_s) with suitable definitions, signs and normalizations. At fixed first order in lambda and away from poles it equals $\lambda +O(\lambda ^{2})$; full resummation is not automatically required just because the channel is spin zero. It is appropriate at leading large N nonperturbatively in lambda, when lambda Pi_s is not small, or near a collective pole. The proposed Pi_s normalization and analytic continuation must be supplied. First-order consistency of spin two still requires the missing complete kernel and Ward checks.

## C8: DISAGREE

**UNVERIFIED** Classification: unresolved simulation prerequisites established for the present proposal; claimed automatic sign problem unproven.

**CITED** No complete lattice action, measure, determinant branch, gauge treatment and regulator prescription has been specified, so there is no well-defined simulation ready to run. That is a missing construction, not a proof that none exists. Negative PV coefficients in loop subtractions do not alone prove a fluctuating-sign/complex Euclidean sampling weight. A product of positive determinants raised to real signed powers remains positive. This does not prove locality, reflection positivity or a convergent auxiliary-boson representation. Those issues differ from a sign problem and require an explicit construction. Li and Vaman, arXiv:2512.06132, formulate PV contributions through weighted determinants/loop prescriptions; that is not automatically a lattice measure. — eq./page UNVERIFIED.

## Missing contributions and checks

- **DERIVED** Full O(lambda), leading-N contractions of the independent potential and mixed derivative/potential six-field vertices, with a completeness/power-counting argument.
- **DERIVED** Constituent self-energy plus matching mass and other required counterterm insertions; source/stress-vertex corrections.
- **DERIVED** Local second metric variations and contact terms, not just separated stress correlators.
- **DERIVED** Vacuum double-bubble contributions, $c1/c2$ shifts and the selected stationary-background condition.
- **DERIVED** A specified interacting regulator and regulated composite-source definition; ordinary matter controls are insufficient.
- **DERIVED** Interacting Ward identities for the complete sum, spin-2 residue sign and the separate scalar spectrum.

## Single decisive next calculation

**DERIVED** Compute the complete regulated O(lambda), leading-N TT inverse response through $O(q^{2})$, including contact and counterterm terms, on a stationary background in one explicit prescription. If a flat reference is used, enforce its interacting tadpole condition. Extract the constant and $q^{2}$ coefficients and compare the latter with the mass-substitution prediction using the same source normalization and parameter definitions. This directly tests C4 and C6. A vanishing projected bridge alone does not.

## Verification and sources

**CITED** The accompanying script passes 55 algebra, quadrature and counterexample checks. They cover general $D$ C1/C2 algebra, the hard-cutoff correction, a specified PV cancellation, normal-ordering combinatorics, the conditional dimensional Planck relation and a positive determinant-ratio counterexample. They do not verify C4 or C6. — eq./page UNVERIFIED.

**CITED** Precisely named primary papers consulted: — eq./page UNVERIFIED.

- **CITED** Christopher D. Carone, Tangereen V. B. Claringbold and Diana Vaman, Composite graviton self-interactions in a model of emergent gravity, Phys. Rev. D 98, 024041 (2018), https://arxiv.org/abs/1710.09367. — eq./page UNVERIFIED.
- **CITED** Christopher D. Carone, Composite Gravitons from Metric-Independent Quantum Field Theories, https://arxiv.org/abs/1911.13166. — eq./page UNVERIFIED.
- **CITED** Chris Li and Diana Vaman, A model of composite gravity with Pauli-Villars regulators, https://arxiv.org/abs/2512.06132. — eq./page UNVERIFIED.

**CITED** No novelty claim. No message was sent to Grok; this is a human-mediated review of the supplied claims. — eq./page UNVERIFIED.
