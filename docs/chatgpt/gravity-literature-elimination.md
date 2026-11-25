# Gravity research by elimination: published results and the remaining bridge

**Conversation date range** 2026-10-29 to 2026-10-30

**ASSUMED** This checkpoint records Gravity research by elimination: published results and the remaining bridge. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Research checkpoint: October 29, 2026, America/Chicago.

**ASSUMED** Method: prioritize peer-reviewed results, inspect the assumptions, and eliminate unsupported inferences before doing new numerical work. This is a literature-backed research assessment and a proposed application, not a proof of quantum/Einstein unification.

## What changed in the research strategy

**DERIVED** The most useful next step is not another parameter scan. It is to use existing auxiliary-field quantization results to decide exactly which quantum formulations may be compared. The terms "same classical action," "same BRST cohomology," "same scattering amplitudes," and "same quantum measure" must be kept distinct.

## Published results that directly address the current gap

**CITED** **[1] Henneaux (1990), Physics Letters B 238, 299-304.** The publisher abstract describes gauge-fixing choices that keep auxiliary variables algebraic and make the gauge-fixed actions agree whether elimination occurs before or after setting up the path integral. Integration-measure contributions are explicitly part of the problem. This supplies an appropriate framework; it does not identify the measure for our model automatically. Access checked: publisher abstract and bibliographic record; full article was not retrieved. — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1016/0370-2693(90)91739-X — eq./page UNVERIFIED.

**CITED** **[2] Barnich, Brandt, Henneaux (1995), Communications in Mathematical Physics 174.** Section 15, Eq. (15.1), characterizes auxiliary variables by a locally solvable algebraic field equation. Theorems 15.1 and 15.2 establish invariance of local BRST cohomology under auxiliary/generalized-auxiliary elimination. The full author manuscript was inspected at printed pages $40-44$. This is a local formal cohomology result, not a theorem of positive spectral density or equality of arbitrary nonperturbative path-integral measures.

**CITED** https://doi.org/10.1007/BF02099464 — eq./page UNVERIFIED.

**CITED** Full manuscript: https://arxiv.org/pdf/hep-th/9405109 — eq./page UNVERIFIED.

**CITED** **[3] Criado and Perez-Victoria (2019), JHEP 03, 038.** Section 2 separates the Jacobian, transformed sources, and on-shell equivalence. The familiar simplification of local perturbative Jacobian ghosts in dimensional regularization does not license ignoring source transformations in off-shell functions. The paper also warns against resumming perturbative redefinition terms into propagators and interpreting the resulting extra poles as physical. Our exact matter-loop response is not automatically such an artifact; that identification would require an actual matching transformation. Full Section 2 and relevant appendix passages were inspected. — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1007/JHEP03(2019)038 — eq./page UNVERIFIED.

**CITED** Full manuscript: https://arxiv.org/pdf/1811.09413 — eq./page UNVERIFIED.

**CITED** **[4] Castela Simao, Cattaneo, Schiavina (2023), Letters in Mathematical Physics 113, 25.** Their treatment distinguishes classical/BV equivalence from stronger compatibility with boundary structures needed for quantization. Classically related examples need not share the same boundary regularity. This prevents treating a bulk cohomology argument as a completed scattering or boundary-quantization argument. Full open-access introduction and stated results were inspected; its examples are not a completed analysis of CEV. — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1007/s11005-023-01646-2 — eq./page UNVERIFIED.

## Applying the known framework to our model

**ASSUMED** The following is our proposed application of those methods, with its unproved steps stated explicitly.

**ASSUMED** Use the clock-free auxiliary classical action in Lorentzian conventions

$$
\begin{gathered}
S0[g,\phi ] = \int \sqrt{|g|} [g^{\mu \nu } P_{\mu\nu} /2 - V(\phi )], \\
P_{\mu\nu} = \sum _a \partial _\mu \phi ^a \partial _\nu \phi ^a.
\end{gathered}
$$

**DERIVED** On a patch with V !$= 0, D$ !$= 2$, and nondegenerate P, its algebraic metric equation has the solution

$$
G_\mu \nu [\phi ] = (D/2-1) P_{\mu\nu} / V.
$$

**CITED** This provides the candidate local auxiliary-elimination map required by [2]. It is branchwise: it does not select a real Lorentzian contour, include singular configurations, or prove that the quantum vacuum lies on this patch. — eq./page UNVERIFIED.

**DERIVED** At the minimal classical BV level, introduce the diffeomorphism ghost c and antifields. Schematically,

$$
\begin{gathered}
S_BV = S0 + \int \phi *_a L_c \phi ^a \\
+ \int g*^{\mu \nu } L_c g_\mu \nu \\
+ \int c*_\mu c^\nu \partial _\nu c^\mu .
\end{gathered}
$$

**ASSUMED** Here g* denotes a BV antifield, while G[phi] denotes the composite metric. Set the metric antifield to zero and solve the algebraic metric equation. The proposed reduced master action is

$$
\begin{gathered}
S_{red},BV = S0[G(\phi ),\phi ] \\
+ \int \phi *_a L_c \phi ^a \\
+ \int c*_\mu c^\nu \partial _\nu c^\mu .
\end{gathered}
$$

**DERIVED** This identifies a concrete classical cohomology comparison rather than an unsupported assertion of full quantum equivalence. The scalar BRST transformation is independent of the auxiliary metric, which makes this map particularly direct. Completing the gauge-fixed quantum construction still requires the measure, regulator, integration domain, boundaries, and transformed observables.

**CITED** A candidate gauge fermion can depend on the composite metric G(phi) rather than on the independent integration variable g. Its gauge condition may be chosen to have harmonic form around a regular flat-metric patch. Because the gauge fermion then contains no independent g, it is designed to preserve algebraic metric elimination, following the objective of [1]. This candidate has not been established as a globally admissible gauge, especially on the actual tuned quantum saddle. Local invertibility at nonzero momentum would not resolve residual zero modes or global gauge copies. — eq./page UNVERIFIED.

**DERIVED** The intended next comparison is therefore a gauge-fixed generating functional with matched measures and matched observables. The target is not equality of two independently chosen flat measures.

## Process-of-elimination ledger

| Proposed shortcut or route | Published test | Current disposition |
|---|---|---|
| **CITED** Eliminate the metric classically and declare quantum equivalence | [1], [3] | Eliminate this shortcut; quantum measure and sources remain necessary. | — eq./page UNVERIFIED.
| **CITED** Use BRST cohomology equality as a proof of unitarity | [2], [4] | Eliminate this inference; these conclusions are different. | — eq./page UNVERIFIED.
| **CITED** Apply an all-orders supersymmetric auxiliary-field identity to the present scalar model | [5] | Not applicable without the required supersymmetric structure. | — eq./page UNVERIFIED.
| **CITED** Obtain a conventional composite graviton while keeping all the hypotheses of the Weinberg-Witten obstruction | [6] | Excluded under those hypotheses; document a real evasion rather than ignore it. | — eq./page UNVERIFIED.
| **CITED** Retain ordinary local bulk observables and locality while demanding nonlinear universal gravitational energy coupling | [7] | Constrained by the theorem's shared-time/shared-locality hypotheses; not an unrestricted matter-only no-go. | — eq./page UNVERIFIED.
| **CITED** Recover only a Newtonian force or one massless pole and declare Einstein gravity | [8], [9] | Insufficient; interactions and universal coupling must also be addressed. | — eq./page UNVERIFIED.
| **CITED** Use a thermodynamic derivation of Einstein's equation as a complete microscopic quantum theory | [10], [11] | Insufficient; it establishes a conditional macroscopic bridge. | — eq./page UNVERIFIED.
| **CITED** Continue the constrained, clock-free matter-only construction using auxiliary-field/BV methods | [1]-[4], [8] | Retain as the most directly relevant next route; quantum matching is open. | — eq./page UNVERIFIED.

**DERIVED** "Eliminate" here means reject a stated inference or an assumption set, not claim to have disproved an entire research field.

## Published constraints on the remaining possibilities

**CITED** **[5] Higashijima and Nitta (2000), Progress of Theoretical Physics 103, 833-846.** They establish cancellation of quantum corrections to the auxiliary Legendre transformation in supersymmetric theories. The current scalar proposal does not specify that structure. This result is therefore a useful applicability control, not a ready-made proof for us. Publisher abstract and publication details checked. — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1143/PTP.103.833 — eq./page UNVERIFIED.

**CITED** **[6] Weinberg and Witten (1980), Physics Letters B 96, 59-62.** The massless-spin obstruction relies on the stipulated Lorentz-covariant conserved energy-momentum/current structures and their matrix elements. A purported composite-gravity construction must identify which hypothesis does not apply. Publisher abstract checked; the full original proof was not retrieved in this batch. The CEV authors discuss their vanishing-total-stress construction as such an evasion, not as a contradiction of the theorem. — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1016/0370-2693(80)90212-9 — eq./page UNVERIFIED.

**CITED** **[7] Marolf (2015), Physical Review Letters 114, 031104.** The published result constrains emergence with universal energy coupling from kinematically local systems. The inspected author manuscript states the shared time evolution and compatible locality assumptions explicitly. Under those conditions bulk local observables freeze. This is a filter on how microscopic observables map to gravitational ones; changing a dispersion relation alone does not address it. — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1103/PhysRevLett.114.031104 — eq./page UNVERIFIED.

**CITED** Full manuscript: https://arxiv.org/pdf/1409.2509 — eq./page UNVERIFIED.

**CITED** **[8] Carone, Erlich, Vaman (2017), JHEP 03, 134; Carone, Claringbold, Vaman (2018), Physical Review D 98, 024041.** The first supplies the constrained matter-only candidate and its leading-large-N graviton; the second tests cubic graviton interactions and universal couplings in its stated regime. These are existing results to reuse, not discoveries from our recent numerical batches. They do not establish global stability of the auxiliary model we tested. Publisher records checked, alongside the earlier full 2018 manuscript reading. — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1007/JHEP03(2017)134 — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1103/PhysRevD.98.024041 — eq./page UNVERIFIED.

**CITED** **[9] Deser (1970), General Relativity and Gravitation 1, 9-18.** His first-order self-coupling construction derives Einstein nonlinearities from a suitable massless spin-2 starting point and consistency requirements. This supplies an established route for the classical low-energy interaction target. It does not construct the microscopic composite state or prove ultraviolet unitarity. Full author-posted original inspected. — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1007/BF00759198 — eq./page UNVERIFIED.

**CITED** Full manuscript: https://arxiv.org/pdf/gr-qc/0411023 — eq./page UNVERIFIED.

**CITED** **[10] Jacobson (1995), Physical Review Letters 75, 1260-1263.** Horizon-area entropy, local Rindler horizons, Unruh temperature, and the Clausius relation yield an Einstein equation of state. This is relevant to the user's thermodynamic/pressure intuition, but does not justify an arbitrary density-pressure force. The required geometric and thermal assumptions must be supplied. Full author manuscript inspected. — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1103/PhysRevLett.75.1260 — eq./page UNVERIFIED.

**CITED** Full manuscript: https://arxiv.org/pdf/gr-qc/9504004 — eq./page UNVERIFIED.

**CITED** **[11] Jacobson (2016), Physical Review Letters 116, 201101.** The entanglement-equilibrium argument is explicitly restricted to first-order vacuum-state variations, with an additional conjectural input for nonconformal fields. Our massive constituent model is not automatically covered by the conformal special case. Published-version abstract and revision notes checked. — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1103/PhysRevLett.116.201101 — eq./page UNVERIFIED.

**CITED** https://arxiv.org/abs/1505.04753 — eq./page UNVERIFIED.

## Consequences for the earlier scalar-pole findings

**DERIVED** The previous numerical results remain conditional statements about their specified auxiliary kernels. They are not yet a demonstrated spectrum of the constrained matter-only quantum theory. Published equivalence methods make the missing comparison more precise; they do not supply a license to transfer every off-shell metric correlator unchanged.

**DERIVED** In particular, the observed pole cancellation in one matter channel cannot prove stability, while a pole in an unmatched auxiliary-field correlator cannot alone establish instability of the original constrained theory. The correct task is a comparison of matched physical observables with an explicit quantum definition. A successful regular reformulation also cannot remove a genuine physical instability merely by changing variables.

**DERIVED** The conditional linearized-curvature response derived in the preceding batch should keep that qualification. A linearized gauge-invariant expression is not automatically a globally defined, diffeomorphism-invariant quantum observable of a different formulation.

## Work queue selected by this elimination

1. **DERIVED** Complete the classical BV elimination map on a stated nonsingular branch and specify boundary conditions. Do not repeat scalar bubble scans to establish this.
2. **DERIVED** Construct an auxiliary-preserving gauge-fixing prescription and the corresponding regulated quantum measure. A general theorem's existence claim must not be substituted for the actual model's measure.
3. **DERIVED** Match matter-scattering insertions and, separately, any composite-curvature observables. Distinguish on-shell amplitude statements from off-shell generating functions.
4. **CITED** If using canonical constraint quantization, classify first- and second-class constraints before inserting any determinant. The fixed-velocity matrix $\delta T/\delta g$ is not, by itself, the full Poisson-bracket constraint matrix. Senjanovic's published construction is the relevant framework, not permission to assume every metric component is second-class. — eq./page UNVERIFIED.
5. **DERIVED** Only after that matching, decide whether the auxiliary scalar pole transfers to the constrained theory. If it does, test remedies within the specified matter-only action rather than adding unexplained geometric terms.

**CITED** Canonical-method reference: P. Senjanovic, Annals of Physics $100 (1976), 227-261$, https://doi.org/10.1016/0003-4916(76)90062-2 . Bibliographic record/abstract checked; full original text not retrieved. — eq./page UNVERIFIED.

## Computation and provenance note

**DERIVED** Immediately before the user redirected this turn toward literature, a standard-library local determinant/source sanity check was completed in constrained_measure_audit.py. It contains twelve nonsingular metric cases, twenty linear harmonic-gauge identities, and six local Gaussian source-contact checks. Those finite-dimensional checks are supplementary project calculations; none is counted as a published theorem or a quantum-equivalence test. No further numerical physics scans were run after the literature-first instruction.

**DERIVED** The local test suggests field-dependent measure factors and local source contacts must be tracked; it cannot select a canonical measure. Full scripts and JSON can be retained as a checkpoint, but the recommended next work is analytic application of the published formalism.

## Current status

**DERIVED** The research has a more precise route and several invalid shortcuts removed. Local classical auxiliary-field equivalence has a known theorem-based framework. Quantum measure/source matching, Lorentzian contours, physical stability, and ultraviolet completion remain open. No novelty, full unification, or experimental confirmation is claimed.

## Further source check: background independence and brane methods

**CITED** **[12] Batz, Erlich, Mrini, Background-Independent Composite Gravity, Classical and Quantum Gravity (2021), DOI below.** The full author manuscript, especially the introduction and conclusions, was inspected. The metric is a composite, not a separate integration variable; its vacuum expectation obeys a self-consistency condition giving Einstein equations with corrections. A physical ultraviolet regulator is still assumed rather than constructed. This supports using the paper's matter-only definition as the primary target, instead of assuming that unrestricted auxiliary-metric integration defines it. It also cautions against treating a coordinate-local curvature expression as a fully defined physical observable. — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1088/1361-6382/abefbd — eq./page UNVERIFIED.

**CITED** Full manuscript: https://arxiv.org/pdf/2011.06541 — eq./page UNVERIFIED.

**CITED** **[13] Kluson (2020), Physics Letters B 810, 135753.** Section 3 and 3.1 present a general-background p-brane volume action and its auxiliary-worldvolume-metric formulation. This is a classical Hamiltonian construction, not a quantum-equivalence theorem. Full author manuscript inspected. It supplies a formula template for the algebraic bridge in the companion brane_formula_bridge.md note. — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1016/j.physletb.2020.135753 — eq./page UNVERIFIED.

**CITED** Full manuscript: https://arxiv.org/pdf/2004.14654 — eq./page UNVERIFIED.

**CITED** **[14] Lee (1998), Physical Review D 57, 5134.** The bosonic D-brane Hamiltonian analysis explicitly separates first- and second-class constraints. The action contains a worldvolume gauge field and the displayed analysis uses a flat target background, so its full constraint list is not automatically our scalar-potential model's list. The full author manuscript was inspected for its setup and primary constraints; the publisher abstract confirms the canonical separation. This is a reusable method, not a direct solution to the present quantum measure. — eq./page UNVERIFIED.

**CITED** https://doi.org/10.1103/PhysRevD.57.5134 — eq./page UNVERIFIED.

**CITED** Full manuscript: https://arxiv.org/pdf/hep-th/9708093 — eq./page UNVERIFIED.

## Regulator filter: what has not become a physical completion

| Existing tool | Useful role | Why it does not complete this program by itself |
|---|---|---|
| **DERIVED** Dimensional regularization | Covariant perturbative integrals and the established large-N calculations | A finite dimensional parameter used as a placeholder does not specify microscopic unitary dynamics. |
| **FITTED** Signed Pauli-Villars determinants | Covariant regulator comparisons and finite response kernels | Their regulator states/weights cannot simply be treated as an ordinary positive-norm physical spectrum. |
| **DERIVED** A fixed local lattice | A finite calculation with a specified measure | Exact gravitational constraints, emergent Lorentz symmetry and the map of observables still need proof. No blanket exclusion of all lattice constructions is claimed. |
| **DERIVED** A classical brane reformulation | A way to reuse established constraint machinery | It does not choose the quantum measure or the real Lorentzian branch. |
| **DERIVED** $A BRST/BV$ reduction | A controlled local gauge-equivalence framework | It must be extended to a matched quantum definition and admissible observables before transferring the spectrum. |

**CITED** This is an applicability filter, not a new regulator no-go theorem. The relevant published limitations are in Carone's review, arXiv:1911.13166, Section V, and the physical-regulator qualification in [12]. A model-specific physical regulator remains an open construction; more comparisons of unphysical regulators cannot establish one. — eq./page UNVERIFIED.
