# Composite-source definitions and the covariant matching dictionary

Date: 2026-10-30 US Central. Baseline: fe85718a13098c30b01be0cd0019dc67a8d0194d.

The source-definition audit does not justify transferring the auxiliary PV pole inventory to a constrained/composite theory. The bounded result is a precise separation of what a vacuum/tadpole condition can determine from the extra finite source and curvature matching data. An explicit covariant dictionary for the already used matching polynomials is derived below. No new physical counterterm value, source prescription, contour or interaction was selected.

## Primary-source reading, with access limits

Batz, Erlich and Mrini: arXiv:2011.06541v1, https://arxiv.org/pdf/2011.06541 , associated DOI 10.1088/1361-6382/abefbd. Read PDF pp.5-16 and corresponding HTML equations. Its clock terms cancel, normal ordering subtracts expectation values, and the displayed effective-action/bubble calculation keeps low adiabatic orders while omitting higher-derivative contributions. It distinguishes free-field stress from total stress including vacuum terms. Section IV.A discusses contact cancellation in its convention. These statements do not specify a common finite-q PV composite-source functional. The arXiv record lists v1 only; the publisher endpoint was blocked, so no final-version equality or correction claim is made. Target locators are recorded in sources.json.

Hollands and Wald: Commun. Math. Phys. 223 (2001) 289-326, DOI 10.1007/s002200100540, author version https://arxiv.org/pdf/gr-qc/0103074 (v2). Read Theorems 5.1/5.2 and equation (81). Under their local covariance, scaling and regularity assumptions, Wick products and time-ordered products retain finite renormalization parameters; the Wick square can shift by curvature and mass terms. Their framework is not a proof that this signed Euclidean PV model satisfies its Lorentzian hypotheses or that a composite-gravity quantum measure exists.

These limited source synopses do not import an unproved finite-momentum equivalence. The following calculations and consequences are derived here in independently specified source conventions.

## 1. First and second source variations contain different information

Fix a source functional W[g,J]. A condition on its first derivative at J=0 does not determine its second derivative there. For example the finite local redefinition

Delta W=lambda/2 * integral sqrt(g) J^2

leaves W[g,0] and its first J derivative unchanged, but shifts the coincident two-source contact by lambda. This is a source-functional example, not a newly adopted microscopic interaction or source value. Specifying a one-point subtraction is therefore insufficient to define the corresponding product of two composite operators.

Likewise, at a flat metric the first variation of integral sqrt(g) R^2 and integral sqrt(g) R_ab R^ab vanishes, whereas their Hessians need not vanish. Even a fixed flat vacuum and a fixed Einstein slope leave two independent curvature-squared directions available. A particular regulator prescription can set values for them; comparing two prescriptions as one physical theory requires those values and source terms to be matched explicitly.

The repository's raw PV prescription remains fixed and its previously computed numbers are unchanged. This freedom concerns missing information in a proposed cross-theory identification, not uncertainty about which formula the saved scripts evaluated.

## 2. Exact covariant dictionary for the existing metric polynomials

Use D=4, Euclidean g=delta+h and q=Q^2. Fix curvature convention

R^(1)(Q)=q tr(h)-Q_a Q_b h_ab,
R_ab^(1)(Q)=[q h_ab+Q_a Q_b tr(h)-Q_a Q_c h_cb-Q_b Q_c h_ca]/2.

Consider the bookkeeping functional

Gamma_ct=integral sqrt(g)[kappa R+alpha R^2+beta R_ab R^ab].

It represents existing local metric-kernel changes; it is not an instruction to add or tune a new action. Work with compact-support or periodic perturbations. Euler and integrated total-derivative terms do not affect this flat quadratic dictionary. K denotes the Hessian in S2=h:K:h/2; the scalar conformal coefficient remains A=3K0/2.

Expanding the integrated Einstein term through the first-order Christoffel products and the quadratic invariants through R^(1),R_ab^(1) gives

Delta K2=-kappa q/2+beta q^2/2,
Delta A=3 kappa q/2+(9 alpha+3 beta)q^2.

Equivalently, if the already specified polynomial increments are Delta A=delta_a1*q+delta_a2*q^2 and Delta K2=delta_k1*q+delta_k2*q^2, covariance requires delta_a1=-3 delta_k1 and

kappa=-2 delta_k1,
beta=2 delta_k2,
alpha=(delta_a2-6 delta_k2)/9.

The map between the two independent q^2 coefficients and alpha,beta is invertible. This makes the independent curvature-squared matching requirement explicit rather than identifying two separate zero-bare-term prescriptions as the same theory.

For the existing free minimal-scalar PV family, let C=(4pi)^-2, Delta T0 be the reference-minus-spectrum change in sum c log z, and Delta T1 the change in sum c z log z. The saved matching increments give

kappa=-C Delta T1/12,
alpha=C Delta T0/240,
beta=C Delta T0/120.

Thus this particular family's regulator differences lie on beta=2alpha. Both tensor coefficients are matched, but these regulator differences span only one direction within the two-dimensional curvature-squared space. No general physical condition forcing every allowed finite source prescription onto that line has been established. This observation does not change the original matching calculation or turn its high-q poles into low-energy predictions.

## 3. Independent exact verification of the dictionary

check.py uses only the standard library and Fraction arithmetic. It computes the Cartesian linearized scalar/Ricci curvatures and the integrated Einstein Gamma-Gamma quadratic form, then compares their mixed Hessians with the transverse spin projectors. Three rational momentum directions, five polarizations (including an arbitrary tensor, conformal, transverse scalar, transverse traceless and gauge components) and five coupling triples give 375 exact bilinear cases. Every difference is exactly zero. The tests establish the finite algebra in this convention; they are not new loop Ward checks.

The action and Hessian factors are explicit. For instance a transverse traceless h has R^(1)=0 and R_ab^(1)=q h_ab/2, so beta contributes Delta K2=beta q^2/2. A conformal h=f delta has R^(1)=3q f and Ricci norm squared 3q^2 f^2, giving Delta A=(9alpha+3beta)q^2. Those independent channels fix both coefficients.

## 4. Canonical contact audit and the limit of a cancellation argument

A separate exact calculation isolates a source-convention issue. Work in our canonical Euclidean free-scalar convention, with fixed-field tadpoles

<partial_a phi partial_b phi>=p delta_ab,
<z phi^2>=-4p,
<T_ab>=p delta_ab.

No response of the state or normal-ordering subtraction is silently differentiated into these fixed-field derivatives. For covariant metric directions h,j define P(h,j)=tr(hj)-tr(h)tr(j)/2. Converting inverse-metric derivatives to covariant metric perturbations produces two local pieces:

inverse-metric insertion contact = -p tr(h)tr(j)/2,
inverse-to-covariant reparametrization contact = p tr(hj).

Their sum is p P(h,j), not an identity equal to zero. Adding lambda integral sqrt(g) contributes -lambda P(h,j)/2. At the stationary choice lambda=p the total fixed-field local term is p P(h,j)/2, consistent with the saved determinant dictionary. Tadpole stationarity alone has not canceled it. The connected bubble must still be kept separately.

To verify this in a way independent of those contractions, the script expands the fixed-field action directly using two nilpotent variables eps,eta. It evaluates det(delta+eps*h+eta*j), its square root and the inverse metric exactly through eps*eta. Fifty mixed second-variation cases at lambda=0,p agree with the separated contacts with exact zero residual.

This supplies a counterexample to importing a contact cancellation as a universal property of stationarity in our convention. It is not a finding of an error in the published composite-gravity theory. A metric-dependent operator subtraction, total-versus-free stress definition, regulator and source counterterms can change the comparison. The accessible source's cancellation must be reproduced inside its own complete source convention before applying it to our Hessian; its unavailable final publication was not adjudicated. No paper formula was substituted for a fixed-field derivative of the present model.

## 5. What has and has not been fixed

The conditional raw source conversion from the preceding checkpoint still requires the same B,p, source map and integration assumptions. The new calculation identifies a covariant representation of its explicit metric-kernel modifications, but it does not fix which composite operator or product prescription accompanies them. Changing the metric functional by curvature terms is different from changing source-only quadratic contacts. Both operations can preserve a flat vacuum while modifying the quantities needed for a finite-momentum comparison.

A source-complete matching statement would need a single definition of: the regulated composite operator; local products and their finite counterterms; metric and mass source variations; the induced measure/Jacobian and allowed configurations; and the Lorentzian integration prescription. None of those missing choices is filled in by declaring two AI derivations consistent or by passing these algebra checks. In particular, the source J^2 example and the curvature functional above are demonstrations of remaining data, not newly adopted choices.

There is useful authorized work left inside the fixed auxiliary model, so this is not an access blocker for the research program. No human model choice is needed to save this bounded audit or to continue that work.

## Reproduction and continuation

Run `python experiments/composite_source_definition_2026_10_30/check.py`. results.json records all 375 covariant bilinears, 50 contact bilinears, exact coefficient maps, conventions and Python version. No previous threshold, pole count, loop Ward or dispersion scan was repeated. No mathematical assertion failed. The literature locators and access limitations are saved separately; no full article was copied into the repository.

Next return to a remaining mathematical gap in the already fixed auxiliary comparison: audit the global first-sheet inventory of the explicitly matched spin-2 kernels, including real-axis zeros and any cut boundary, using the established raw Stieltjes representation. Existing isolated roots and finite contours alone should not be called a global count. Keep the incomplete composite-source definition and physical Lorentzian contour as separate unresolved problems rather than silently selecting them.
