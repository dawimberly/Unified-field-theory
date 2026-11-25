# Applying a known brane formula to the scalar determinant action

**Conversation date range** 2026-10-29 to 2026-10-30

**DERIVED** This checkpoint records Applying a known brane formula to the scalar determinant action. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Companion to gravity_literature_elimination.md. October 29, 2026.

**DERIVED** This is an algebraic identification and an applicability audit. It is not a new quantum theory, an independent unification proof, or a novelty claim. The composite-gravity authors already identify the resemblance to scalar Dirac-Born-Infeld dynamics.

## The existing formula to reuse

**DERIVED** The induced-volume action for a D-dimensional worldvolume has the form

    S_volume = integral d^D x sqrt(|det h_mu nu|),
    h_mu nu = H_ab(phi) partial_mu phi^a partial_nu phi^b.

**DERIVED** The target metric H_ab is allowed to depend on the embedding coordinates phi. For a regular conformal target metric H_ab=Omega(phi)^2 delta_ab,

    h_mu nu = Omega^2 P_mu nu,
    sqrt(|det h|) = |Omega|^D sqrt(|det P|).

**DERIVED** These determinant identities are exact at the classical algebraic level.

## Matching our action

**DERIVED** Write a=D/2-1. The clock-free scalar action used in the CEV family is

    S_CEV = integral d^D x (a/V(phi))^a sqrt(|det P|),
    P_mu nu = sum_a partial_mu phi^a partial_nu phi^a.

**DERIVED** On a real patch with D>2 and V>0, choose

    Omega(phi) = (a/V(phi))^(a/D),
    H_ab(phi) = (a/V(phi))^((D-2)/D) delta_ab.

**DERIVED** Then |Omega|^D=(a/V)^a, so the two volume densities coincide. In D=4 the target metric is simply

    H_ab(phi) = V(phi)^(-1/2) delta_ab.

**DERIVED** The potential has been absorbed into a curved field-space metric. No new spacetime field has been introduced. A constant tension or reference scale can be restored to fix units; it does not alter this local determinant identity. At V=0 this mapping is singular. Nonpositive or complex branches require an explicitly different branch prescription and are not included in the displayed real formula.

## What this identification buys us

**DERIVED** It identifies an existing family of classical constrained systems whose auxiliary-metric and Hamiltonian techniques can be reused. The microscopic scalar fields act as embedding coordinates into field space in this representation; the identity does not assert that extra physical spatial dimensions have been discovered.

**DERIVED** The concrete next application is to import the general curved-target brane variational and constraint machinery, specialize its target metric to the expression above, and derive the first/second-class split before constructing a constrained phase-space measure. This avoids guessing a determinant from the configuration-space equation T=0.

## Four transfers that are not justified

1. **DERIVED** **Signature:** real fields with delta_ab give P as a positive-semidefinite Gram matrix. A positive conformal factor cannot turn this into a classical Lorentzian induced metric. The real V>0 algebraic identification is therefore a volume identity, not a construction of a real timelike brane. The composite-gravity quantum vacuum involves regulated composite expectations; that separate construction must not be replaced with an unjustified signature assignment. This observation is not a proof excluding every quantum realization.
2. **DERIVED** **Quantum measure:** equal classical volume densities do not identify a flat scalar measure with a brane measure, an induced surface measure, or the measure obtained after auxiliary-metric integration.
3. **DERIVED** **Worldvolume dimension:** D=4 corresponds to a 3-brane. Special results relying on the Weyl symmetry of a two-dimensional string worldsheet do not transfer automatically. A four-dimensional auxiliary metric has different trace constraints.
4. **DERIVED** **Matter content/background:** flat-target D-brane constraints with a gauge field are a methodological reference, not the final constraint list for a curved-target scalar-only model. The field-dependent potential is now target curvature and has to remain in the calculation.

**DERIVED** Thus the process-of-elimination verdict is: retain the classical curved-target volume map; reject automatic transfer of timelike signature, string quantization, or the complete D-brane quantum measure.

## Where the published input ends

**DERIVED** Kluson, *Note about covariant Hamiltonian formalism for strings, p-branes and unstable Dp-branes*, Physics Letters B 810 (2020), 135753, Section 3, Eqs. (12)-(21), supplies the general brane and auxiliary-metric templates. The conformal specialization above is our direct algebraic application. The source does not claim a quantum completion of CEV.

**DERIVED** https://doi.org/10.1016/j.physletb.2020.135753

**DERIVED** https://arxiv.org/pdf/2004.14654

**DERIVED** Lee, *Hamiltonian formulation of the D-brane action and the light-cone Hamiltonian*, Physical Review D 57 (1998), 5134, illustrates why lapse/shift and spatial auxiliary data require a genuine constraint classification. Its setting differs from ours as stated above.

**DERIVED** https://doi.org/10.1103/PhysRevD.57.5134

**DERIVED** https://arxiv.org/pdf/hep-th/9708093

**UNVERIFIED** Carone, Erlich, Vaman, JHEP 03 (2017), 134, and Batz, Erlich, Mrini, *Background-Independent Composite Gravity*, provide the composite-gravity starting action. Their prior brane analogy is acknowledged; no originality is claimed for identifying a determinant-volume structure.

**DERIVED** https://doi.org/10.1007/JHEP03(2017)134

**DERIVED** https://arxiv.org/pdf/2011.06541
