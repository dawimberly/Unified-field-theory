# Auxiliary integration, reduced measures and source matching
Research batch: 2026-10-29 US Central; directory date is UTC.
Baseline: aaa13ebfcfff2284a9d62e964bab1d1fd4c48a33.
This batch applies published quantization principles and checks one previously uncomputed model-specific Hessian. No old pole scans were repeated.

## Published inputs and access limits
1. Henneaux, Elimination of the auxiliary fields in the antifield formalism, Physics Letters B 238 (1990), 299-304, DOI https://doi.org/10.1016/0370-2693(90)91739-X . The publisher abstract establishes the existence of gauge fermions preserving auxiliary elimination and discusses induced measures. Only that abstract and bibliographic record were accessible; the full proof was not obtained. It cannot supply this model's actual measure or contour by assertion.
2. Criado and Perez-Victoria, Field redefinitions in effective theories at higher orders, JHEP 03 (2019) 038, DOI https://doi.org/10.1007/JHEP03(2019)038 . Published open-access PDF inspected at Section 2, equations (2.1)-(2.6), printed pages 5-7: https://inspirehep.net/files/84c3b1b40fac77e9429e45c460412403 . Exact generating-functional changes require transformed sources and Jacobians. Their local perturbative dimensional-regulator simplification and on-shell equivalence have stated restrictions. They do not identify an integration-out operation with an invertible redefinition.
3. Carone, Composite Gravitons from Metric-Independent Quantum Field Theories, author manuscript https://arxiv.org/pdf/1911.13166 . Footnote 4 explicitly leaves unrestricted-auxiliary/constrained quantum equivalence unestablished. Section IV, equations (4.1)-(4.2), removes the clock fields and organizes the leading-N perturbation series using cancelling background/counterterm contributions. Those statements are not a specification of a nonperturbative real Lorentzian metric-integration cycle.

The following finite-regulator formulas and algebraic limitations are our application. No source above is credited with proving their applicability to a complete CEV quantum measure.

## 1. An invertible shift is distinct from elimination
On a regular classical branch, take
S0[g,phi]=integral sqrt(|g|)[g^(mu nu)P_mu nu/2 -V(phi)],
P_mu nu=sum_a partial_mu phi^a partial_nu phi^a,
G_mu nu[phi]=(D-2)P_mu nu/(2V).
Assume D!=2, V!=0, nonsingular P, a specified finite regulator and a defined integration domain/cycle. Set g=G(phi)+h while leaving phi unchanged.

At finite dimension the Jacobian matrix of (phi,h)->(phi,g) is block triangular with identity diagonal blocks, hence its determinant is one, even when G depends on finitely regulated derivatives of phi. This says nothing about the value of the subsequent h integral, the transformed cycle, the original measure density or the admissibility of gauge fixing. Eliminating h is not itself an invertible field redefinition.

Let Sred(phi)=S0[G(phi),phi]. For sources j for phi and J for an independent g insertion, define
I_phi[J]=integral_C_h(phi) dh mu_g(G+h,phi)
 exp{i[S0(G+h,phi)-Sred(phi)+J.h]/hbar}.
Then, with the same regulator and normalization,
Zaux[j,J] proportional to integral dphi mu_phi(phi)
 exp{i[Sred+j.phi+J.G(phi)]/hbar} I_phi[J].

The zero-source induced reduced measure is mu_eff=mu_phi I_phi[0]. A normalization Z[0]=1 removes a single overall constant, not this field-dependent weight. Matching mu_red to mu_eff is required for equality of matter-source generating functionals. Starting instead with a prescribed mu_red and the determinant-weighted delta constraint is an alternative exact representation, as in the imported matched_constraint_identity.md; it does not show that mu_red=mu_eff.

For composite sources, equality additionally requires I_phi[J]/I_phi[0]=1 or compensating matched insertions/source terms. If the auxiliary source couples directly to G(phi), rather than g, no J.h is introduced and the matter/composite generating functional with mu_eff matches by construction. That observable is different from an unmatched independent-g insertion.

In the determinant-weighted constraint representation, the delta constraint is the source-free T(g,phi)=0. An insertion J.O does not authorize replacing it by the stationary equation of a different source-dependent action.

## 2. Source contact term and a useful covariance identity
For a local saddle expansion with nondegenerate algebraic Hessian H, a specified Gaussian cycle gives formally
I_phi[J]/I_phi[0]=exp[-i J.H^(-1).J/(2hbar)]
at Gaussian order, and the conditional connected covariance is i hbar H^(-1). There can be additional local source cumulants and a quantum shift of the conditional mean beyond that order.

If an auxiliary-preserving gauge, regulator, measure and cycle factorize the h integrations between cells at fixed phi, let m_x(phi)=conditional_mean(g_x|phi) and C_xy(phi)=conditional_connected_covariance(g_x,g_y|phi). The exact algebraic law of total covariance is
Cov(g_x,g_y)=Cov_phi(m_x,m_y)+E_phi[C_xy].
It applies formally to normalized complex weights where the indicated conditional integrals exist; it does not require probabilistic positivity. Cell factorization makes C_xy a contact in x,y. However m(phi) need not equal the classical G(phi), and the phi weight is mu_eff. Those two changes can affect nonlocal correlators. A contact polynomial alone cannot cancel an isolated nonlocal pole if the measure and nonlocal insertions are otherwise already matched.

No equality of the previously computed independent-metric response with a constrained composite-metric correlator follows from the shift or Gaussian formula alone.

## 3. Actual model Hessian: a positive Gaussian is not automatic
At the branch solution, P=kG with k=2V/(D-2). Expanding covariant g=G+h, M=G^(-1)h, gives the quadratic density
S0^(2)=sqrt(|G|) V/[2(D-2)] {tr(M^2)-(tr M)^2/2}.
Thus
H[h,j]=sqrt(|G|) V/(D-2)
 {tr(G^(-1)h G^(-1)j)-tr(G^(-1)h)tr(G^(-1)j)/2}.

On a definite D=4 patch normalized to G=I, in the orthonormal symmetric basis,
H=(V/2)[I-(v v^T)/2], v=(1,1,1,1,0,...,0).
There are nine traceless eigenvalues V/2 and one trace eigenvalue -V/2. The Hessian is nondegenerate for V!=0 but indefinite. Therefore a convergent positive real Gaussian exp(-h.H.h/2) is not supplied by this algebraic saddle. One must state a cycle/prescription before using that integral. This auxiliary trace direction is not, by itself, a physical propagating ghost diagnosis.

For T_mu nu=P_mu nu-g_mu nu[g^(rho sigma)P_rho sigma/2 -V], at G=I and P=VI,
J_T[h]=-V[h-I tr(h)/2]=-2H[h].
The orthonormal-basis determinants are det H=-(V/2)^10 and det J_T=-V^10.
The leading stationary-phase magnitude of the metric integral is proportional to |det H|^(-1/2)=(2/|V|)^5, times the specified measure density and universal hbar factors. If V depends on phi, this is not generically an ignorable constant. By contrast the determinant-weighted delta constraint cancels its own constraint-change Jacobian exactly. These are different operations.

These H and J_T calculations are classical fixed-phi algebraic quantities; H is not the loop-resummed tensor kernel K2 or A from the earlier determinant calculations, and J_T is not the complete canonical constraint Poisson matrix.

## 4. Correction to the classical branch qualification
For real fields with positive flavor metric,
u^mu P_mu nu u^nu=sum_a (u^mu partial_mu phi^a)^2>=0.
Therefore raw P is positive semidefinite and, on a nondegenerate patch, positive definite. Multiplication by the nonzero scalar (D-2)/(2V) makes G definite with one overall sign. In D>2 it cannot produce a real Lorentzian inertia with both timelike and spacelike directions.

Accordingly, the imported note's clock-free algebraic map cannot simply be assumed to define a real Lorentzian classical branch for this positive-flavor raw Gram matrix. Its nonsingular real branch is definite; a Lorentzian interpretation needs an independently justified quantum/renormalized-composite construction or complex continuation. Clock fields with an indefinite target metric alter the premise and are not being introduced here.

This is not a no-go theorem for a Lorentzian quantum expectation, a normal-ordered operator, the clock-containing CEV action, or Carone's perturbative construction. A renormalized composite is not the same pointwise positive raw matrix, and a Lorentzian expectation with an oscillatory/continued measure cannot be analyzed as an ordinary positive probability average. No such quantum construction is completed in this batch.

## 5. Conditions still needed before physical spectrum transfer
A concrete reduced theory must supply its regulator, measure, gauge quotient, cycle/state/boundary data and renormalized observables. A comparison must then specify the induced measure and source operators just derived, with a common admissible gauge and regulator. A composite-dependent gauge fermion can preserve algebraic g elimination because its BRST variation can avoid independent g, but invertibility, residual modes, anomalies and the quantum master equation are still model-specific obligations.

These are technical research targets rather than a request to select unexplained new physical fields or couplings. The current auxiliary determinant remains a conditional comparison calculation. Its pole results are not discarded, but they are not automatically the constrained spectrum.

## 6. Reproducible new calculations
Run python experiments/auxiliary_source_matching_2026_10_30/check.py (standard library only).
Nine 70-digit directional second derivatives of the actual D=4 algebraic S0 at G=I, for V=1,4,9 and trace/traceless/off-diagonal directions, agree with the displayed Hessian within 6.75e-24. Exact fraction determinants verify det H and det J_T at all three V values. Assertions passed on Python 3.12.14. The code records the limited definite-patch scope; it constructs no quantum cycle or vacuum.

Next: inspect the precise renormalized composite and gauge/measure definition in the background-independent clock-free literature, rather than infer it from the raw classical Gram matrix or substitute the old unrestricted determinant. A source-matched large-N comparison remains a missing calculation; no fundamental model change or human choice has been made.
