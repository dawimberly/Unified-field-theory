# Scalar source dispersion, subtractions and the conditional contact dictionary

Date: 2026-10-30 US Central. Baseline: 5a6a498fb88c02e837e755d58a5ca38ce625bc47.

This batch audits the scalar component of the specified constituent connected-source functional. It does not calculate the full tensor scattering amplitude. The previous checkpoint's shorthand 'full constituent source exchange W' should be read with this qualification: the spin-2 contribution is separate, and the two-distinct-source cross coefficient has a factor of two relative to a single-source quadratic coefficient.

The source dispersion formulas agree with new numerical evaluations: fifteen raw/reference second-subtraction checks, twelve non-reference matched absolute checks, eight source moment checks and nine Cauchy-derivative checks. Exact inherited-and-transformed UV inequalities bound the omitted source continuum. A separate exact block-Gaussian calculation tests why specifying a metric Hessian alone does not specify its composite source functional. No physical model, matter coupling, contour or source prescription was selected or changed.

## 1. Normalization from a conserved source

Use D=4, q=Q_E^2, M=1 and the completed stationary transverse metric Hessian. Coupling -h:T/2 gives the formal Gaussian connected quadratic coefficient

W_h[T]=T:K^-1:T/8.

This is the established Euclidean algebraic convention, not a prescription for an indefinite Lorentzian or Euclidean integration cycle. Restrict inverses to the conserved source sector, or fix gauge consistently before taking ordinary inverses. For a conserved source, let tau=tr T, theta=delta-Q Q/q and P0=theta theta/3. The saved scalar normalization K0=2A/3 gives

W0=tau^2/(16A).

For the specified identical minimally coupled constituent, tau=2-q, hence W0=(q-2)^2/(16A). Throughout this report W means this scalar W0. The full conserved quadratic coefficient also includes T:P2:T/(8K2). A scalar-channel sign by itself is not a diagnosis of a full Einstein exchange ghost. For two distinct sources the scalar cross term is tau1*tau2/(8A), not the single-source quadratic coefficient above.

This derivation was checked against the original saved script scripts/scalar_source_projection_check.py through the connected repository. Its old projector tests and finite-D cancellation scans were not rerun. The PV model considered here and the finite-D DR identity remain distinct.

## 2. Poles, cut and UV behavior

Let F=1/A, rho=Im F(-s-i0)/pi, and (q_j,r_j) be its complete first-sheet inventory including q=0. Then

r_j^W=(q_j-2)^2 r_j/16,
rho_W(s)=(s+2)^2 rho(s)/16.

All regulator cuts and both members of each complex pair are retained. The prior first-threshold disk theorem excludes a coefficient zero at q=2. Therefore W has a double zero there, rather than a pole cancellation at that point in these PV spectra. The known coefficient poles have nonzero source factors and survive in this specified scalar channel. This says nothing about a different external source or a constrained theory.

For non-reference matched kernels, write the added UV polynomial as a2*q^2+a1*q, so

A(q)=a2*q^2+a1*q+a_infinity+O(Log q/q).

These a1,a2 are the added polynomial coefficients, not the total matched low-energy coefficients. The density rho_W decays as s^-3. Polynomial division of the complete F dispersion, using its finite M0=0 and M1=1/a2 identities, gives

W(q)=1/(16a2)+sum_j r_j^W/(q-q_j)+integral_4^infinity rho_W(s)/(s+q)ds.

The constant is fixed by the specified kernel and source projection. It was not fitted to a reconstruction. It must not be confused with the distinct auxiliary/composite source-conversion contact below. The q=2 checks test its cancellation against the pole and continuum terms.

For raw/reference kernels, F_infinity=1/a_infinity and

W(q)=F_infinity*(q-2)^2/16+O(q Log q).

Here rho_W grows as s. A naive absolute or once-subtracted source continuum diverges. Define V=W-F_infinity*(q-2)^2/16. With q0=1 and h=q-q0, two subtractions give

W(q)-W(q0)-h W'(q0)
=F_infinity*h^2/16
+sum_j r_j^W*h^2/[(q-q_j)(q0-q_j)^2]
+integral_4^infinity rho_W(s)*h^2/[(s+q)(s+q0)^2]ds.

This follows from Cauchy's theorem applied to V=O(q Log q); its second-subtracted outer boundary vanishes. The known quadratic polynomial has not vanished under two subtractions and supplies the displayed h^2 term. The implementation uses these factored kernels to avoid subtractive cancellation in the UV. Subtraction constants W(q0),W'(q0) are computed from the unchanged Gaussian response, not chosen as new source counterterms.

## 3. Exact source tail estimates and finite source moments

For each distinct prescription use the preceding exact certificate 0<rho(s)<=K/s^p, s>=R. There are nine distinct cases: five raw and four non-reference matched; matched L=10 duplicates the raw L=10 case. All rational parent certificates are preserved inside the new results. Multiplication by the known source factor gives

0<rho_W(s)<=K_W/s^(p-2),
K_W=K*(1+2/R)^2/16.

Thus p-2=3 for the non-reference matched source, and -1 for raw/reference sources. K_W is computed exactly from the rational parent bound. No floating-point sample decides the inequality.

For the matched absolute source integral and |q|<=H<R,

|tail|<=K_W/[3(1-H/R) R^3].

For the raw/reference second-subtraction integral, with |q|,|q0|<=H<R,

|tail|<=K_W*|h|^2/[(1-H/R)^3 R].

The program uses rational component-sum upper bounds on complex moduli. These certify omitted continuum tails only; they do not certify every digit of quadrature, derivatives or saved numerical residues.

Only two ordinary non-reference matched source moments are finite. With

N_n=sum_j r_j^W*q_j^n+(-1)^n integral s^n rho_W(s)ds,

polynomial multiplication of the four already derived inverse moments gives

N0=(-4/a2-a1/a2^2)/16,
N1=(4/a2+4a1/a2^2+a1^2/a2^3-a_infinity/a2^2)/16.

Their tail bounds are K_W/[(2-n)R^(2-n)] for n=0,1. The n=2 source continuum moment diverges logarithmically. Raw/reference source zeroth and higher moments diverge more strongly; their convergent twice-subtracted representation does not authorize those moment sum rules. No positive global spectral measure follows from these signed identities.

## 4. The raw subtraction polynomial matches the conditional source contact

The saved stress-contact dictionary defines

p=-C*T2/4>0, a_infinity=3C*T2/8=-3p/2.

The earlier formal Schur bridge, with the same B and p and D=4, gave

D_h=C_O+2P/p, P=I_sym-vv^T/2.

Therefore the connected single-source quadratic coefficients differ by

W_h[T]-W_O[T]=T:P:T/(4p).

For the pure transverse scalar source T0=tau*theta/3, tr(T0^2)=tau^2/3 and (tr T0)^2=tau^2, giving

W_h,0-W_O,0=-tau^2/(24p).

The raw UV subtraction polynomial is exactly

F_infinity*tau^2/16=-tau^2/(24p).

Thus the raw scalar polynomial removed before the two-subtraction dispersion agrees with the scalar projection of the conditional auxiliary-minus-composite contact, including its sign and factor of four/eight normalization. This is an algebraic consequence of the audited determinant normalization and the conditional Schur definitions. It is not proof that a constrained CEV functional uses those definitions, the same regulator, or the same source counterterms. The spin-2 source-conversion contact is separate and contributes to the complete tensor contact.

The scalar source zero at q=2 is kinematic here; the above polynomial is a local source-conversion term. Neither should be conflated with the non-reference matched constant 1/(16a2), which instead comes from its changed UV asymptote.

## 5. Why matched metric kernels do not determine a unique composite source

The original formal bridge has metric Hessian K=H-B/4, with H=pP/2. Adding an explicitly specified metric-kernel polynomial Delta can be allocated, algebraically, in either of two ways:

A: H_A=H+Delta, B_A=B;
B: H_B=H, B_B=B-4Delta.

Both give exactly the same K+Delta and auxiliary metric covariance D_h. They need not give the same composite operator or generating functional. In the following Schur formulas X is the internal stress-channel variable, distinct from the prescribed external test tensor T above. In the general quadratic bridge,

O_H=H^-1*X/2,
C_O=H^-1*C_T*H^-1/4,
D_h-C_O=H^-1.

Allocation A changes the composite map and contact to (H+Delta)^-1; allocation B changes the stress covariance/contact prescription while retaining H^-1. In particular, a momentum polynomial in H generally has a momentum-dependent rational inverse and is not the original ultralocal source-conversion contact. This does not generate a new physical degree of freedom by itself: singularities in such a representation can cancel between its composite covariance and contact.

The report does not select either allocation as a physical composite theory. Actual local stress-source counterterms and operator renormalization must be derived from one common regulated generating functional. The equality of total metric Hessians alone is insufficient information. Full nonlinear measure and contour equivalence remain unresolved. At stationary gauge zero modes the ordinary inverse Schur algebra must first be restricted or gauge fixed; no unrestricted path-integral inverse is asserted.

counterterm_allocation.py makes this non-uniqueness concrete with a two-channel, noncommuting finite-dimensional example. It directly inverts the full 4x4 coupled block covariance, rather than merely substituting the expected Schur expression. At q=0,1,2,3 it checks both allocations: eight block inversions and twelve rational source contractions. The metric covariance agrees exactly, the composite covariances differ at nonzero q, and their difference equals the difference of conversion contacts. Maximum residual is exactly zero in Fraction arithmetic. This is a formal matrix illustration, not the actual PV counterterm tensor, a positive measure or a physical path integral.

## 6. Primary-literature convention check

Bzowski, McFadden and Skenderis, JHEP11(2018)153, DOI 10.1007/JHEP11(2018)153, published PDF https://inspirehep.net/files/a9e0845c8d7d1acb1809e3504cd3c00e . Inspected Appendix C, C.1-C.6: changing the placement of source-dependent metric factors relative to functional differentiation changes contact contributions. This supports specifying the entire source convention; it does not provide the massive PV formulas or establish composite-gravity equivalence. All dispersion and finite-matrix formulas above are derived here. The established mathematical input remains Cauchy's residue theorem; it is applied with the stated asymptotic subtractions.

## 7. Execution, limitations and continuation

Run:

1. `python experiments/source_dispersion_audit_2026_10_30/check.py`
2. `python experiments/source_dispersion_audit_2026_10_30/counterterm_allocation.py`

New execution results:

- Fifteen raw/reference second-subtraction reconstructions: maximum normalized closure 8.61e-14. Omitting the known quadratic remainder produces error up to 0.0235 in the tested normalization.
- Twelve non-reference matched absolute source reconstructions: maximum normalized closure 1.69e-14. Omitting the source constant produces error up to 265; this includes q=2 where the direct source coefficient is zero and errors are absolute.
- Eight finite source moment checks: maximum component-normalized discrepancy 2.44e-16; largest absolute discrepancy about 3.23e-9. These are empirical numerical agreements, not interval-certified digits.
- Nine newly required derivatives at q0=1, checked by independent 32/64-node Cauchy integrals: maximum normalized disagreement 4.70e-16.
- Eight exact block-matrix cases and twelve source forms: exact zero residual, with no positivity or convergence assertion.

For reconstruction and derivative comparisons normalization is max(1,|direct|); moment normalization is the sum of the magnitudes of pole, cut, tail and target terms, bounded below by one. Every threshold and the stable-UV switch is split in quadrature; tails use x=R/s. No quadrature warnings or failed mathematical assertions occurred. The recorded long-double mantissa is 63 bits. Exact tail certificates and finite Gaussian algebra are separate from numerical continuum and residue accuracy. Prior passing Ward, threshold, root-count and inverse dispersion calculations were not repeated.

No discovery, physical UV completion or constrained-measure equivalence follows. Next examine the published composite-source and operator-renormalization prescriptions directly, especially their clock/metric variations and local counterterms, to determine what finite-momentum source matching is actually specified and what remains an additional assumption. Do not choose a new composite measure or identify either illustrative counterterm allocation with the physical model without justification. This literature/source-definition work can continue without a new physical choice.
