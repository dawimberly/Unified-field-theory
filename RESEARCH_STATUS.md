# Research checkpoint

## Goal
Derive a healthy universal gravitational response from quantum dynamics, recovering Einstein gravity at low energy. No solution or discovery is established.

## Current calculational scope
Free, minimally coupled, leading-N auxiliary-metric determinant with stationary cosmological contacts. The raw baseline has no bare curvature terms; explicitly labeled regulator comparisons add the previously specified linear/quadratic scalar and spin-2 matching polynomials, without adopting a physical UV action. This is not yet demonstrated to be equivalent to the constrained CEV path integral with its clock fields and measure. Constituent M=1. Four-species PV spectra z=[1,L,2L,5L+1], L=3,5,10,30,100, with first three signed moments vanishing.

## Completed checks
- Full momentum scalar kernel in finite D has an exact zero at Q_E^2=4M^2/(D-2), not a derivative-expansion artifact.
- Identical on-shell constituent scalar source traces vanish at that same zero, cancelling it in that exchange channel. Unequal-mass probes need not cancel it; vacuum stability does not follow.
- Three raw PV spectra have subthreshold scalar poles with positive connected Euclidean source residues. Exact coefficient/sign arguments now give complete extra-zero counts within |q|<4: one for raw L=3,5,10, none for raw L=30,100; one for all five explicitly Lref=10 matched spectra. The later cut-boundary argument now gives complete global counts for these five raw and explicitly matched spectra; see the latest batch.
- Full 4D 10x10 metric Hessian, with bubble and stationary contacts, passes longitudinal Ward checks in 40 cases (maximum absolute residual about 1.9e-13). Scalar projection K0=2A/3 agrees.
- Sixteen finite-D longitudinal/trace tensor checks agree with the corresponding scalar kernel.
- Physical-only continuum scalar weights are positive analytically below the first PV threshold, confirmed in 750 samples. The spin-2 inverse cut sign is positive there too.
- Above regulator thresholds each tested PV spectrum develops negative scalar spectral weight. PV is not thereby a physical unitary UV completion.

- Fixed four-species no-bare-curvature spin-2 determinant: the Genocchi-Hermite representation makes K2/q a positive Stieltjes transform and excludes additional first-sheet spin-2 poles. This does not extend automatically to modified kernels or constrained CEV.

## Open problems
1. Determine what finite-momentum composite-source and operator-renormalization prescription is actually fixed by the primary composite-gravity literature, before transferring auxiliary poles or matched kernels. The scalar constituent connected-source dispersion/contact audit is complete within the stated auxiliary approximation; it is not the full tensor scattering amplitude. Other spectra, adjacent sheets and physical interpretation remain separate. The fixed-prescription spin-2 coefficient has an analytic positive-measure zero-exclusion result.
2. Lorentzian contour from the microscopic metric path integral, rather than only Euclidean sign calibration.
3. Constrained metric/clock-field measure equivalence.
4. Matching Einstein and independent curvature-squared coefficients across regulator prescriptions. Separately setting bare terms to zero is not a matched comparison of one physical theory.
5. Interacting scalar mixing, 1/N corrections, and eventual nonlinear gravity/universality checks for the extended model.

## Latest bounded batch (2026-10-30 US Central)
See experiments/source_dispersion_audit_2026_10_30/REPORT.md, code and saved results.
- Clarified the observable: W=(q-2)^2/(16A) is the scalar component of the identical-constituent connected single-source quadratic coefficient. Spin-2 exchange is separate; distinct-source cross terms carry a factor of two. No S-matrix normalization or isolated scalar ghost conclusion is implied.
- Audited transformed pole residues and cut weights. Nine distinct exact source-tail certificates follow from the saved coefficient certificates. Raw/reference source densities grow as s and require two subtractions after retaining/removing the known quadratic polynomial; four non-reference matched source densities fall as s^-3 and admit an absolute formula with constant 1/(16a2).
- Fifteen raw/reference second-subtraction reconstructions, twelve matched absolute reconstructions and eight source moments agree numerically. Maximum normalized reconstruction errors are below 8.61e-14 and 1.69e-14. Nine Cauchy-derivative checks agree below 4.70e-16. No quadrature warnings occurred. Tail certification remains separate from numerical residue/quadrature accuracy.
- The raw polynomial subtraction equals the scalar projection of the previously conditional auxiliary/composite source contact: a_infinity=-3p/2 implies F_infinity*(q-2)^2/16=-(q-2)^2/(24p). This closes its algebraic normalization dictionary under the same B,p/source assumptions, not microscopic measure equivalence.
- A new exact 4x4 block-Gaussian example shows that allocating a kernel modification to H or to B in K=H-B/4 leaves the metric response identical while changing the composite covariance/contact. Eight block inversions and twelve source forms have exactly zero residual. This illustrates missing source-functional information; neither allocation was adopted as a physical composite model.
- Primary-literature source conventions were checked directly against published JHEP11(2018)153 Appendix C. No matter coupling, physical regulator, contour, source prescription or model assumption changed.

## Next bounded work batch
Inspect the published composite-gravity source and operator-renormalization prescriptions directly, focusing on clock/metric variations, normal ordering, local counterterms and the discarded higher-derivative terms. Determine which finite-momentum source matching is specified and which requires additional physical input. Do not select a new composite measure or counterterm allocation merely to transfer auxiliary poles. Preserve the distinction between the scalar connected-source component and full tensor exchange. Constrained measure, Lorentzian contour, interactions and 1/N remain distinct open problems.

## Working policy
Keep established literature, reproduced calculations, assumptions, and speculation distinct. Agreement between AI assistants is not independent verification. Preserve unsuccessful tests and corrections. Do not infer ghost freedom from one residue, or unification from consistency checks. Do not add a material/spin interaction solely to fit a desired conclusion.

## Automated continuation
Hourly task: Advance gravity research. Resume from this checkpoint and latest experiment reports. Avoid routine user notifications; request input only for consequential model choices, genuine access blockers, or major results requiring review. Work in bounded runs, not an implied continuously running process. No paid compute or unrelated projects.
