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
1. Audit full source-projected dispersion, contact terms and required subtractions using the now complete scalar inventories and infinite-boundary coefficient-inverse reconstruction. Keep the different raw/reference and non-reference matched asymptotics explicit. Other spectra, adjacent sheets and physical interpretation remain separate. The fixed-prescription spin-2 coefficient has an analytic positive-measure zero-exclusion result.
2. Lorentzian contour from the microscopic metric path integral, rather than only Euclidean sign calibration.
3. Constrained metric/clock-field measure equivalence.
4. Matching Einstein and independent curvature-squared coefficients across regulator prescriptions. Separately setting bare terms to zero is not a matched comparison of one physical theory.
5. Interacting scalar mixing, 1/N corrections, and eventual nonlinear gravity/universality checks for the extended model.

## Latest bounded batch (2026-10-30 US Central)
See experiments/infinite_scalar_dispersion_2026_10_30/REPORT.md and executed code/results.
- The preceding exact cut-boundary argument supplies complete additional scalar counts: raw 1,1,1,0,0; explicitly matched 3,3,1,4,2 for L=3,5,10,30,100. This batch uses those inventories rather than inferring completeness from reconstruction residuals.
- Derived infinite-cut absolute and once-subtracted coefficient-inverse dispersion relations. Raw/reference F=1/A retains F_infinity=1/a_infinity; four non-reference matched inverses vanish at infinity. Every isolated pole, complex conjugate and regulator cut is retained.
- Ten exact rational UV certificates bound the inverse density by K/s for raw/reference cases and K/s^5 for non-reference matched cases. Explicit tail bounds cover reconstructions and all convergent matched moments. They certify continuum truncation only, not every numerical quadrature/residue digit.
- Thirty subtracted and twenty-seven absolute numerical reconstructions close to normalized errors below 6.02e-15 and 2.44e-14 respectively. Sixteen non-reference matched moment checks agree with the derived coefficients; the largest absolute moment discrepancy is about 2.91e-8 and the largest component-normalized discrepancy is below 4.29e-16. No quadrature warnings occurred.
- The four finite inverse moments are M0=0, M1=1/a2, M2=-a1/a2^2, M3=a1^2/a2^3-a_infinity/a2^2, with a1,a2 the added UV polynomial coefficients. The fourth continuum moment diverges logarithmically. Even the zeroth raw/reference density moment diverges; no matched superconvergence rule is transferred to L=10.
- Derived the algebraic source-projection implications and identified the different source subtractions requiring the next audit. No physical source prescription, constrained measure, Lorentzian contour or action was changed, and no health or UV-completion claim follows.

## Next bounded work batch
Independently audit the full constituent source exchange W=(q-2)^2/(16A): transformed pole residues and cut weights, the non-reference matched constant contact 1/(16a2), and two subtractions after removing the known quadratic polynomial in raw/reference cases. Check the relationship to the previously saved stress-bubble/composite-source contact dictionary without assuming microscopic equivalence. Preserve rigorous tail control separately from numerical error estimates. Constrained measure, the common source prescription, Lorentzian contour, interactions and 1/N remain distinct open problems.

## Working policy
Keep established literature, reproduced calculations, assumptions, and speculation distinct. Agreement between AI assistants is not independent verification. Preserve unsuccessful tests and corrections. Do not infer ghost freedom from one residue, or unification from consistency checks. Do not add a material/spin interaction solely to fit a desired conclusion.

## Automated continuation
Hourly task: Advance gravity research. Resume from this checkpoint and latest experiment reports. Avoid routine user notifications; request input only for consequential model choices, genuine access blockers, or major results requiring review. Work in bounded runs, not an implied continuously running process. No paid compute or unrelated projects.
