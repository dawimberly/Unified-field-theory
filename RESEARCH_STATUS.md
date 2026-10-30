# Research checkpoint

## Goal
Derive a healthy universal gravitational response from quantum dynamics, recovering Einstein gravity at low energy. No solution or discovery is established.

## Current calculational scope
Free, minimally coupled, leading-N auxiliary-metric determinant; stationary cosmological contacts; no bare curvature terms. This is not yet demonstrated to be equivalent to the constrained CEV path integral with its clock fields and measure. Constituent M=1. Four-species PV spectra z=[1,L,2L,5L+1], L=3,5,10,30,100, with first three signed moments vanishing.

## Completed checks
- Full momentum scalar kernel in finite D has an exact zero at Q_E^2=4M^2/(D-2), not a derivative-expansion artifact.
- Identical on-shell constituent scalar source traces vanish at that same zero, cancelling it in that exchange channel. Unequal-mass probes need not cancel it; vacuum stability does not follow.
- Three raw PV spectra have subthreshold scalar poles with positive connected Euclidean source residues. Exact coefficient/sign arguments now give complete extra-zero counts within |q|<4: one for raw L=3,5,10, none for raw L=30,100; one for all five explicitly Lref=10 matched spectra. This is not global first-sheet pole control.
- Full 4D 10x10 metric Hessian, with bubble and stationary contacts, passes longitudinal Ward checks in 40 cases (maximum absolute residual about 1.9e-13). Scalar projection K0=2A/3 agrees.
- Sixteen finite-D longitudinal/trace tensor checks agree with the corresponding scalar kernel.
- Physical-only continuum scalar weights are positive analytically below the first PV threshold, confirmed in 750 samples. The spin-2 inverse cut sign is positive there too.
- Above regulator thresholds each tested PV spectrum develops negative scalar spectral weight. PV is not thereby a physical unitary UV completion.

- Fixed four-species no-bare-curvature spin-2 determinant: the Genocchi-Hermite representation makes K2/q a positive Stieltjes transform and excludes additional first-sheet spin-2 poles. This does not extend automatically to modified kernels or constrained CEV.

## Open problems
1. Global scalar first-sheet zero control outside |q|<4. The scalar first-threshold disk now has an analytic zero classification; prior bounded contours and finite-radius inverse dispersion checks remain. The fixed-prescription spin-2 coefficient has an analytic positive-measure zero-exclusion result.
2. Lorentzian contour from the microscopic metric path integral, rather than only Euclidean sign calibration.
3. Constrained metric/clock-field measure equivalence.
4. Matching Einstein and independent curvature-squared coefficients across regulator prescriptions. Separately setting bare terms to zero is not a matched comparison of one physical theory.
5. Interacting scalar mixing, 1/N corrections, and eventual nonlinear gravity/universality checks for the extended model.

## Latest bounded batch (2026-10-30 US Central)
See experiments/scalar_conformal_count_2026_10_30/REPORT.md, check.py, verify.py and saved results.
- Derived a whole-first-sheet conformal representation q=16w/(1-w)^2 retaining every regulator cut on |w|=1. Counted only |w|<4/5, not the full unit disk.
- Exact rational Chebyshev coefficients, bounded logarithmic constants, an explicit analytic tail and directed-Decimal arc enclosures give validated additional-zero counts. Raw L=3,5,10,30,100: 1,1,1,0,0. Explicitly Lref=10 matched: 3,3,1,2,2.
- Newly isolated matched scalar poles: L=3 q=-17.4662104652 +/-12.7654098460i; L=5 q=-25.3822329364 +/-30.4480072636i; L=30 q=+34.2768433995; L=100 q=+27.5993608103. The earlier subthreshold poles remain.
- Four new representative isolating counts and eight known subthreshold isolating counts each equal one. Conjugation supplies the complex partners. Independent principal-log parameter quadrature at three node counts gives residuals below 2.81e-14 for the new poles.
- Focused N=160/80-digit refinement preserves matched counts for L=3,100. Twelve conformal/compact overlap checks and exact-Fraction arithmetic enclosure audits pass.
- An initial rounding-context issue and conservative interval-wrapping attempts were discarded before accepted results; corrections are documented. Final counts use explicit tail margins and exact rational winding, not floating-point phase agreement alone.
- All new poles are outside |q|<4 (|q|/4=5.41--9.91). They concern the full explicitly modified kernels, not controlled low-energy predictions or an adopted physical UV action. No global first-sheet count, constrained equivalence or physical health is established.

## Next bounded work batch
Audit large-q and cut asymptotics in the conformal representation before extending beyond |w|<4/5. Seek explicit tail and boundary control for the unresolved collar 4/5<=|w|<1; do not infer a global count from convergence of bounded counts. Keep the common regulated composite-source/measure dictionary, Lorentzian contour, interactions and 1/N as distinct open problems. No new physical action or regulator is chosen.

## Working policy
Keep established literature, reproduced calculations, assumptions, and speculation distinct. Agreement between AI assistants is not independent verification. Preserve unsuccessful tests and corrections. Do not infer ghost freedom from one residue, or unification from consistency checks. Do not add a material/spin interaction solely to fit a desired conclusion.

## Automated continuation
Hourly task: Advance gravity research. Resume from this checkpoint and latest experiment reports. Avoid routine user notifications; request input only for consequential model choices, genuine access blockers, or major results requiring review. Work in bounded runs, not an implied continuously running process. No paid compute or unrelated projects.
