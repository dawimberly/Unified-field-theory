# Research checkpoint

## Goal
Derive a healthy universal gravitational response from quantum dynamics, recovering Einstein gravity at low energy. No solution or discovery is established.

## Current calculational scope
Free, minimally coupled, leading-N auxiliary-metric determinant; stationary cosmological contacts; no bare curvature terms. This is not yet demonstrated to be equivalent to the constrained CEV path integral with its clock fields and measure. Constituent M=1. Four-species PV spectra z=[1,L,2L,5L+1], L=3,5,10,30,100, with first three signed moments vanishing.

## Completed checks
- Full momentum scalar kernel in finite D has an exact zero at Q_E^2=4M^2/(D-2), not a derivative-expansion artifact.
- Identical on-shell constituent scalar source traces vanish at that same zero, cancelling it in that exchange channel. Unequal-mass probes need not cancel it; vacuum stability does not follow.
- Three PV spectra have subthreshold scalar poles with positive connected Euclidean source residues. Two have none in the sampled subthreshold interval. Finite scans do not prove global pole absence.
- Full 4D 10x10 metric Hessian, with bubble and stationary contacts, passes longitudinal Ward checks in 40 cases (maximum absolute residual about 1.9e-13). Scalar projection K0=2A/3 agrees.
- Sixteen finite-D longitudinal/trace tensor checks agree with the corresponding scalar kernel.
- Physical-only continuum scalar weights are positive analytically below the first PV threshold, confirmed in 750 samples. The spin-2 inverse cut sign is positive there too.
- Above regulator thresholds each tested PV spectrum develops negative scalar spectral weight. PV is not thereby a physical unitary UV completion.

- Fixed four-species no-bare-curvature spin-2 determinant: the Genocchi-Hermite representation makes K2/q a positive Stieltjes transform and excludes additional first-sheet spin-2 poles. This does not extend automatically to modified kernels or constrained CEV.

## Open problems
1. Global scalar first-sheet zero control. Bounded scalar contours and finite-radius inverse dispersion checks exist; the fixed-prescription spin-2 coefficient has an analytic positive-measure zero-exclusion result.
2. Lorentzian contour from the microscopic metric path integral, rather than only Euclidean sign calibration.
3. Constrained metric/clock-field measure equivalence.
4. Matching Einstein and independent curvature-squared coefficients across regulator prescriptions. Separately setting bare terms to zero is not a matched comparison of one physical theory.
5. Interacting scalar mixing, 1/N corrections, and eventual nonlinear gravity/universality checks for the extended model.

## Latest bounded batch (2026-10-30 US Central)
See experiments/stress_contact_dictionary_2026_10_30/ and experiments/low_energy_remainders_2026_10_30/ (each has REPORT.md, check.py, results.json).
- Derived the fixed-field seagull, connected stress bubble and stationary cosmological dictionary for the specified free Euclidean PV determinant: seagull=pP, cosmological Hessian=-pP/2, Pi(0)=2pP, full K=pP/2-Pi/4. Metric-dependent operator subtractions still require explicit local source terms.
- New independent radial propagator/Wick integrals for all five spectra agree with closed logarithmic tadpoles and the zero-momentum stress normalization to relative errors below 4.82e-16. No old finite-q Ward scan was repeated.
- Derived geometric absolute Taylor-tail bounds for both explicitly Lref=10 matched channels inside |q|<4, retaining that spectrum's own higher coefficients. Propagated the bound to inverse responses with a nonvanishing-denominator condition.
- All 320 new kernel comparisons and 320 inverse comparisons satisfy their bounds. Forty disk inequalities satisfy the sufficient no-extra-zero condition, including the full disk |q|<=2 for both matched channels and all five spectra. Numerical constants are not outward-rounded interval certificates.
- Uniform quadratic inverse relative-error bounds on |q|<=0.5 are below 0.964% for A and 0.0541% for K2 across these spectra. They bound the specified kernel approximation, not unknown UV, interaction or measure corrections.
- Primary literature contact-convention context checked in the published Bzowski/McFadden/Skenderis JHEP11(2018)153 introduction and Appendix C; no conformal three-point formula was transplanted to massive PV.
- The high-q resummed matched spin-2 poles, scalar threshold continuation and conditional quadratic composite source bridge remain in earlier folders. No contour, physical bare action, composite measure equality or healthy UV completion was adopted.

## Next bounded work batch
Extend scalar first-sheet zero control beyond the present local disk using subtracted analytic representations, not another unmotivated broad scan. Audit whether signed scalar divided differences admit a useful integral sign or variation bound; preserve negative-continuum obstruction to finite-index Nevanlinna methods. Before any composite-spectrum transfer, supply one common regulated generating functional including local operator/source subtractions and gauge treatment. Lorentzian contour, interactions and 1/N remain distinct open problems.

## Working policy
Keep established literature, reproduced calculations, assumptions, and speculation distinct. Agreement between AI assistants is not independent verification. Preserve unsuccessful tests and corrections. Do not infer ghost freedom from one residue, or unification from consistency checks. Do not add a material/spin interaction solely to fit a desired conclusion.

## Automated continuation
Hourly task: Advance gravity research. Resume from this checkpoint and latest experiment reports. Avoid routine user notifications; request input only for consequential model choices, genuine access blockers, or major results requiring review. Work in bounded runs, not an implied continuously running process. No paid compute or unrelated projects.
