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

## Open problems
1. First-sheet complex poles and dispersion-relation reconstruction.
2. Lorentzian contour from the microscopic metric path integral, rather than only Euclidean sign calibration.
3. Constrained metric/clock-field measure equivalence.
4. Matching Einstein and independent curvature-squared coefficients across regulator prescriptions. Separately setting bare terms to zero is not a matched comparison of one physical theory.
5. Interacting scalar mixing, 1/N corrections, and eventual nonlinear gravity/universality checks for the extended model.

## Latest bounded batch (2026-10-29 US Central)
See experiments/inverse_scalar_dispersion_2026_10_30/REPORT.md and check.py/results.json.
- Fifteen cut-avoiding subthreshold contours reproduce A/q zero counts [1,1,1,0,0] for L=[3,5,10,30,100], stable across three resolutions. Maximum phase step 0.04184 radians.
- Isolated roots and residues agree with the earlier scalar calculation. The massless residue is included as tensor-structure bookkeeping, not a standalone physical scalar interpretation.
- Forty-five once-subtracted reconstructions of 1/A include massless and isolated poles, all regulator cuts and an explicit finite outer circle. Maximum relative closure error is 2.00e-13.
- Omitting the outer circle leaves maximum relative errors 2.164e-4,1.508e-6,5.486e-8 at R/L=100,1000,5000. The outer term quantitatively closes these residuals.
- A tends to (3C/8) sum c_i z_i^2 log(z_i); 1/A therefore has a constant asymptote. An unsubtracted inverse spectral representation cannot silently drop it.
- Preceding bounded scalar complex-contour/branch audit remains at experiments/complex_scalar_first_sheet_2026_10_30/. Global pole absence, physical stability, regulator matching and constrained equivalence remain unproved.

## Next bounded work batch
Derive and audit the spin-2 first-sheet kernel from the saved full tensor Hessian, then apply complex-root/dispersion checks with explicit branches and boundaries. Scalar near-axis strips and rigorous global zero/asymptotic control remain open. Match Einstein and independent curvature-squared terms before interpreting regulator differences physically. Preserve all prior experiments.

## Working policy
Keep established literature, reproduced calculations, assumptions, and speculation distinct. Agreement between AI assistants is not independent verification. Preserve unsuccessful tests and corrections. Do not infer ghost freedom from one residue, or unification from consistency checks. Do not add a material/spin interaction solely to fit a desired conclusion.

## Automated continuation
Hourly task: Advance gravity research. Resume from this checkpoint and latest experiment reports. Avoid routine user notifications; request input only for consequential model choices, genuine access blockers, or major results requiring review. Work in bounded runs, not an implied continuously running process. No paid compute or unrelated projects.
