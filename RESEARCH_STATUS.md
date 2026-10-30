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
See experiments/complex_scalar_first_sheet_2026_10_30/REPORT.md and its executable check.py/results.json.
- Scalar first-sheet principal-branch formula agrees with independent 1024-node parameter integration at 25 off-axis points (max relative error 4.58e-12).
- Ninety adaptively resolved upper-half-plane contour evaluations have zero winding for A/q. R=4L,20L,100L; imaginary floors 0.05,0.005; three initial resolutions. Max final phase step 0.233 radians. This excludes no untested near-axis region, real-axis root, second-sheet pole, or remote complex root; numerical contour sampling is not a certified proof.
- Fifteen twice-subtracted dispersion reconstructions of A agree to max relative error 3.19e-12. This is the kernel, not 1/A: no isolated-pole terms enter this check. Subtraction data are supplied from the exact kernel.
- Existing subthreshold scalar poles and all earlier limitations remain.

## Next bounded work batch
Count known real subthreshold roots by cut-avoiding contours as a nonzero-count calibration. Reconstruct 1/A including massless and isolated poles, regulator cuts, and large-circle/asymptotic control. Improve near-axis error control; spin-2 complex poles remain untested. Keep regulator matching and constrained measure separate. See latest report for bounds and reproduction.

## Working policy
Keep established literature, reproduced calculations, assumptions, and speculation distinct. Agreement between AI assistants is not independent verification. Preserve unsuccessful tests and corrections. Do not infer ghost freedom from one residue, or unification from consistency checks. Do not add a material/spin interaction solely to fit a desired conclusion.

## Automated continuation
Hourly task: Advance gravity research. Resume from this checkpoint and latest experiment reports. Avoid routine user notifications; request input only for consequential model choices, genuine access blockers, or major results requiring review. Work in bounded runs, not an implied continuously running process. No paid compute or unrelated projects.
