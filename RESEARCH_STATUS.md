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
See experiments/spin2_stieltjes_2026_10_30/REPORT.md and check.py/results.json.
- Derived K2 from the saved transverse-traceless full tensor Hessian including stationary contacts; 25 tensor comparisons agree within 9.37e-14 absolute.
- Applying the standard Genocchi-Hermite divided-difference identity gives K2(q)/q=H integral dx d lambda [1/Z]/[q+Z/t], H=-C/(4w0)>0, Z=sum lambda_i z_i, t=x(1-x). This positive Stieltjes representation analytically excludes all additional first-sheet spin-2 zeros in the fixed four-species no-bare-curvature model.
- The massless spin-2 residue is 24/G>0. Its continuum sign remains positive above PV thresholds as well. This does not repair the scalar negative spectral regions or establish physical UV unitarity.
- 20 off-axis quadrature checks agree within 2.67e-14 relative. Ten simplex normalization checks converge to at worst 1.15e-6; 30 bounded contour windings are zero. Fifteen kernel dispersion checks agree within 3.81e-11 relative.
- This is a conditional analytic result ready for independent review, not external validation or a novelty claim. Adding bare curvature terms can invalidate the representation.
- Prior scalar contour and inverse-dispersion audits remain under their existing experiment directories. Scalar global pole control, physical stability, regulator matching and constrained equivalence remain unproved.

## Next bounded work batch
Apply divided-difference/measure methods to the scalar kernel and attempt analytic control of its remaining zeros. Then match the induced Einstein and two independent curvature-squared coefficients before physical regulator comparisons. Keep bare-curvature changes distinct from the proven fixed-prescription spin-2 statement. Preserve all previous experiments.

## Working policy
Keep established literature, reproduced calculations, assumptions, and speculation distinct. Agreement between AI assistants is not independent verification. Preserve unsuccessful tests and corrections. Do not infer ghost freedom from one residue, or unification from consistency checks. Do not add a material/spin interaction solely to fit a desired conclusion.

## Automated continuation
Hourly task: Advance gravity research. Resume from this checkpoint and latest experiment reports. Avoid routine user notifications; request input only for consequential model choices, genuine access blockers, or major results requiring review. Work in bounded runs, not an implied continuously running process. No paid compute or unrelated projects.
