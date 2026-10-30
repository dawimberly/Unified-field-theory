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
See experiments/scalar_threshold_disk_2026_10_30/REPORT.md, check.py and results.json.
- Proved positive scalar inverse moments and coefficient-ratio bound d_(k+1)<d_k/4 for the whole PV family L>=3, raw or explicitly matched to Lref=10. The first-ratio certificates use exact rational Riemann sums; the higher ratios use a finite exact prefix and an infinite induction argument.
- A deflated convergent-series argument gives at most one additional scalar zero in the entire open disk |q|<4; if present it is simple and on the negative real q axis.
- Exact rational logarithm/arctangent enclosures certify all ten tested threshold signs: raw extra-zero counts 1,1,1,0,0 at L=3,5,10,30,100; all five matched counts equal one. The known massless zero is separate.
- Eight independent 512-term positive-series brackets confirm archived real pole locations within +/-1e-8; these location checks use long double and are not interval-certified. Maximum saved tail bound is 5.09e-14.
- Eighty inverse-moment identities were checked exactly. All new assertions passed. This improves the prior disk-|q|<=2 bound without repeating a broad contour or Ward scan.
- The signed scalar continuum outside regulator thresholds remains negative in places; positive local coefficients do not give a global Stieltjes representation. No global first-sheet, constrained-equivalence or physical-health conclusion follows.

## Next bounded work batch
Seek scalar zero control outside |q|<4 using a conformal map or cut-subtracted representation retaining every regulator cut. Do not assume finite negative-index methods apply or transfer the disk theorem to a global stability claim. Keep the common regulated composite-source/measure dictionary, Lorentzian contour, interactions and 1/N as distinct open problems.

## Working policy
Keep established literature, reproduced calculations, assumptions, and speculation distinct. Agreement between AI assistants is not independent verification. Preserve unsuccessful tests and corrections. Do not infer ghost freedom from one residue, or unification from consistency checks. Do not add a material/spin interaction solely to fit a desired conclusion.

## Automated continuation
Hourly task: Advance gravity research. Resume from this checkpoint and latest experiment reports. Avoid routine user notifications; request input only for consequential model choices, genuine access blockers, or major results requiring review. Work in bounded runs, not an implied continuously running process. No paid compute or unrelated projects.
