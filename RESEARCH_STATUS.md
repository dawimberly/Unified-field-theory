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
See experiments/matched_spin2_2026_10_30/REPORT.md and check.py/results.json.
- Audited the explicitly Lref=10 coefficient-matched spin-2 kernels; no physical bare action was adopted.
- Additional principal-sheet zeros: L=3 q=60.3330308219; L=5 q=82.4469390832; L=30 q=8.72600156453 +/-78.0620223980i; L=100 q=10.0532001673 +/-66.1017934580i. L=10 retains its original no-extra-pole theorem.
- Independent parameter integrals at 256/512/1024 nodes agree at the roots (maximum absolute matched-kernel residual 4.56e-14). Thirty-eight adaptive isolating/upper-rectangle contour checks pass.
- Derived positive-measure restrictions for F=K2_matched/q=S+a+bq: b<=0 excludes off-axis zeros and gives a unique positive real zero here; b>0 permits at most one simple conjugate pair. A Cauchy-Schwarz real-axis lower bound excludes real zeros for tested L=30,100, with floating-point coefficient signs rather than interval certification.
- Extra zeros have |q|/4=15.08--20.62 and quadratic/Einstein ratios 1.38--1.88. These are full modified-kernel poles, not controlled low-energy expansion predictions or constrained-CEV stability conclusions.
- Earlier independently verified scalar threshold/sheet continuation and matched scalar poles remain in experiments/threshold_independent_2026_10_30/. The original no-bare spin-2 theorem remains valid under its original hypotheses.
- All new assertions passed. No previous full scans were repeated. No UV completion or healthy theory is established.

## Next bounded work batch
Inspect the pending published background-independent composite/source definition (arXiv:2011.06541) before transferring any auxiliary spectra to constrained CEV. Then formulate a controlled low-energy matched response comparison with explicit truncation/error assumptions, preserving the distinction between that expansion and fully resummed polynomial-modified kernels. Global scalar zero control, measure, Lorentzian contour and interactions remain open.

## Working policy
Keep established literature, reproduced calculations, assumptions, and speculation distinct. Agreement between AI assistants is not independent verification. Preserve unsuccessful tests and corrections. Do not infer ghost freedom from one residue, or unification from consistency checks. Do not add a material/spin interaction solely to fit a desired conclusion.

## Automated continuation
Hourly task: Advance gravity research. Resume from this checkpoint and latest experiment reports. Avoid routine user notifications; request input only for consequential model choices, genuine access blockers, or major results requiring review. Work in bounded runs, not an implied continuously running process. No paid compute or unrelated projects.
