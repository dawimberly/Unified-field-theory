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
See experiments/threshold_independent_2026_10_30/REPORT.md and check.py/results.json.
- Read the newly authorized provisional threshold archive (Library libfile_fa52b432f24c8191ac948aabcca17d35) and verified it independently with u=sqrt(4/s-1), not its root solver.
- The physical lower q-lip has u=-i beta and Im J=-pi beta. The below-threshold adjacent sheet changes J by -2pi u, confirming the archived minus sign.
- Decimal threshold root: L_c=14.5572064556843705..., agreeing with the provisional value within 1.73e-13. Eight offsets confirm the local sheet-I to sheet-II real pole trajectory and vanishing square-root residue. No global uniqueness or complete second-sheet search is established.
- Independent principal-log quadrature initially missed narrow complex-log layers; retained failure data and corrected root-neighborhood splits reduce the maximum comparison error to 1.93e-15. All corrected assertions passed.
- Explicitly matched the Einstein slope and both tensor q^2 coefficients to reference L=10, recording both A and K2 polynomial changes. All five tested original spectra then have first-sheet subthreshold scalar poles, including L=30,100. The unmatched crossing is not a fixed-UV-theory stability prediction; higher-order coefficients differ.
- No bare model, physical regulator, contour or constrained-equivalence assumption was adopted. Prior source/measure and fixed-prescription pole results remain in their existing folders.

## Next bounded work batch
Audit zeros of the explicitly coefficient-matched spin-2 kernel under its changed polynomial terms, without transferring the no-bare-curvature positive-measure theorem to it. Keep such resummed kernel behavior distinct from controlled low-energy EFT predictions. Continue pending background-independent composite/source-matching literature (arXiv:2011.06541) before any auxiliary-to-constrained spectrum transfer. Preserve all prior data.

## Working policy
Keep established literature, reproduced calculations, assumptions, and speculation distinct. Agreement between AI assistants is not independent verification. Preserve unsuccessful tests and corrections. Do not infer ghost freedom from one residue, or unification from consistency checks. Do not add a material/spin interaction solely to fit a desired conclusion.

## Automated continuation
Hourly task: Advance gravity research. Resume from this checkpoint and latest experiment reports. Avoid routine user notifications; request input only for consequential model choices, genuine access blockers, or major results requiring review. Work in bounded runs, not an implied continuously running process. No paid compute or unrelated projects.
