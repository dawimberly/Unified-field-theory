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
1. Global scalar first-sheet zero control in the remaining finite collar near the cuts. The first-threshold disk and |w|<4/5 conformal region are controlled, and explicit UV exterior regions are zero-free. These regions do not yet cover the whole sheet. The fixed-prescription spin-2 coefficient has an analytic positive-measure zero-exclusion result.
2. Lorentzian contour from the microscopic metric path integral, rather than only Euclidean sign calibration.
3. Constrained metric/clock-field measure equivalence.
4. Matching Einstein and independent curvature-squared coefficients across regulator prescriptions. Separately setting bare terms to zero is not a matched comparison of one physical theory.
5. Interacting scalar mixing, 1/N corrections, and eventual nonlinear gravity/universality checks for the extended model.

## Latest bounded batch (2026-10-30 US Central)
See experiments/uv_asymptotic_audit_2026_10_30/REPORT.md, check.py, precision_check.py and saved results/failure diagnostics.
- Derived a convergent uniform exterior expansion retaining every PV species and the explicit principal Log q branch. Confirmed A_infinity=3C T2/8 and K2_infinity=-C T2/8, with exact leading 1/q logarithmic correction coefficients.
- Exact rational remainder/log bounds exclude raw scalar zeros for all |q| beyond R=1024,1664,3264,38656,128256 at L=3,5,10,30,100. Explicitly matched scalar exterior radii are 256,416,3264,2416,8016. These are sufficient, nonoptimal bounds.
- The scalar and raw spin-2 lower-cut imaginary parts are strictly negative for all s>=16 z_max in each tested spectrum. The inverse scalar cut weight is therefore eventually positive again; intervening negative intervals remain and global positivity/physical unitarity do not follow.
- Sixty direct 160-digit comparisons satisfy analytic UV tails plus tracked log-coefficient uncertainty. Forty exterior/cut candidate certificates retain successes and insufficient bounds separately.
- An initial comparison omitted finite log-midpoint uncertainty when the UV tail became smaller; the failure is retained and the error budget corrected before extension.
- Twenty conditioning checks demonstrate catastrophic cancellation in naive compact long-double evaluation at extreme q; use the stable UV series for future exterior work. Fifteen archived moderate-q asymptotic rows differ from high-precision references by at most 2.37e-9 relative, preserving their qualitative interpretation.
- No original files or full passing Ward/dispersion scans were rerun or overwritten. No new physical regulator, action, contour, constrained equivalence or healthy UV completion is established.

## Next bounded work batch
Audit inverse response reconstruction for the explicitly matched scalar kernels with real and complex isolated poles, an explicit finite boundary and stable UV evaluations. Keep pole-inventory completeness separate from a dispersion residual. A whole-first-sheet count still requires the finite collar between the conformal interior and UV exterior, including every cut. Preserve the common regulated composite-source/measure, Lorentzian contour, interactions and 1/N as distinct open problems.

## Working policy
Keep established literature, reproduced calculations, assumptions, and speculation distinct. Agreement between AI assistants is not independent verification. Preserve unsuccessful tests and corrections. Do not infer ghost freedom from one residue, or unification from consistency checks. Do not add a material/spin interaction solely to fit a desired conclusion.

## Automated continuation
Hourly task: Advance gravity research. Resume from this checkpoint and latest experiment reports. Avoid routine user notifications; request input only for consequential model choices, genuine access blockers, or major results requiring review. Work in bounded runs, not an implied continuously running process. No paid compute or unrelated projects.
