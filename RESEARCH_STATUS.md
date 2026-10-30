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

## Latest bounded batch (2026-10-29 US Central)
See experiments/scalar_identity_matching_2026_10_30/REPORT.md and check.py/results.json.
- Read and imported three current conditional notes under bots/gpt/research_notes/2026-10-29/: literature elimination, matched constraint identity, and stochastic regulator screening. They do not establish quantum equivalence or adopt a stochastic microscopic theory.
- Applying the established divided-difference identity yields an exact scalar simplex formula, but its integrand has no fixed sign. The spin-2 positive-measure zero theorem cannot be copied to the scalar channel. Global scalar zero control is still open.
- Six scalar simplex checks converge to maximum relative error 3.806e-8. A failed 48-node resolution check (1.436e-5 error) is retained; the wide L=100 cases were refined to 96 nodes.
- Computed a formal reference-L=10 map matching the Einstein slope and both independent tensor q^2 coefficients. This is bookkeeping, not an adopted bare action or physical regulator.
- Thirty exact-kernel matched differences agree with independently derived Taylor differences through q^8 within 5.81e-17 absolute. First unmatched cubic coefficients are a3=(17C/1120)sum c/z and k3=(C/3360)sum c/z; Delta a3=51 Delta k3 in this family. Higher-order mismatch is not a fixed-theory prediction.
- The prior spin-2 analytic proof and earlier scalar contour/dispersion data remain in their existing experiment folders. Matching polynomials can invalidate the fixed-prescription spin-2 theorem and were not used to claim a pole cure.

## Next bounded work batch
Use established signed-measure/generalized-Nevanlinna results only where their hypotheses can be verified for the scalar kernel; otherwise retain the unresolved global-zero question. Prioritize the imported literature's matched-observable and reduced-measure construction before transferring auxiliary poles to constrained CEV. Inspect the stochastic comparison's source/state-space assumptions without adopting them. Preserve all earlier experiments; do not repeat passing scans without a new reason.

## Working policy
Keep established literature, reproduced calculations, assumptions, and speculation distinct. Agreement between AI assistants is not independent verification. Preserve unsuccessful tests and corrections. Do not infer ghost freedom from one residue, or unification from consistency checks. Do not add a material/spin interaction solely to fit a desired conclusion.

## Automated continuation
Hourly task: Advance gravity research. Resume from this checkpoint and latest experiment reports. Avoid routine user notifications; request input only for consequential model choices, genuine access blockers, or major results requiring review. Work in bounded runs, not an implied continuously running process. No paid compute or unrelated projects.
