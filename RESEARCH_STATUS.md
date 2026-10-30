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
See experiments/auxiliary_source_matching_2026_10_30/REPORT.md, sources.json and check.py/results.json.
- Inspected published source/Jacobian conditions in Criado-Perez-Victoria Section 2 and Carone's auxiliary/constrained distinction (footnote 4, Section IV). Henneaux's auxiliary-preserving gauge result was accessible only as an abstract; no model-specific quantum measure was inferred.
- Derived the finite-regulator triangular shift g=G(phi)+h (unit change-of-variable Jacobian), followed by the distinct induced weight I_phi[0] and source factor I_phi[J]/I_phi[0]. Matching matter sources requires the induced reduced measure; matching independent-g insertions additionally requires their source cumulants/conditional-mean map.
- Computed the actual algebraic metric Hessian on a definite D=4 branch: nine eigenvalues V/2 and one -V/2; det H=-(V/2)^10. The constraint Jacobian has det J_T=-V^10 and is not the canonical Poisson matrix.
- Nine new high-precision directional checks and exact fraction determinants passed. The leading stationary-phase magnitude scales as |V|^-5, so it is not generically a field-independent normalization. No quantum integration cycle has been defined.
- Qualified the imported classical bridge: a raw real positive-flavor gradient Gram matrix is definite when nonsingular and cannot itself yield a Lorentzian G by scalar multiplication. This does not rule out a renormalized quantum composite, complex continuation, clock-containing CEV, or the published perturbative construction.
- All previous auxiliary pole/dispersion results remain conditional; no spectrum transfer or physical ghost diagnosis follows from the algebraic trace direction.

## Next bounded work batch
Inspect the renormalized composite, self-consistency, and gauge/measure definition in background-independent clock-free composite gravity (arXiv:2011.06541 and published version) before trying to transfer the auxiliary spectrum. Specify which limits/observables are actually controlled. A source-matched large-N comparison and Lorentzian cycle remain missing; preserve all earlier data and avoid repeated passing scans.

## Working policy
Keep established literature, reproduced calculations, assumptions, and speculation distinct. Agreement between AI assistants is not independent verification. Preserve unsuccessful tests and corrections. Do not infer ghost freedom from one residue, or unification from consistency checks. Do not add a material/spin interaction solely to fit a desired conclusion.

## Automated continuation
Hourly task: Advance gravity research. Resume from this checkpoint and latest experiment reports. Avoid routine user notifications; request input only for consequential model choices, genuine access blockers, or major results requiring review. Work in bounded runs, not an implied continuously running process. No paid compute or unrelated projects.
