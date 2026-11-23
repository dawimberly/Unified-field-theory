# Full scalar-response and regulator audit

## Outcome

The finite-D minimal auxiliary-metric determinant has an exact spacelike zero at Q^2=4 M^2/(D-2), not merely a zero of a derivative expansion. For five PV spectra, that same spacelike zero is absent on the sampled spacelike intervals. The subthreshold timelike pole is also not regulator-independent. These comparisons change the UV prescription without matching the low-energy curvature coefficients, so they do not establish scheme dependence of a fixed physical theory.

## Equations and assumptions

Set M=1 and q=Q_E^2. In dimensional regularization: U=1-(D-2)q/4; G=(4pi)^(-D/2) Gamma(1-D/2); P(q)=(4pi)^(-D/2) Gamma(2-D/2) integral_0^1 [1+x(1-x)q]^(D/2-2) dx; A(q)=-U(PU+(D-2)G/2)/4.

Thus U=0 is an exact zero and A′(q*)=(D-2)^2 G/32>0 for 4<D<6. This is a conformal inverse-kernel zero, not yet a verified physical scattering pole. A positive slope alone does not establish a positive-norm particle: source overlaps, the full tensor constraints, Lorentzian continuation and sign normalization must be established.

For PV, each species has its own U_i=z_i-q/2. The summed kernel is C/4 sum_i c_i[(z_i-q/2)^2 integral log(z_i+x(1-x)q) dx -(z_i-q/2)z_i log(z_i)], C=1/(4pi)^2. Consequently the shared U factor of the single-mass dimensional-regulator expression is lost. The first three PV moments vanish.

The exact low-q coefficients are a1=-G_PV/8 and a2=C sum c_i log(z_i)/16. The nonzero root of a1*q+a2*q^2 is compared with the full nonlocal kernel. No expansion is used to locate the full roots.

## PV results

| Regulator L | Quadratic root q | Full subthreshold roots q | No sampled spacelike zero |
|---|---:|---|---|
| 3 | -3.743659 | -2.09002986 | True |
| 5 | -5.354691 | -2.77851951 | True |
| 10 | -8.764793 | -3.74730196 | True |
| 30 | -19.582842 | none found | True |
| 100 | -49.137673 | none found | True |

## What this establishes

Seven finite-D root locations and slopes agree with the analytic factorization. Five PV spectra satisfy moment cancellation and positive induced Einstein coefficient. Each PV spectrum is scanned at 400 spacelike and 500 subthreshold timelike points, with bracketing root refinement. This finite sampling does not prove global root absence. Signed PV weights are not a unitary UV completion.

A pole near the pair-production threshold is sensitive to higher derivative terms; extrapolating an R+R^2 truncation there is not reliable. Conversely, the exact finite-D zero cannot be dismissed solely as a truncation artifact.

## Biggest remaining weakness and next calculation

The conformal kernel has not been converted into a conserved-source, gauge-consistent physical scattering amplitude. Therefore none of these numbers is a final ghost verdict. Next derive the full tensor projector decomposition and the scalar source overlap using the same contact prescription, then analytically continue below threshold and compare the pole residue with a canonical positive-energy scalar exchange. Only after that should interactions or a claimed healthy scalar be added.

For regulator comparisons of the same physical theory, match the Einstein term and independent curvature-squared coefficients; imposing zero bare terms separately in different prescriptions does not perform this matching.

## References and reproduction

Carone, Claringbold and Vaman, Composite graviton self-interactions in a model of emergent gravity, arXiv:1710.09367, Phys. Rev. D 98, 024041 (2018): https://arxiv.org/abs/1710.09367 . Carone, Composite Gravitons from Metric-Independent Quantum Field Theories, arXiv:1911.13166: https://arxiv.org/abs/1911.13166 . These references motivate the program; the conditional kernels and scans above are our calculation, not claims attributed to those papers.

Run python scalar_full_response_scan.py with NumPy and SciPy. Dimensionless masses are in constituent-mass units. No claim of novelty or experimental confirmation.