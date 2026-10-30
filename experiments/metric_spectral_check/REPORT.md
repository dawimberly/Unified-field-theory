# Continuum spectral check

## Main result

All 750 sampled scalar source spectral weights between the physical pair threshold and first PV pair threshold are positive. The spin-2 inverse-response cut has the positive spectral sign at the same 750 points. Explicit tensor cut integrals reproduce the scalar and spin-2 discontinuities and satisfy the longitudinal Ward check. This closes a real-axis continuum consistency check in the specified free auxiliary-metric model, not a proof of full unitarity.

## Conventions and derivation

Constituent M=1, s is timelike invariant mass squared, q=Q_E^2=-s-i0. A canonical positive Euclidean propagator F(q)=1/(q+m^2) has rho(s)=Im F(-s-i0)/pi positive at its pole. This fixes the spectral sign used here.

The exact Feynman-parameter logarithm has Im integral log(z-s*x*(1-x)-i0) dx = -pi*beta, beta=sqrt(1-4z/s), if s>4z, and zero otherwise. Hence Im A=-C*pi/4 sum c_i(z_i+s/2)^2 beta_i. For the transverse traceless tensor coefficient, Im K2=-C*pi*s^2/240 sum c_i beta_i^5. C=1/(4pi)^2.

Below the first PV threshold only the positive-weight physical scalar contributes. Thus Im A<0 and Im(1/A)>0 analytically throughout this interval; similarly the spin-2 inverse has nonnegative continuum weight wherever finite. This sign conclusion is stronger than the finite scan. The identical-constituent source contributes (2+s)^2/16 times Im(1/A)/pi, which is positive.

The full 10 by 10 tensor cut is independently obtained by integrating the shifted quadratic-vertex polynomial over each open pair interval. Contacts are real and add no discontinuity. Three energies per PV spectrum check the longitudinal rows and the scalar and spin-2 projections. The analytic real logarithm is also checked against independent adaptive integration at five points.

## Results by spectrum

| L | Physical-only s interval | Positive scalar samples | Negative scalar samples above PV threshold | Negative spin-2 cut samples above PV threshold |
|---|---|---:|---:|---:|
| 3 | 4.0001 to 11.988 | 150/150 | 88/300 | 0/300 |
| 5 | 4.0001 to 19.98 | 150/150 | 88/300 | 0/300 |
| 10 | 4.0001 to 39.96 | 150/150 | 88/300 | 0/300 |
| 30 | 4.0001 to 119.88 | 150/150 | 88/300 | 0/300 |
| 100 | 4.0001 to 399.6 | 150/150 | 88/300 | 0/300 |

## What remains open

Positive real-axis continuum weight and isolated-pole residues do not exclude complex poles, poles on other sheets, a bad contour, or failures of a dispersion relation. The scalar massless term cannot be interpreted separately from the Einstein spin-2 tensor structure. Signed PV determinants do not themselves specify a physical unitary UV completion. Above-PV negative weights, if found, concern the regulator-open model; they must not be called low-energy experimental particles.

The next decisive test is a first-sheet complex-pole search and a subtracted dispersion-relation reconstruction below the regulator scale. Agreement there would be stronger than more real-axis parameter scans. The constrained-CEV measure and interactions remain separate gaps.

Reproduce: python metric_spectral_check.py (NumPy and SciPy). No experimental claim or novelty claim. Full spectra, sampled response, integration checks and cut-tensor checks are in results.json.