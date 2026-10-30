# Inverse scalar dispersion with isolated poles
Research batch: 2026-10-29 US Central; directory date is UTC.
Parent checkpoint: a7ebff203b992df51f84ea60a60efa3101b5b977.

## Scope and outcome
Free minimally coupled leading-N auxiliary-metric determinant with stationary contacts and the same five fixed PV prescriptions. No bare curvature terms are added. No constrained-CEV measure equivalence or Lorentzian stability result is inferred.

The real subthreshold contour counts reproduce one zero of A/q for L=3,5,10, and zero for L=30,100. The inverse kernel F=1/A is reconstructed from its massless pole, those isolated poles, the regulator cuts and an explicitly evaluated finite outer boundary. Forty-five finite-radius Cauchy reconstructions close to a maximum relative error of 2.00e-13.

This checks a dispersion identity of the specified auxiliary model. It does not prove a complete global pole inventory, positivity, a physical UV completion, or a discovery.

## Branch and precision
Continue the principal-log parameter integral from positive Euclidean q. The q-plane cut starts at -4 and the regulator branch points are -4z_i. The off-axis formula is unchanged from the preceding branch audit.

This script computes PV coefficients with long-double barycentric weights:
w_i=1/product_(j!=i)(z_i-z_j), c_i=w_i/w_0.
The vanishing first three moments then follow from the four-point interpolation identity. Extended-precision arithmetic reduces cancellation in the large-q signed sum. The recorded platform has a 63-bit long-double mantissa; this is not arbitrary-precision certification. No extra package or compute service was purchased.

## Subthreshold contour and residue checks
A bracket scan locates roots in (-3.999999,-0.001), then Brent refinement and complex-step slopes determine residues. Three slope steps (1e-4,1e-5,1e-6) check derivative convergence. The largest relative change between the last two slopes is 2.61e-10.

| L | Isolated zero q_p | Residue of 1/A | Contour count of A/q |
|---|---:|---:|---:|
| 3 | -2.0900298566858115 | 738.2971085399244 | 1 |
| 5 | -2.7785195111637395 | 265.78946081336886 | 1 |
| 10 | -3.7473019640017764 | 57.220118706499406 | 1 |
| 30 | none in the tested interval | — | 0 |
| 100 | none in the tested interval | — | 0 |

For every L, a counterclockwise rectangle with real bounds [-3.9999,-0.001] and imaginary bounds [-0.2,0.2] avoids the cut and the massless origin. Counts agree at 512,1024,2048 points per edge (15 contour evaluations). The maximum phase increment is 0.04184 radians. These contours provide a nonzero-count calibration missing from the previous all-zero contour batch. Sampling still supplies no interval-certified bound between nodes or in the tiny excluded threshold/origin gaps.

The massless residue is r_0=1/a_1, a_1=-C sum c_i z_i log(z_i)/8. It is negative in this scalar coefficient. Its physical interpretation requires the Einstein tensor structure; this isolated scalar sign is not a standalone negative-norm-particle diagnosis.

## Once-subtracted finite-radius reconstruction
Use q_0=1 and q=0.5,2,10. Let Delta F=F(q)-F(q_0). For a slit disk of radius R:
Delta F = sum_p r_p[1/(q-q_p)-1/(q_0-q_p)]
+ integral_4^R rho(s)[1/(s+q)-1/(s+q_0)] ds + C_R,
rho(s)=Im[1/A(-s-i0)]/pi,
C_R=(1/(2pi i)) integral_|z|=R F(z)[1/(z-q)-1/(z-q_0)] dz.

The pole sum includes q_p=0. The continuum includes all four species, separated at all pair thresholds, including negative regulator spectral weight. The outer circle is parameterized on the principal sheet from theta=-pi to pi; its endpoint limits lie on opposite lips of the cut. Gauss-Legendre evaluations use 128,256,512 nodes. Contacts and subtraction data are the original model's; nothing is fitted to force agreement.

At R/L=100,1000,5000, the maximum relative errors when C_R is omitted are respectively 2.164e-4,1.508e-6,5.486e-8. Including C_R reduces the maximum relative closure error over all 45 cases to 1.999e-13. The finite outer boundary therefore quantitatively accounts for the apparent truncation residual. This is stronger than just assuming that contribution is negligible.

## Asymptotic caution
Write t=z_i/q. The large-q expansion of the parameter integral is
I_i(q)=log(q)-2+2t(1-log t)+t^2(1+2log t)+O(t^3 log t).
Substitution in the kernel, followed by sum c_i z_i^k=0 for k=0,1,2, gives
A(infinity)=(3C/8) sum c_i z_i^2 log(z_i).
Hence F approaches a nonzero constant for these spectra. An unsubtracted inverse spectral formula that silently omits that constant would be wrong. The difference used here removes the constant, while retaining the exact finite-R boundary explicitly. Numerical large-positive-q ratios approach this limit; these sampled asymptotics are not a uniform complex-plane bound.

## Reproduction and next work
Run python experiments/inverse_scalar_dispersion_2026_10_30/check.py.
Dependencies: NumPy and SciPy. Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 were recorded. Assertions executed successfully. Full pole rows, 15 contour counts, 45 reconstructions, outer-boundary convergence and 15 asymptotic samples are in results.json. Earlier scripts/results were preserved.

Next priority: derive and audit the spin-2 first-sheet kernel and apply comparable complex-root and dispersion checks. Scalar near-axis strips, global asymptotic control and a rigorous complete zero count remain open. Then match Einstein and independent curvature-squared coefficients before physical regulator comparisons. Constrained measure, Lorentzian contour, interactions and 1/N remain distinct questions. No creative model decision is needed for the next technical batch.
