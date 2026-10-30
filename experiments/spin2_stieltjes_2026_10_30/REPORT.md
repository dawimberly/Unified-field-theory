# Spin-2 positive-measure representation and first-sheet pole exclusion
Research batch: 2026-10-29 US Central. Directory date is UTC.
Baseline: 42f08ab25b3bb1d547bab9c65e06f86b772984b7.

## Result and precise scope
An analytic application of a standard divided-difference identity gives K2(q)/q a positive Stieltjes integral representation. Within the specified free, minimally coupled, leading-N auxiliary-metric determinant, with four distinct ordered positive PV masses, the three cancelling moments and the stationary contacts, K2 has no first-sheet zeros away from the massless origin. This applies to all five spectra already studied, including beyond their PV thresholds.

This is a conditional mathematical statement about one tensor coefficient of the saved model. It is not an independent review, a novelty claim, constrained-CEV equivalence, a microscopic Lorentzian contour construction, or a healthy interacting theory. The scalar channel's negative regulator-open spectral regions remain. Adding bare curvature terms can invalidate the representation and its zero-exclusion conclusion.

## 1. Derive the spin-2 coefficient from the tensor Hessian
Choose unit transverse traceless h01=h10=1/sqrt(2), with momentum along axis 3. In the saved full tensor calculation, tr h=0, hQ=0, tr(h^2)=1. Its fourth-order angular factor is 1/12; its second-order and constant vertex factors vanish. The stationary contact contribution is p/2=-C sum c_i z_i^2 log(z_i)/8. Thus
K2(q)=C/8 sum_i c_i {integral_0^1 (z_i+tq)^2 log(z_i+tq) dx - z_i^2 log(z_i)},
where t=x(1-x), C=(4pi)^(-2).
This uses the actual contact prescription, not a projection imposed to force transversality.

An integration by parts gives the compact evaluation:
K2(q)=C/15 sum_i c_i {(z_i+q/4)^2 I_i(q) -(z_i^2-z_i q/8)log(z_i)},
I_i(q)=integral_0^1 log(z_i+tq) dx.
The principal-log continuation from q>0 has cut q<=-4z_0.

Low energy:
K2(q)=(G/24)q + (C/240 sum c_i log z_i)q^2 + O(q^3),
G=C sum c_i z_i log z_i>0.
The existing K0=2A/3 gives K0=-Gq/12+O(q^2), so K0=-2K2 at leading order, as required by the linear Einstein projector relation in these conventions. The independent q^2 coefficients are not matched between PV spectra here.

## 2. Standard identity, applied to this kernel
Set 0<z_0<z_1<z_2<z_3 and
w_i=1/product_(j!=i)(z_i-z_j), c_i=w_i/w_0.
The unique weights with c_0=1 and sum c_i z_i^k=0 for k=0,1,2 have this form; w_0<0.
Let lambda range over the three-simplex, Z=sum lambda_i z_i, with its ordinary coordinate volume d lambda (total volume 1/6).

The Genocchi-Hermite divided-difference identity, in Carl de Boor, Divided Differences, Surveys in Approximation Theory 1 (2005), 46-69, equation (52), represents a third divided difference by a simplex integral of the third derivative:
sum_i c_i f(z_i)=(1/w_0) integral_simplex f'''(Z) d lambda.
Primary author-hosted reference: https://ftp.cs.wisc.edu/Approx/dvdsurvey.pdf .
This established identity is the imported mathematical ingredient; the following application is our derivation.

For f_a(z)=(z+a)^2 log(z+a), f_a'''(z)=2/(z+a). Apply the identity at a=tq and at a=0, subtract, and integrate over x:
K2(q)=C/(4w_0) integral_x integral_simplex [1/(Z+tq)-1/Z] d lambda dx
=H q integral_x integral_simplex t/[Z(Z+tq)] d lambda dx,
H=-C/(4w_0)>0.

Consequently
B2(q):=K2(q)/q
=H integral_x integral_simplex [1/Z]/[q+Z/t] d lambda dx.
Endpoints t=0 have zero x measure and are understood by limits. Every weight is positive and s=Z/t>=4z_0. This is a positive Stieltjes representation with support in [4z_0,infinity). It follows first for positive q; analyticity and the identity theorem extend it to the first sheet. The simplex identity also applies directly for complex a away from its log cut.

The Stieltjes/complete-Bernstein terminology is standard; see Schilling, Song and Vondracek, Bernstein Functions: Theory and Applications (2010), chapters 2,6,7, https://doi.org/10.1515/9783110215311 . The zero argument below is explicit and does not depend on an uninspected theorem from that book.

## 3. Why this excludes additional first-sheet zeros
For Im q>0,
Im B2(q)=-H Im(q) integral_x integral_simplex (1/Z)/|q+Z/t|^2 d lambda dx <0.
Therefore B2 cannot vanish in the upper half-plane. Reflection excludes zeros in the lower half-plane. On the remaining real first-sheet interval q>-4z_0, the integrand in B2 is positive. Thus K2 has only its simple zero at q=0 on the open first sheet, and the reciprocal has only the massless pole there. Its residue is 24/G>0 in the fixed Euclidean tensor convention.

This proof excludes off-axis zeros globally for this exact kernel; the finite numerical contours below serve as checks, not as the reason for the global statement. The cut itself is not an open analytic domain. The following cut-sign argument excludes hidden real zeros on its lips above threshold.

## 4. Regulator-open spin-2 discontinuity
The independent cut formula is
Im K2(-s-i0)=-C pi s^2/240 sum_i c_i (1-4z_i/s)_+^(5/2).
For fixed s, g_s(z)=(1-4z/s)_+^(5/2) has nonpositive third derivative, strictly negative below its threshold:
g_s'''(z)=-(15/8)(4/s)^3(1-4z/s)^(-1/2), for z<s/4, and zero above.
The first two derivatives are continuous at the threshold and the third derivative is locally integrable, so the simplex identity extends by smoothing/absolute-continuity limits. With w_0<0, its PV sum is positive for every s>4z_0. Thus Im K2 is strictly negative and the inverse spin-2 continuum has positive weight on both physical and regulator-open intervals. This conclusion does not transfer to the scalar coefficient.

The same measure gives
K2(infinity)=-C/8 sum c_i z_i^2 log z_i>0
along the positive axis, by dominated convergence. Its finite limit does not justify dropping contact constants from an unsubtracted inverse representation.

## 5. Execution evidence
Run: python experiments/spin2_stieltjes_2026_10_30/check.py.
Dependencies: NumPy/SciPy plus the unchanged repository script scripts/metric_tensor_ward_check.py. The new script loads only its definitions, leaving historical report generation unexecuted. All original repository files are preserved.
Recorded versions: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0; long-double mantissa 63 bits. Assertions executed successfully.

- 25 compact-kernel comparisons with the full tensor Hessian: maximum absolute difference 9.37e-14.
- 20 off-axis principal-log parameter integrals at 1024 nodes: maximum relative difference 2.67e-14.
- 10 positive-measure simplex comparisons at 12,24,48 nodes per dimension: maximum final relative error 1.15e-6, with improving convergence. This is the least accurate numerical check, especially for L=100. It checks the representation's normalization; it is not a numerical proof of positivity.
- 30 adaptively resolved upper-half-plane contour checks of K2/q: zero winding throughout; maximum final phase step 0.12144 radians. Bounds R=4L,100L and imaginary floor 0.005; three initial resolutions.
- 15 once-subtracted kernel dispersion reconstructions, including all PV cuts to infinity: maximum relative difference 3.81e-11. The formula is K2(q)=-(q/pi) integral_4^infinity Im K2(-s-i0)/[s(s+q)] ds.
- 2,000 cut samples over s=4.00001 through 5000L: no negative inverse-continuum sign. The analytic representation supplies the sign result beyond sampling.

## 6. Review and next work
This derivation is ready for independent mathematical/physics review. Main audit points: contact normalization; the divided-difference factor w_0 and simplex volume; extension of the log identity to the slit plane; distinction between this fixed prescription and matched physical curvature terms.

Next apply the same mathematical tools to the scalar kernel, looking for a signed-measure decomposition or a finite-rank perturbation that controls its zeros. Do not infer scalar health from this spin-2 result. Then match G and both independent quadratic-curvature coefficients before comparing regulator prescriptions as one physical theory. Constrained measure, microscopic Lorentzian contour, interactions and 1/N corrections remain separate open problems.
