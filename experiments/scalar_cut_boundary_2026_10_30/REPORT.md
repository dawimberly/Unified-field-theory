# Global first-sheet scalar count from the cut boundary

Date: 2026-10-30 US Central. Baseline: b1d6e710fe71ef0e549fc42f31bba20530a4194f.

The remaining first-sheet counting gap closes for the five specified raw PV spectra and their previously defined Lref=10 matched scalar kernels. This is a mathematical result for the free leading-N auxiliary determinant, not a statement about a constrained CEV path integral, a Lorentzian prescription, unitarity or a physical UV completion. No model parameter or matching prescription changed.

| L | Raw additional zeros | Explicitly matched additional zeros |
|---|---:|---:|
| 3 | 1 | 3 |
| 5 | 1 | 3 |
| 10 | 1 | 1 |
| 30 | 0 | 4 |
| 100 | 0 | 2 |

Counts include multiplicity and exclude the simple massless zero at q=0. The known isolated zeros exhaust these counts: raw subthreshold roots where present; matched subthreshold roots, previously isolated complex pairs or positive-axis roots, and the newly saved L=30 collar pair. Thus there are no further first-sheet zeros for these ten coefficient kernels. The statement is specific to the five spectra, not all L or other regulator families, and concerns A, not an unqualified arbitrary source-contracted observable.

## Branch and cut formula

Use q=Q_E^2, M=1, C=(4pi)^-2, z=[1,L,2L,5L+1] and c_i=w_i/w_0, w_i=1/product_(j!=i)(z_i-z_j). In particular w_0<0. The defining logarithm integral is analytic on C minus (-infinity,-4], with the principal Log of z+q x(1-x). On the lower lip q=-s-i0,

Im I_z=-pi sqrt(1-4z/s) for s>4z, and zero otherwise.

Consequently

Im A(-s-i0)=-C pi H(s)/4,
H(s)=sum_(4z_i<=s) c_i (z_i+s/2)^2 sqrt(1-4z_i/s).

The previously specified real matching polynomial leaves H unchanged. It changes Re A on the cut, so the crossing classification is computed separately for raw and matched kernels.

For each open species contribution, put beta=sqrt(1-4z/s). Its real integral is log z-2+beta log[(1+beta)/(1-beta)]. For a closed species put u=sqrt(4z/s-1); it is log z-2+2u atan(1/u). At threshold the added piece is zero. These formulas, exact rational interval logs/atan/sqrt, and the existing kernel formula enclose Re(A/C) throughout each candidate interval. Matching increments are included explicitly, with the signs from q=-s.

## Complete classification of imaginary-part zeros

1. For 4<s<4L only the physical species is open and H>0.
2. On 4L<s<8L there are two species. Define polynomial X_i=c_i^2 (s+2z_i)^4(s-4z_i). Any H zero is a zero of P_2=X_0-X_1, degree five.
3. On 8L<s<4(5L+1) there are three species. Any H zero is a zero of P_3=(X_2-X_0-X_1)^2-4X_0X_1, degree ten.
4. Above the last threshold all four species are open. For f_s(z)=(z+s/2)^2 sqrt(1-4z/s), y=1-4z/s>0,

f_s'''(z)=-(27 y^(-5/2)+18 y^(-3/2)+15 y^(-1/2))/(2s)<0.

The signed sum is the third divided difference divided by w_0. The Genocchi-Hermite formula therefore gives H>0 for every s>4z_max, not merely the sufficiently large UV tail. The last threshold itself is checked separately.

The established divided-difference input is Carl de Boor, *Divided Differences*, Surveys in Approximation Theory 1 (2005), 46-69, section 9, equation (52), inspected at https://ftp.cs.wisc.edu/Approx/dvdsurvey.pdf . Only the simplex-integral sign property is imported; the cut derivative and PV application are derived here.

Exact rational Sturm sequences show exactly one polynomial root on each of the two finite intermediate intervals, for every spectrum. Rational bisection isolates each in width below 1e-14. The saved polynomial coefficients, sequence degrees and rational interval endpoints make the root counts reproducible. Polynomial squaring could introduce extraneous roots; this is not assumed away. Exact cut bounds give H signs +,-,+ at the three regulator thresholds. The intermediate-value theorem therefore supplies an actual H zero in each interval. Since the candidate polynomial has exactly one root there, these are the only H zeros and both are sign crossings. No zero is missed by squaring.

At all ten isolated crossings, Re A is separated from zero for both prescriptions. All twenty thresholds are also nonzero, including the physical threshold where H=0. Thus neither lip of the entire cut has a kernel zero. Fractions with outward rounding onto the rational 10^-30 grid are used; no numerical phase or approximate root is needed for this exclusion.

Useful approximate crossing summaries (not replacements for the rational data):

| L | First s crossing | Second s crossing | Re A signs, raw | Re A signs, matched |
|---|---:|---:|---|---|
| 3 | 12.9351260 | 41.5612022 | +,- | +,- |
| 5 | 21.1379165 | 67.6056109 | +,- | +,- |
| 10 | 41.7740086 | 132.7581532 | +,- | +,- |
| 30 | 124.4783593 | 393.4290282 | +,- | +,- |
| 100 | 414.0349321 | 1305.8158825 | +,- | +,+ |

## Argument-principle count, with orientation explicit

A has no first-sheet poles. Its zero at q=0 is simple: A'(0)/C=-sum c z log z/8<0, enclosed separately in the code (reference value after matching). Existing UV analysis proves uniform asymptotics on both first-sheet half-planes and excludes zeros outside explicit finite radii. For raw kernels,

A/C -> 3 sum c z^2 log z/8 < 0.

For non-reference matched kernels A/C ~ d_2 q^2, with d_2 the nonzero scaled matching coefficient; its sign is negative for L=3,5 and positive for L=30,100. Matched L=10 is unchanged. The limiting outer-circle phase contribution is therefore 2pi*m, m=0 for raw/reference and m=2 otherwise. This uses the previously bounded asymptotic formula, not a finite number of outer samples.

Orient the slit disk positively: outer circle counterclockwise from its lower negative-real endpoint to its upper endpoint; upper lip from s=R down to s=4; lower lip from 4 up to R. If phi(s) is a continuous lower-lip phase, reflection makes the two lips contribute 2[phi(R)-phi(4)]. Small threshold detours have vanishing contributions because the kernel is continuous and nonzero there; integrable square-root derivative singularities do not require A to extend analytically across a threshold. Nonvanishing on a compact cut segment allows contours approaching its lips without zeros crossing. The existing exterior bound prevents zeros escaping during the large-R limit.

Initially Im A<0. Hence phi(4)/pi is -1 if A(-4)<0, and 0 if A(-4)>0. Exact enclosures give -1 for raw L=3,5,10; 0 for raw L=30,100; and -1 for every matched case.

At the first H crossing, Im A changes from negative to positive. A crossing of the negative real axis would decrease the phase branch by one full turn; here every first crossing has positive Re A and does not change that branch index. At the second crossing Im A changes from positive to negative. Negative Re A adds one full turn; positive Re A adds none. Thus the branch-turn index k is 1 in every case except matched L=100, where k=0.

At infinity the lower lip still has Im A<0. Its principal phase tends to -pi if the leading real coefficient is negative, and to zero from below if positive. Adding 2pi*k yields the continuous phi(infinity). Therefore

N_A = m + [phi(infinity)-phi(4)]/pi.

The exact sign data give total counts 2,2,2,1,1 raw and 4,4,2,5,3 matched. Subtracting the independently checked simple massless zero gives the table above. Counts are derived from cut signs, asymptotics and the argument principle; a small dispersion residual or a converged sampling grid is not used as a completeness proof.

The standard argument principle is NIST DLMF 1.10.9, https://dlmf.nist.gov/1.10#E9 , inspected directly. The orientation, threshold limiting argument, phase bookkeeping and model-specific signs are given explicitly here.

## Independent diagnostics and corrections

Thirty numerical slit windings at three resolutions (512,2048,8192 outer nodes) agree with the exact counts. They use radii twice the previously certified exterior radii, retain every threshold, and use the stable UV evaluator on the outer boundary. The maximum sampled phase step is about 0.804 rad, below pi/2. These are diagnostic discretizations, not replacements for the analytic boundary proof.

Thirty independent real parameter-log integrals check the new cut real-part expressions at the candidate centers and every threshold. Maximum relative discrepancy is below 1.88e-12 in the final evaluation. Reported quadrature warnings and error estimates are retained in numerical_results.json; they do not enter the exact proof. An initial direct integrand log|z-s*x*(1-x)| suffered cancellation near threshold and issued a slow-convergence warning. Its results are preserved in initial_numerical_results.json. The final evaluation factors the quadratic into its two real roots on the open cut before taking real logarithms. Residual quadrature warnings remain explicit rather than being declared certified accuracy.

The initial exact implementation also stopped at an assertion when rounded reciprocal operations enclosed an exactly zero threshold beta^2 by a small interval straddling zero. The corrected branch recognizes s=4z exactly and sets beta^2=0. This avoids choosing a wrong branch from rounding; no sign tolerance or physical condition was relaxed. audit_notes.json records this failure and correction.

Reproduce with:

1. `python experiments/scalar_cut_boundary_2026_10_30/check.py`
2. `python experiments/scalar_cut_boundary_2026_10_30/numerical_audit.py`

The exact proof uses Fraction arithmetic and existing rational elementary-function enclosures. Existing helper modules import the repository's numerical dependencies, but floating-point values do not decide its inequalities. No original scripts, hypotheses or historical results were overwritten.

## Continuation

Use the now complete inventories for an infinite-boundary subtracted dispersion reconstruction. For the four non-reference matched cases, investigate the convergent superconvergence sum rule of the coefficient inverse, with every isolated pole and regulator cut included, and derive a quantitative tail bound before claiming a numerical precision. Reference L=10 and all raw kernels have a different constant asymptote and cannot be assigned the same unsubtracted moment rule. Preserve the distinction between 1/A and the complete source exchange. The microscopic constrained measure, common source prescription, Lorentzian contour, interactions and 1/N remain open and unchanged.
