# Matched inverse reconstruction and an additional collar pole pair

Date: 2026-10-30 US Central. Baseline: fd6102fd0420d6ac7ffdabb22317f3a11707a9b0.

## Outcome

The inverse-reconstruction audit exposed an additional first-sheet pole pair in the explicitly Lref=10 matched L=30 scalar kernel:

q=-407.9744557683 +/- 12.2684914375 i.

An independent exact-rational local Rouché certificate isolates the upper root within radius 1e-10; reflection gives its conjugate. Its conformal modulus is |w|=0.9970138267..., outside the previously validated |w|<0.8 region. Its |q| is below the previously established matched UV exterior radius 2416. Thus it lies precisely in the still-open finite collar and contradicts no earlier bounded count.

The pair was not fitted to the reconstruction residual. Its locations come from the kernel, residues from an analytic derivative independently checked by parameter quadrature, and existence from a separate rational interval calculation. Adding those residues repairs the failed slit reconstruction without changing any model parameter. The resulting closure is a consistency check, not a proof of a globally complete pole inventory.

This remains the free leading-N auxiliary determinant with the already explicit polynomial matching. No physical UV action, contour, constrained-CEV equivalence, discovery or healthy theory is established. The new poles lie well outside the low-energy derivative expansion.

## 1. Two reconstruction domains must be distinguished

Use F(q)=1/A_matched(q), with q0=1. These are scalar coefficient inverses, not the complete source-contracted amplitude W; source factors and their contact terms must be applied separately.

Inside the prior validated conformal domain q(w)=16w/(1-w)^2, |w|<0.8, define H(w)=F(q(w)). If F has residue r_j at q_j, then H has residue

r_j^(w)=r_j/q'(w_j), q'(w)=16(1+w)/(1-w)^3.

The massless residue transforms by q'(0)=16 as well. Cauchy's residue theorem gives

H(w)-H(w0)=sum_(inside |w_j|<0.8) r_j^(w)[1/(w-w_j)-1/(w0-w_j)]
+ (1/(2pi i)) integral_(|z|=0.8) H(z)[1/(z-w)-1/(z-w0)]dz.

The inventory inside that contour is complete by the prior validated count. Poles outside it belong in the finite boundary's analytic contribution, not in the explicit residue sum. The newly identified pair must therefore be excluded from this bounded residue sum. Including it there would double-count its effect.

Fifteen such checks (five spectra, q=0.5,2,1+2i) agree to maximum relative closure error 3.86e-12, with 128/256/512-node boundary convergence data saved. Omitting the earlier extra matched poles changes this reconstruction by up to 0.00149 relative error. These checks verify coordinate/residue normalization, including complex pairs.

## 2. Slit-cut reconstruction, with finite boundary retained

For the slit q disk of radius R, define the actual matched density

rho_matched(s)=Im[1/A_matched(-s-i0)]/pi.

The matching polynomial is real on the cut, so it leaves Im A unchanged, but it changes |A| and hence the inverse spectral density. Reusing the raw inverse density would be incorrect.

With all relevant poles inside the slit disk, the once-subtracted relation is

F(q)-F(q0)=sum_j r_j[1/(q-q_j)-1/(q0-q_j)]
+ integral_4^R rho_matched(s)[1/(s+q)-1/(s+q0)]ds
+ (1/(2pi i)) integral_(|z|=R) F(z)[1/(z-q)-1/(z-q0)]dz.

Every physical and regulator threshold is a separate cut interval. A sin-squared parameter change smooths both endpoint square roots. The lower-lip sign is imposed explicitly, not via a small untracked imaginary regulator. Outer radii are 16z_max and 64z_max. For those outer boundaries and the far-cut portion, the new stable UV series is used. Thirty wrapper comparisons against the previous 160-digit reference values validate that software path; this is not a repeat of the old physical calculation.

Initially, using only the previously known poles, the L=30 residual reached 2.53503e-7 relative and did not decrease when R increased. That failed assertion is preserved in initial_results.json. Independent fixed-node 128/256/512/1024 quadrature agreed with adaptive integration on every cut segment; quadrature_diagnostic.json retains the values. The residual was not removed by loosening the acceptance threshold.

A targeted 48-seed search of the previously uncounted L=30 collar located the new pair. This search alone was not treated as proof of existence or completeness. See collar_search_diagnostic.json and the independent certificate below.

## 3. Independent local existence certificate

local_pole_certificate.py performs the certificate entirely with Fraction arithmetic, integer square roots, and elementary series enclosures. It does not use NumPy complex evaluation to certify the root. All interval endpoints are rounded outward onto the rational grid 10^-30 after each operation.

For q in the upper half-plane, beta_z=sqrt(1+4z/q) has positive real part and negative imaginary part. The quotient (beta_z+1)/(beta_z-1) lies strictly in the upper half-plane. Therefore

I_z=log z-2+beta_z Log[(beta_z+1)/(beta_z-1)]

is evaluated with an explicit principal logarithm. This expression is analytic in the upper half-plane and agrees with the defining integral on the positive-q boundary. The root disk never crosses a cut. Log-modulus bounds use rational positive-log series; arguments use a Machin enclosure for pi, reciprocal/arctangent identities and alternating series at arguments of magnitude at most 1/2. Square-root bounds use exact integer square roots. No phase is guessed from a floating-point argument routine.

Differentiate this same expression analytically:

I'_z=1/q-2z Log[(beta_z+1)/(beta_z-1)]/(q^2 beta_z).

The script encloses A/C and A'/C at a refined rational center. Midpoint Newton steps help choose the center; they are not the certificate itself. A separate analytic second-derivative bound supplies the certificate. For a circle radius r, let Q bound |q| and y>0 bound Im q from below. With t=x(1-x),

|I'_z|<=1/y, |I''_z|<=1/y^2,
|z+qt|>=z y/Q, |z+qt|<=z+Q/4.

The latter two inequalities bound |I_z| by log-modulus bounds plus pi<4. Since

(A/C)''=(1/4) sum_i c_i[I_i/2-2U_i I'_i+U_i^2 I''_i]+2 Delta a2/C,

where U_i=z_i-q/2, these estimates give an explicit rational bound M2 for the circle. Rouché's inequality is

|A(center)/C|+M2 r^2/2 < |A'(center)/C| r.

The modulus upper bounds use the sum of absolute component bounds; the derivative lower bound uses one separated component. At r=1e-10, the left side is at most 1.152e-17 and the right side at least 5.84765e-9. The inequality is checked exactly as a rational comparison. Hence the disk contains exactly one zero, counted with multiplicity, so the root is simple. Schwarz reflection supplies the lower partner. A coarser radius-0.001 certificate is retained separately as part of the refinement history.

All displayed decimal endpoints are summaries of rational enclosures, not assertions that every written center digit is certified. The radius-1e-10 localization is the explicit claim. The method relies on the stated model formula and elementary series inequalities, not on AI agreement or a small floating-point residual.

## 4. Independent parameter checks, residues and repaired reconstruction

Defining principal-log parameter integrals at both conjugate locations have numerical kernel residual about 1.10e-13. Their conservative accumulated quadrature error estimate is about 1.25e-8; that numerical check alone is therefore not an existence proof. It supports the branch assignment, while the exact local certificate establishes the root. Root-neighborhood split points resolve narrow endpoint layers. SciPy split points are converted explicitly to float64; an intermediate runtime type mismatch was corrected without changing the mathematics.

Analytic residues for F at the new pair are approximately

r_upper=-1.6761537244+1.3103048575 i,
r_lower=conjugate(r_upper).

The independent derivative quadrature agrees at about 2.05e-15 relative. The residue values are numerical approximations, not positivity or norm statements.

Adding this pair to the slit-domain residue sum produces thirty closures with maximum relative error 6.02e-15. Omitting the pair restores the original maximum residual 2.53503e-7. Omitting the finite outer boundary causes errors up to 6.11e-6 in the tested cases, so it is not silently dropped even when a particular boundary happens to be small. No cutoff, counterterm, source or parameter was adjusted to achieve closure.

The successful reconstructions do not prove that no further poles lie elsewhere in the finite collar: only finitely many probes were tested, and omitted contributions could be small or partially cancel. In particular, a hypothetical extrapolation from the previously stable |w|<0.8 counts would have missed this pair at |w| approximately 0.997. It is a concrete reason to retain the global-completeness qualification.

## 5. Reproduction and continuation

Run:

1. `python experiments/matched_inverse_reconstruction_2026_10_30/local_pole_certificate.py`
2. `python experiments/matched_inverse_reconstruction_2026_10_30/check.py`
3. `python experiments/matched_inverse_reconstruction_2026_10_30/collar_parameter_check.py`

The certificate's seed diagnostic is saved and is only a locator; the Rouché comparison proves local existence independently. Artifacts include code, the initial failed reconstruction, independent cut quadrature diagnostics, targeted search seeds, coarse/fine certificates, final pole residues, all boundary convergence values and parameter checks. The executed numerical reconstruction used a 63-bit long-double mantissa, recorded in results.json; the rational local certificate does not depend on that floating-point precision. Existing helpers and original model files are unchanged. NumPy/SciPy are used for reconstruction and cross-checks; the local certificate calculations use exact rational arithmetic.

Established mathematical inputs are Cauchy's residue theorem and Schwarz reflection, checked in NIST DLMF sections 1.10(iv) and 1.10(ii), https://dlmf.nist.gov/1.10. The application, subtractions, residue transformation, spectral density and certificate bounds are specified here. No external result is cited as proving this model physically viable.

Next prioritize boundary/cut-based global scalar counting with the established UV exclusion, rather than a prohibitively slow uniform expansion toward |w|=1. Continue auditing the remaining finite collar, including possible cut-boundary zeros and every regulator threshold. Keep constrained measure, the common source prescription, Lorentzian contour, interactions and 1/N distinct. No new physical model choice is needed for these mathematical audits.
