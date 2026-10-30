# Uniform UV expansion, exterior zero exclusion and cut-sign audit

Date: 2026-10-30 US Central. Baseline: bb5a0266bdf9ed03fc68a2b19af390cae0f30e3d.

## Outcome and scope

The raw scalar kernel is nonzero outside a finite radius in every tested spectrum, with an explicit uniform first-sheet remainder bound. The scalar kernels with the previously specified reference-L=10 matching also have certified exterior regions without zeros. These exterior results complement the saved conformal count; a finite intermediate region near the cuts remains unresolved, so they do not constitute a whole-first-sheet pole inventory.

The lower-cut scalar imaginary part is strictly negative sufficiently far above all regulator thresholds. Thus the inverse scalar spectral weight is eventually positive again. This does not repair the known intervening negative intervals, establish unitarity, or select a physical regulator. The words 'negative weight above PV thresholds' in earlier checkpoints describe the existence of such intervals, not an assertion that the entire infinite tail is negative.

A numerical-conditioning audit also finds that direct compact long-double formulas become unusable at sufficiently large q because of cancellations. The new convergent UV series avoids that problem. No original script or result was overwritten.

## 1. Explicit branch convention and independent expansion

Work with the existing free minimally coupled four-species PV determinant. Set C=(4pi)^-2, S_n=sum_i c_i z_i^n, T_n=sum_i c_i z_i^n log z_i. The PV conditions give S_0=S_1=S_2=0, while S_3=1/w_0<0 exactly.

For a species of mass squared z, x=z/q and beta=sqrt(1+4x), choose beta near one at infinity and the principal Log q. On the cut exterior q plane,

I_z(q)=Log q-2+(beta-1)Log(q/z)+2 beta log[(1+beta)/2].

This agrees with the defining parameter logarithm on positive q and by analytic continuation on both exterior half-planes. On the lower cut q=-s-i0, Log q=log s-i pi; beta is real positive for s>4z, yielding Im I=-pi beta. Regulator species are all retained.

Expand beta=sum b_j x^j and 2 beta log[(1+beta)/2]=sum a_j x^j. The script builds both coefficient lists exactly with rational polynomial recurrences. Initial coefficients are

b=[1,2,-2,4,-10,28,...], a=[0,2,1,-10/3,59/6,-449/15,...].

After canceling the three PV moments algebraically, the leading terms are

A/C=3 T_2/8 + [5(S_3 Log q-T_3)+S_3/6]/(4q)+O(Log q/q^2),
K2/C=-T_2/8 + [(5/4)(S_3 Log q-T_3)+(55/24)S_3]/(15q)+O(Log q/q^2).

In particular the previously used constants A_infinity=3C T_2/8 and K2_infinity=-C T_2/8 are confirmed, including their relative factor -3. The new calculation adds explicit correction coefficients and a uniform remainder, rather than inferring a limit from a few large-q samples.

The elementary logarithm and binomial series are established results: NIST DLMF 4.6.1, 4.6.4 and 4.6.7, https://dlmf.nist.gov/4.6, inspected directly. The regulator sums and zero/cut bounds below are derived here; they are not attributed to a new gravity theorem in that reference.

## 2. Uniform analytic remainder

Take the x-plane comparison circle rho=1/8. On |x|=rho, the principal beta has positive real part, and

|beta|<5/4, |beta-1|<=1/2,
|log[(1+beta)/2]|<=1/3,
|2 beta log[(1+beta)/2]|<=5/6.

Cauchy coefficient bounds therefore control the tail after N terms. For |q|>=R>=16 z_max, u_i=8z_i/R<=1/2. The species logarithm remainder is bounded by

[(5/4)(log(R/z_i)+4)+5/6] u_i^(N+1)/(1-u_i).

Here pi<4 bounds the principal argument. Multiply by (R/2+z_i)^2 |c_i|/4 for A/C, or (R/4+z_i)^2 |c_i|/15 for K2/C, then sum. All bound constants, mass moments and logarithmic enclosures are exact Fractions. The log enclosures use 128 positive-series terms with a geometric tail.

Multiplying the truncated species series produces terms through q^(-N), including its last two partial edge coefficients; these are retained explicitly rather than silently dropped. Let B_k,D_k be the rational coefficient convolutions with (1/4,-1,1) for A, or (1/16,1/2,1) for K2. The complete truncated correction is

normalization * sum_(k=3)^(N+2) q^(2-k)[B_k(S_k Log q-T_k)+D_k S_k],

with normalization 1/4 or 1/15. Bounding each moment term and adding the species remainder yields D(R)>=|K(q)-K_infinity| in scaled kernel units, uniformly over both exterior first-sheet half-planes. The terms decrease beyond R: powers q^(2-k) dominate their logarithms, and the species remainder has N=64, decreasing much faster than its quadratic prefactor. Thus a successful radius gives the entire exterior, not only a circle sample.

If D(R)<|K_infinity|, no raw zero exists at |q|>=R. The constant's lower modulus comes from its exact rational log interval. For the explicitly matched scalar kernel, compare the quadratic coefficient instead:

|Delta a2/C| > |Delta a1/C|/R + (|A_infinity/C|+D(R))/R^2.

This suffices for no matched scalar zeros beyond R. For L=10 the matching increments vanish and the raw criterion applies. Non-reference matched spin-2 kernels were not audited by this dominance test in this batch.

Smallest successful radii among the four candidates tested are:

| L | Raw scalar exterior | Explicitly matched scalar exterior |
|---|---:|---:|
| 3 | 1024 | 256 |
| 5 | 1664 | 416 |
| 10 | 3264 | 3264 |
| 30 | 38656 | 2416 |
| 100 | 128256 | 8016 |

Every entry means all |q| at or above that radius, off the cut. These are sufficient bounds, not optimal radii. Failed candidate inequalities are saved as such; their failure is not a claim that a zero exists there. Across both raw channels and all candidate radii, 33 of 40 exterior tests succeed. This does not supersede the stronger existing global raw spin-2 theorem.

## 3. Eventual cut sign

On q=-s-i0 the logarithmic corrections give

Im A/(C pi)=5 S_3/(4s)+O(s^-2),
Im K2/(C pi)=S_3/(12s)+O(s^-2).

Both leading signs are negative because S_3<0. Unlike the real corrections, the cut imaginary remainder has no log s factor: only beta's series contributes to it. Bounding every subleading imaginary term and the beta tail proves Im A<0 and Im K2<0 for all s>=16 z_max in each tested spectrum. All 40 candidate cut-sign inequalities succeed; the smallest radius is already sufficient.

Since Im(1/A)=-Im A/|A|^2, the inverse scalar cut weight is positive on this sufficiently high-energy tail. Real matching polynomials do not change Im A on a real cut, so the eventual sign statement also holds for those modified scalar inverses. The previously exhibited finite negative-weight intervals remain, and still obstruct a positive global spectral measure and the finite-negative-index method considered earlier.

## 4. Numerical checks and retained correction

Sixty comparisons use direct Decimal evaluation of the compact species formula at 160 digits, independently of the rational UV coefficient recurrence: five spectra, both channels, three radii, positive q and lower cut. All satisfy the analytic tail plus the separately tracked logarithmic-coefficient midpoint uncertainty.

The initial comparison code used only the analytic UV truncation tail when its logarithmic constants were themselves finite-accuracy midpoints. It failed when that tail became smaller than the uncounted log uncertainty (an error of order 1e-63 was being compared to a tail of order 1e-92). initial_failure.json retains twelve diagnostic cases. The corrected code includes the exact log-coefficient uncertainty and increases the reference precision/log truncation. This is a correction to the comparison error budget, not a changed physical assumption or evidence that a failed formula was declared healthy.

A separate twenty-point conditioning audit compares the compact long-double scalar expression with the high-precision and stable UV evaluations at q=(100,10000,1000000,100000000) z_max. The compact formula's largest relative error exceeds 2.5 million at these extreme points; the UV series agrees with the reference to better than 7.1e-91 relative error in the same sample. Those observed digits are numerical cross-checks, not interval-certified assertions about the Decimal reference implementation. The exterior/cut certificates themselves use the exact analytic inequalities, not these empirical comparisons.

The fifteen saved asymptotic-value rows from experiments/inverse_scalar_dispersion_2026_10_30/results.json were also checked specifically because of this new conditioning concern. Their largest relative discrepancy is 2.37e-9. The previously stated slow approach to the limiting constant remains qualitatively correct at those moderate radii; microscopic last-digit accuracy should not be inferred from those archived asymptotic samples. No previous full dispersion or Ward test was rerun.

## 5. Reproduction and continuation

Run:

1. `python experiments/uv_asymptotic_audit_2026_10_30/check.py`
2. `python experiments/uv_asymptotic_audit_2026_10_30/precision_check.py`

Files include code, all successful/unsuccessful bound candidates, sixty high-precision comparisons, twenty conditioning comparisons, fifteen legacy-value audits and retained initial failure diagnostics. Existing helper imports require the repository's NumPy/SciPy environment; the bound algebra itself uses the standard-library Fraction and Decimal modules.

Use the stable UV representation rather than unchecked compact formulas for future outer-boundary work. The exterior bounds do not fill the finite gap between the validated |w|<4/5 region and the cut boundary. Next audit the matched inverse dispersion with isolated real/complex poles and an explicit finite boundary, keeping completeness assumptions separate; a vanishing residual alone cannot prove an omitted pole inventory is complete. A whole-first-sheet count still needs control of the remaining finite collar and every regulator cut. Composite measure, Lorentzian contour, interactions and 1/N remain distinct open problems. No new action, regulator or physical coupling was selected.
