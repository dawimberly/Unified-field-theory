# Scalar zero count throughout the first-threshold disk

Date: 2026-10-30 US Central. Baseline: bc6c373e6f9909b00c89b98738f9e6c2c5831606.

## Result and scope

For the specified four-species PV family with L>=3, the raw scalar kernel, or its previously specified Lref=10 polynomial-matched version, has at most one additional zero in the entire open disk |q|<4. Any such zero is simple and lies at negative real q. This is stronger than the previous disk-|q|<=2 bound and the bounded numerical contours. It is not a theorem about the whole cut first sheet.

Exact rational enclosures for the threshold signs give the following complete additional-zero counts (excluding the known massless zero at q=0):

| L | Raw kernel | Explicitly matched kernel |
|---|---:|---:|
| 3 | 1 | 1 |
| 5 | 1 | 1 |
| 10 | 1 | 1 |
| 30 | 0 | 1 |
| 100 | 0 | 1 |

All listed extra zeros are the previously located subthreshold real poles. This batch adds uniqueness and off-axis exclusion throughout the disk, with exact threshold-sign certificates, rather than another broad root scan. Eight independent positive-series brackets check the archived locations. Zeros outside |q|<4, regulator cuts, the adjacent sheet, nonlinear measure equivalence, Lorentzian contour and physical health remain outside this result.

## 1. Positive inverse moments despite signed PV weights

Let z=[1,L,2L,5L+1], c_i=w_i/w_0, w_i=product_(j!=i)(z_i-z_j)^-1. For L>=3, w_0<0. Define

f(u)=sum_i c_i/(u+z_i)=-1/[w_0 product_i(u+z_i)]>0, u>=0.

This rational identity follows by partial fractions. At u=0,

f_0=f(0)=product_(i>0)(1-1/z_i), 0<f_0<1.

Let h_k(x_0,...,x_3) be the complete homogeneous polynomial of degree k with all monomial coefficients one, h_0=1. Expanding the product resolvent at zero gives

M_m=sum_i c_i z_i^-m=f_0 h_(m-1)(1,1/L,1/(2L),1/(5L+1))>0.

Thus positivity of all these inverse moments follows algebraically; it is not inferred from finitely many samples. The script verifies 80 instances with exact rational arithmetic.

The same resolvent yields the convergent integral identities

J=-sum_i c_i log z_i=integral_0^infinity f(u)du>0,
G/C=sum_i c_i z_i log z_i=integral_0^infinity u f(u)du>0,
C=(4pi)^-2.

## 2. Subtracted scalar series

Set s=-q and F(s)=A(-s)/s, filling in the removable value at s=0. In |s|<4,

F(s)=c-sum_(k>=1) d_k s^k,
c=-a_1>0, d_1=-a_2>0,
d_(n-1)=C rho_n M_(n-2)>0 for n>=3,

where the coefficients audited in the preceding checkpoint simplify to

rho_n=3 mu_(n-2)(3n^2-3n-1)/[16(n-2)(4n^2-1)],
mu_j=(j!)^2/(2j+1)!.

For the raw kernel, c=C integral u f/8 and d_1=CJ/16. For explicit matching, replace only c and d_1 by their reference-L=10 values. Higher d_k still belong to that spectrum. No bare action is adopted here.

We prove the coefficient ratio bound d_(k+1)<d_k/4 for every k>=1 and all L>=3, for both prescriptions.

## 3. First ratio: exact rational certificates

Since rho_3=17/1120,

d_2/d_1=(17/70) f_0/J

for a raw spectrum. The normalized positive integral J/f_0 increases as any regulator mass increases: its integrand is product_i z_i/(u+z_i). Therefore its minimum over L>=3 is at L=3. A right-endpoint Riemann sum on [0,16], using 256 intervals of width 1/16, is an exact rational lower bound. The script proves

(J/f_0)_(L=3)>0.9859582453...>34/35.

The decimal is illustrative; the comparison against the rational 34/35 is exact. It follows that d_2/d_1<1/4 for every raw L>=3.

For reference matching, f_0<1 and a separate exact rational Riemann bound gives

J_(L=10)>1.4101250841...>34/35.

Hence d_2/d_1=(17/70) f_0/J_ref<1/4 for every matched L>=3 as well. These certificates do not use numerical logarithms or interval-library rounding.

## 4. All higher ratios: finite prefix and induction

For n>=3,

rho_(n+1)/rho_n=(n-2)(3n^2+3n-1)/[2(2n+3)(3n^2-3n-1)].

Write k=n-3 and y=(1/L,1/(2L),1/(5L+1)). Because the remaining inverse mass is one,

M_(n-1)/M_(n-2)=h_(k+1)(1,y)/h_k(1,y)
=1+h_(k+1)(y)/h_k(1,y)<=1+h_(k+1)(y).

For L>=3, y is componentwise bounded by (1/3,1/6,1/16). At n=3,4,5, exact upper bounds on the product of these two ratios are respectively

875/4896, 167147/887040, 3555817/18849792,

all strictly below 1/4.

For n>=6, k>=3. Since every component of y is at most 1/3,

h_(k+1)(y)<=binomial(k+3,2)/3^(k+1)<=1/(k+1).

The second inequality starts with 60<=81 at k=3; the ratio of successive numerator polynomials (k+1)(k+2)(k+3)/2 is (k+4)/(k+1)<=3, proving the induction. Also

rho_(n+1)/rho_n < (n-2)/[4(n-1)],

because the reduced numerator difference is 3n^2-3n-5>0. Multiplying the two bounds gives d_n/d_(n-1)<1/4 for n>=6, covering the remaining indices. The infinite statement rests on these inequalities, not on the saved finite coefficient samples.

## 5. Complete disk zero count by deflation

The series F(s)=c-sum d_k s^k decreases strictly on 0<s<4. If F(4)>=0, triangle inequality shows it has no zero anywhere in |s|<4: for every r<4, sum d_k r^k<c. The finite threshold value is the continuous limiting value of the kernel.

If F(4)<0, continuity and monotonicity give a unique real r in (0,4) with F(r)=0. Factor it exactly:

F(s)=(r-s)E(s),
E(s)=sum_(k>=0) e_k s^k,
e_k=sum_(n>=k+1) d_n r^(n-1-k)>0.

These sums converge; the ratio bound gives e_(k+1)<=e_k/4 and ensures E is analytic in |s|<4. Set b_k=4^k e_k, a positive nonincreasing sequence. For |z|<1,

(1-z)E(4z)=b_0-sum_(k>=1)(b_(k-1)-b_k)z^k.

The absolute modulus of the sum on the right is at most |z| b_0<b_0. Therefore E has no zero in |s|<4. The real zero is simple since E(r)>0, and no complex zeros remain anywhere in that disk. This proves the classification using only the threshold sign. The massless zero of A is separate, since A(q)=-q F(-q).

This is the coefficient-monotonicity mechanism associated with Enestrom-Kakeya zero bounds, extended here directly to the convergent deflated power series. Primary mathematical background located: N. Anderson, E. B. Saff and R. S. Varga, *On the Enestrom-Kakeya theorem and its sharpness* (1979), author-hosted scan https://www.math.kent.edu/~varga/pub/paper_116.pdf. The scan has no extractable text in the retrieval service; no uninspected infinite-series corollary is imported. The short proof needed here is supplied explicitly above.

## 6. Exact threshold signs and new numerical cross-check

The script encloses A(-4)/C with rational arithmetic. It uses log z=k log2+log(z/2^k), evaluating each logarithm by the positive atanh series with a geometric tail, with argument at most 1/3. For regulator z>=3,

sqrt(z-1) atan[1/sqrt(z-1)]=sum_(n>=0)(-1)^n/[(2n+1)(z-1)^n].

Its alternating partial sum and next term enclose the value. The physical species at threshold is inserted exactly. With 64 logarithm terms and 128 arctangent terms, every raw/matched threshold interval has a strict sign. Maximum enclosure width for A/C is below 2.75e-40. Results store rounded endpoints, actual interval widths and classification labels; exact Fraction endpoints and comparisons are reproduced by check.py. The rounded displayed endpoints can coincide even though the underlying rational interval has positive width.

The existing eight real pole locations are independently enclosed numerically at plus/minus 1e-8 using 512 positive coefficients and their geometric tail; all signs bracket the zero. These location checks use long double and are not interval-certified. They are distinct from the exact threshold-sign/count proof. Maximum saved series tail is 5.09e-14. All assertions passed; no old broad contour or Ward scan was repeated.

Run `python experiments/scalar_threshold_disk_2026_10_30/check.py` from the repository root. Dependencies: existing scalar/threshold scripts and NumPy/SciPy. The coefficient and endpoint certificates themselves use only Fraction arithmetic and elementary series inequalities. Version information, 80 exact moment checks, two Riemann certificates, finite ratio-prefix bounds, ten disk classifications and eight new location checks are in results.json.

## Limits and continuation

Positive inverse moments and local Taylor coefficients do not imply a globally positive scalar spectral measure. The previously demonstrated negative continuum above regulator thresholds and infinite negative-index obstruction remain intact. No root count has been established for the rest of the first sheet outside |q|<4. The theorem also depends on the stated PV family and the stated low-energy matching; arbitrary extra local terms can invalidate its coefficient inequalities.

Next examine whether a conformal map or a cut-subtracted representation can extend scalar control beyond this disk while retaining all regulator cuts. Do not silently promote this disk theorem to global stability, constrained equivalence or UV completion. Continue the common regulated source/measure problem separately.
