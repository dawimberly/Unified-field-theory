# Scalar conformal continuation and validated bounded pole count

Date: 2026-10-30 US Central. Baseline: d766b9696dd7e0f5fd900eff799ec7fe9a6d4051.

## Outcome and interpretation

A first-sheet conformal representation, with exact rational coefficients and directed Decimal contour enclosures, extends the previous |q|<4 analysis. The explicit reference-L=10 matched scalar kernels have additional poles beyond their known subthreshold poles:

| L | Newly isolated matched scalar q poles |
|---|---|
| 3 | -17.4662104652 +/- 12.7654098460 i |
| 5 | -25.3822329364 +/- 30.4480072636 i |
| 30 | +34.2768433995 |
| 100 | +27.5993608103 |

These are zeros of the full polynomial-modified scalar kernel on the principal sheet. A positive Euclidean q pole would be tachyonic if that particular full kernel were adopted with a physical Lorentzian prescription. No such microscopic action or prescription is adopted here. The extra roots have |q|/4=5.41--9.91, outside the derivative expansion's first-threshold radius. They are not controlled low-energy EFT predictions or evidence about constrained CEV. Together with the earlier matched spin-2 findings they reinforce why finitely matched coefficients do not identify full resummed UV theories.

The result is bounded: it covers |w|<4/5 under the map below, not the entire first sheet. No physical health, discovery, unification or UV completion follows.

## 1. Explicit principal-sheet map and all cuts

Define

w=(sqrt(1+q/4)-1)/(sqrt(1+q/4)+1), q=16w/(1-w)^2.

The square root is positive at q=0 and has positive real part off the cut q<=-4. The open cut q plane maps one-to-one to |w|<1. The first threshold is w=-1, infinity is w=1. Every regulator threshold q=-4z_i maps to

w_i=(i sqrt(z_i-1)-1)/(i sqrt(z_i-1)+1)

on the upper semicircle and its conjugate on the lower semicircle. No regulator cut is dropped. Only the bounded region |w|<4/5 is counted in this run; on the real q axis its endpoints are -320/81=-3.950617... and +320. Complex points beyond the Taylor radius are included, but sufficiently close approaches to any cut and to infinity remain outside.

## 2. Exact conformal coefficients and a rigorous tail

Let t=x(1-x), a_i=1-8t/z_i in [-1,1], and

P_i(w)=z_i(1-w)^2-8w,
b_(i,n)=integral_0^1 T_n(1-8x(1-x)/z_i)dx,

where T_n is a Chebyshev polynomial. Factorizing the parameter logarithm gives, throughout |w|<1,

I_i=log z_i-2 log(1-w)-2 sum_(n>=1) b_(i,n) w^n/n.

The potentially large common log(1-w) cancels exactly: sum_i c_i P_i(w)^2=0 by the three PV moments. Define the analytic, massless-deflated function

G(w)=4(1-w)^4 A(q(w))/(C w), C=(4pi)^-2.

Its removable value at zero is nonzero. The remaining representation is

G(w)=-8 G_L (1-w)^2+64 S_L w
-2 sum_i c_i P_i(w)^2 sum_(n>=1) b_(i,n) w^(n-1)/n,

G_L=sum_i c_i z_i log z_i, S_L=sum_i c_i log z_i.

Explicit Einstein/q^2 matching replaces only G_L,S_L in this low-degree logarithmic polynomial by their reference-L=10 values. This is the same specified matching as earlier checkpoints, not a new model choice.

All b_(i,n) are evaluated exactly by rational Chebyshev recurrence and the moments integral t^k dx=(k!)^2/(2k+1)!. For z_0=1 the independent simplification b_(0,n)=1/(1-4n^2) is used. Logarithmic constants are enclosed by the exact rational log series already audited in the previous checkpoint.

Since |T_n(a_i)|<=1, truncation after N logarithmic coefficients gives, uniformly for |w|<=r<1,

|G-G_N|<=E_N(r),
E_N(r)=2 r^N/[(N+1)(1-r)] sum_i |c_i|[z_i(1+r)^2+8r]^2.

This tail is evaluated as an exact rational number then rounded upward. It retains every PV species. The function normalization removes C, so no numerical enclosure of pi is required for the counts.

## 3. Validated contour implementation

For N=128 and r=4/5, the coefficients of G_N are real rational intervals. Interval arithmetic rounds each elementary operation outward at 60 decimal digits. Negation uses exact copy_negate; midpoints for the winding polygon are exact rational averages of Decimal endpoints.

The circle is parameterized by rational points r(1-t^2+2it)/(1+t^2), t in [0,1], rotated through four quadrants. On each adaptively bisected arc, coordinates are monotone. Hence the distance from its first endpoint is at most the exact endpoint L1 coordinate change delta.

At endpoint a, evaluate G_N(a) and G_N'(a) with directed interval arithmetic. For the whole radius-r disk, the exact coefficient norm bounds

sup |G_N''|<=M2=sum_(n>=2) n(n-1) max|coefficient_n| r^(n-2).

Every polynomial value on the arc lies in the endpoint image rectangle expanded by

sup|G_N'(a)| delta+M2 delta^2/2.

The endpoint derivative uses the conservative sum of absolute real/imaginary component bounds. The rectangle is enlarged to contain the next endpoint's evaluation interval as well. An arc is accepted only when one real or imaginary component of this complete rectangle is separated from zero by more than E_N. Thus both the polynomial and full analytic function avoid zero along the arc, and the polynomial path can be deformed to the segment joining the rational nodal midpoints without crossing zero. The full-function and polynomial counts agree by the Rouché inequality. The nodal polygon's winding is computed by exact rational crossing/orientation tests, not floating-point phase unwrapping.

The same routine counts small shifted circles around the candidate poles, while retaining the global r=4/5 tail/derivative bounds. Each of four newly identified representatives has isolating count one; conjugation supplies the lower partner for each complex pair. Eight known subthreshold poles are also isolated. Those inventories exhaust the bounded counts below. Software execution is not a formal proof-assistant verification; the analytic enclosures and arithmetic code are exposed for independent review.

## 4. Complete bounded counts and independent checks

Additional zero counts inside |w|<4/5, excluding q=0:

| L | Raw | Matched to L=10 |
|---|---:|---:|
| 3 | 1 | 3 |
| 5 | 1 | 3 |
| 10 | 1 | 1 |
| 30 | 0 | 2 |
| 100 | 0 | 2 |

Every accepted arc has a strictly positive saved separation after subtracting the analytic tail. The largest N=128 tail is about 4.31e-8 in G units; the smallest saved contour margin is about 2.99e-4. Focused recomputation at N=160 and 80 Decimal digits preserves the matched counts 3 and 2 for L=3 and L=100, respectively. Increasing order/precision supplies a new reason for these focused repeats; no previous broad contour or Ward scan was rerun.

For each new root, independently evaluate the defining principal-log parameter integral at 128,256,512 Gauss-Legendre nodes, then add the explicit matching polynomial. Maximum absolute scalar-kernel residual is 2.81e-14. The compact-kernel roots themselves have residuals below 9.30e-16; these decimal root locations are numerical approximations, not certified last digits. The isolating w circles have rational radius 1/200. Three derivative steps support nonzero, convergent residues; no positive-norm interpretation is inferred from these residues.

Twelve additional overlap points compare the independently constructed conformal representation against the compact kernel. Maximum relative overlap error is 2.64e-15. Fifteen real arithmetic pairs check sum/product/negation enclosure with exact Fractions; a separate fifteen-case audit checks complex products and midpoint containment. All final assertions pass.

## 5. Development corrections retained as part of the audit

An early interval draft used ordinary unary Decimal negation and Decimal midpoint arithmetic outside the directed contexts, risking ambient 28-digit rounding. It was rejected before any accepted count was reported; exact copy_negate and Fraction midpoints replaced these operations. The saved arithmetic audits cover the corrected path. Direct Horner propagation of a whole complex arc rectangle was also too conservative due to interval wrapping; those attempts were interrupted without a count claim. A global first-derivative-only enclosure was then refined to the rigorous endpoint-derivative/global-second-derivative bound above for usable performance. Only final accepted results are used.

These corrections are mathematical/software verification steps, not evidence of a changed physical model. No failed numerical result was reclassified as a passing physical test.

## 6. Reproduction, sources and continuation

Run, in order:

1. `python experiments/scalar_conformal_count_2026_10_30/check.py`
2. `python experiments/scalar_conformal_count_2026_10_30/verify.py`

check.py uses the standard library plus existing helper imports; verify.py additionally uses NumPy/SciPy for independent root/quadrature checks. All model definitions, raw/matched counts, strict margins, refinement counts, isolated roots, derivatives and overlaps are in results.json and verification_results.json; arithmetic_audit.json records the separate exact complex arithmetic check.

Established numerical-method context: S. M. Rump, *Verification methods: Rigorous results using floating-point arithmetic*, Acta Numerica 19 (2010), 287--449. Publisher extract inspected at https://www.cambridge.org/core/journals/acta-numerica/article/abs/verification-methods-rigorous-results-using-floatingpoint-arithmetic/770FE58E5293985CCAB770AF09C4F3FF. The author PDF retrieval timed out. This is context for verification rather than a borrowed implementation or a gravity result; the specific map, tail and contour arguments used here are supplied above.

Next: extend the radius only with explicit tail/control estimates, and audit large-q/cut asymptotics before attempting a whole-first-sheet count. The remaining collar 4/5<=|w|<1 is unresolved. The prior signed scalar continuum problem is unchanged. Constrained measure, source matching, a physical Lorentzian contour, interactions and 1/N remain distinct problems.
