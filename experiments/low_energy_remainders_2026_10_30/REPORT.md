# Controlled low-energy matched-kernel remainder bounds

Date: 2026-10-30 US Central. Baseline: d21e4d562a7924a053d49475433e414192fb9076.

## Fixed comparison and result

Continue the explicit reference-L=10 matching of Einstein slopes and both independent tensor q^2 coefficients, without adopting the polynomials as a microscopic action. Work with A and K2; the scalar tensor coefficient is K0=2A/3. The physical first branch point is q=-4, and the following bounds apply only inside |q|<4.

This batch derives an absolute Taylor-tail inequality, propagates it to inverses, and checks it at 320 new complex/real points for both channels, L=3,5,10,30,100 and truncations N=2,4. All comparisons satisfy the bound, as do 320 inverse comparisons. No numerical tolerance allowance was added to these assertions. These are bounds on the specified model functions, not uncertainty bounds for an unknown UV completion or for omitted interactions.

## Analytic coefficients and a geometric bound

Let mu_n=integral_0^1[x(1-x)]^n dx=(n!)^2/(2n+1)!. The principal log expansion converges uniformly on every closed |q|<=r<4z_i:

I_i(q)=log z_i+sum_{n>=1}(-1)^(n+1) mu_n (q/z_i)^n/n.

For n>=3 write the channel coefficient as a_n=C r_n sum_i c_i z_i^(2-n). Directly multiplying the polynomial prefactors in the existing integral kernels gives

r_n^A=(-1)^(n+1)/4 [mu_n/n+mu_(n-1)/(n-1)+mu_(n-2)/(4(n-2))],
r_n^K2=(-1)^(n-3) mu_n/[4n(n-1)(n-2)].

These agree with the earlier exact rational coefficient formulas; this batch adds a remainder estimate rather than repeating the earlier small-q series scan. The identity

mu_(n+1)/mu_n=(n+1)/[2(2n+3)]<1/4

implies |r_(n+1)|<=|r_n|/4 for n>=3 in both channels, term by term in the scalar expression. Therefore for N>=2 and |q|<=r<4,

|K_matched(q)-T_N(q)| <= B_N(r,L),
B_N=C |r_(N+1)| r^(N+1) sum_i [|c_i| z_i^(1-N)/(1-r/(4z_i))].

T_N retains the matched coefficients at n=1,2 and that spectrum's own n>=3 coefficients. Real local matching changes no higher coefficients or tail. Absolute PV weights make the bound conservative and independent of signed cancellations. For two matched spectra, the difference from a common quadratic kernel is bounded by the sum of their B_2 values. The same reasoning supplies a pairwise bound after explicitly retaining their different cubic and higher terms.

This is an exact-arithmetic analytic inequality; the script evaluates its constants with floating-point arithmetic, not outward-rounded interval arithmetic. It is not appropriate to label all displayed decimal digits certified.

## Inverse response and zero exclusion

Where |T_N(q)|>B_N,

|1/K-1/T_N| <= B_N/[|T_N|(|T_N|-B_N)],
|(1/K-1/T_N)/(1/K)| <= B_N/|T_N|.

The second relation follows exactly from (1/K-1/T_N)/(1/K)=1-K/T_N. A common finite nonzero source factor preserves this relative error; cancellation of tensor channels in a full amplitude can invalidate a relative amplitude bound, requiring an absolute bound instead. q=0 is treated by factoring out the known massless zero, not by evaluating 1/K there.

For the quadratic truncation, T_2/q=a1_ref+a2_ref q. On a disk radius r,

|T_2/q|>=|a1_ref|-|a2_ref|r,
|(K-T_2)/q|<=B_2(r,L)/r.

If the lower bound exceeds the tail, K/q cannot vanish anywhere on that disk. All 40 evaluated disk inequalities (two channels, five spectra, radii 0.25,0.5,1,2) satisfy this condition with substantial margins. In particular the formulas exclude extra zeros of each explicitly matched channel throughout |q|<=2, conditional on the stated kernels and coefficients. Numerical evaluation of the margins is not an independent interval certificate; the analytic sufficient condition is explicit and reproducible.

The worst uniform quadratic inverse relative-error bounds over the five spectra are:

| Disk radius | A channel | K2 channel |
|---|---:|---:|
| 0.25 | 0.223% | 0.0129% |
| 0.5 | 0.964% | 0.0541% |
| 1 | 4.57% | 0.244% |
| 2 | 28.1% | 1.34% |

Values in the table are rounded upward. The scalar approximation is appreciably less controlled as the threshold is approached; these bounds do not extend through its known subthreshold poles near s=3.7. The previously found high-q spin-2 poles have |q| much greater than 4 and are outside the present analysis.

## Execution and next work

Run `python experiments/low_energy_remainders_2026_10_30/check.py`. Dependencies are the existing portable series and compact-kernel scripts, NumPy and SciPy. Tested radii 0.25,0.5,1,2 with phases 0,pi/3,2pi/3,pi and truncations 2,4. Maximum actual kernel error divided by its absolute bound is 0.980299 (rounded up); no bound failure occurred. Results retain all values, inverse bounds and disk margins. This does not independently verify every compact-kernel formula; those comparisons were already performed in previous batches.

Next extend scalar first-sheet zero control with subtracted analytic representations, using the present local bounds as a controlled region rather than another broad scan. Before interpreting composite poles, finish the common regulated source/measure dictionary. Preserve the distinction between exact model mathematics and physical applicability; no new coupling or contour is chosen.
