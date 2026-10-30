# Complete matched spin-2 first-sheet zero inventory

Date: 2026-10-30 US Central. Baseline: 649a2d4b287407342ab5fb9e7fd2b4c8a4591d3c.

The additional-zero counts of the five explicitly Lref=10 matched spin-2 kernels are now complete on the open first sheet, with existence and real-axis signs certified by rational inequalities. They are 1,1,0,2,2 for L=3,5,10,30,100. These are the previously located roots; no new pole was found. The gap closed here is the distinction between numerical evidence and complete existence/exclusion control.

| L | Additional spin-2 zeros | Location/type |
|---|---:|---|
| 3 | 1 | q approximately 60.3330308219, positive real |
| 5 | 1 | q approximately 82.4469390832, positive real |
| 10 | 0 | Unmodified reference, by positive Stieltjes representation |
| 30 | 2 | q approximately 8.72600156453 +/-78.0620223980 i |
| 100 | 2 | q approximately 10.0532001673 +/-66.1017934580 i |

Every listed zero is simple. The common massless zero is excluded from this table and remains simple. There are no cut-lip zeros, including at the first or regulator thresholds. The counts concern coefficient functions, not a count of independent physical polarizations or particles.

Together with the preceding scalar result, both conserved tensor coefficient inventories are complete for the stated raw and explicitly matched free leading-N auxiliary prescriptions. This does not establish constrained-CEV equivalence, a physical UV completion, a Lorentzian integration contour or a healthy interacting theory. The matched roots remain far outside the low-energy derivative expansion radius. The raw fixed-prescription spin-2 theorem is unchanged; scalar negative continuum intervals remain.

## 1. Reuse the established measure, do not repeat scans

Write S_L(q)=K2_raw(q)/q and

F_L(q)=K2_matched(q)/q=S_L(q)+a_L+b_L q,
a_L=k1_ref-k1_L, b_L=k2_ref-k2_L.

The raw divided-difference derivation supplies

S(q)=integral dmu(t)/(q+t), dmu>0, support t>=4,
S(0)=k1>0, S'(0)=k2<0.

The measure has nontrivial continuous support. This is the representation established in experiments/spin2_stieltjes_2026_10_30/REPORT.md, from the actual stationary determinant contacts. The matched function is not assumed to retain the Stieltjes sign property. The earlier matched audit derived the analytic restrictions used below, but its real-axis discriminants and complex-root existence were numerical. This batch fills precisely those gaps.

The first sheet is q in C minus (-infinity,-4]. Branches are fixed by the principal parameter logarithm I_z(q)=integral_0^1 Log[z+q x(1-x)]dx. Reflection holds because masses, weights and matching coefficients are real.

## 2. Global off-axis restriction

For Im q=y>0,

Im F(q)/y=b-integral dmu(t)/|q+t|^2.

If b<=0 this is strictly negative, excluding all nonreal zeros globally. If b>0 and z,w are upper-half-plane zeros, subtract F(z)-conjugate(F(w)). Since z-conjugate(w) is nonzero,

integral dmu(t)/[(z+t)(conjugate(w)+t)]=b.

The same identity at z=w gives the two diagonal squared norms b. Therefore the squared norm of (z+t)^-1-(w+t)^-1 is zero. Positivity of the measure forces z=w. There is at most one distinct upper-half-plane zero. At a zero,

|integral dmu/(z+t)^2| < integral dmu/|z+t|^2=b,

because its phase varies over continuous support. Hence F'(z)=b-integral dmu/(z+t)^2 is nonzero. Any upper zero is simple; reflection supplies its simple lower partner. This proof covers the entire half-plane, not only rectangles sampled previously.

## 3. Complete real-axis exclusion and existence

All coefficient signs are recomputed using outward rational log enclosures on the grid 10^-30, factoring out the common positive C=(4pi)^-2. The code checks k1>0, k2<0, k1_ref>0, k2_ref<0, b<0 for L=3,5 and b>0 for L=30,100. For L=10 the matching vanishes exactly.

For b<0,

F'(q)=b-integral dmu/(q+t)^2<0 for q>-4.

Since F(0)=k1_ref>0 and F(q) tends to negative infinity at positive infinity, exactly one real zero exists, it is positive, and it is simple. No negative real zero is possible. Exact kernel endpoint signs isolate the L=3 and L=5 roots inside the saved brackets of radius 1e-8 about their old numerical locations.

For b>0, F''(q)=2 integral dmu/(q+t)^3>0 on q>-4, while F'(0)=k2_ref<0. Thus F(q)>F(0)>0 on -4<q<0. The same inequality holds in the finite threshold limit q=-4.

For q>=0, the positive measure and Cauchy-Schwarz give

S(q)>=k1^2/(k1-k2 q).

For example, set dnu=dmu/t and apply Cauchy-Schwarz to integral dnu/(1+q/t); its total mass is k1 and its inverse-t moment is -k2. The lower bound for F has positive denominator and numerator

(-b k2)q^2+(b k1-a k2)q+k1*k1_ref.

Exact interval arithmetic proves a positive leading coefficient and negative discriminant for L=30,100. The scaled discriminants, with their common C^4 removed, are strictly below -3.93e-5 and -2.12e-4, respectively. Rational endpoints, not rounded displayed decimals, decide the strict inequalities. Therefore neither prescription has a nonnegative real zero. This replaces the earlier floating-point discriminant sign checks.

## 4. Cut boundary and thresholds

The saved raw cut theorem is

Im K2_raw(-s-i0)=-C*pi*s^2/240 * sum_i c_i(1-4z_i/s)_+^(5/2)<0, s>4.

Its strict sign follows from the third divided difference, including the absolutely continuous extension at each species threshold. The real polynomial a q+b q^2 changes no cut imaginary part. Thus no matched kernel can vanish on either lip anywhere with s>4, including regulator thresholds. At s=4 the imaginary part vanishes, but the real threshold value of F is positive: monotonicity for b<=0 or the convexity argument above gives F(-4)>=F(0)>0. The raw S(-4) limit is finite because the first density vanishes as (s-4)^(5/2).

The cut is not an extra analytic real interval on the first sheet; adjacent-sheet resonances are a separate problem. No real roots on another sheet are included or excluded by this count.

## 5. Exact complex existence certificates

For L=30,100 the old numerical root is only a locator. The new evaluator encloses K2/C and its derivative using Fraction arithmetic, integer square roots, principal-log argument bounds and elementary log/atan series. It uses the existing independently audited rational interval helpers, not NumPy complex values, to decide existence.

For Im q>0 choose beta_z=sqrt(1+4z/q) with positive real and negative imaginary part. The ratio (beta_z+1)/(beta_z-1) lies in the upper half-plane; evaluate its principal logarithm. Then

I_z=log z-2+beta_z Log[(beta_z+1)/(beta_z-1)],
I'_z=1/q-2z Log[(beta_z+1)/(beta_z-1)]/(q^2 beta_z).

With U=z+q/4, the exact compact tensor expression and derivative are

K2_raw/C=(1/15) sum c[U^2 I_z-(z^2-zq/8)log z],
K2_raw'/C=(1/15) sum c[U I_z/2+U^2 I'_z+z log z/8].

Add (a/C)q+(b/C)q^2 and the corresponding derivative explicitly. Midpoint Newton steps refine the rational center but do not themselves prove a zero.

For a disk of radius r, choose rational bounds |q|<=Q and Im q>=y>0. Parameter integrals give

|I'|<=1/y, |I''|<=1/y^2,
z*y/Q <= |z+q*x(1-x)| <= z+Q/4.

The modulus bounds plus pi<4 enclose |I|. Since

(U^2 I)''=I/8+U I'+U^2 I'',

these estimates supply an explicit rational upper bound M2 for |(K2/C)''| over the disk, including 2|b/C|. At each center c the script checks exactly

|K2(c)/C|+M2*r^2/2 < |K2'(c)/C|*r.

Component absolute bounds overestimate the left side; one separated derivative component underestimates the right side. Rouché's theorem gives exactly one zero counting multiplicity inside each circle, so that zero is simple. The disks have radius 1e-9 and remain strictly in the upper half-plane.

| L | Error upper bound | Strict Rouché margin |
|---|---:|---:|
| 30 | below 4.30e-19 | above 2.37e-10 |
| 100 | below 1.80e-18 | above 2.85e-10 |

Saved rational centers and intervals, rather than every printed decimal digit, specify the localization claim. Combining certified existence with the global at-most-one theorem proves exactly one conjugate pair and no other off-axis zeros for each b>0 case.

The standard complex-analysis input is Rouché's theorem, reference NIST DLMF section 1.10(iv), https://dlmf.nist.gov/1.10 . This reference was inspected in prior checkpoints; the current retrieval returned a service error. The inequality and analytic assumptions used here are written explicitly. No inaccessible physics result is imported.

## 6. Uniform exterior check

The prior exact UV expansion also supplies a bound |K2_raw/C-k_infinity/C|<=D(R), uniformly on the exterior first sheet, with k_infinity=-C*T2/8>0. For non-reference matching, the exact inequality

|b/C| > |a/C|/R+(|k_infinity/C|+D(R))/R^2

proves quadratic dominance for every |q|>=R. Five certificates are saved, using R=4096,6656,13056,38656,128256 for L=3,5,10,30,100. The unchanged reference uses constant dominance instead. These sufficient radii are not optimized. The global zero exclusion above does not rely on tracing a sampled outer circle; this separate bound makes the UV behavior explicit.

## 7. Execution, result and continuation

Run `python experiments/spin2_global_count_2026_10_30/check.py`.

Execution completed with all assertions passing. results.json retains exact log-derived coefficient intervals, the two negative-discriminant certificates, two signed real brackets, two complex disk certificates and five exterior certificates. Root-refinement traces are labeled as locators. No broad scan, previous loop Ward test, contour mesh or successful parameter quadrature was repeated. No mathematical assertion failed and no tolerance was relaxed. The proof relies on the stated kernel and measure, not on passing scripts or AI agreement.

The raw spin-2 counts remain zero additional zeros for all five spectra. The matched spin-2 counts are now complete at 1,1,0,2,2; the existing matched scalar counts are 3,3,1,4,2. These functions now have complete first-sheet coefficient inventories within their stated prescriptions. Source overlap, tensor multiplicity, nonlinear interactions and physical interpretation are different questions.

With these inventories closed, avoid further routine pole scans of the same kernels. Next use established causal/analytic-continuation results to audit what follows from the fixed Euclidean response functions without choosing a microscopic metric integration contour. Keep a formal frequency continuation, a retarded Green function, and an actual constrained quantum source functional distinct. No new Lorentzian prescription, matter coupling or regulator is selected by this mathematical inventory.
