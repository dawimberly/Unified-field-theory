# Infinite-boundary scalar inverse dispersion and moment identities

Date: 2026-10-30 US Central. Baseline: ab6f85ff8e1d150e39c8f10782615df8b140213f.

This batch uses the complete scalar first-sheet inventories established by the preceding cut-boundary argument. It derives infinite-cut reconstructions and explicit UV tail bounds, then evaluates the reconstructions and the applicable inverse moment identities. The raw and explicitly matched cases have different UV falloffs and must be treated separately. No physical action, source, contour or regulator prescription was changed.

Ten exact tail inequalities, thirty once-subtracted numerical reconstructions, twenty-seven absolute reconstructions and sixteen matched inverse moment checks are saved. The mathematical conclusions concern the specified free leading-N auxiliary coefficient A and its inverse F=1/A. They establish neither a positive global spectral measure nor constrained-CEV equivalence, Lorentzian health or a UV completion.

## 1. Formula and branch audit

The first sheet is the q=Q_E^2 plane cut along (-infinity,-4], M=1. Define

rho(s)=Im F(-s-i0)/pi, s>=4.

All physical and PV species are retained, and the logarithms have the previously independently verified lower-lip sign. Matching adds a real polynomial on the cut, which leaves Im A unchanged but changes rho through |A|^2. Using the raw inverse density for a matched response would be wrong.

Let S_n=sum c_i z_i^n and T_n=sum c_i z_i^n log z_i, C=(4pi)^-2. The audited uniform UV expansion is

A_raw(q)=a_infinity+[b Log q+d]/q+O(Log q/q^2),
a_infinity=3C T_2/8<0, b=5C S_3/4<0,
d=C[-5T_3+S_3/6]/4.

The existing explicit Lref=10 matching gives

A_matched(q)=a_2 q^2+a_1 q+a_infinity+O(Log q/q),
a_1=a_1,low(ref)-a_1,low(L),
a_2=a_2,low(ref)-a_2,low(L).

Here a_1 and a_2 denote the added UV polynomial coefficients, not the matched total low-energy slope and curvature coefficient. The raw contributions vanish polynomially at infinity and instead approach a_infinity. Confusing these two definitions would give incorrect inverse moment targets. For L=10 both increments vanish exactly; this reference belongs to the raw asymptotic class.

It follows that

rho_raw(s)=alpha_raw/s+O(log s/s^2), alpha_raw=-b/a_infinity^2>0;
rho_matched(s)=alpha_matched/s^5+O(1/s^6), alpha_matched=-b/a_2^2>0

for the four non-reference matched cases. The positive eventual tail does not remove the known finite intervals of negative density. These rho functions are coefficient-inverse densities, not complete source-projected weights.

## 2. Infinite-boundary reconstruction

Write the complete pole inventory, including the simple massless pole, as (q_j,r_j), where r_j=1/A'(q_j). Complex conjugates are both included. The previous global count, rather than a small residual, supplies completeness. The formulas apply off the cut and away from the poles.

Cauchy's residue theorem gives the absolute relation

F(q)=F_infinity+sum_j r_j/(q-q_j)+integral_4^infinity rho(s)/(s+q) ds,

where F_infinity=1/a_infinity for raw/reference kernels, and zero for the four non-reference matched kernels. For raw kernels apply the vanishing-boundary argument to F-F_infinity=O(Log q/q), not to F itself. The continuum integral converges even though the raw zeroth unweighted density moment diverges.

Subtracting at q0=1 gives

F(q)-F(q0)=sum_j r_j[1/(q-q_j)-1/(q0-q_j)]
+integral_4^infinity rho(s)[1/(s+q)-1/(s+q0)]ds.

This has no residual finite outer boundary. The uniform UV asymptotics justify its removal; an outer term was retained in the preceding finite-radius experiment and is not simply forgotten here.

The established input is Cauchy's residue theorem, inspected directly in NIST DLMF 1.10.8, https://dlmf.nist.gov/1.10#E8 . The model-specific cut sign, convergence, constant term, residues and limiting estimates are derived explicitly here. No external quantum-field-theory positivity or superconvergence theorem is assumed.

## 3. Exact quantitative UV tails

Use R=256 z_max: R=4096,6656,13056,38656,128256 for L=3,5,10,30,100. Let the previous exact UV bound be

|A_raw(q)/C-a_infinity/C|<=D(R), for |q|>=R.

Let ell=-5 S_3/4>0, and let E(R) bound the subleading cut imaginary terms, so

|Im A(-s-i0)/(C pi)| <= B/s, B=ell+R E(R), s>=R.

Multiplication by s makes every subleading power decrease; the species-series remainder times s decreases too. Hence this bound covers the entire tail, not just its starting point. The saved preceding cut-sign bound also proves strict positivity of rho on this tail.

For raw/reference cases use the positive lower bound

m_raw=|a_infinity/C|-D(R), |A(-s-i0)|>=C m_raw.

For non-reference matched cases use

m_matched=|a_2/C|-|a_1/C|/R-(|a_infinity/C|+D(R))/R^2,
|A_matched(-s-i0)|>=C m_matched s^2.

All ten lower bounds are positive in exact Fraction arithmetic. The script rounds them downward and B upward onto the rational grid 10^-30, then uses an outward rational Machin enclosure of pi. Therefore

0<rho(s)<=K/s^p,
K=(16 pi_upper^2) B_upper/m_lower^2,
p=1 (raw/reference), p=5 (non-reference matched).

K is rounded upward. Saved rational endpoints, exponents, radii and code reproduce the certificates; floating-point evaluations do not decide these inequalities.

For |q|,|q0|<=h<R, the once-subtracted omitted continuum satisfies

|tail| <= K |q-q0|/[(p+1)(1-h/R)^2 R^(p+1)].

The absolute reconstruction tail satisfies

|tail| <= K/[p(1-h/R)R^p].

The code uses rational component-sum bounds on the complex moduli of the specified probes. For a non-reference matched density moment n=0,1,2,3,

0<integral_R^infinity s^n rho(s) ds <= K/[(4-n)R^(4-n)].

These are rigorous truncation bounds for the continuum only. They do not certify finite-interval quadrature errors or all digits of numerical pole residues. The numerical program also integrates each tail using x=R/s and checks that its value is consistent with the certified bound.

## 4. Which inverse moments exist

For non-reference matched kernels,

F(q)=1/(a_2 q^2)-a_1/(a_2^2 q^3)
+[a_1^2/a_2^3-a_infinity/a_2^2]/q^4+O(Log q/q^5).

Define M_n=sum_j r_j q_j^n+(-1)^n integral_4^infinity s^n rho(s)ds. The density tail makes these integrals absolutely convergent for n=0,1,2,3. Expansion of the convergent dispersion representation, or its finite expansion with a dominated remainder along positive q, gives

M_0=0,
M_1=1/a_2,
M_2=-a_1/a_2^2,
M_3=a_1^2/a_2^3-a_infinity/a_2^2.

The massless pole contributes to M_0 and must not be dropped. Both members of each complex pair contribute to all moments. These are signed continuum and pole identities, not positive-measure probability normalizations.

For n=4 the continuum moment diverges logarithmically because alpha_matched>0. The next inverse term contains Log q/q^5, so an ordinary finite M_4 must not be inferred from a formal power-series manipulation. For every raw/reference case even the zeroth density moment diverges logarithmically, despite the convergent absolute and once-subtracted dispersions. The four non-reference sum rules cannot be imposed on L=10.

## 5. Numerical execution and its limits

Run:

`python experiments/infinite_scalar_dispersion_2026_10_30/check.py`

The program reads the saved isolated locations, refines them using the unchanged kernel and analytic derivative, and checks the inventory length against the exact global count. It uses long-double complex arithmetic with a recorded 63-bit mantissa. Existing model helpers and numerical dependencies are imported from the repository; all paths are relative to the repository. Original scripts and historical files are unchanged.

Every finite continuum is split at all four thresholds and at the stable-UV switch. A sine-squared change of variable resolves endpoint square roots. The tail is integrated after x=R/s. Adaptive quadrature error estimates are saved separately from the exact tail certificates. No quadrature warnings occurred.

- Thirty once-subtracted reconstructions: five spectra, raw/matched, q=0.5,2,1+2i, q0=1. Maximum normalized closure error 6.02e-15.
- Twenty-seven absolute reconstructions: the same probes, with the duplicate matched-L=10 case omitted. The raw constant 1/a_infinity is fixed by the UV formula, not fitted. Maximum normalized error 2.44e-14. This independently checks the additive constant that subtraction cannot determine.
- Sixteen inverse moment checks: four non-reference matched cases and n=0,1,2,3. Maximum discrepancy divided by the sum of absolute pole, cut, tail and target terms is 4.29e-16. Largest absolute discrepancy is about 2.91e-8 in the large n=3 moment. These observed errors are not interval-certified digits.
- Thirty density-asymptotic ratios at R,100R,10000R show the expected leading falloffs. They are diagnostics; the exact inequalities establish tail control.

For the first two categories normalization is max(1,|direct|). Omitting the numerical tail from the once-subtracted calculation gives an error up to 1.64e-7 in the tested probes. The tail cannot generally be set to zero just because its upper bound is known. In matched n=3 moments, tails at these radii are approximately 1937,19383,316843,797421 for L=3,5,30,100; their certified bounds are approximately 1953,19530,318842,802186. Reporting high precision for those moments without accounting for the tail would be misleading.

There were no failed mathematical assertions in this batch. No acceptance tolerance was relaxed to obtain these results.

## 6. Source projection: implications and a remaining audit

The constituent scalar source factor saved in earlier checkpoints is

W(q)=(q-2)^2 F(q)/16.

It changes residues to r_j^W=(q_j-2)^2 r_j/16 and cut weight to rho_W=(s+2)^2 rho/16. One must not identify the coefficient inverse moment rules with source-exchange moment rules.

For a non-reference matched case, the identities above algebraically imply

W(q)=1/(16a_2)+sum_j r_j^W/(q-q_j)+integral_4^infinity rho_W(s)/(s+q)ds.

The constant follows by polynomial division: the linear contact coefficient is M_0/16=0 and the remaining constant is (M_1-4M_0)/16=1/(16a_2). This is a mathematical source-projection contact, not an identification with the separate composite/auxiliary contact dictionary or its microscopic measure.

For raw/reference kernels W grows with a known quadratic polynomial plus q Log q terms. Its source density grows as s, so a naive unsubtracted or once-subtracted continuum formula diverges. Subtract the known F_infinity(q-2)^2/16 polynomial first, then use two subtractions for the remainder; alternatively retain that polynomial explicitly in an appropriately subtracted formula. Do not silently apply the coefficient-inverse reconstruction to W.

Next audit these source-projected formulas and their subtraction/contact normalizations independently, including their relation to the previously saved stress-bubble/composite-source dictionary. The algebraic implications above are not a new microscopic source model. Constrained measure, Lorentzian contour, interactions, 1/N and physical universality remain separate unresolved issues.
