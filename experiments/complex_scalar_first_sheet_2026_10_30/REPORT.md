# Bounded scalar first-sheet and kernel dispersion audit
Research batch: 2026-10-29 (US Central); directory date is UTC. Baseline: dae0591e60a931c000a41328c40dc962ea3c077a.

## Scope
Free minimally coupled leading-N auxiliary-metric PV determinant, stationary contacts, no bare curvature terms. Scalar A only. No constrained-CEV equivalence, Lorentzian stability verdict, spin-2 pole search, matched-regulator prediction, or unification claim.

## Branch and audited formula
q is Euclidean momentum squared. The first sheet is the q plane cut along (-infinity,-4]. Each species has its cut at -4 z_i. Continue from positive real q using principal log in I_i(q)=integral_0^1 log(z_i+q x(1-x)) dx.
The evaluated closed form is I_i=log(z_i)-2+2 atanh(r_i)/r_i, r_i=sqrt(q/(q+4z_i)), with principal sqrt and atanh. Values at q=0 are taken by limits, never by that singular-looking expression.
A=C/4 sum c_i[(z_i-q/2)^2 I_i-(z_i-q/2)z_i log z_i].
All five spectra are unchanged from baseline; coefficients have not been matched across regulators.

Independent principal-log Gauss-Legendre integration at 25 moderately off-axis points checks this expression. At 1024 nodes, maximum relative error is 4.58e-12; 512-to-1024 change is at most 2.00e-9. Schwarz reflection is also checked. These are numerical checks, not a formal proof of branch correctness on every point.

## Bounded contour results
Use B(q)=A(q)/q to deflate the stationary massless zero. Counterclockwise rectangles have real range [-R,R] and imaginary range [h,R]. For each L=3,5,10,30,100, take R=4L,20L,100L, h=0.05,0.005, and initial 512,1024,2048 points per edge: 90 contours.

Every sampled winding count is zero and is stable across these refinements. Midpoints and phase variation trigger adaptive edge subdivision. The maximum final phase increment is 0.232971 radians. The initial fixed grid had increments near pi; those unresolved counts were not accepted as evidence. The script now asserts phase resolution and integer consistency.

Since B is analytic inside these upper-half-plane rectangles and has no poles there, the argument principle associates a resolved winding with its zero count. The numerical evidence finds no zeros inside these tested rectangles. Reflection gives the corresponding lower-half-plane evidence.

Limitations: no rigorous interval bounds between samples; h<0.005, |Re q|>100L, Im q>100L, the real axis, and second sheets are not excluded. Existing real subthreshold zeros are outside these rectangles and remain untouched. This does not establish global pole absence or a healthy theory.

## Two-subtracted kernel dispersion cross-check
The kernel, rather than its inverse, is reconstructed first:
A(q)=a1 q + (q^2/pi) integral_4^infinity Im A(-s-i0)/[s^2(s+q)] ds,
a1=-C sum c_i z_i log(z_i)/8.
The discontinuity is -C pi/4 sum c_i(z_i+s/2)^2 sqrt(1-4z_i/s), with each term zero below its threshold.
All four thresholds are integration breakpoints. A(0)=0 and a1 are subtraction data from the baseline kernel. This cross-check does not independently derive a1.

At q=0.5,2,10 for all five spectra (15 reconstructions), maximum relative difference from the exact kernel is 3.19e-12. No isolated-pole terms belong in this dispersion relation for A itself. A dispersion relation for 1/A, including its massless and isolated poles, remains to be done. Its scalar massless contribution cannot alone be assigned an independent positive-particle interpretation.

## Reproduction and next work
Run python experiments/complex_scalar_first_sheet_2026_10_30/check.py with NumPy and SciPy. results.json contains versions, all 25 quadrature audits, 90 contour rows, and 15 reconstructions. Assertions ran successfully on Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0.

Next: independently count the known real subthreshold zeros with contours that avoid the physical cut; reconstruct the inverse scalar response with isolated poles, massless subtraction, regulator cuts, and large-circle control; extend the near-axis search with explicit error control. Spin-2 complex roots remain a separate missing test. Match independent curvature-squared and Einstein coefficients before comparing regulator prescriptions physically. No new model choice is needed for these calculations.
