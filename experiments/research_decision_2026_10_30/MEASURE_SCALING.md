# Conditional large-N sensitivity of measure and source matching

Date: 2026-10-30 US Central. Baseline: 2345a4ecc8eb6e7b96b8b1018b54943fbcfb1cdd.

This is a new analytic bookkeeping derivation, not a selected quantum measure, a numerical experiment, or a completed 1/N calculation. It identifies when a proposed change could affect the leading-N conclusions and when it could not. Every hypothesis below matters.

## Hypotheses

1. Start at a fixed finite regulator, with a fixed number m of algebraic auxiliary components per site, independent of N.
2. Work on one nonsingular branch with an isolated saddle and invertible auxiliary Hessian after any required consistent gauge restriction. The metric trace direction may require a complex steepest-descent cycle; this derivation does not choose one.
3. The leading collective action and the definition of the observable are the same in the compared formulations. Their regulator, subtraction and background tunings are matched at order N.
4. Additional measure weights have log density of order one at fixed collective fields, with regular derivatives, rather than an order-N exponential or a singular branch restriction.
5. The relevant large-N limit and stationary expansion are uniform on the region where they are used. Continuum/cutoff limits are not interchanged without an estimate.

These hypotheses have NOT been established for the constrained CEV functional. The result must not be used to assert they hold there.

## 1. Ordinary auxiliary saddle factors are subleading in the action

At fixed matter configuration or fixed collective fields, write a formal local saddle expansion

S_N(g,phi)=N S_red(phi)+N h:H0(phi):h/2+O(N h^3),
g=G(phi)+h.

A regular auxiliary saddle integration gives, up to the phase of its as-yet unspecified cycle,

integral dg mu_g(g,phi) exp[-S_N]
~ exp[-N S_red] mu_g(G,phi) (2pi/N)^(m/2) |det H0|^(-1/2) [1+O(1/N)].

The field-dependent log measure contribution is

Gamma1_measure=-log mu_g(G,phi)+(1/2)log|det H0(phi)|.

The term (m/2)log N is field-independent at fixed regulator and N. It does not change a source derivative. The log determinant is order one under the stated hypotheses, whereas N S_red is order N. This is a formal saddle-power-counting statement, not convergence of a real integral with an indefinite metric Hessian.

In D=4 there are m=10 symmetric metric components per site, so the explicit Gaussian N power is N^-5 per site. A configuration-space constraint T_N=N T0 has det(delta T_N/delta g)=N^m det J0. A plain delta[T_N] produces the inverse-Jacobian factor; a determinant-weighted normalized delta cancels it. These definitions are not the same, even though their ordinary finite-component field-dependent log factors are subleading relative to a common regular order-N action.

An N-independent power of det g for the fixed auxiliary components has the same order-one classification. A matter-measure factor repeated independently for N species can instead contribute order N, even when its per-species factor looks ordinary; it must be counted and regulated before applying this classification. This is not a license to guess its exponent. Nor does a configuration-space constraint determinant become a canonical second-class constraint determinant by this argument.

## 2. Source conversion is a different order-counting question

Include a small ordinary source j for g. At fixed phi, formal completion of the same local Gaussian square produces the factor

exp[j:G(phi)+j:H0(phi)^-1:j/(2N)]

inside the remaining matter integral, with the sign/phase appropriate to the chosen Gaussian convention. It cannot generally be pulled outside that integral: H0 depends on phi. At quadratic source order the auxiliary-to-composite covariance difference includes (1/N)<H0^-1> in the common zero-source measure. A regular collective saddle permits replacing that expectation by its leading background value for this order; otherwise its dependence must be retained.

The displayed Euclidean-sign algebra does not fix Lorentzian i factors. The source-dependent local contribution is order 1/N, the same order as the leading connected covariance of normalized collective metric fluctuations. It cannot be discarded merely because the source-independent log determinant is order one in an effective action of order N.

This explains why two separate statements coexist:

- ordinary smooth measure corrections are subleading in the collective action;
- the auxiliary-to-composite source conversion contact can be leading in the metric two-point function.

The previous quadratic source bridge and its fixed normalization remain conditional. This scaling argument does not reproduce a general source-dependent nonlinear integral or identify a microscopic generating functional. Macroscopic sources scaled as N require a different source expansion and are not covered by the small-source formula.

## 3. Conditional sensitivity of a simple pole

On a fixed physical/gauge-restricted channel suppose

K_N(q)=N K0(q)+K1(q)+O(1/N),

where K1 is regular near a simple zero q_star of K0 and all relevant coefficients remain bounded as N grows. The implicit-function expansion yields

q_star(N)=q_star-K1(q_star)/(N K0'(q_star))+O(1/N^2).

Away from zeros the inverse expands as

K_N^-1=(1/N)K0^-1-(1/N^2)K0^-1 K1 K0^-1+... .

The latter inverse expansion is not uniform at the zero; it is not used there as a proof. The root shift instead follows from expanding the root equation under its stated regularity assumptions. In a matrix problem the same formula uses the simple eigenchannel's left/right null vectors and its nonzero projected derivative. Multiple zeros and mixing require their own analysis.

Consequently an ordinary order-one regular measure correction cannot generically shift a fixed simple leading-N pole by order one or remove the leading-N branch of poles. It can shift it by order 1/N and alter residues. This is a necessary-order filter, not a proof that the poles of the present auxiliary model belong to the constrained theory.

Possible failures of the filter include: an order-N measure term; a change in leading action or regulator; singular Jacobian/saddle/observable maps; a different contributing contour or branch; a cutoff scaling with N; nonuniform large-N limits; and nonanalytic self-energy effects near the pole. None has been ruled out merely by the saved finite-dimensional algebra.

A finite regulator volume can have many sites, but its number is fixed in this counting. If the cutoff or number of regulator degrees of freedom grows with N, a nominal per-site order-one term may compete with the leading action. That joint limit needs an actual bound, not an interchange by notation.

## 4. Local metric corrections cannot repair the saved cut sign

For the already specified scalar coefficient, let A_new=A+P(q), where P is real on the real cut and has no discontinuity there. At a nonzero finite denominator,

Im[1/A_new(-s-i0)] = -Im A(-s-i0)/|A(-s-i0)+P(-s)|^2.

The sign of the specified inverse coefficient density is therefore unchanged. The fixed identical-constituent scalar source multiplies it by (s+2)^2/16>0, also preserving its sign. A real local polynomial may change isolated poles and the density magnitude, but it cannot turn a strictly negative continuum interval of this coefficient/source into a positive one. New singular cut points are not a repair of that statement.

This observation applies even beyond order-one corrections whenever only such real local terms change. It does NOT state that arbitrary nonlocal measures, different operator insertions, different matter content or a physical regulator leave the response invariant. It also is not a standalone full-gravity unitarity theorem: the observable, tensor structure and quantum source functional still matter.

## 5. Research consequence

'Changing the measure might fix it' is too imprecise to be a calculation. A proposed repair must identify whether it changes the leading action, only an ordinary subleading measure, the observable/source overlap, or the contour/branch, and must satisfy the assumptions appropriate to that claim. A local source contact alone changes neither a nonlocal pole nor the continuum discontinuity of the complete unchanged response.

The most discriminating next step is a common regulated leading-N source functional, not another scan of the same kernels. If its leading action and source map genuinely agree with the audited auxiliary approximation up to the stated local contacts and ordinary subleading measure changes, then the nonlocal leading response must be compared rather than dismissed. If they do not agree, the leading discrepancy itself must be calculated. No physical equivalence or repair has been assumed to make either outcome occur.

No new numerical checks were run for this derivation. Its status is a conditional analytic sensitivity argument awaiting application to a fully specified regulated source definition, not a verified property of a chosen constrained measure or a 1/N correction calculation.
