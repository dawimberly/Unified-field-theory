# Scalar divided-difference identity and low-energy matching
Research batch: 2026-10-29 US Central. Directory date is UTC.
Baseline: d824c4edc05b4b2e3a1d6865635b7bbe191772df.
This is a bounded algebraic/matching batch, not a new global scalar spectrum claim.

## Continuity and literature
The current repository status and preceding spin-2 proof were read before selecting work. Three newer Library notes were also read in full and imported, unchanged in substance, under bots/gpt/research_notes/2026-10-29/. Their bibliography and claimed past inspections are provenance from those notes, not a new verification of all cited papers in this batch.

Their conclusions qualify all the determinant calculations: no spectrum transfer to constrained CEV has been established; a determinant-weighted delta constraint can represent a specified reduced measure but does not select that measure; stochastic composite gravity is an existing comparison proposal, not an adopted regulator for our model. No unconstrained-to-constrained equivalence, stochastic quantization assumption, or physical regulator is introduced here.

The applied mathematical result is the same established Genocchi-Hermite identity used for spin-2: Carl de Boor, Divided Differences, Surveys in Approximation Theory 1 (2005), 46-69, equation (52). Primary author-hosted manuscript checked: https://ftp.cs.wisc.edu/Approx/dvdsurvey.pdf . It converts a third divided difference into an integral of the third derivative over a simplex. The following formula is our application of that identity, not a claim that the reference solved the gravitational model.

## Scalar identity: valid, but no positivity shortcut
Use the unchanged scalar kernel
A(q)=C/4 sum_i c_i {integral_0^1 (z_i-q/2)^2 log(z_i+tq) dx -(z_i-q/2)z_i log(z_i)},
t=x(1-x), C=(4pi)^(-2).
Let w_i=1/product_(j!=i)(z_i-z_j), c_i=w_i/w_0, w_0<0, Z=sum lambda_i z_i. The simplex has volume 1/6.

At fixed q,x define
f(z)=(z-q/2)^2 log(z+tq) -(z-q/2)z log z.
Writing y=Z+tq and d=q(t+1/2), its third derivative is
f'''(Z)=2/y+2d/y^2+2d^2/y^3-2/Z-q/(2Z^2).
Therefore
A(q)=C/(4w_0) integral_0^1 dx integral_simplex f'''(Z) d lambda.

The branch is the principal-log first sheet continued from positive q. For q outside (-infinity,-4z_0], every z+tq on the relevant convex mass interval avoids the logarithmic singularity/cut. The simplex identity first holds on the positive axis and extends analytically to that sheet.

Unlike the spin-2 formula, this integrand has no pointwise fixed sign: for Z=1,t=1/4 it is positive at small positive q but negative at q=100. The identity does not justify classifying A/q as a positive Stieltjes transform, importing the spin-2 zero theorem, or declaring that all scalar complex poles are absent. Prior scalar contour bounds and known real poles remain the available evidence. This failed shortcut is retained rather than hidden.

## Match the three low-energy tensor data explicitly
Define G_L=C sum c_i z_i log z_i and S_L=sum c_i log z_i. In the existing conventions:
A_L(q)=-G_L q/8 + C S_L q^2/16 + a3_L q^3+O(q^4),
K2_L(q)=G_L q/24 + C S_L q^2/240 + k3_L q^3+O(q^4).

The single Einstein coefficient fixes both q slopes. The independent q^2 coefficients in the spin-0 and spin-2 projector sectors are both specified; the ratio happens to be locked to 15 for this particular minimal-scalar regulator family. That ratio is not a general condition on arbitrary gravitational EFTs.

For a bookkeeping comparison choose L_ref=10, without choosing it as a physical regulator or adding terms to the adopted microscopic theory:
A_L^match=A_L -(G_ref-G_L)q/8 + C(S_ref-S_L)q^2/16,
K2_L^match=K2_L +(G_ref-G_L)q/24 + C(S_ref-S_L)q^2/240.
These are formal kernel counterterm maps matching the Einstein slope and both independent curvature-squared projector data to the reference. They are not an argument that the remaining full nonlocal functions define the same UV theory. No numerical pole search of these modified kernels is performed or interpreted physically.

The first unmatched terms are
a3_L=(17C/1120) sum c_i/z_i,
k3_L=(C/3360) sum c_i/z_i.
Consequently the cubic differences obey Delta a3=51 Delta k3 for this family. Matching through q^2 leaves higher derivative Wilson coefficients unequal. Those residuals cannot be called regulator-independent predictions without additional matching data or a specified UV theory. The no-bare-curvature spin-2 proof must not be silently transferred to the formally matched polynomial-modified kernels.

## Reproducible checks, restricted to the new gaps
Run python experiments/scalar_identity_matching_2026_10_30/check.py.
It imports definitions from the two preceding portable experiment scripts; no older batch's full test suite or scans are rerun.

1. Six off-axis scalar identities compare simplex integration with the exact kernel at L=3,10,100 and q=1+2i,-2+0.3i. Initial 24/48-node integration failed the proposed 1e-5 relative tolerance at L=100,q=-2+0.3i: error 1.43584e-5. This is retained in initial_48_node_results.json. Both L=100 cases were then refined to 96 nodes per integration dimension; the maximum final relative error is 3.806e-8. Final figures are in results.json.
2. Thirty formal matching comparisons cover all five spectra, both tensor coefficients and q=0.05,0.1,0.2. Exact-kernel matched differences agree with Taylor differences through q^8 to maximum absolute error 5.81e-17. These low momenta are inside the analytic radius set by the physical threshold.
3. The new q^3 through q^8 coefficients are derived using exact rational beta integrals, not fitted to the comparison values. With m_n=integral_0^1 t^n dx=(n!)^2/(2n+1)!,
a_n/C=(-1)^(n+1)[m_n/n+m_(n-1)/(n-1)+m_(n-2)/(4(n-2))] sum c_i z_i^(2-n)/4,
k_n/C=(-1)^(n-3) m_n sum c_i z_i^(2-n)/[4n(n-1)(n-2)], n>=3.
The code verifies the displayed cubic rational factors independently of numerical fits.

The final command completed successfully with all assertions passing on Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0. All figures are consistency evidence for the formulas, not a proof of physical equivalence. Code, final results, retained failed-resolution results and the imported-note provenance are preserved.

## Next work
The established spin-2 theorem remains restricted to its fixed bare prescription. Scalar global zero control is still unresolved; investigate applicable signed-measure/generalized-Nevanlinna results rather than extend a positive-measure theorem beyond its hypotheses. Separately prioritize the literature notes' matched-observable/reduced-measure question before claiming a scalar pathology or cure in constrained CEV. A stochastic comparison can be investigated without adopting its microscopic axioms. No foundational user choice is required for those investigations.
