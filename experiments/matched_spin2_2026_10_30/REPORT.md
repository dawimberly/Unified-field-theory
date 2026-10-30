# Explicitly matched spin-2 kernel audit

Baseline: d153f018725678865a31c4d1bf48664c6809e4a1. Research date: 2026-10-30 US Central.

## Scope and outcome

Continue the previous reference-L=10 bookkeeping comparison, without adopting the added polynomials as a microscopic physical action. Write S_L(q)=K2_L(q)/q and

F_L(q)=K2_matched(q)/q=S_L(q)+a_L+b_L q,

a_L=k1_ref-k1_L, b_L=k2_ref-k2_L.

The scalar matching is still the explicit Delta A in the previous checkpoint; neither tensor coefficient has been silently removed. The matched spin-2 inverse has additional poles on the open principal sheet:

| L | Extra zero q of K2_matched |
|---|---|
| 3 | 60.3330308218942 |
| 5 | 82.4469390832120 |
| 10 | none, by the existing unmodified Stieltjes theorem |
| 30 | 8.72600156453 +/- 78.0620223980 i |
| 100 | 10.0532001673 +/- 66.1017934580 i |

The zeros at positive Euclidean q would be tachyonic poles if this particular resummed kernel were adopted literally with an appropriate Lorentzian interpretation. Complex first-sheet poles likewise require attention for that modified action. Neither conclusion transfers to the constrained CEV theory or to a controlled low-energy effective expansion. No model change, healthy theory, discovery, or UV completion is established.

## Branches and independent execution evidence

The principal parameter logarithm defines the q-plane cut (-infinity,-4]. The compact atanh/sqrt kernel is reused from the previously checked spin-2 determinant. The new roots all have positive real part and therefore avoid cut-lip ambiguity. Their complex conjugates satisfy conjugation exactly at stored precision.

For every root, an independent Gauss-Legendre parameter integral of (z+qx(1-x))^2 log(z+qx(1-x)) at 256,512,1024 nodes is evaluated before adding exactly a q+b q^2. This does not reuse the compact atanh expression. Maximum absolute residual of the independently evaluated matched K2 is 4.56e-14; the compact F residual is at most 4.63e-18. Finite-difference derivative steps 1e-3,1e-4,1e-5 support simple poles and stable residues (individual data retained). A spin-2 channel residue alone is not a physical source-amplitude positivity proof.

Thirty-eight adaptively refined argument-principle contours pass: eight isolating squares (two initial meshes per root) each count one zero; thirty upper-half-plane rectangles, R=200,1000,10000 and initial 256/512 points per edge with bottom Im q=0.001, count one upper zero for L=30,100 and none for L=3,5,10. Each mesh refines midpoint phase increments and variation in magnitude; maximum accepted phase increment is below 0.21. The real roots are intentionally outside those upper rectangles and counted by their own squares. This supplies convergence evidence, not an interval-certified contour proof. The analytic restriction below provides control beyond bounded scans.

Run `python experiments/matched_spin2_2026_10_30/check.py` from the repository root. It writes results.json and asserts independent quadrature agreement, contour counts, and the real-axis bound signs. The run completed with all assertions passing; versions are recorded. No previous full tensor scan was repeated.

## Analytic restrictions from the already established positive measure

The previous no-bare-curvature proof supplies

S(q)=integral dmu(t)/(q+t), dmu positive, support t>=4,

with nontrivial continuous support. Its finite inverse moments yield S(0)=k1>0 and S'(0)=k2<0. The following deductions apply to F=S+a+bq, not to the scalar response.

For Im q=y>0,

Im F(q)/y=b-integral dmu(t)/|q+t|^2.

If b<=0 there are no upper or lower off-axis zeros. For L=3,5, b<0 and F'(q)=b-integral dmu/(q+t)^2<0 on q>-4. Since F(0)=k1_ref>0 and F(q) tends to negative infinity as q increases, there is exactly one positive real zero and none on (-4,0]. Thus the two computed real zeros exhaust the additional open-sheet zeros for these prescriptions. The cut is a boundary, not part of the open sheet.

If b>0, two distinct upper zeros z,w would imply

integral dmu(t)/[(z+t)(conjugate(w)+t)]=b,

as follows by subtracting F(z)-conjugate(F(w)). Diagonal entries also equal b. Equality in Cauchy-Schwarz forces (z+t)^-1=(w+t)^-1 almost everywhere, hence z=w. There is at most one distinct upper zero. A multiple zero is also excluded: at a zero, |integral dmu/(z+t)^2| is strictly less than integral dmu/|z+t|^2=b because its phase varies across continuous support, so F'(z) is nonzero. The conjugate half-plane follows by reflection. This argument is a direct application of the positive-measure identity, not an assumption that the modified F remains Stieltjes.

For L=30,100, no real zero occurs on (-4,0]: F''(q)>0 and F'(0)=k2_ref<0 imply F(q)>F(0) there. For q>=0, Cauchy-Schwarz gives

S(q)>=k1^2/(k1-k2 q).

For completeness take dnu=dmu/t, so S(q)=integral dnu/(1+q/t); Jensen or Cauchy-Schwarz gives this bound using integral dnu=k1 and integral dnu/t=-k2.

The lower bound for F has positive denominator and numerator

(-b k2) q^2+(b k1-a k2) q+k1 k1_ref.

Its discriminants are -6.3345540e-14 and -3.4245908e-13, respectively, with positive leading and constant coefficients. These numerical coefficient signs, also accompanied by direct positive minima in results.json, verify the real-axis exclusion for the tested parameter values. They are floating-point sign evaluations, not rigorous interval arithmetic. Together with the positive-measure upper-zero restriction and the numerically isolated roots, this strongly supports a single conjugate pair for each tested b>0 prescription; existence is numerically checked rather than interval certified. No analogous theorem for the scalar channel is asserted.

## Why these poles are not low-energy predictions

All extra roots have |q|/4 between 15.08 and 20.62. The first physical threshold fixes the derivative expansion radius at 4; these points lie far outside it. At the roots the ratio |k2_ref q/k1_ref| is 1.38 to 1.88, also outside a hierarchy where the quadratic correction is small compared with the Einstein term. The quadratic-only kernel would predict a common positive root q=-k1_ref/k2_ref about 43.824; the full modified determinants instead produce the different real/complex values above. This demonstrates the sensitivity of nonperturbative inversion to the higher coefficients left unmatched. It does not justify discarding exact poles of a separately adopted full action.

Primary literature check: Jonathan Z. Simon, Physical Review D 41 (1990) 3720, DOI 10.1103/PhysRevD.41.3720, https://journals.aps.org/prd/abstract/10.1103/PhysRevD.41.3720. The publisher abstract was inspected; full text is subscription gated. It discusses perturbative constraints for higher-derivative effective/nonlocal theories. We use that as limited methodological context only; no detailed theorem from inaccessible text is imported. The scale and expansion diagnostics above are calculated here, not attributed to Simon.

## Continuation

The no-bare spin-2 theorem remains correct under its original assumptions. Matching a finite set of low-energy coefficients changes the full resummed functions and cannot retain that theorem automatically. Next prioritize the pending background-independent composite/source definition (arXiv:2011.06541), then formulate a controlled low-energy response comparison with explicit truncation/error assumptions. Further global scalar analysis remains open. Preserve all previous checkpoints and do not choose a microscopic contour or additional curvature action without distinguishing that physical assumption.
