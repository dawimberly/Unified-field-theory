# Independent threshold continuation and matched-regulator comparison
Research batch: 2026-10-30 US Central.
Repository baseline: d148b9142f57586084f8e7761078f13fc9865335.

## Scope and main findings
The provisional scalar_threshold_transition.zip was read before this independent calculation. Its local adjacent-sheet sign and threshold scaling are supported. An independent Decimal calculation places the unmatched threshold root at
L_c=14.5572064556843705...,
agreeing with the provisional 14.557206455684199 within 1.73e-13.

The comparison changes after explicit low-energy matching. With the Einstein and both quadratic-curvature tensor data matched to the existing bookkeeping reference L=10, every tested spectrum L=3,5,10,30,100 has a first-sheet subthreshold scalar pole. The raw threshold crossing near L_c therefore cannot be treated as a scheme-independent stability improvement of one fully matched UV theory.

This remains a fixed free leading-N auxiliary-metric determinant calculation. No physical UV completion, constrained-CEV equivalence, global pole inventory, healthy vacuum or microscopic Lorentzian contour is established. The local pole event occurs near timelike s=4, not at a spacelike q>0.

## Provisional input and provenance
Input Library identifier: libfile_fa52b432f24c8191ac948aabcca17d35.
Filename: scalar_threshold_transition.zip.
SHA256: f06f5ccd590fcc18d930552f9ad5c8be0e65197bf795fc6c8ddeb99704416615.
The original provisional report/results are preserved under provisional_input/. The original archive remains the input source. The new script does not import or execute its threshold/root solver; it uses a new uniformizing coordinate and Decimal threshold evaluation with independently constructed barycentric weights. It reuses only the repository's already-audited scalar kernel for overlap comparisons.

The terminology bound versus virtual poles is standard in scattering theory; the primary paper Matuschek, Baru, Guo and Hanhart, European Physical Journal A 57 (2021), 101, https://link.springer.com/article/10.1140/epja/s10050-021-00413-y distinguishes physical-sheet subthreshold poles from real subthreshold unphysical-sheet poles. We use only that location terminology, with a response-pole qualification. No hadronic compositeness criterion or physical positive-norm S-matrix is assumed for this PV model.

## 1. Local uniformizer fixes the branch sign independently
Write q=-s and introduce
u=sqrt(4/s-1), s(u)=4/(1+u^2).
For real subthreshold s, the physical sheet has u>0 and the sheet adjacent across the physical constituent cut has u<0. In this local chart the physical-species parameter logarithm is the analytic function
J(u)=-2+pi u-2u atan(u).
It follows from J=-2+2u atan(1/u) on the positive real u axis and atan(1/u)=pi/2-atan(u) there. Analytic continuation then fixes the nearby branches without imposing a discontinuity sign by hand.

For s>4 below every PV pair threshold, let beta=sqrt(1-4/s). The lower q lip, q=-s-i0, corresponds to u=-i beta, yielding
J=-2+2 beta atanh(beta)-i pi beta.
This agrees with a separate parameter-space argument: the principal log is negative-imaginary on the interval where 1-s x(1-x)<0; that interval has length beta. Hence Im J=-pi beta.
The upper q lip is its conjugate. Continuing the lower physical lip through the cut keeps u near -i beta on the adjacent upper-lip chart, instead of replacing it with the physical upper-lip +i beta.

Below threshold,
J(-u)-J(u)=-2pi u.
Only the physical species changes in this local continuation, so
A_II(s,L)=A_I(s,L)-(C pi/2)(1+s/2)^2 sqrt(4/s-1),
C=(4pi)^(-2).
This independently reproduces the provisional formula's minus sign. No regulator cut is crossed in this test.

## 2. New branch and quadrature checks
Eight complex overlap points compare the uniformizer calculation with the repository principal-sheet kernel around q=-3.9 +/-0.01i and -4.1 +/-0.01i, for L=3 and L_c. Maximum recorded difference is zero at the stored precision. These evaluations share the regulator log formula, so the independent parameter quadrature is also needed.

Six cut cases (s=4.01,4.3,6 at both L values) verify the analytic scalar discontinuity and conjugate lips; maximum discontinuity error is 8.68e-19. Six below-threshold u-sign flips independently check the sheet-change formula.

Twelve finite-imaginary-part parameter integrations directly use principal complex log at eta=1e-4,1e-6. The initial generic adaptive integration missed narrow log layers and failed the 1e-8 comparison tolerance, with a worst error about 2.45e-7 despite optimistic quadrature estimates. Those results are preserved in initial_unsplit_quadrature_results.json.
Adding split points at each physical root and around its local complex-log width eta*z/(s^2 beta) resolves the layers. The corrected maximum difference from the uniformizer kernel is 1.93e-15. The actual eta-dependent approach to the boundary is retained separately from the integration error; it must not be mistaken for a branch-sign mismatch.

## 3. Threshold root and scaling
The Decimal computation uses 70-digit working precision, native log/sqrt and an alternating arctangent series for the regulator terms at s=4. It brackets L in [14,15] and performs 190 bisections. The long decimal output in results.json is working-precision evidence, not an interval-certified assertion of every displayed digit or proof of global uniqueness over all L.

Near u=0,
J=-2+pi u-2u^2+...,
On the physical u>0 side, sqrt(delta)=2u+O(u^3), delta=4-s,
A=A(4,L)+k sqrt(delta)+O(delta), k=9 C pi/8.
The neighboring root therefore obeys u_p approximately -A(4,L)/(2k) and delta approximately [A(4,L)/k]^2.
For the existing identical-constituent source response W=(2+s)^2/(16A), the physical-sheet Euclidean residue obeys
Res W approximately 4.5 sqrt(delta)/k.
It vanishes at the threshold. On the adjacent real sheet the corresponding q-residue has the opposite sign because dq/du changes sign; this is not a physical-sheet ghost diagnosis.

Eight signed offsets (L-L_c)/L_c=+/-0.03,+/-0.003,+/-0.0003,+/-0.00003 were solved directly in u. Negative offsets give u_p>0 on sheet I; positive offsets give u_p<0 on sheet II. At the closest offsets, the gap/prediction ratios are 1.0000139545 and 0.9999860459. The signed source-residue ratios approach +1 and -1 respectively. Complex-step derivatives at three steps check residue convergence. No finite-difference derivative is taken at the branch point.

This is local response-pole tracking. It supplies neither an exhaustive second-sheet search nor a proof that no other off-axis poles exist.

## 4. Explicit matching and what it changes
Retain the previous reference convention L_ref=10. Define a1=-G/8, a2=C sum c_i log z_i/16 for A, and k1=G/24, k2=C sum c_i log z_i/240 for K2.
For each spectrum add the explicit bookkeeping polynomials
Delta A=(a1_ref-a1_L)q+(a2_ref-a2_L)q^2,
Delta K2=(k1_ref-k1_L)q+(k2_ref-k2_L)q^2.
Both tensor sectors' quadratic coefficients are specified, as is the shared Einstein slope. The numerical increments are recorded in matched_kernel_map. These are not silently installed as a microscopic bare action; the fixed-prescription spin-2 theorem does not automatically apply to the modified kernels.

The scalar poles below were found in the same real-u bracket for the explicitly matched A. No spin-2 pole search of the modified K2 is claimed.

| L | Matched first-sheet pole s | Matched A(s=4) |
|---|---:|---:|
| 3 | 3.8559544985 | -0.00774076952 |
| 5 | 3.7957424507 | -0.00923311578 |
| 10 | 3.7473019640 | -0.01031031376 |
| L_c | 3.7316391515 | -0.01064237093 |
| 30 | 3.7136858056 | -0.01101511971 |
| 100 | 3.7016942467 | -0.01125990297 |

In particular, the raw A(4,L_c)=0 becomes negative after matching, and the unmatched absence of a physical subthreshold root at L=30,100 does not persist under this comparison. This shows the relevance of specifying the low-energy action before interpreting regulator variation. It is not a proof that any particular matched pole is physical or that the full functions represent the same UV completion: higher-order coefficients still differ. The previously proven scalar negative continuum sign is also not repaired by real polynomial matching.

## 5. Reproduction, evidence and next work
Run python experiments/threshold_independent_2026_10_30/check.py.
Dependencies: NumPy, SciPy and the earlier portable repository scalar-kernel script. The corrected final command passed all assertions on Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0. No provisional helper script was executed. No older full spectrum scans were repeated.

Code, final results, retained failed quadrature, provisional input report/results and provenance are saved. Next: verify how the explicitly matched spin-2 kernel's zeros change, since its no-bare-curvature positivity theorem has changed hypotheses; keep that a kernel/EFT comparison rather than a UV-health verdict. Independently continue the pending background-independent composite/source-matching literature work before transferring any auxiliary spectrum to constrained CEV. No new physical model assumption was adopted.
