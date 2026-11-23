# Independent review of the FFZ three-loop Newton poles

Commit reviewed: `196c696` (FFZ Newton matching review packet).
Prescription: fixed \(\xi=1/6\), \(d=4-2\varepsilon\), equal masses, units \(h^2 z/K\), \(L=\ln(z/\mu^2)\).

Runnable checks in this directory (they do not call the checkpoint scripts):

- `ffz_independent_pole_review.py` — QCc and CcFf poles
- `ffz_ffff_A_poles.py` — FFff(A) soft-corner Taylor averages
- `ffz_ffff_A_regulated_conversion.py` — projective \(\to\) graph conversion and B targets
- `ffz_ffff_B_gaussian.py` — FFff(B) parity-even Gaussian algebra
- `ffz_ffff_B_chatgpt_audit.py` — independent FFff(B) poles vs parallel ChatGPT claims (`FFFF_B_CHATGPT_AUDIT.md`)
- `ffz_ffff_B_contraction_contamination.py` — Wick multiplicity and ε×triple→double shift
- `ffz_ffff_B_slope16_hunt.py` — locate the −16 slope term (vertex on double-log overlap)
- `ffz_ffff_B_overlap_fix.py` — apply overlap assembly fix; verify projected double cancels (`FFFF_B_OVERLAP_FIX.md`)
- `ffz_ffff_B_simple_ledger.py` — projected simple-pole ledger vs \(1/24\); vertex-ll \(\Delta c_0=64/3\) (`FFFF_B_SIMPLE_LEDGER.md`)
- `ffz_ffff_B_soft_linear_audit.py` — soft-linear of extra vertex-ll copy; excluded from overlap \(c_0\) (`FFFF_B_SOFT_LINEAR_AUDIT.md`)
- `ffz_ffff_B_overlap_rigor.py` — sector inclusion-exclusion: overlap = soft-soft only; coupling \(\equiv\lambda_F\) (`FFFF_B_OVERLAP_RIGOR.md`)
- `ffz_ffff_A_simple_reextract.py` — FFff(A) projective \(c_0\) with cut/seed error (`FFFF_A_SIMPLE_REEXTRACT.md`)
- `ffz_simple_gap_audit.py` — full simple-pole ledger + Y/CT/QFf re-derivation; gap stands (`FFFF_SIMPLE_GAP_AUDIT.md`)
- `ffz_gap_physical_audit.py` — is the gap physical? topology coverage; CcFf±1e-4; \(a_Y\) scan (`FFFF_GAP_PHYSICAL.md`)
- `ffz_ffff_B_nonoverlap_finite.py` — forest non-overlap \(\Delta c_0=0\) (not a second \(64/3\)); gap stands (`FFFF_B_NONOVERLAP_FINITE.md`)
- `ffz_finite_eps_crosscheck.py` — independent finite-ε Laurent cross-check of every ledger class (`FFFF_FINITE_EPS_CROSSCHECK.md`)
- `ffz_ffff_B_literal_coupling.py` — bare-Lagrangian coupling weight; literal global vs overlap-fixed poles (`FFFF_B_LITERAL_COUPLING.md`)
- `ffz_dred_yukawa_ledger.py` — DRED + global Yukawa weight + BPHZ/R; gap stands (`FFFF_DRED_YUKAWA.md`)
- `ffz_heat_kernel_a1.py` — Seeley–DeWitt / HK lower-loop inputs for \(1/24\) and CT; 3-loop HK not done (`FFFF_HEAT_KERNEL_A1.md`)
- `ffz_superspace_ward.py` — SUSY non-renormalization + partner Ward audit; gap stands (`FFFF_SUPERSPACE_WARD.md`)
- Ladder: `docs/PROVEN_METHODS_LADDER.md`

SymPy and mpmath own the exact identities. The CcFf simple pole and the FFff(A) regulated remainder are separate quasi-Monte Carlo estimates of integrands written for this review.

## Verdicts

| Claim | Verdict |
| --- | --- |
| QCc double pole \(1/(72\varepsilon^2)\) | Agree |
| QCc simple pole \(-0.01649742330-L/24\) | Agree (constant recomputed as \(-0.0164974232969\)) |
| CcFf double pole \(1/(144\varepsilon^2)\) | Agree |
| CcFf simple pole about \(-0.0247-L/48\) | Agree at the checkpoint's stated precision. This review's extrapolation is \(-0.02480\) |
| FFff(A) short-short triple pole | None (earlier \(1/7776\) retracted; see correction) |
| FFff(A) corner double pole \(5/(1728\varepsilon^2)\) | Agree with corrected cycle Wick check |
| FFff(A) long-pair triple add-back \(1/(576\varepsilon^3)\) | Agree |
| FFff(A) full graph higher poles | Numerical inference only (`FFFF_A_REGULATED.md`); conversion exact |
| FFff(B) triple \(1/48\), bare double \(-1/96\) | Agree (residue / regulated) |
| FFff(B) overlap-fixed double \(-1/288\) | Applied (vertex-ll \(\Delta a_1=-16\); see `FFFF_B_OVERLAP_FIX.md`) |
| FFff(B) literal coupling \(\Delta a_1=-16\), \(\Delta c_0=64/3\) | Soft-soft residues confirmed; global doubling refuted; overlap-restricted reading confirmed (`FFFF_B_LITERAL_COUPLING.md`) |
| DRED + global Yukawa + BPHZ | Higher poles cancel only after R (soft-soft remainder \(-16,64/3\)); simple gap \(\sim+0.01\) stands (`FFFF_DRED_YUKAWA.md`) |
| Heat-kernel / Seeley–DeWitt \(a_1\) | Lower-loop inputs for \(c_Y=1/24\) and CT verified; full 3-loop HK infeasible; gap unchanged (`FFFF_HEAT_KERNEL_A1.md`) |
| Superspace / supergraph Ward | SUSY forces \(c_Y\) not baseline; inventory/partner Wards PASS; gap stands (`FFFF_SUPERSPACE_WARD.md`) |
| Projected double vs Y at \(L=0\) | Cancels after overlap fix |
| Projected finish line at \(a_X=0\), \(a_Y=-1\) | Higher poles OK; simple gap \(\approx+0.0097\) stands (`FFFF_SIMPLE_GAP_AUDIT.md`) |
| First corrected equation in disagreement | Earlier FFff(A) bilinear \(-4/27\) (now \(0\)) |

The superseded corner numerator \(2d\varepsilon/3\) is the error the latest QCc note corrects. The numerator before any \(\varepsilon\to 0\) truncation is
\[
A+(k+1)B=-\frac{d(4-d)}{3}+3\varepsilon\cdot\frac{2d}{3}.
\]
Substituting \(d=4-2\varepsilon\) gives \(4d\varepsilon/3\), identically. Integrated on the unit square against \(p^{\varepsilon-1}q^{\varepsilon-1}\), the change from \(2d\varepsilon/3\) to \(4d\varepsilon/3\) adds
\[
\frac{8}{3\varepsilon}-\frac{4}{3}
\]
to the projective integral. The \(1/\varepsilon\) coefficient moves from \(16/3\) to \(8\). With \(\Gamma(3\varepsilon)\sim 1/(3\varepsilon)\), Wick weight \(+1/4\), and \(1/[4d(d-1)]=1/48\) at \(d=4\),
\[
\frac{8}{3\cdot 4\cdot 48}=\frac{1}{72}.
\]
The same \(z^{-3\varepsilon}\) factor turns the double pole into \(-L/(24\varepsilon)\).

The long-pair tail, kept at \(d=4-2\varepsilon\), is
\[
\frac{8d\varepsilon}{3}\,\mathrm{B}(1-2\varepsilon,1+\varepsilon)\,\mathrm{B}(2\varepsilon,2\varepsilon)\,\mathrm{B}(\varepsilon,1)\,\mathrm{B}(1-\varepsilon,1-\varepsilon)
=\frac{32}{3\varepsilon}+\frac{80}{3}+O(\varepsilon).
\]
The geometric Dirichlet integral is \(-8/\varepsilon-32+O(\varepsilon)\). The once-subtracted evanescent face integrates to
\[
\frac{2}{3}-\frac{4(r+2)\ln r}{3(r-1)}
\]
before the fraction average. Those pieces reproduce the corrected finite projective value \(-18.83584915236\). The factor \(1/[d(d-1)]=(1/12)(1+(7/6)\varepsilon+\cdots)\) then gives the simple pole above. No \(\zeta(2)\) term enters that pole: the constant term in \(\Gamma(3\varepsilon)e^{3\gamma_E\varepsilon}\) cancels, and the \(O(\varepsilon)\) piece multiplies the double pole into a finite remainder.

## CcFf

The \(p=0\) triangle was rebuilt from the conductance Laplacian, not from the checkpoint's covariance formulae. With the spin weight \(qN=d/(2S)\) and one factor \(S^{-d/2}\), the fraction averages are
\[
H_0=-\frac{d^2(d-4)}{6},\qquad
H_1=-\frac{d(d-4)(17d^2-72d+4)}{72}.
\]
The second expands as \(-4\varepsilon/3+O(\varepsilon^2)\). Its unit-square integral is the overlap \(-4/(3\varepsilon)\) plus a finite piece. The \(d=4\) cubic-vertex face integrates to \(+2\) against \(dP/P\), and that \(2\) is the coefficient of \(1/\varepsilon\) coming from the short-fermion measure. The long-scalar-pair model, continued in the fermion-fraction beta function before \(\varepsilon\to 0\), is
\[
-\frac{8}{3\varepsilon}-\frac{20}{3}+O(\varepsilon).
\]
The projective residue is \(2-4/3-8/3=-2\). The flat weight \(-1/2\) and the same curvature denominator produce
\[
(-2)\cdot\frac{-1/2}{3\cdot 48}=\frac{1}{144},
\]
and \(z^{-3\varepsilon}\) supplies \(-L/(48\varepsilon)\).

An independent quasi-Monte Carlo of this integrand (long-pair model subtracted and restored analytically, short faces subtracted by a finite difference) extrapolates to
\[
\varepsilon I\to -1.99997,\qquad \text{simple pole}\to -0.02480
\]
at \(L=0\). That matches the checkpoint value \(-0.0247\) inside the tolerance written on that check.

## FFff(A)

Class A is the multiplicity-8 matching with Wick weight \(-1/8\). The corrected soft-corner calculation is in `FFFF_A_CORNER_CORRECTION.md` (connected cycle, full two-difference Wick covariance). That check was re-run here and passes.

At \(d=4\) the angular averages of the projective density \(H\) are
\[
\bigl\langle H,\partial_q H,\partial_r H,\partial_q\partial_r H\bigr\rangle
=\bigl(-\tfrac{32}{3},\,0,\,-8,\,0\bigr).
\]
The earlier bilinear average \(-4/27\) is wrong; the mixed derivative vanishes. There is **no triple pole from the short-short corner**. With unit physical cutoffs the finite-\(d\) corner Laurent series is \(10/(3\varepsilon)-169/3\). Multiplication by \(\Gamma(3\varepsilon-1)\sim-1/(3\varepsilon)\), the Wick factor, and \(1/[4d(d-1)]\) gives the corner double pole
\[
+\frac{5}{1728\,\varepsilon^2}.
\]
The once-soft faces are rational with fifth-order denominators in the spectator effective lengths. At fixed spectator fraction they decay at large pair length; the nonuniform \(g\to 0\) region with \(Q=rg\) fixed is the long-pair endpoint.

The long-pair add-back models (`FFFF_A_LONG_PAIR.md`) were re-checked: each long-pair/short-spectator residue is \(+1\) at projective order \(1/\varepsilon^2\), and the two models together contribute
\[
+\frac{1}{576\,\varepsilon^3}
\]
to the graph. The analytic two-tail Dirichlet form has leading coefficient \(2/\varepsilon^2\) in the projective integral.

### Regulated full-graph inference

`FFFF_A_REGULATED.md` subtracts both long-pair models and short Taylor faces at nonzero \(\varepsilon\). Independently re-run at \(\varepsilon=0.00375\), power \(2^{17}\), the three cuts \(0.5,1,2\) give \(\varepsilon^2 I\approx 1.9877\) and \(\varepsilon(I-2/\varepsilon^2)\approx -3.28\), consistent with the note's finer extrapolation to projective
\[
I=\frac{2}{\varepsilon^2}-\frac{10}{3\varepsilon}+c_0,\qquad c_0\approx 14.93.
\]
The conversion through \(-\Gamma(3\varepsilon-1)e^{3\gamma_E\varepsilon-3\varepsilon L}/[32d(d-1)]\) is exact (`ffz_ffff_A_regulated_conversion.py`) and yields, in units \(h^2 z/K\),
\[
\frac{1}{576\varepsilon^3}+\Bigl(\frac{5}{1152}-\frac{L}{192}\Bigr)\frac{1}{\varepsilon^2}
+\Bigl(\frac{c_0}{1152}+\frac{79}{6912}+\frac{\pi^2}{768}-\frac{5L}{384}+\frac{L^2}{128}\Bigr)\frac{1}{\varepsilon}.
\]
With \(c_0\approx 14.93\) the constant simple piece is about \(0.03724\). The higher poles match the long-pair triple with no extra triple from the remainder. This is numerical evidence for the full A graph, not an analytic remainder evaluation.

## FFff(B)

`FFFF_B_GAUSSIAN.md` reduces the parity-even four-kinetic contraction and the Weyl covariance tilt; the runnable Gaussian check passes. Later notes evaluate the graph poles: bare triple \(1/48\), bare double \(-1/96\) at \(L=0\), and after the two-log overlap assembly fix (`FFFF_B_OVERLAP_FIX.md`) double \(-1/288-L/16\), matching the budget target that used to read
\[
\frac{5}{192\varepsilon^3}+\Bigl(\frac{11}{1152}-\frac{5L}{64}\Bigr)\frac{1}{\varepsilon^2}
\]
before the A \(\times 4\) trace correction (post-correction B target is \(1/48\) and \(-1/288-L/16\)).

## Counterterm budget (from `COUNTERTERM_BUDGET.md`)

The packet's two-loop insertion into the two-loop bare baseline, combined with the one-loop mass insertion, cancels the triple pole and leaves double pole \((L/12-5/48)/\varepsilon^2\). Against the Y target and the verified bare subtotal \(5/144\), the remaining **bare FFff(A)+FFff(B)** budget is
\[
\frac{1}{36\varepsilon^3}+\Bigl(\frac{1}{72}-\frac{L}{12}\Bigr)\frac{1}{\varepsilon^2}.
\]
The runnable asserts in that note were re-run and pass. The simple-pole counterterm formula there cancels the explicit \(L\) and \(L^2\) against the implied FFff higher-pole logarithms; the constant simple pole remains open.

## Finish line

Not reached. Independently verified / inferred pieces in \(h^2 z/K\):

| Piece | Triple | Double | Simple |
| --- | --- | --- | --- |
| QCc (doubled) | \(0\) | \(1/72\) | \(-0.01649742330-L/24\) |
| CcFf | \(0\) | \(1/144\) | \(\approx-0.02480-L/48\) |
| FFff(A) (regulated inference) | \(1/576\) | \(5/1152-L/192\) | \(\approx 0.03724-5L/384+L^2/128\) |
| Combined one- and two-loop CT | \(0\) | \(L/12-5/48\) | see budget note |
| FFff(B) (overlap-fixed) | \(1/48\) | \(-1/288-L/16\) | simple constant open |

Gap audit (`FFFF_SIMPLE_GAP_AUDIT.md`): every simple-pole piece re-checked (Y \(1/24\), CT, QQ, CCcc, QCc tadpole pole, QFf derived; QCc/CcFf/A/B within budget). **No piece is off by \(\sim0.01\).** Projected total \(\approx+0.0513\) vs \(1/24\) leaves gap \(\approx+0.0097\pm0.0024\). Physical follow-up (`FFFF_GAP_PHYSICAL.md`): \(c_Y=1/24\) is RG-fixed, but baseline\(=1/24\) is a matching test (not a t Hooft theorem) — residual is scheme-absorbable (\(\sim+23\%\) of the three-loop \(1/G_{\mathrm{ind}}\) piece); no missing Wick topology; no simultaneous \((a_X,a_Y)\) solution. Non-overlap finite of the B forest copy: \(\Delta c_0=0\) exact — not a second \(64/3\) (`FFFF_B_NONOVERLAP_FINITE.md`); gap stands. Independent finite-ε Laurent cross-check (`FFFF_FINITE_EPS_CROSSCHECK.md`): no graph-by-graph disagreement; new gap \(\approx+0.0095\pm0.0036\) still consistent with \(\sim+0.01\).

## What was derived, and what was only run

Derived in this review, with new SymPy/mpmath/NumPy code:

- QCc corner \(4d\varepsilon/3\), overlap, tail, faces, double pole \(1/72\), simple pole constant
- CcFf corner, overlap, geometric face, tail, double pole \(1/144\), QMC simple pole \(\approx-0.02480\)
- Re-run of the corrected FFff(A) cycle Wick corner (averages, no triple pole, corner double \(5/1728\))
- Re-run of the counterterm-budget asserts
- Re-run of the FFff(A) long-pair add-back asserts (triple \(1/576\))
- Exact conversion of the regulated projective series to graph poles and B targets
- Re-run of the FFff(B) Gaussian algebra check
- Smoke re-run of the regulated FFff(A) QMC (consistent with \(2/\varepsilon^2\) and \(-10/(3\varepsilon)\))

Retracted: the earlier FFff(A) bilinear \(-4/27\) and triple pole \(1/7776\) (incorrect Wick/assembly for the mixed term).

Not done: full analytic FFff(A) remainder; independent operator form of the B coupling-weight overlap (assembly fix applied). Gap audit re-derived QQ/CCcc/QFf/CT/Y; the \(\approx+0.01\) residual versus \(1/24\) is now a standing result under the present graph+CT list.
