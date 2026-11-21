# Independent review of the FFZ three-loop Newton poles

Commit reviewed: `196c696` (FFZ Newton matching review packet).
Prescription: fixed \(\xi=1/6\), \(d=4-2\varepsilon\), equal masses, units \(h^2 z/K\), \(L=\ln(z/\mu^2)\).

Runnable checks in this directory (they do not call the checkpoint scripts):

- `ffz_independent_pole_review.py` — QCc and CcFf poles
- `ffz_ffff_A_poles.py` — FFff(A) soft-corner Taylor averages and triple pole

SymPy and mpmath own the exact identities. The CcFf simple pole is a separate quasi-Monte Carlo of an integrand written for this review.

## Verdicts

| Claim | Verdict |
| --- | --- |
| QCc double pole \(1/(72\varepsilon^2)\) | Agree |
| QCc simple pole \(-0.01649742330-L/24\) | Agree (constant recomputed as \(-0.0164974232969\)) |
| CcFf double pole \(1/(144\varepsilon^2)\) | Agree |
| CcFf simple pole about \(-0.0247-L/48\) | Agree at the checkpoint's stated precision. This review's extrapolation is \(-0.02480\) |
| FFff(A) short-short triple pole | None (earlier \(1/7776\) retracted; see correction) |
| FFff(A) corner double pole \(5/(1728\varepsilon^2)\) | Agree with corrected cycle Wick check |
| FFff(A) full double / simple / graph triple | Open (faces, long overlaps, FFff(B)) |
| Projected finish line at \(a_X=0\), \(a_Y=-1\) | Not reached |
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
The once-soft faces are rational with fifth-order denominators in the spectator effective lengths (printed by the correction script). Their \(r\)-independent \(q\)-linear piece is nonuniform at the spectator fraction endpoint. Short/long overlaps and FFff(B) are still required for the full graph poles.

## Counterterm budget (from `COUNTERTERM_BUDGET.md`)

The packet's two-loop insertion into the two-loop bare baseline, combined with the one-loop mass insertion, cancels the triple pole and leaves double pole \((L/12-5/48)/\varepsilon^2\). Against the Y target and the verified bare subtotal \(5/144\), the remaining **bare FFff(A)+FFff(B)** budget is
\[
\frac{1}{36\varepsilon^3}+\Bigl(\frac{1}{72}-\frac{L}{12}\Bigr)\frac{1}{\varepsilon^2}.
\]
The runnable asserts in that note were re-run and pass. The simple-pole counterterm formula there cancels the explicit \(L\) and \(L^2\) against the implied FFff higher-pole logarithms; the constant simple pole remains open.

## Finish line

Not reached. Independently verified pieces in \(h^2 z/K\):

| Piece | Triple | Double | Simple |
| --- | --- | --- | --- |
| QCc (doubled) | \(0\) | \(1/72\) | \(-0.01649742330-L/24\) |
| CcFf | \(0\) | \(1/144\) | \(\approx-0.02480-L/48\) |
| FFff(A) short-short corner | \(0\) | \(+5/1728\) (corner only) | open |
| Combined one- and two-loop CT | \(0\) | \(L/12-5/48\) | see budget note |

Still open: FFff(A) faces and long overlaps, FFff(B), and the projected simple pole versus \(1/24\).

## What was derived, and what was only run

Derived in this review, with new SymPy/mpmath/NumPy code:

- QCc corner \(4d\varepsilon/3\), overlap, tail, faces, double pole \(1/72\), simple pole constant
- CcFf corner, overlap, geometric face, tail, double pole \(1/144\), QMC simple pole \(\approx-0.02480\)
- Re-run of the corrected FFff(A) cycle Wick corner (averages, no triple pole, corner double \(5/1728\))
- Re-run of the counterterm-budget asserts

Retracted: the earlier FFff(A) bilinear \(-4/27\) and triple pole \(1/7776\) (incorrect Wick/assembly for the mixed term).

Not done: checkpoint QMC scripts not executed; QQ, QCc tadpole, CCcc, QFf, FFff(B), and full FFff(A) faces not re-derived as complete Laurent series.
