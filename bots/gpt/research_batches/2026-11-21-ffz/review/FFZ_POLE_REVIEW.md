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
| FFff(A) triple pole \(1/7776\) | Derived (soft-corner bilinear only) |
| FFff(A) double / simple poles | Open (once-soft faces unfinished) |
| Projected finish line at \(a_X=0\), \(a_Y=-1\) | Not reached |
| First corrected equation in disagreement | None |

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

Class A is the multiplicity-8 matching: fermion pairs \((0,2)\) and \((1,3)\), scalars \((0,3)\) and \((1,2)\), vertex 0 fixed at the origin. The flat Wick weight is \(-1/8\).

Each parallel pair has spin factor \(N_a=d\,c_a/(4 s_{a1}s_{a2})\) and endpoint difference \(\Delta_a\). The curvature uses the two-constraint moment
\[
M_{ij}=\frac{\mathbb{E}[(x_i\cdot x_j)\,\Delta_1^2\Delta_2^2]}{\mathbb{E}[\Delta_1^2\Delta_2^2]}.
\]
Vertex weights are \(-\lambda_F M_{ii}\) with \(\lambda_F=(4-d)/2\) (the root is omitted). Scalar edges carry the evanescent mass \(-d(4-d)s/3\). The mass profile is
\[
B=\frac{2}{3}\sum_e s_e\bigl(M_{ii}+M_{jj}+M_{ij}+d_e s_e\bigr),
\]
with \(d_e=d-1\) on fermion edges and \(d_e=d\) on scalar edges. The Symanzik polynomial is exactly \(U=qr(v+w+Q+R)\). Each \(N_a\) scales as \(\rho^{-1}\). The radial factor is
\[
\Gamma(3\varepsilon-1)\,z^{1-3\varepsilon}\bigl[A+(3\varepsilon-1)B\bigr],
\]
equivalently \(\Gamma(3\varepsilon)\,z^{1-3\varepsilon} B\) when \(A=0\) at \(d=4\). Divide by \(4d(d-1)\) and multiply by \(-1/8\).

On the gauge \(v+w=1\), the density multiplying \(q^{\varepsilon-2}r^{\varepsilon-2}\,dq\,dr\,dg\,dh\,dv\) is
\[
H=(qN_1)(rN_2)\,S^{-2}\,B,\qquad S=1+Q+R.
\]
At the double-soft corner, \(H\to 32/3\). The exact flat angular averages of the Taylor coefficients are
\[
\langle a_{00}\rangle=\frac{32}{3},\qquad
\langle\partial_q H\rangle=0,\qquad
\langle\partial_r H\rangle=8,\qquad
\langle\partial_q\partial_r H\rangle=-\frac{4}{27}.
\]
(The root/non-root asymmetry makes \(\partial_q H\) vanish while \(\partial_r H\) does not.) The constant term integrates to \(1/(\varepsilon-1)^2\) and is non-singular. The bilinear coefficient supplies the projective \(1/\varepsilon^2\). With \(\Gamma(3\varepsilon)\sim 1/(3\varepsilon)\) and weight \(-1/384\),
\[
\text{FFff(A) triple pole}=\frac{1}{7776\,\varepsilon^3}.
\]
The corner linear piece alone would give a double pole \(1/(144\varepsilon^2)\), but that is incomplete: the once-soft faces (integrated \(\partial_r H\) at \(r=0\) with the short-\(q\) overlap subtracted, and the symmetric \(q\)-face) still have to be evaluated. Those faces converge at large pair length without an extra growing tail model for \(H\) itself, but the derivative face has not been reduced to a stable number in this review.

## Finish line

Not reached. The Y-matching targets for the baseline at \(\xi=1/6\) are
\[
\frac{1}{36\varepsilon^3}-\frac{1}{18\varepsilon^2}+\frac{1}{24\varepsilon}
\]
(so that at \(a_Y=-1\) the higher poles cancel and the simple pole cancels the known Y residue).

Independently verified pieces in \(h^2 z/K\):

| Piece | Triple | Double | Simple |
| --- | --- | --- | --- |
| QCc (doubled) | \(0\) | \(1/72\) | \(-0.01649742330-L/24\) |
| CcFf | \(0\) | \(1/144\) | \(\approx-0.02480-L/48\) |
| FFff(A) soft corner | \(1/7776\) | faces open | open |
| Mass-only \(Z_z\) CT (packet algebra, rechecked) | \(1/12\) | \(-1/16-L/12\) | \(\zeta(2)/24+L/16+L^2/24\) |

Bare doubles verified without FFff: \(1/72+1/72+1/144=5/144\) (matches the packet's scalar+QFf+CcFf double-pole subtotal; QFf and CCcc add no double pole). After the mass CT constant \(-1/16\), that subtotal is \(5/144-1/16=-1/36\). The target double is \(-1/18\), so FFff(A) faces, FFff(B), and the remaining two-loop counterterms must still supply \(-1/18-(-1/36)=-1/36\) in the double pole.

For triples: mass CT \(1/12\) plus FFff(A) \(1/7776\) overshoots the target \(1/36\) by \(433/7776\). Other counterterms must cancel that excess.

Still open:

1. FFff(A) once-soft faces (completes its double and feeds its simple pole).
2. FFff(B), the four-vertex spin-trace graph.
3. The remaining two-loop curvature counterterms beyond mass-only \(Z_z\).
4. The full projected sum versus \(1/24\).

## What was derived, and what was only run

Derived in this review, with new SymPy, mpmath, or NumPy code:

- QCc corner \(4d\varepsilon/3\), overlap shift \(8/(3\varepsilon)-4/3\), residue \(8\), double pole \(1/72\), logarithm \(-L/24\)
- QCc tail \(32/(3\varepsilon)+80/3\), geometric pole \(-8/\varepsilon-32\), evanescent face antiderivative, and the simple-pole constant \(-0.0164974232969\)
- CcFf corner \(H_0\) and \(H_1\) from the Laplacian, overlap \(-4/(3\varepsilon)\), geometric residue \(+2\), tail \(-8/(3\varepsilon)-20/3\), double pole \(1/144\), logarithm \(-L/48\)
- CcFf simple pole \(-0.02480\) from a new quasi-Monte Carlo of that integrand
- FFff(A) polynomial \(U\), radial factor, corner value \(32/3\), exact averages \(\langle\partial_r H\rangle=8\) and \(\langle\partial_q\partial_r H\rangle=-4/27\), triple pole \(1/7776\)

Checked by re-running algebra already in the packet (not a new curved integral): mass-only \(Z_z\) Laurent coefficients \(1/12\), \(-1/16-L/12\), and the simple-pole formula with \(\zeta(2)\).

Not done: the checkpoint QMC scripts were not executed. QQ, the QCc tadpole, CCcc, QFf, and FFff(B) were not re-derived. FFff(A) double and simple poles are not claimed.
