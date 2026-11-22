# FFff(B) audit of parallel ChatGPT claims

Parallel ChatGPT notes are treated as **claims to test**, not as inputs. This note records an independent FFff(B) derivation that does **not** use endpoint/region-splitting of the full regulated integrand for the poles, then checks the quoted finite-ε numbers and the 59/648 ledger claim against this branch’s bookkeeping.

Runnable check: `ffz_ffff_B_chatgpt_audit.py` (OPENBLAS_NUM_THREADS=1).

## Kernel (derived)

Four-vertex Weyl fermion loop with two internal scalar lines parallel to fermion edges 02 and 13. Edge order `(02,12,13,03,02,13)`. The signed flat Wick / symmetry / vacuum-cumulant weight after extracting \(z|y|^4\) is \(+1\); four kinetic denominators contribute \(1/16\); curvature conversion with that factor is \(1/[64\,d(d-1)]\). The trace keeps the full Weyl contraction (the factor of two is inside \(\mathrm{tr}\)). Evanescent pieces kept at finite \(d\):

- vertex weights \(\lambda_F=(4-d)/2\) in \(A\);
- scalar curvature mass \(d(4-d)/3\) on the two scalar edges in \(A\);
- radial / mass-profile weight \(B\), including \((d-1_{\mathrm{fermion}})\) on lines.

Radial reduction uses \(\Gamma(k-1)\) with \(k=6-3d/2\) and numerator \(WA+(k-1)B\).

## Poles (independent scheme)

| Quantity | Result | Method |
| --- | --- | --- |
| Projective double-log residue | \(-48\) (families \(-8,-16,-16,-8\)) | Exact ordered-sector integration of the two-log families only |
| Bare triple in \(h^2z/K\) | \(1/48\) | \(\Gamma(3\varepsilon-1)/(64d(d-1))\) on that residue |
| Numerator \(\varepsilon\)-derivative residue | \(+32\) | Double-projective residue on \(\partial_\varepsilon J\) at \(d=4\) |
| Measure contribution to slope | \(+192\) | Frozen-\(J_4\) diagnostics (reproduced, not re-derived here) |
| Full projective slope | \(224\) | \(192+32\) |
| Bare double at \(L=0\) | \(-1/96\) | Same conversion on \(-48/\varepsilon^2+224/\varepsilon\) |

Power/log families integrate to zero at the triple order (existing ordered-sector check). This is sector analysis of the leading singular families and a residue calculus for the subleading numerator, not a Taylor endpoint subtraction of the full \(\varepsilon\neq 0\) integrand.

### Internal splits (not full-graph doubles)

| Piece | Graph double at \(L=0\) | Status |
| --- | --- | --- |
| Scalar curvature-mass part of \(B\) | \(+1/18\) | Derived conversion from residue \(-128\) |
| Homogeneous large-length E-split | \(-1/54\) | Derived: \((128/3)/(-3\cdot 64\cdot 12)\) |
| ChatGPT homogeneous simple | \(-11/324+\ln c/54\) | Not re-derived; used only to decode \(59/648\) |

The E-split \(-1/54\) is already inside the regulated graph that yields \(-1/96\). It must not replace the full double.

## Finite-\(\varepsilon\) checks

### \(\varepsilon=0.7\) (two fixed-vertex anchors)

Under this kernel’s Dirichlet/simplex projective normalization the script reports \(\Gamma(k-1)\) times the curvature density for anchors at vertices 0 and 2. The two anchors agree at the percent level; the flat anchors agree to sampling noise. The absolute ChatGPT values \(1.00333538\) and \(1.00333622\) are **not** reproduced in these units. Likely a different absolute normalization of a convergent test integral. Anchor independence of the finite-\(\varepsilon\) kernel is supported; the quoted absolute number is not.

### \(\varepsilon=1/2\) and \((2\varepsilon-1)\)

Distributional identity (derived):
\[
(2\varepsilon-1)\int_0^1 x^{2\varepsilon-2}\,\mathrm{d}x = 1.
\]
Setting \(\varepsilon=1/2\) before integrating drops that finite part. If a particular endpoint’s residue density integrates to \(4\pi\), the missed piece is \(4\pi\). The **full** projective curvature integral at \(\varepsilon=1/2\) in the simplex chart is finite and anchor-stable but is **not** equal to \(4\pi\) in these units—so the qualitative cancellation claim is accepted; the identification of the whole integral with \(4\pi\) is not.

## Bookkeeping vs \(59/648\) and \(L/\varepsilon\)

With corrected FFff(A) double \(5/288\) (at \(L=0\)), this B double \(-1/96\), bare non-FFff subtotal \(5/144\), and CT \(-5/48\):

\[
\frac{5}{144}+\frac{5}{288}-\frac{1}{96}-\frac{5}{48}=-\frac{1}{16},\qquad
-\frac{1}{16}-\Bigl(-\frac{1}{18}\Bigr)=-\frac{1}{144}.
\]

- Triple pole: cancels in the projected sum (\(A\) \(1/144\) + \(B\) \(1/48\) = \(1/36\)).
- Double residual vs Y: **\(-1/144\)**, not \(\pm 59/648\).
- First term responsible: **FFff(B) double constant \(-1/96\) vs budget \(-1/288\)** (gap exactly \(-1/144\)).
- ChatGPT arithmetic \(5/144-1/16=-1/36\) (mass CT only) is correct but incomplete for the projected double.
- Identity \(59/648=1/8-11/324\): the quoted gap equals the homogeneous-region **simple** \(-11/324\) mixed with \(1/8\). That is an assembly error if that simple is entered into the double-pole ledger.

### Scale logarithms

Double-pole \(L\) coefficients: \(A\) supplies \(-L/48\), \(B\) must supply \(-L/16\) from its triple, CT supplies \(+L/12\). Sum vanishes (cutoff / scale logs cancel in the double). The residual double \(\Delta=-1/144\) still produces
\[
e^{-3\varepsilon L}\ \Longrightarrow\ \frac{+L}{48\varepsilon}
\]
in the simple pole. So **\(\ln(z/\mu^2)/\varepsilon\) does not cancel** while the double residual remains—agreeing with the ChatGPT qualitative claim, with coefficient \(+L/48\) rather than a \(59/648\) constant gap.

## What is derived vs reproduced / assumed

**Derived in the runnable script**

- B kernel conventions (edges, Wick \(+1\), \(1/16\), finite-\(d\) \(A,B\))
- Exact double-log triple \(\to 1/48\)
- Numerator \(\varepsilon\)-derivative residue \(+32\)
- Conversions \(+1/18\) (curvature mass) and \(-1/54\) (E-split)
- Full slope \(224\to\) double \(-1/96\) once the measure \(192\) is included
- Finite-\(\varepsilon\) simplex anchors at \(0.7\) and \(1/2\)
- \((2\varepsilon-1)\) endpoint identity
- Ledger residual \(-1/144\), decoding of \(59/648\), \(+L/48\) survival

**Reproduced / assumed (not re-derived here)**

- Measure slope \(+192\) (frozen-\(J_4\) regulated / residue diagnostics on this branch)
- Corrected A double \(5/288-L/48\) and A curvature-mass identity \(5/216\) under the A pref
- ChatGPT homogeneous simple \(-11/324+\ln c/54\) (used only algebraically)
- Absolute normalization that would make \(\varepsilon=0.7\) evaluate to \(\sim 1.00333\)

## Relation to prior review notes

Consistent with `FFFF_B_REGULATED.md` (intercept \(-48\), slope \(\sim 224\), double \(\sim -1/96\)), `FFFF_TRIPLE_NORMALIZATION.md` (triples \(1/144+1/48=1/36\)), `FFFF_B_MIXED_PAIR_IBP.md` (E-split \(-1/54\)), and `FFFF_A_ANALYTIC_DOUBLE.md` (projected residual \(\sim -1/144\), uncanceled \(+L/48\)).

## Follow-up: contractions and ε×triple contamination

Runnable: `ffz_ffff_B_contraction_contamination.py`.

### Why a loop (PSE / Tong framing)

The Yukawa vertex is one scalar + two Weyl fermions (no \(\phi^2\psi\) operator). Four scalar legs therefore need four vertices. Pairing those legs into two internal scalars and wiring the fermions as a single 4-cycle gives \(V=4\), \(I=6\), \(L=3\). The 2+2 fermion wiring is FFff(A); the single 4-cycle is FFff(B). Signed weight used in the poles: \(32\times(1/2)^4/2!=+1\). The packet’s “\(2!^2\)” wording would give \(+1/2\); that global rescaling is incompatible with the passing triple sum \(1/36\), so it is flagged as wording ambiguity only.

### Why the double fails while the triple passes

Independent B extraction gives projective slope \(224\) and graph double \(-1/96\). The budget needs double \(-1/288\), i.e. slope \(208\). Difference \(\Delta a_1=-16\).

A relative factor \((1+\delta\varepsilon)\) on the triple density shifts the projective \(1/\varepsilon\) coefficient by \(-48\delta\). Taking \(\delta=1/3\) gives exactly \(\Delta a_1=-16\), and in graph units
\[
-\frac1{96}-\Bigl(-\frac1{288}\Bigr)=-\frac1{144}=-\delta\times\Bigl(\frac1{48}\Bigr).
\]
So ChatGPT’s qualitative suspect (“\(\varepsilon\) times a triple pole can change the double while leaving the triple intact”) matches this review’s residual **quantitatively** for \(\delta=1/3\).

The numerator \(\varepsilon\)-derivative residue splits as explicit \((k-1)\) piece \(+144\) plus \(A,B\) \(d\)-dependence, netting \(+32\). Dropping the explicit \((k-1)\) piece removes \(144\), not \(16\), and is the wrong fix: \(A+(k-1)B\) with \(\Gamma(k-1)\) already equals \(\Gamma(k-1)A+\Gamma(k)B\), so that factor is required once. The note in `FFFF_A_ANALYTIC_DOUBLE.md` that the budget “needs \(+16\) not \(+32\)” is the cancellation requirement, not an independent evaluation that the residue is wrong by half.

**Still open:** identify a concrete missing/double-counted overlap term that removes exactly \(16\) from the slope (or an equivalent \(\delta=1/3\) relative error). Scalar-mass and E-split internals are not that term.
