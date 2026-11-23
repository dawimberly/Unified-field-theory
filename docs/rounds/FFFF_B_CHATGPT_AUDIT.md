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

Double-pole \(L\) coefficients: \(A\) supplies \(-L/48\), \(B\) supplies \(-L/16\) from its triple, CT supplies \(+L/12\). Sum vanishes. After the overlap fix the double residual is \(0\), so the residual-driven \(+L/48\) in the simple pole cancels as well. (Bare field-weight-only assembly had residual \(-1/144\) and uncanceled \(+L/48\); ChatGPT’s \(59/648\) remains an absent/assembly-error claim.)

## What is derived vs reproduced / assumed

**Derived**

- B kernel conventions; triple \(1/48\); bare slope \(224\to\) double \(-1/96\)
- Vertex-ll residue exactly \(-16\); overlap fix \(\to\) double \(-1/288\); projected residual \(0\)
- Conversions \(+1/18\), \(-1/54\); \((2\varepsilon-1)\) identity; \(59/648\) decoding

**Assumed in the fix**

- The missing soft-soft factor is one extra copy of the vertex-ll piece (coupling-weight / overlap forest), not a different pl-active operator

**Reproduced**

- Measure slope \(+192\); corrected A double \(5/288-L/48\)

## Relation to prior review notes

See `FFFF_B_OVERLAP_FIX.md`. Bare regulated slope \(\sim 224\) is unchanged as a field-weight-only diagnostic; the ledger now uses the overlap-corrected B double \(-1/288\).

## Follow-up: contractions and ε×triple contamination

Runnable: `ffz_ffff_B_contraction_contamination.py`.

### Why a loop (PSE / Tong framing)

The Yukawa vertex is one scalar + two Weyl fermions (no \(\phi^2\psi\) operator). Four scalar legs therefore need four vertices. Pairing those legs into two internal scalars and wiring the fermions as a single 4-cycle gives \(V=4\), \(I=6\), \(L=3\). The 2+2 fermion wiring is FFff(A); the single 4-cycle is FFff(B). Signed weight used in the poles: \(32\times(1/2)^4/2!=+1\). The packet’s “\(2!^2\)” wording would give \(+1/2\); that global rescaling is incompatible with the passing triple sum \(1/36\), so it is flagged as wording ambiguity only.

### Why the bare double failed while the triple passed

Field-weight-only extraction gives slope \(224\) and double \(-1/96\). Budget needs \(-1/288\) (slope \(208\)): \(\Delta a_1=-16\). A relative \((1+\delta\varepsilon)\) on the triple density with \(\delta=1/3\) produces that shift and equals \(-\delta\times(1/48)\) in graph units. The overlap fix implements that \(\Delta a_1\) as one extra vertex soft-soft copy (see above). Dropping the explicit \((k-1)\) piece (\(+144\)) is the wrong fix.

## First concrete \(-16\) term: vertex on double-log overlap

Runnable: `ffz_ffff_B_slope16_hunt.py`. Sector split of the numerator \(\varepsilon\)-derivative residue \(+32\):

| Piece | double-log (40 sectors) | power/log (72) | sum |
| --- | ---: | ---: | ---: |
| vertex \(\lambda_F\) | **\-16** | \(-48\) | \(-64\) |
| scalar curvature mass | \(0\) | \(-128\) | \(-128\) |
| explicit \((k-1)\) | \(+144\) | \(0\) | \(+144\) |
| \(B_d\) | \(+77.\overline{3}\) | \(+2.\overline{6}\) | \(+80\) |
| **full** | \(+205.\overline{3}\) | \(-173.\overline{3}\) | **\(+32\)** |

The fermion-vertex piece on the two overlapping logarithmic endpoints is **exactly \(-16\)**.

## Fix applied (overlap assembly)

Runnable: `ffz_ffff_B_overlap_fix.py`. Details: `FFFF_B_OVERLAP_FIX.md`.

Global \(\lambda_F\) doubling also shifts power/log by \(-48\) and overshoots. The applied fix adds one extra copy of the **double-log** vertex piece only (\(\Delta a_1=-16\)):

\[
\text{slope }224\to 208,\qquad
\text{B double }-\frac1{96}\to -\frac1{288}\quad(L=0).
\]

Projected double residual vs Y is then **0**; the residual-driven \(+L/48\) cancels. Triples unchanged.

Reading: Yukawa coupling Weyl weight \(\varepsilon=\lambda_F\) at this order, missed on the two-log overlap assembly. Field-weight-only finite-\(\varepsilon\)/Gaussian kernels are unchanged; the correction sits in the pole assembly step.
