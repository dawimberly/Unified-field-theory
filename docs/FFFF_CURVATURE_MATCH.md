# UV→IR curvature matching for hadronic \(a_1\)

Runnable: `ffz_curvature_match.py` + test. Stacked on #39.  
**No tuning.** Narrows [#36] band \([1/6,\ 7/6]\).

## Methods

| Piece | Role | Label |
| --- | --- | --- |
| Adler / Sakharov \(1/G\sim\int\mathrm{d}s\,\rho(s)/s\) | spectral structure | **derived** structure |
| QCD \(\rho_{TT}(s)\) piecewise model | hadronic + continuum | **assumed** |
| Soft/hard-wall AdS/QCD glueballs | residual after confinement | **assumed** (\(n=0\)–\(2\), central \(1\)) |
| ChPT / Donoghue–Leutwyler GFFs | pion \(\xi\) preference \(\approx 1/6\) | **assumed** soft preference |
| Dispersive continuum \(\sim 2\,\mathrm{GeV}\) | cancels free-gluon IR Kabats | model consistency, **not** a theorem |

## Narrowed band

\[
a_1^{\mathrm{matched}} = \tfrac13 \pm 0.29
\quad\in\quad \bigl[\tfrac16,\ \tfrac34\bigr]
\]

| | [#36] | matched |
| --- | ---: | ---: |
| low | \(1/6\) | \(1/6\) |
| high | \(7/6\) | \(3/4\) |
| central | — | \(1/3\) |
| sign | all \(+\) | all \(+\) |

Proxy \(+2/3\) still inside. Width shrinks \(1.0\to 0.583\).

## Coincidence flag

Continuum threshold \(\sim 2\,\mathrm{GeV}\sim\) few \(\times\Lambda_{\mathrm{QCD}}\) — structural, not miraculous.

## Plain

Stays positive. Narrower, but still model-dependent — do not treat \(1/3\) as a lattice theorem.
