# Sharp phase-lock boundary vs smooth (AI)

Runnable: `analysis/ffz/ffz_phase_lock.py` + test. Stacked on #61 (AH).  
**No tuning** as a *claim*; when a free BC/contact parameter must be set to GAP, that is stated plainly.

## Results

| Boundary | Derived \(\delta C\) into \(h^2z/K\) | Free parameter? |
| --- | --- | --- |
| Smooth (no \(\partial M\)) | **0** | none |
| Dirichlet | **0** (surface \(a_1\sim K/3\) is \(\Lambda^2\) channel) | none |
| Robin (any \(S\)) | **0** (same channel split) | \(S\) free but irrelevant to \(C\) |
| Contact ansatz \(\delta C=\alpha\) | \(\alpha\) | **yes** — \(\alpha=\mathrm{GAP}\) tuned |

## Tuning statement

To hit GAP \(+0.009651\) via a boundary/contact ansatz requires setting free \(\alpha=\mathrm{GAP}\) (identical to scheme \(\kappa=-\mathrm{GAP}\) from #23). **That is tuning / scheme, not a derived phase-lock prediction.**

## Coincidences (not derived)

| Name | Value | Note |
| --- | --- | --- |
| \(\ln 2/72\) | \(0.009627\) | within \(1\sigma\) of GAP |
| \(\alpha=\mathrm{GAP}\) | \(0.009651\) | tautology if chosen |
| \(S=\mathrm{GAP}\) sample | — | identifying Robin \(S\) with pole units is not derived |

## Plain

Dan’s sharp interaction surface can be modeled; it does not move the three-loop pole without a free coefficient dialed to the gap.
