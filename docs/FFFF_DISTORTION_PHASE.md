# Distortion phase and the Diósi–Penrose \(R_0\)

Runnable: `bots/gpt/research_batches/2026-11-21-ffz/review/ffz_distortion_phase.py`  
Test: `.../test_ffz_distortion_phase.py` (hooked in `scripts/run_checks.py`).  
Stacked on induced-collapse (#24).  
**No tuning.** DERIVED vs ASSUMED labeled. Honest negatives.

## (1) Three regimes (GL / phase-field)

Relative phase \(\theta\) between superposition branches. Standard energy density:

\[
f=\frac{\hbar^2 n}{2m}\,|\nabla\theta|^2+V(\theta).
\]

| Regime | Meaning |
| --- | --- |
| **Locked** | \(\theta\) near a minimum of \(V\); coherent (wave-like) |
| **Distorted** | Finite core where \(\theta\) slips across the barrier; width \(\xi_\theta\); **recoverable** (cf. Langer–Ambegaokar / McCumber–Halperin TAPS, quantum phase slips, Minev et al. 2019 reversible jump window) |
| **Collapsed** | One branch selected; relative-phase coherence lost (DP completed) |

**Derived:** regime classification from this free-energy structure.  
**Assumed:** barrier \(V(\theta)\) is **not** fixed by the WZ+Yukawa Newton ledger / \(C_{\mathrm{ledger}}\).

## (2) Natural width of the distortion region

GL healing / slip-core scaling (**derived** formula):

\[
\xi=\frac{\hbar}{\sqrt{2m|\alpha|}},\qquad
\xi_\theta\sim\hbar\sqrt{\frac{n}{2m V_0}},\qquad
\xi_{\mathrm{GP}}=\frac1{\sqrt{8\pi n a_s}}.
\]

\(|\alpha|\), \(V_0\), \(a_s\) need inputs (**assumed** when numbers are quoted).

## (3) Candidates vs Gran Sasso \(R_0\ge5.4\times10^{-11}\,\mathrm{m}\)

| Candidate | Width | Meets bound? | From induced-G / WZ \(C_{\mathrm{ledger}}\)? |
| --- | ---: | --- | --- |
| \(\bar\lambda_e\) | \(3.86\times10^{-13}\,\mathrm{m}\) | no (\(\times0.007\)) | mass only |
| \(\bar\lambda_p\) | \(2.10\times10^{-16}\,\mathrm{m}\) | no | mass only |
| **Bohr \(a_0\)** | \(5.2918\times10^{-11}\,\mathrm{m}\) | **no** (ratio \(0.980\) — just under) | **No** — atomic EM; **COINCIDENCE** (\(\sim2\%\) from the bound) |
| BEC GP healing (\(10^6\) Rb, \(R=10\,\mu\mathrm{m}\), \(a_s=5.3\,\mathrm{nm}\)) | \(1.77\times10^{-7}\,\mathrm{m}\) | yes (\(\times3.3\times10^3\)) | condensate case; **not** universal \(R_0\) from \(C_{\mathrm{ledger}}\) |
| Cooper \(\xi\) (Al-like) | \(1.0\times10^{-7}\,\mathrm{m}\) | yes | **assumed** material scale |
| \(\ell_P(G_N)\) | \(1.62\times10^{-35}\,\mathrm{m}\) | no | \(G_N\) |
| \(\ell_P(G_{\mathrm{ind}})\) | \(1.46\times10^{-35}\,\mathrm{m}\) | no | \(C_{\mathrm{ledger}}\) via \(G_{\mathrm{ind}}\propto1/C\) |

**Plain verdict:** **No** — nothing in the induced-gravity / WZ+\(C_{\mathrm{ledger}}\) content derives a universal Gran-Sasso-viable \(R_0\) without tuning. Native lengths are far too small; condensate healing lengths are system-specific and much larger. \(a_0\) is \(\sim2\%\) below \(5.4\times10^{-11}\,\mathrm{m}\) — a **coincidence**, not a derivation (same honesty rule as \(\ln2/72\)).

## (4) Modified collapse times (\(R_0=\) derived/assumed width)

Induced factor \(\tau_{\mathrm{ind}}/\tau(G_N)=C_{\mathrm{ledger}}/(1/24)\approx1.232\) (**derived** scaling; **assumed** \(G_N\leftrightarrow1/24\)).

| Case | \(R_0\) choice | \(\tau(\mathrm{ind})\) | vs coherence | Gravity dominates? |
| --- | --- | ---: | --- | --- |
| BEC \(10^6\) Rb, \(d=1\,\mathrm{mm}\) | GP healing \(1.77\times10^{-7}\,\mathrm{m}\) | \(1.38\times10^{7}\,\mathrm{s}\) | \(\sim1\,\mathrm{s}\) | **No** |
| Flux \(N_{\mathrm{eff}}=10^9\) | Cooper \(\xi=10^{-7}\,\mathrm{m}\) | \(5.33\times10^{10}\,\mathrm{s}\) | \(\sim50\,\mu\mathrm{s}\) | **No** |
| Nanosphere \(10^{-14}\,\mathrm{kg}\) | geometric \(R=1.06\,\mu\mathrm{m}\) | \(6.48\times10^{-2}\,\mathrm{s}\) | \(\sim1\,\mathrm{ms}\) | **No** |
| Nanosphere (probe) | \(R_0=a_0\) **coincidence only** | \(8.58\times10^{-7}\,\mathrm{s}\) (not a prediction) | \(\sim1\,\mathrm{ms}\) | probe only (would dominate — **not** claimed) |

## (5) Distinct testable prediction?

**Yes, qualitatively:** a **recoverable distortion window** (phase-slip / incomplete jump) with duration and barrier/mass scaling from GL dynamics — distinct from irreversible plain DP at \(\tau=\hbar/E_G\). Concrete window times need \(V(\theta)\) (**assumed**).

**No for \(R_0\):** the distortion phase does **not** fix DP’s smearing problem from \(C_{\mathrm{ledger}}\) alone.

## Plain statement

The distortion phase is a useful intermediate between locked and collapsed and suggests Minev-like reversible-window tests, but it does **not** supply a tuning-free Diósi–Penrose \(R_0\) from induced gravity. Treating \(a_0\approx\) Gran Sasso bound as the answer would be repeating the \(\ln2/72\) coincidence mistake.
