# Simple-pole gap audit vs \(1/24\)

Runnable: `ffz_simple_gap_audit.py`. **No piece is tuned to hit \(1/24\).**

Prescription: fixed \(\xi=1/6\), \(a_X=0\), \(a_Y=-1\), \(L=0\), units \(h^2 z/K\).

## Ledger (every simple-pole contribution)

| Piece | Value | Error | Source | \(L\) tracked |
| --- | ---: | ---: | --- | --- |
| QQ | \(+1/16=+0.0625\) | exact | **derived** (basketball) | yes \(-L/24\) |
| CCcc | \(-5/72\approx-0.069444\) | exact | **derived** (projective \(16/3\), Wick \(-5/8\)) | none (simple only) |
| QCc doubled | \(-0.01649742330\) | \(\sim10^{-11}\) | reproduced + review analytic faces | yes \(-L/24\) |
| QCc tadpole | \(-1/36\approx-0.027778\) | exact pole | **derived**; \(O(1)\) is finite matching, **not** a \(1/\varepsilon\) term | none |
| QFf | \(-0.07466151326\) | \(\sim10^{-11}\) | **derived** (packet \(j_1\) integral) | none (cancels) |
| CcFf | \(\approx-0.02480\) | \(\sim0.001\) | review QMC | yes \(-L/48\) |
| FFff(A) \(\times4\) | \(+0.1496\) | \(\pm0.002\) | re-extracted \(c_0=15.11\pm0.58\) | yes from higher poles |
| FFff(B) overlap-only | \(-0.032709\) | \(\pm0.001\) | \(a_1=208\), \(\Delta c_0=64/3\); field \(c_0\approx192\pm2\) | yes from higher poles |
| CT (1+2 loop) | \(+0.085120\) | exact | **derived** \(\equiv(C+6)/18+(1-\zeta_2)/24\) | yes \(7L/48-L^2/8\) |
| **Total** | \(\approx+0.05137\) | \(\sim0.002\) | | |
| **Target \(1/24\)** | \(+0.041667\) | **derived** \(c_Y\) | Martin \(\beta=3/2\), \(-\beta/(2\cdot3)\) | |
| **GAP** | \(\approx+0.0097\) | \(\sim0.002\) | | |

## Priority re-derivations (this audit)

1. **Y target \(1/24\)** — verified. Contractions give \(\beta_\Lambda=(3/2)h^2 z s_Y/K\); with \(s_Y=-a_Y R/6\) and pole \(=-\beta/(2\cdot3)\), \(c_Y=1/24\). Lower-loop checks \(1/6\) and \(-1/24\) reproduce.
2. **CT constant** — verified independently from \(Z_z\) one-loop mass insertion plus two-loop sunset with \((1-2\varepsilon)/36\). Algebraically identical to the packet formula (\(\zeta_2=\pi^2/6\)).
3. **QQ, CCcc, QCc tadpole pole, QFf** — verified / exact as in the table.
4. **QCc tadpole \(O(1)\)** — confirmed **not** a simple-pole (\(1/\varepsilon\)) contribution; cannot move the gap.
5. **FFff(B) soft-linear** — already excluded by sector inclusion-exclusion (prior rigor PR); including it would worsen the gap and spoil the pl\(=0\) overlap story.

## First piece off?

**None.** Every audited piece agrees with the ledger inside an error budget \(\sim0.002\), while the gap is \(\approx+0.0097\) (\(\gtrsim 3\sigma\) in that budget).

## Status

The \(\sim+0.01\) discrepancy versus \(1/24\) **stands** under the present graph+CT list and overlap-only B finite content. That is a real residual: it needs a new idea (missing operator/graph, a scheme convention beyond the checks above, or a derived change to B’s non-overlap finite content that does not reintroduce a pl shift), not a retune of an existing number.
