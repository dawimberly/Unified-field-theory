# Background-field + R* audit (ladder rung #6)

Runnable: `analysis/ffz/ffz_bgfield_rstar.py`.  
**No tuning to \(1/24\).** Feasible subset — not a full diagrammatic R* re-integration of every QMC graph.

## Scope

| Done | Not done |
| --- | --- |
| Identify packet metric split as BG-field R extraction | Dynamical graviton / quantum \(h_{\mu\nu}\) loops |
| Identify soft-sector Taylor IE as \(R^*\) IR dual \(\tilde R\) | Scratch R* recomputation of all CcFf/A QMC integrands |
| Soft-linear \(1/\varepsilon\) test | |
| Per-piece IR/UV checklist; re-derive closed-form poles | |

## Background field

\[
g=\bar g+h,\quad \sigma=-\kappa x^2,\quad \bar R=4d(d-1)\kappa,
\]
with matter on \(\bar g\) and \(R\) read via \(1/[4d(d-1)]\). Already the packet. No dynamical graviton (**assumed** continuation).

## R* (Chetyrkin–Tkachov / Smirnov)

\[
R^*=R_{\mathrm{UV}}\circ\tilde R_{\mathrm{IR}}.
\]

On soft projective axes, \(\tilde R\) is the sector Taylor forest
\(G-G(x,0)-G(0,y)+G(0,0)\) — the same IR operation used for FFff(B) overlap (soft-soft only). **Label:** identification derived; equivalence sector-IE \(\equiv\tilde R\) assumed as the feasible dictionary.

### Soft-linear

Single-soft with \(O(\varepsilon)\) numerator: linear term \(\to\) finite as \(\varepsilon\to0\), **no** \(1/\varepsilon\). Excluding soft-linear \(\approx-105.304\) from overlap \(c_0\) is R*-correct, not a missed IR CT.

### Tadpole

QCc tadpole pole \(-1/36\) is UV mass; \(O(1)\) is finite matching. At \(\xi=1/6\), \(a_1=0\) — no extra linear-\(R\) IR. R* does not promote \(O(1)\to1/\varepsilon\).

## Per-piece L=0 simple (all \(\Delta=0\))

| Piece | Value | Err | Missed IR CT? |
| --- | ---: | ---: | --- |
| QQ | \(+0.062500\) | exact | no |
| CCcc | \(-0.069444\) | exact | no |
| QCc doubled | \(-0.016497\) | \(\sim10^{-11}\) | no |
| QCc tadpole | \(-0.027778\) | exact | no |
| QFf | \(-0.074662\) | \(\sim10^{-11}\) | no |
| CcFf | \(-0.024800\) | \(0.001\) | no (within QMC) |
| FFff(A)×4 | \(+0.149587\) | \(0.002\) | no (within err) |
| FFff(B) | \(-0.032709\) | \(0.00087\) | no (R̃=soft-soft) |
| CT | \(+0.085120\) | exact | no |
| **TOTAL** | \(\boldsymbol{+0.051317\pm0.0024}\) | | |
| target \(1/24\) | \(+0.041667\) | | |
| **GAP** | \(\boldsymbol{+0.009651\pm0.0024}\) | | |

Higher poles: triple/double still cancel. **Pieces changed: none.**

## Plain statement

R* / background-field feasible subset **confirms** the existing ledger (soft-soft kept, soft-linear out, tadpole pole-only). It does **not** close the \(\sim+0.01\) gap. Full R* re-integration of every QMC graph remains out of scope.

With rung #6 done (partial/negative on the gap), the established peer-reviewed ladder is exhausted without closing the residual — a genuinely new idea is required beyond the methods tried.
