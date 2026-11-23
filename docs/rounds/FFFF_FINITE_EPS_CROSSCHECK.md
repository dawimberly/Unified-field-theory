# Finite-ε Laurent cross-check of the projected simple pole

Runnable: `analysis/ffz/ffz_finite_eps_crosscheck.py`. **No tuning to \(1/24\).**  
pySecDec is not installed; this is a pure-Python finite-ε / sector implementation.

## Method (different from the PR #10 pipeline)

Previous ledger path: endpoint/region-split → projective \((a_0,a_1,c_0)\) → convert.

**This check:** for each class, build the graph contribution \(G(\varepsilon)\) in units \(h^2 z/K\) at finite \(\varepsilon>0\) (including \(\Gamma\) / measure prefactors at \(L=0\)), then extract

\[
G(\varepsilon)=\frac{t}{\varepsilon^3}+\frac{d}{\varepsilon^2}+\frac{s}{\varepsilon}+O(1)
\]

by staged polynomial fits / high-precision series. **No reuse of stored ledger \(c_0/a_1\) values** — every coefficient is recomputed in the runnable.

| Class | Finite-ε implementation |
| --- | --- |
| CT | closed sunset/\(Z_z\) form |
| QCc tadpole | \(T\times J'\) short-distance |
| QQ | Wick + curvature normalization |
| CCcc | Wick×residue/\([4d(d-1)]\) expanded at finite \(d=4-2\varepsilon\) |
| QFf | Schwinger \(J(\varepsilon)\to J_{\mathrm{fin}}\to K^2 D_x\) (no \(j_1\)@0) |
| QCc doubled | full Schwinger integrand QMC → graph(\(\varepsilon\)) → Laurent |
| CcFf | QMC graph(\(\varepsilon\)) → Laurent (no \(1/144\) subtract-first) |
| FFff(A)×4 | regulated QMC → fresh \(c_0\); mpmath \(\mathrm{pref}(\varepsilon)I(\varepsilon)\) |
| FFff(B) | recomputed sector \(a_0/a_1/c_0\)+forest; \(\mathrm{pref}(\varepsilon)I(\varepsilon)\) series |

## Graph-by-graph simple vs PR #10 ledger

| Piece | New simple | Err | Ledger | \(\Delta\) | Flag |
| --- | ---: | ---: | ---: | ---: | --- |
| QQ | \(+0.062500\) | exact | \(+0.062500\) | \(0\) | ok |
| CCcc | \(-0.069444\) | \(8\times10^{-4}\) | \(-0.069444\) | \(0\) | ok |
| QCc doubled | \(-0.016483\) | \(8\times10^{-5}\) | \(-0.016497\) | \(+1.5\times10^{-5}\) | ok |
| QCc tadpole | \(-0.027778\) | exact | \(-0.027778\) | \(0\) | ok |
| QFf | \(-0.074643\) | \(1\times10^{-5}\) | \(-0.074662\) | \(+2\times10^{-5}\) | ok |
| CcFf | \(-0.024605\) | \(2\times10^{-4}\) | \(-0.024786\) | \(+1.8\times10^{-4}\) | ok |
| FFff(A)×4 | \(+0.149170\) | \(0.003\) | \(+0.149587\) | \(-4\times10^{-4}\) | ok |
| FFff(B) | \(-0.032626\) | \(0.0018\) | \(-0.032709\) | \(+8\times10^{-5}\) | ok |
| CT | \(+0.085120\) | exact | \(+0.085120\) | \(0\) | ok |

**Disagreements \(>\) combined error: none.**

Higher-pole sanity (same run): QCc double \(\approx 1/72\), CcFf double \(\approx 1/144\), A×4 triple \(\approx 1/144\), B triple \(\approx 1/48\), B double \(\approx -1/288\).

## New total vs \(1/24\)

| Quantity | Value |
| --- | ---: |
| New projected total | \(+0.05121\pm0.0036\) |
| Target \(1/24\) | \(+0.041667\) |
| **New GAP** | \(\boldsymbol{+0.00954\pm0.0036}\) |
| PR #10 gap | \(+0.0098\pm0.0021\) |

**Gap survives:** yes — central value still \(\sim+0.01\), consistent with PR #10 within combined error. The independent finite-ε method does **not** remove the residual.
