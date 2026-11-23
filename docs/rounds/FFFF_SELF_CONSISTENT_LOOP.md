# Self-consistent loop: \(v\leftrightarrow G_{\mathrm{ind}}\leftrightarrow m_s\) (new-idea rung)

Runnable: `analysis/ffz/ffz_self_consistent_loop.py`.  
Stacked on BG-field/R* (`FFFF_BGFIELD_RSTAR.md`) and condensate rungs #16–#17.  
**No tuning to \(1/24\).** Derived vs assumed labeled.

## Loop dictionary

| Step | Equation | Label |
| --- | --- | --- |
| Contact | \(c_{\mathrm{eff}}=c_{\mathrm{ledger}}+K v^2/(6 h^2 z)\) | #16; **assumed** as the finite \(R\) shift |
| Induced Planck | \(M_P^2=16\pi\,c_{\mathrm{eff}}\,h^2 z/K\) | **assumed** normalization from \(|C_R|\) |
| Gravity mediation | \(m_s^2=-|F|^2/M_P^2\) (tachyonic \(X\) direction) | **assumed** hidden-sector / sign |
| Soft gap | \(v^2/z=-2(1+m_s^2/z)/(h K)\) if \(r<-1\) else \(0\) | #17 T2; **derived** given soft \(r\) |
| Map algebra | fixed-point iteration of the above | **derived** |

\(c_{\mathrm{ledger}}=+0.051317\) (R* total). Gap target \(=+0.009651\).

### F-closure (no external scale)

| Scheme | \(F\) | Free param? |
| --- | --- | --- |
| **F_A** | \(\varphi=F/z\) free | yes — dimensionless \(\varphi\) |
| **F_B** | \(F\sim z\) \(\Rightarrow\varphi=1\) | no continuous free param |
| **F_C** | \(F=\|y\|v^2\) (condensate closes \(F\)) | none |

## (1) Existence / uniqueness

| Scheme | Nontrivial solution? | Unique / stable? |
| --- | --- | --- |
| F_A | yes when \(\varphi\) large enough that \(r(0)<-1\) | attractor unique for scanned \(\varphi\) at \(h=1\) |
| F_B \(\varphi=1\) | **yes** | unique attractor at \(h=1\) (trivial unstable) |
| F_C | **no** — only \(v=0\) ( \(\varphi\propto x\) never drives \(r<-1\) ) | trivial stable |

## (2) Free parameters

| Scheme | Remaining free param | Natural vs tuned |
| --- | --- | --- |
| F_A | \(\varphi=F/z\) (dimensionless) | O(1) natural for \(F\sim z\); **ratio\(\approx1\) needs \(\varphi\approx0.136\) at \(h=1\)** — tuned |
| F_B | none (continuous) | \(\varphi=1\) is a natural discrete choice |
| F_C | none | only \(v=0\) |

## (3) Predicted \(v^2/z\) and \(\delta_R\) vs gap \(+0.00965\pm0.0024\)

### F_B (\(\varphi=1\), natural closed scale)

| \(h\) | \(v^2/z\) | \(\delta_R\) | **ratio \(\delta_R/\mathrm{GAP}\)** | \(v^2/z\) / target |
| ---: | ---: | ---: | ---: | ---: |
| 1 | \(3.194\times10^{-2}\) | \(0.841\) | \(\boldsymbol{87.1}\) | \(86.3\) |
| 0.1 | \(7.496\times10^{-2}\) | \(197.3\) | \(\boldsymbol{2.04\times10^4}\) | \(2.02\times10^4\) |

Target absorber: \(v^2/z=3.703\times10^{-4}h^2\).

### F_C
\(v^2/z=0\), \(\delta_R=0\), **ratio \(=0\)**.

### F_A
Can approach ratio \(\sim1\) only by tuning \(\varphi\approx0.136\) (\(h=1\): best ratio \(\approx0.83\) in the scan).

## (4) Higher poles

Triple/double cancelation **unaffected** (finite \(\xi v^2\) + IR soft mass). **Derived** from #16/#18.

## Plain statement

A gravity-mediation self-consistent loop with \(F\sim z\) **does** admit a nontrivial fixed point, but it **overshoots** the \(\sim+0.01\) Newton gap by \(\sim10^2\) (\(h=1\)) to \(\sim10^4\) (\(h=0.1\)). Fully closing \(F\) through the condensate yields only \(v=0\). Leaving \(\varphi\) free restores tuning. This whole-system dynamics rung does **not** explain the gap without tuning or a huge mismatch.
