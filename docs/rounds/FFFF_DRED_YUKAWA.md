# DRED redo of Yukawa/fermion Newton graphs

Runnable: `analysis/ffz/ffz_dred_yukawa_ledger.py`.  
**No tuning to \(1/24\) or \(-1/288\). No sector-restricted overlap rule as an input.**

## Setup

| Ingredient | Choice |
| --- | --- |
| Regulator | DRED: momenta/Schwinger/orbital Laplacian in \(d=4-2\varepsilon\); Dirac/spinor algebra at \(D=4\) |
| ε-scalars | orbital↔spin mismatch on fermion lines (\(d_{\mathrm{spin}}=4\) vs \(d_{\mathrm{orb}}\)) |
| Yukawa weight | global \(w_f=w_c=1\) from \([y]=\varepsilon=\lambda_F\) (literal Lagrangian) |
| Subdivergences | BPHZ / R-operation: \(R=G-K_1G-K_2G+K_{12}G\) on soft UV forests |
| Two-loop CTs | DRED-consistent; finite DRED↔DREG \(Z\) absorbed into coupling redefinition at fixed renormalized \(h,z\) |

Compare to DREG with the same global+BPHZ assembly (and to raw global without BPHZ).

## (1) Triple / double poles

Coupling-only soft residues (DRED = DREG at this order — see §3):

| | ll | pl | \(W/U\) |
| --- | ---: | ---: | ---: |
| coup_only | \(-16\) | \(-48\) | \(+64/3\) |
| global (field+coup) | \(-32\) | \(-96\) | \(+128/3\) |

BPHZ on the **extra** global coupling copy keeps soft-soft only (\(\Delta a_1=\mathrm{ll}=-16\), \(\Delta c_0=64/3\)) and subtracts single-soft (\(\Delta\mathrm{pl}=0\), raw pl \(-48\) removed by \(K_i\)).

| Assembly | triple | double | simple | proj. double resid |
| --- | ---: | ---: | ---: | ---: |
| DRED global RAW (no BPHZ) | \(+0.020833\) | \(+0.017352\) | \(+0.054179\) | \(+2.08\times10^{-2}\) **fail** |
| **DRED global + BPHZ** | \(+0.020833\) | \(\boldsymbol{-0.003481}\) | \(\boldsymbol{-0.032626}\) | \(\boldsymbol{-8.9\times10^{-6}}\) **OK** |
| DREG global + BPHZ | \(+0.020833\) | \(-0.003481\) | \(-0.032626\) | \(-8.9\times10^{-6}\) |
| DREG field-only | \(+0.020833\) | \(-0.010426\) | \(-0.052302\) | \(-7.0\times10^{-3}\) |

A+B triple residual vs \(1/36\): \(0\).  
Projected double residual (DRED+BPHZ): \(\approx0\).

**Cancelation without a sector-restricted rule as input: YES**, but only after BPHZ/R. Global weight alone (RAW) does **not** cancel the double pole. After R, the remainder is numerically the old soft-soft \((-16,\,64/3)\) content — now derived as forest output, not assumed as an overlap patch.

## (2) Simple pole vs \(1/24\)

| Piece | Value |
| --- | ---: |
| QQ | \(+0.062500\) |
| CCcc | \(-0.069444\) |
| QCc doubled | \(-0.016497\) |
| QCc tadpole | \(-0.027778\) |
| QFf | \(-0.074662\) |
| CcFf | \(-0.024800\) |
| FFff(A)×4 | \(+0.1496\pm0.002\) |
| FFff(B) DRED+BPHZ | \(-0.032626\) |
| CT | \(+0.085120\) |
| **TOTAL** | \(\boldsymbol{+0.05141\pm0.0024}\) |
| target \(1/24\) | \(+0.041667\) |
| **GAP** | \(\boldsymbol{+0.00975\pm0.0024}\) |

**DRED does not close the gap.** The residual remains \(\sim+0.01\).

## (3) Evanescent / ε-scalar vs DREG (graph by graph)

| Graph | DRED vs DREG at kept order |
| --- | --- |
| **FFff(B)** | Soft-soft/pl residues of the vertex/coupling \(\varepsilon\)-numerator: \(\Delta\mathrm{ll}=\Delta\mathrm{pl}=\Delta(W/U)=0\). Reason: Isserlis \(T'(d)\) is multiplied by \((4-d)/2\to0\) at \(d=4\) in \(\partial_d A_{\mathrm{vert}}\), so orbital↔spin mismatch drops out of these leading residues. Finite (beyond \(W/U\)) content can still differ; not needed for the pole test. |
| **FFff(A)** | ×4 Wick already matched 4D Dirac traces; no further leading-pole shift. Same global+BPHZ policy as B. |
| **CcFf** | \(d_{\mathrm{spin}}\to4\) vs \(d\to4\) agree on double \(1/144\); ε-scalar \(O(\varepsilon)\) faces inside existing \(\pm0.001\) simple QMC error. |
| **QFf** | Finite \(d\to4\) limit of \(D_x\) identical; no ε-scalar in this topology’s \(1/\varepsilon\). |
| **CT (fermion/Yukawa)** | One-loop SUSY poles universal. Two-loop finite DRED↔DREG \(Z\) absorbed into coupling redefinition; pole CT at fixed renormalized \(h,z\) unchanged (`COUNTERTERM_BUDGET.md` algebra). |
| **QQ / CCcc / QCc** | Pure scalar — DRED=DREG. |

## Plain statement

DRED + literal global Yukawa weight does **not** by itself fix the Newton double-pole cancelation. BPHZ/R is required; after R the soft-soft remainder is \((-16,\,64/3)\). The projected simple gap \(\sim+0.01\) vs \(1/24\) **stands** — DRED does not remove it.
