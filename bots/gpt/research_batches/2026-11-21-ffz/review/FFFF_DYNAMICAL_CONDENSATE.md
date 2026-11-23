# Dynamical condensate / BCS–NJL gap equation (new-idea rung)

Runnable: `ffz_dynamical_condensate.py`.  
Stacked on the ξv² contact analysis (`FFFF_CONDENSATE_R.md`).  
**No tuning to \(1/24\).** Derived vs assumed labeled below.

Target from tree \(\xi v^2\) absorber: \(v^2/z = 3.703\times10^{-4}\,h^2\) (\(|v|/\sqrt{z}=0.0192\) at \(h=1\)).

## (1) Gap-equation setup

Model: \(W=mXY+yX^3/6\), \(h=|y|^2/K\), \(z=|m|^2\).

| Truncation | Content | SUSY breaking | Label |
| --- | --- | --- | --- |
| **T1** | Exact F-term stationarity \(F_X=F_Y=0\) | none | **derived** |
| **T2** | Tree \(V_F+m_s^2|X|^2\) on \(Y=0\): \(V=(z+m_s^2)|X|^2+(|y|^2/4)|X|^4\) | soft \(m_s^2\) **assumed** | min **derived** |
| **T3** | Mean-field NJL/BCS: \(G=\kappa y^2/z\) (\(\kappa=1\)), cutoff \(\Lambda\); \(v=\Delta/|y|\) | cutoff breaks SUSY **assumed** | **assumed** |

NJL gap (T3):
\[
1=\frac{G}{4\pi^2}\Bigl[\Lambda^2-\Delta^2\ln\frac{\Lambda^2+\Delta^2}{\Delta^2}\Bigr],\qquad G=\frac{hK}{z}.
\]

## (2) Does a nontrivial solution exist?

| Truncation | Nontrivial \(v\)? | Critical condition |
| --- | --- | --- |
| T1 (exact SUSY) | **No** — \(F_Y=mX=0\Rightarrow X=0\); non-renormalization forbids loop shift of the F-flat locus | none (SUSY forbids) |
| T2 (soft tree) | **Yes** iff \(r=m_s^2/z<-1\) (tachyonic soft) | soft ratio is an **input** |
| T3 (NJL) | **Yes** iff \(\Lambda/\sqrt{z}>2\pi/\sqrt{hK}\) | cutoff ratio is an **input** |

SUSY does **forbid** a dynamical condensate in the unbroken theory. Soft/cutoff scales must be supplied.

## (3) Predicted \(v^2/z\) vs target \(3.703\times10^{-4}\,h^2\)

### T1
\(v^2/z=0\) ⇒ **ratio predicted/target \(=0\)**.

### T2
\[
\frac{v^2}{z}=\frac{-2(1+r)}{hK}\quad(r<-1).
\]
No unique prediction — depends on free \(r\).

| \(r\) | \(h\) | \(v^2/z\) | target | **ratio** |
| ---: | ---: | ---: | ---: | ---: |
| \(-1.5\) | 1 | \(6.33\times10^{-3}\) | \(3.70\times10^{-4}\) | \(\boldsymbol{17.1}\) |
| \(-2.0\) | 1 | \(1.27\times10^{-2}\) | \(3.70\times10^{-4}\) | \(\boldsymbol{34.2}\) |
| \(-1.5\) | 0.1 | \(6.33\times10^{-2}\) | \(3.70\times10^{-6}\) | \(\boldsymbol{1.71\times10^4}\) |

To **hit** the target: \(r_*=-1.029\) at \(h=1\) (\(-1.000029\) at \(h=0.1\)) — tuned soft input.

### T3 (sample)
| \(h\) | \(\Lambda/\sqrt{z}\) | \(v^2/z\) | target | **ratio** |
| ---: | ---: | ---: | ---: | ---: |
| 1 | \(0.505\) (\(\approx\lambda_c\times1.01\)) | \(5.62\times10^{-6}\) | \(3.70\times10^{-4}\) | \(\boldsymbol{0.0152}\) |
| 1 | \(0.75\) | \(1.83\times10^{-3}\) | \(3.70\times10^{-4}\) | \(\boldsymbol{4.95}\) |
| 1 | \(1.5\) | \(1.42\times10^{-2}\) | \(3.70\times10^{-4}\) | \(\boldsymbol{38.4}\) |
| 0.1 | \(1.60\) | \(5.62\times10^{-4}\) | \(3.70\times10^{-6}\) | \(\boldsymbol{152}\) |

Again: ratio spans orders of magnitude with free \(\Lambda/\sqrt{z}\); no unique prediction.

## (4) Free parameters — no tuning claim

| Truncation | Remaining free parameter |
| --- | --- |
| T1 | — (only \(v=0\)) |
| T2 | soft ratio \(m_s^2/z\) |
| T3 | cutoff \(\Lambda/\sqrt{z}\) and \(\kappa\) |

Making the condensate “dynamical” **relocates** the tuning into the SUSY-breaking / UV input; it does **not** predict the Newton gap-absorbing \(v\).

## Plain statement

Exact SUSY forbids a dynamical \(\langle X\rangle\). Soft or NJL truncations can produce nontrivial \(v\), but always with a free soft/cutoff parameter. Predicted/target ratios are \(0\) (SUSY) or \(O(10^{-2}\)–\(10^{6})\) depending on that input — not fixed to 1. This rung does not close the \(\sim+0.01\) gap without tuning.
