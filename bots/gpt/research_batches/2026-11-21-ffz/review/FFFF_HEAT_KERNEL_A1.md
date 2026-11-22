# Heat-kernel / Seeley–DeWitt rung for induced \(R\)

Runnable: `ffz_heat_kernel_a1.py`.  
Ladder entry: `docs/PROVEN_METHODS_LADDER.md` (rung 4).  
**No tuning to \(1/24\).** Independent of the three-loop Schwinger/projective integrand.

## Naming clash

| Symbol | Meaning here |
| --- | --- |
| Seeley–DeWitt \(a_1\) | First heat-kernel coefficient of a Laplace-type operator |
| Projective \(a_1\) (FFff(B) ledger) | \(1/\varepsilon\) coefficient of the stripped projective integral — **different object** |

## What a full three-loop HK would require (not done)

An interacting three-loop Seeley–DeWitt expansion of the Wess–Zumino theory on a curved background — all \(a_n\) with Yukawa vertices, overlapping subdivergences, and metric insertions — does **not** exist in the packet and is not a feasible review cross-check. The designed three-loop path remains the Wick/Schwinger ledger. This rung therefore does the **feasible subset**: verify the lower-loop HK inputs that fix the \(1/24\) target and the CT, then compare to the standing gap.

## Feasible subset — results

| Check | Result |
| --- | --- |
| Scalar \(a_1=(1/6-\xi)R\) at \(\xi=1/6\) | \(0\) (no tadpole linear \(R\)) |
| Dirac coincident \(R/6-R/4\) | \(-1/12\) (matches packet spinor check) |
| 1-loop Y-linear \(F_R\) | \(+1/6\) in units \(z a_Y/(K\varepsilon)\) |
| 2-loop Y-linear \(F_R\) | \(-1/24\) in units \(h z a_Y/(K\varepsilon)\) |
| CT higher poles | triple \(0\), double \(L/12-5/48\) |
| CT simple \(L=0\) | \(+0.085120\) |
| \(c_Y\) from Martin \(\beta_\Lambda=(3/2)\) + \(s_Y=-a_Y R/6\) | \(1/24\) |
| Y-linear higher poles | \(1/36\), \(-1/18\), \(+1/24\) |

All verified in the runnable. Soft-mass dictionary \(s_Y=-a_Y R/6\) is the same \(R/6\) that appears in the conformal heat-kernel coefficient.

## Comparison to \(1/24\) and the gap

This rung does **not** recompute FFff/CcFf/QFf. It confirms the targets those graphs are tested against. Against the standing ledger total:

| Quantity | Value |
| --- | ---: |
| Ledger simple total (DRED+BPHZ / overlap-fixed B) | \(\approx +0.05141\) |
| Target \(1/24\) (HK+β confirmed) | \(+0.041667\) |
| **GAP** | \(\boldsymbol{+0.00975}\) (unchanged) |

Prior methods: DREG \(\sim+0.0097\), finite-\(\varepsilon\) \(+0.00954\pm0.0036\), DRED+BPHZ \(+0.00975\pm0.0024\).

## Plain statement

The heat-kernel / Seeley–DeWitt route, restricted to what can be done without a full three-loop interacting HK computation, **confirms** \(c_Y=1/24\) and the CT and **does not close** the \(\sim+0.01\) simple-pole gap. Next ladder rung: superspace / supergraph method.
