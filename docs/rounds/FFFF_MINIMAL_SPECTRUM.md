# Minimal broken-gauge spectrum for \(S>0\)

Runnable: `analysis/ffz/ffz_minimal_spectrum.py` + test. Stacked on #51.  
**No tuning.** Fields-only: massive vector **field modes** (Proca); IR massless photon/gluon modes.

## Setup

- Above \(M_*\): \(12+N_V\) massless gauge modes (Kabats).
- Below \(M_*\): 12 massless SM + \(N_V\) massive vector modes.
- \(\Lambda_{\mathrm{UV}}=M_{\mathrm{Pl}}\), \(M_*=10^{16}\,\mathrm{GeV}\) (**assumed** scales).

## Result

| \(N_V\) | \(S\) | Note |
| ---: | ---: | --- |
| 0 | \(-\) | SM-like |
| **15** | \(+\) | **minimal** (arith **derived** at fixed scales) |
| 12 (SU(5)-like) | \(\sim-54\) | **not** enough at these scales |
| 33 (SO(10)-like) | \(\sim+405\) | **allowed**, not selected |
| 66 (E6-like) | \(\sim+1127\) | **allowed**, not selected |

## Constraints

| Constraint | Status |
| --- | --- |
| IR massless γ/g modes | requires residual 12; \(M_*\gg\Lambda_{\mathrm{QCD}}\) (**derived**) |
| Proton decay | \(M_*\gtrsim10^{15}\,\mathrm{GeV}\) (**literature**) |
| SM coupling unification | SM alone **does not** unify (1-loop); extra thresholds **assumed** |

## Plain

\(N_V^{\min}=15\) is derived arithmetic under scale assumptions. SO(10)/E6 are **allowed**, not uniquely **selected**. No simple group is derived from \(S>0\) alone.
