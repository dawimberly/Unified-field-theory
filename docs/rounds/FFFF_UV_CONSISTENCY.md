# Consistency constraints on V survivors

Runnable: `analysis/ffz/ffz_uv_consistency.py` + test. Stacked on #49.  
**No tuning.**

## Checks

| Constraint | SO(10)/E6 | Complete Proca→EW | Label |
| --- | --- | --- | --- |
| Unitarity near \(M_{\mathrm{Pl}}\) | GUT Proca OK-ish | Planck Proca **marginal/fail** | sketch **derived**; \(\alpha\) **assumed** |
| Proton decay \(\tau\gtrsim 2\times10^{34}\) y | rough OK / marginal | N/A classic X/Y | formula **assumed**; bound **literature** |
| Gauge unification | **motivating** | **spoiled** | lore |
| IR massless γ/g | OK (residual SM below \(M_{\mathrm{GUT}}\)) | **ruled out** (QCD/QED) | IR requirement **derived** |

## IR crossover (quantified)

- Gluons must be massless/confining for \(\mu\lesssim\Lambda_{\mathrm{QCD}}\sim0.3\,\mathrm{GeV}\).
- Photon: \(m_\gamma\lesssim10^{-18}\,\mathrm{eV}\) (**literature**).
- Complete Proca→EW sets crossover \(\sim v_{\mathrm{EW}}\) → **conflicts**.

## Scorecard

| Mechanism | Flips \(S\) (V)? | Passes W? | Overall |
| --- | --- | --- | --- |
| SO(10)/E6 long broken Proca | yes | tentatively | **best remaining (ASSUMED GUT)** |
| Complete Proca→EW | yes | **no** (IR) | **ruled out** |
| Short Proca / KK / TC-4 | no | — | dead from V |

## Plain

Arithmetic Proca→EW is not a physical option. GUT-style heavy Proca remains the only tentatively consistent flip path — still not derived from SM alone.
