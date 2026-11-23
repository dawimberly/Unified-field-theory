# GUT spectra SU(5) / SO(10): UV Higgsing vs \(a_1\) sign

Runnable: `analysis/ffz/ffz_gut_spectrum.py` + test. Stacked on #45.  
**No tuning.**

## Question

Can a UV “confinement-like” mechanism (all extra GUT gauge bosons massive / Proca above \(M_{\mathrm{GUT}}\)) flip the UV \(a_1\) / integrated \(S\)?

## Local \(a_1\)

| Spectrum | Gauge | Fermions | Scalars (assumed) | Net | Sign |
| --- | ---: | ---: | ---: | ---: | --- |
| SU(5) unbroken (24 massless) | \(-16\) | \(+15/4\) | \(+5/3\) (\(5_H\)) | \(-16+15/4+5/3\) | **−** |
| SU(5)→SM (12 Proca X/Y) | \(-2\) | \(+15/4\) | \(+5/3\) | \(\boldsymbol{+41/12}\) | **+** |
| SO(10) unbroken (45 massless) | \(-30\) | \(+4\) | \(+10/3\) (\(10_H\)) | large − | **−** |
| SO(10)→SM (33 Proca) | \(+17/2\) | \(+4\) | \(+10/3\) | **+** | **+** |

Gauge arithmetic **derived**; Higgs multiplets **assumed** minimal.

## Spectral integral (conservative: gauge+fermion only)

| Case | SU(5) flip? | SO(10) flip? |
| --- | --- | --- |
| \(M_{\mathrm{GUT}}=\mu_{\mathrm{UV}}\) | **yes** (\(S\sim+69\)) | **yes** |
| Planck \(\to M_{\mathrm{GUT}}\) | **no** (\(S\sim-105\)) | **yes** (more Proca) |

With minimal Higgs ASSUMED, SU(5) can flip even the wide bin — do not treat as robust.

**Coincidence flag:** \(M_{\mathrm{GUT}}\sim\mu_{\mathrm{UV}}\) often assumed — structural.

## Plain

GUT X/Y Proca help locally. Unbroken GUT makes UV worse. Flip is **not** scale-independent — SU(5) fails the Planck→GUT stress test without optimistic Higgs.
