# Why \(S\sim-265\) is negative while IR \(a_1\) is positive

Runnable: `analysis/ffz/ffz_uv_decompose.py` + test. Stacked on #44.  
**No tuning.**

## Decomposition (#41 integral)

| Bin | \(a_1\) | weight \(2\ln(\mu_hi/\mu_lo)\) | contrib |
| --- | ---: | ---: | ---: |
| UV SM (\(10^{16}\to v\)) | \(-43/12\) | \(\sim63\) | \(\sim-225\) (**dominant**) |
| EW→top | \(-1/12\) | small | \(\sim0\) |
| top→QCD | \(-13/4\) | \(\sim13\) | \(\sim-41\) |
| Hadronic IR | \(+1/3\) | \(\sim2\) | \(\sim+0.7\) |

UV SM split: gauge Kabats \(12\times(-2/3)=-8\) dominate fermions \(+15/4\) and scalars \(+2/3\).

## Minimal flip (physically motivated, not dial-tuning)

| Option | \(\Delta a_1\) UV | Flips \(S\)? | Label |
| --- | ---: | --- | --- |
| A. All 12 vectors → Proca | \(+14\) | **Yes** | arithmetic **derived**; spectrum **assumed** |
| B. Conformal Higgs \(\xi=1/6\) | \(-2/3\) | **No** (worse) | **derived** |
| C. Adler \(\langle TT\rangle\) sum rule | ? | unknown | **assumed** UV completion |
| D. Asymptotic safety screening | ? | unknown | **assumed** |
| E. Extra minimal scalars | need \(\sim O(10^2)\) | yes if added | count **derived**; content **assumed** BSM |

Need \(\Delta a_1^{\mathrm{UV}}\gtrsim +4.2\) in the UV bin alone.

## Plain

IR positivity is negligible in the log-weighted integral. Gauge Kabats above \(\Lambda_{\mathrm{QCD}}\) own the minus sign. Cleanest arithmetic flip: complete vector Higgsing.
