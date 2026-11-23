# Spectral total \(1/G_{\mathrm{ind}}\) (heavy + hadronic)

Runnable: `analysis/ffz/ffz_spectral_gind.py` + test. Stacked on #40.  
**No tuning.**

## Question

Is the **integrated** (not IR-only) induced \(1/G\) positive once the matched hadronic band is included? Near \(M_{\mathrm{Pl}}^2\)?

## Method

Adler / Sakharov spectral sketch:

\[
\frac1{16\pi G_{\mathrm{ind}}}\ \propto\ \sum_i a_1(\mu_i)\,\Delta\log\mu_i^2
\]

| Piece | \(a_1\) | Label |
| --- | ---: | --- |
| UV SM (\(\xi_H=0\)) | \(-43/12\) | **derived** content |
| EW Proca | \(-1/12\) | **derived** |
| top→QCD free | \(-13/4\) | **derived** |
| Hadronic IR | \(1/3\pm0.29\) | matched #40 (**assumed** ρ) |

## Result

| Quantity | Value |
| --- | --- |
| Log-weighted \(S\) (central) | **negative** (\(\sim -300\)) |
| Entire hadronic band flips sign? | **No** |
| Near \(M_{\mathrm{Pl}}^2\) without cutoff? | **No** |
| \(\Lambda_{\mathrm{UV}}\) for \(M_{\mathrm{Pl}}\) | \(\sim M_{\mathrm{Pl}}\) (**assumed**) |

## Coincidence flag

\(\Lambda_{\mathrm{UV}}\sim M_{\mathrm{Pl}}\) is the usual Sakharov assumption — not derived here.

## Plain

Honest negative: confinement-band positivity does **not** make the integrated \(1/G_{\mathrm{ind}}\) positive or Planckian. UV Kabats still win the spectral integral.
