# Stage 52: audit of the 504.5 zs slope (Grundmann et al. 2020)

**Why.** R. Dörner (Frankfurt) asked how we got 504.5 zs from the data in the paper. We re-fetched everything on 2026-11-25 and checked every step for convention errors.

**Sources fetched this session** (quotes in `results/stage52_sources.csv`):
- the arXiv:2010.08298v1 PDF, the author's version of Science 370, 339; supplementary material included;
- Zenodo 10.5281/zenodo.3899920 `Fig3_data.txt`, which is byte-identical to the Stage 31 copy;
- NIST CODATA 2022.

Script: `scripts/stage52_run.py`. Output: `results/stage52.json`.

## What the paper defines (CITED, arXiv v1)
- **Eq. 1–2, p. 5:** "An electron wave is emitted and propagates with the phase velocity vph = Ee/pe for the time τb before a second electron wave is born at a distance R … τb = cos(α0) R/vph (2)". Here R is the full internuclear distance, not R/2.
- **Fig. 1 caption:** "the angle α is enclosed by the electron momentum vector and the molecular axis". α0 is the position of the central fringe.
- **Eq. 3, p. 6–7:** τb′ = cos(β) R/c, "the expected birth time delay was ±247 zs for R = 0.74 Å … linear dependence on cos(β)". So 247 zs is the value at cos β = ±1, which is also the slope per unit cos β. The span from −1 to +1 is 494 zs.
- **p. 6:** cos(α0) is found "for each row of the histogram in Fig. 3 B via a Gaussian fit". The paper does not state the fit window or whether there is a background term.
- **Fig. 3C caption:** "Left axis: cos(α0), right axis: birth time delay calculated using Eq. 2." The paper prints no numerical slope.

## Conversion audit (DERIVED)

| step | check | result |
|---|---|---|
| Eq. 2 factor | R/v_ph, with v_ph = E_e/p_e, E_e = 735 eV, R = 0.74 Å | **9204.3 zs per unit cos α0** (relativistic correction 0.04 %; E_e ±15 eV gives ±1.0 %) |
| Paper's own axes | Fig. 3C ticks digitized at 600 dpi: the 200 zs tick on the right sits at the same height as the 0.022 tick on the left | 200/9204 = 0.0217, which rounds to the printed label 0.022. **Our factor matches their right axis.** |
| Phase vs group velocity | Group velocity v_g = p_e/m = 2 v_ph would give 4602 zs per unit, i.e. half | Eq. 1 is a phase condition, so v_ph is the correct choice. **No error.** |
| R vs R/2 | R/2 would also halve the factor | Eq. 1 uses the centre-to-centre path difference R cos α. **No error.** |
| Slope vs span | We fit τ against cos β, which is a slope per unit cos β | This is directly comparable to 247 zs. The paper's dashed line reaches about 445 zs at cos β = 1. **No error.** |
| Sign | Positive at cos β > 0 in both our fit and the paper's figure | **Consistent.** |

**Coincidence flag.** Our 504.5 zs is 2.04 × R/c. That is numerically close to two unrelated factors of 2: v_g/v_ph = 2, and span/slope = 2 (494 zs). Neither is the cause, because the paper's own points, converted with the paper's own axis, already sit at **1.80 × R/c**.

## Where 504.5 vs ~445 comes from: the fit model, not the conversion

- **Paper's dashed fit**, digitized: slope **445.2 zs**, intercept −0.4 zs. An unweighted fit through the origin to the 15 digitized points gives 445.3 zs.
- **Refits of the public histogram** (weighted slope; ± is the statistical error):

| per-row fit model | slope (zs) | our row centres ÷ paper's |
|---|---|---|
| Gaussian + constant, \|cos α\| ≤ 0.40 (Stage 31; Poisson version gives 504.5 ± 36.2, reproduced) | 500.7 ± 35.7 | 1.10 |
| Gaussian + constant, \|cos α\| ≤ 0.30 | 479.6 ± 34.0 | 1.07 |
| **pure Gaussian, \|cos α\| ≤ 0.25** | **441.0 ± 35.8** (χ² 12.7/13) | **0.99** (rms 29 zs) |
| pure Gaussian, \|cos α\| ≤ 0.30 | 426.7 ± 32.2 | 0.94 |

- The pure Gaussian at |cos α| ≤ 0.25 was picked because it reproduces the paper's points. That makes it a diagnostic of what they likely did, not a physics fit.
- With a constant background term and a wider window, the centres move about 10 % further from 90°. That is the whole 504.5 vs 445 difference.
- **Fit-model systematic** (half the difference): **±32 zs**. Stage 31 left it out; it only counted ±24 zs from the window.

## Verdict
- **504.5 zs is not a convention error.** It is the paper's own conversion applied to the public histogram with a different per-row fit model: Gaussian + background, |cos α| ≤ 0.4.
- **Corrected headline:** with a fit that reproduces the paper's points, the slope is **441 ± 36 (stat) ± 32 (fit model) zs**. That is 0.99 × the paper's dashed fit (445 zs) and 1.79 × R/c (246.8 zs). Stage 31's 504.5 is 1.13 × their fit and 2.04 × R/c.
- **Residual over the paper's nondipole + travel-time model** (red line, about 270 zs, DIGITIZED in Stage 31):
  - **+171 ± 130 zs (1.3σ)** with the matched fit;
  - +235 ± 130 zs (1.8σ) with the Stage 31 fit.
- **Bound on d₁.** The 95 % bound is |d₁| < 425 zs with the matched fit, or **< 489 zs** as a conservative envelope that uses 504.5 plus the fit-model systematic. This replaces the earlier 485 zs (a change of 1 %).
- **Not resolvable from public data:** which per-row model is "right". Choosing between them needs the authors' fit settings or the event-level data.
- **Side note (DERIVED).** If the interference phase is set by the photon frequency (800 eV) rather than E_e, the light-travel prediction becomes 268.7 zs-equivalent. The paper's p. 7 remark about Ref. 20 "if one neglects the ionization potential" points at the same 8 % difference. **Coincidence flag:** 268.7 is close to the red line's ~270, which comes from a different physical effect, the atomic nondipole tilt.
