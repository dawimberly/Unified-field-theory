# Stage 54 — The 8.2 kHz – 4.6 MHz resonant-release gap: cantilevers, trampolines, levitated nanospheres, nanostrings

**Question.** The 8.2 kHz – 4.6 MHz window was flagged open in Stages 33/36/40. Find published experiments that bound a resonant-release add-on rate there, convert them to the standing parametrization (CSL-equivalent, mass-proportional λ per nucleon, m₀ = 1 u, r_C = 10⁻⁷ m), and map exactly what is still uncovered.

**Rules.** Every input was fetched this session (arXiv PDFs → `pdftotext`, NIST CODATA pages, Wikipedia raw wikitext for two densities). URL, location and verbatim quote for each are in `results/stage54_sources.csv`. Labels: **CITED** (from a fetched source), **DERIVED** (computed here from CITED inputs), **ASSUMED** (our stated choice). No parameter was tuned. All measured heating is attributed to the effect (conservative upper limits).

**Files:** `scripts/stage54_run.py`, `results/stage54.json` (rows, validations, all bands Stages 31–54, coverage), `results/stage54_bounds.csv`, `results/stage54_coverage_8k_4p6M.csv`, `results/stage54_sources.csv`.

## 54.1 Formulas (DERIVED; same linearised machinery as Stage 43)
η = λ(4πr_C²)^{3/2}/(m₀²(2π)³)∫d³k e^{−k²r_C²}|F(k)|², heating Γ = ħη/(2M_eff ω) (phonons/s).

1. **Thin plate (Stage 36), unchanged:** Γ/λ = 4πħρr_C²(1−e^{−h²/4r_C²})/(m₀² h ω), valid when every lateral dimension ≫ r_C.
2. **Beam of finite width (new here, needed because nanoladder rails and string defects are only 100–500 nm wide):** the k-integral separates, giving
   Γ/λ = Γ_plate · g(w), g(w) = erf(w/2r_C) − (2r_C/(√π w))(1 − e^{−w²/4r_C²}).
   - g → 1 for w ≫ r_C, so it reduces **exactly** to the Stage 36 plate formula (checked analytically).
   - Closed-form g checked against direct numerical integration: 100 nm 0.27090 vs 0.27090; 500 nm 0.77435 vs 0.77435; 2 µm 0.94358 vs 0.94357.
   - For a structure made of several beams, Γ/λ is a mass-weighted average of local values, and in-phase cross terms between neighbouring elements are ≥ 0. Taking the smallest-g element is therefore conservative.
3. **Rigid sphere (levitated particles):** Γ/λ = ħ m f(R/r_C)/(4 m₀² r_C² ω), f(x) = (12/x⁶)[e^{−x²} − 1 + (x²/2)(e^{−x²} + 1)]. f → 1 gives the Stage 33 trapped-ion (point-mass) formula. This is Zheng et al. 2020 Eq. (3) with S_CSL = ħ²η (two-sided) and Γ = S_CSL/(2mħω).
   - Closed form vs direct k-space integration of the sphere form factor: agreement to ≤ 2×10⁻⁶ at R = 20, 50, 71.5, 1000 nm.

## 54.2 Validation against a published bound (the stated cross-check)
**Zheng et al. 2020 (diamagnetically levitated 1 µm PEG sphere, 4.7 pg).** Inputs CITED: R = 1.0 µm, m = 4.7 pg, 95 % CL excess force noise σ = 3.3×10⁻²⁰ N/√Hz, criterion ħ²η ≤ δS_total.
- Our sphere formula gives **λ ≤ 4.16×10⁻⁷ = 10^−6.38**, vs the published **10^−6.4** (ratio **1.04**): **pass**.
- This also settles the Stage 33 flag. The review's 6.7×10⁻⁷ for Zheng is *not* reproduced from Zheng's own inputs; the paper's 10^−6.4 is. The older stages used the weaker review value, which was conservative and does not need changing.
- Not used as a check: Goldwater et al. 2016 (arXiv:1506.08782). Their Eq. (1) diffusion constant is 2× ours in the point limit, and their Fig. 2 (≈2500 vs ≈350 phonons/s at λ = 10⁻⁸) lies between the two conventions. It is a proposal with figure-read numbers, so it was not used either way.

## 54.3 New bounds (r_C = 10⁻⁷ m, ±0.1 % narrow lines as in Stage 31)
**Tier A** means the total heating, mode temperature or force noise was *measured*. **Tier B** means it is bath-implied (k_BT/ħQ with T not measured), as Goryachev was in Stage 42.

| f | λ_max (s⁻¹) | device / source | heating used (s⁻¹) | Γ/λ | tier, labels |
|---|---|---|---|---|---|
| **25.22 kHz** | **3.6×10⁻⁶** | diamond nanoladder cantilever, Héritier 2018 | 2.3×10⁵ | 6.4×10¹⁰ | A. Inputs CITED: f, Q 131k–162k, mode T ≈ 140 mK (measured by PSD integration), F_th = 190 (+42) zN, m = 4.1 ± 0.6 pg, t ≈ w = 100–300 nm. Heating DERIVED: max of k_BT/ħQ_lo and F_th(+1σ)²/(4m(−1σ)ωħ). Beam formula with t = w = 100 nm, the minimum of the range (ASSUMED conservative). ρ = 3500 (lower end of the fetched 3.5–3.53 g/cm³) |
| 41.4 kHz | 2.4×10⁻⁵ | SiN trampoline 44 nm, Reinhardt 2016 | 1.2×10⁶ | 5.1×10¹⁰ | A. Q = πfτ_m = 3.7×10⁷ (DERIVED from τ_m = 285 s amplitude ringdown, CITED). T = 293 K. Measured x_rms 165 ± 5 pm vs thermal 161 pm, so heating = thermal × ((165+10)/161)². Tether 2.1 µm, g = 0.95. ρ = 2700 (their value) |
| 77.6 kHz | 2.2×10⁻⁵ | levitated SiO₂ sphere, 4 K, Tebbenjohanns 2021 | 3.8×10⁴ | 1.7×10⁹ | A. Γ_tot/2π = 5.5 ± 0.3 kHz, +2σ. d = 106 ± 5 nm (−1σ used). m = min(“∼1 fg”, ρV) |
| **100.5 kHz** | **5.6×10⁻⁶** | levitated SiO₂ sphere, Jain 2016 | 1.96×10⁴ | 3.5×10⁹ | A. Reheating 19.4 ± 0.1 ×10³ s⁻¹ (+2σ), R = 68.0 nm, ρ = 2200 |
| 104 kHz | 3.8×10⁻⁵ | levitated SiO₂ sphere, Magrini 2021 | 1.2×10⁵ | 3.2×10⁹ | A. Γ_ba + Γ_th = 2π·19.4 kHz (no errors quoted). R = 71.5 nm, m = 2.8 fg |
| 140 kHz | 5.5×10⁻⁵ | SiN tethered membrane 20 nm, Norte 2016 | 4.0×10⁵ | 7.3×10⁹ | B. Q = 9.8×10⁷ CITED. T = 300 K ASSUMED (“room temperature”, not stated). Plate limit ASSUMED (tether width not in fetched text) |
| 148.8 kHz | 1.8×10⁻⁵ | levitated SiO₂ sphere, Jain 2016 | 2.1×10⁴ | 1.2×10⁹ | A. R = 52.7 nm |
| 151.2 kHz | 1.1×10⁻⁵ | levitated SiO₂ sphere, Jain 2016 | 3.0×10⁴ | 2.7×10⁹ | A. R = 71.6 nm |
| 158.8 kHz | 1.7×10⁻⁵ | levitated SiO₂ sphere, Jain 2016 | 3.9×10⁴ | 2.2×10⁹ | A. R = 68.0 nm |
| 305 kHz | 1.4×10⁻⁴ | levitated SiO₂ sphere, Delić 2020 | 1.6×10⁵ | 1.1×10⁹ | A. Γ_x/2π = 20.6 ± 2.3 kHz (+2σ), m = 2.83 fg, d = 143 nm |
| **1.406 MHz** | **4.7×10⁻⁹** | soft-clamped SiN string 20 nm, 46 mK, Gisler 2022 | 2.6 | 5.7×10⁸ | B. Heating = 1/(0.38 s), the authors' thermal decoherence time (= k_BT/ħQ with Q = 2.3×10⁹; our recomputation gives 2.6 s⁻¹). Defect w = 500 nm, g = 0.77. ρ = 2700, the lower of the two fetched SiN values |

**Reading.**
- **First bounds in 8.2–121 kHz.** Before this stage nothing was converted between 8174 Hz and 121 kHz. The best new points are 3.6×10⁻⁶ at 25.2 kHz and 5.6×10⁻⁶ at 100.5 kHz.
- **They are weak.** They sit 2–3 orders above the cantilever (10⁻⁸–10⁻¹⁰) and Yuan-membrane (1.8×10⁻⁸ at 242 kHz) points, and all lie above Adler's 10⁻⁸±² central value. Levitated nanospheres are limited by photon-recoil heating (10⁴–10⁵ s⁻¹). Room-temperature trampolines are limited by k_BT/ħQ ≈ 10⁶ s⁻¹.
- **Gisler 1.406 MHz (4.7×10⁻⁹)** matches Seis 1.486 MHz (2.9×10⁻⁹), but it is tier B: the mode temperature was not measured. Treat it as indicative.

## 54.4 Checked and not converted (honest negatives)
- **Tao 2014** single-crystal diamond cantilever (32.14 kHz, Q = 1.51×10⁶ at 93 mK) would give ≈10⁻⁸ if the mode sat at bath temperature. It was **dropped**: the authors state the thermomechanical temperature measurement “was inaccurate in our experiments due to slight mechanical vibrations introduced by the dilution circuit”. Excess noise was present but unquantified, so no valid upper limit exists.
- **Beccari 2022** strained-Si string (1.46 MHz, “60 quanta-per-second”) was dropped: the string width is not in the fetched text, and g(w) can change Γ/λ by ×3.7 between 100 nm and ≫r_C.
- **Fogliano 2021** SiC nanowire at 11.6 kHz (the only 8–25 kHz candidate found) was dropped: the noise temperature is reported only for the 3.94 kHz wire.
- **Héritier 2018** devices at 9.50 kHz and 13.9 kHz: no mode temperature or heating is reported for them.
- **Dania 2024** Paul-trap nanoparticle: 1.28 kHz, out of band.
- **Quartz BAW / phononic-crystal membranes / AURIGA-NAUTILUS harmonics:** no fetched paper gave a heating number in 8.2 kHz – 4.6 MHz. Bars resonate near 900 Hz, and no measured harmonic in band was found.
- **CSL papers with an in-band bound:** none found. Vinante 2016/2017/2020 are at 3.1, 8.2 and 3.5 kHz. Zheng 2020 is about 10 Hz. Pontin 2020 has no frequency in the fetched text (Stage 33).

## 54.5 Coverage map, 8.2 kHz – 4.6 MHz (2.75 decades)
- **25 narrow lines** now sit in the window: 11 new, plus Stages 31/33/36/40/41. At ±0.1 % they cover **0.019 decades (0.7 %)**.
- **Open stretches with no line at all (> 0.1 decade):**

| gap | width (dec) | comment |
|---|---|---|
| **8.18 – 25.2 kHz** | **0.49** | the largest; no converted experiment |
| 25.2 – 41.4 kHz | 0.21 | |
| 41.4 – 77.5 kHz | 0.27 | |
| 77.7 – 100.4 kHz | 0.11 | |
| **305 – 599 kHz** | **0.29** | |
| 601 – 799 kHz | 0.12 | ion lines only (λ ~ 0.7) on either side |
| **1.80 – 4.60 MHz** | **0.41** | only the ⁴⁰Ca⁺ line (λ = 1.4) at the top edge |

- **At a meaningful level (λ ≤ 10⁻⁷, near Adler's band)**, the points are: 8174 Hz (Vinante, likely instrumental), 242 kHz, 1.406 MHz (tier B), 1.486 MHz and 1.8 MHz. Between them:
  - **8.2 – 242 kHz (1.47 decades) has no point at ≤ 10⁻⁷.**
  - 242 kHz – 1.406 MHz (0.76 dec) and 1.8 – 4.6 MHz (0.41 dec) have none either.
- **At ≤ 10⁻⁵**, 8.2 – 25.2 kHz, 25.2 – 121 kHz and 242 kHz – 1.14 MHz are still empty.

## 54.6 Flags, anomalies, coincidences
- **No excess candidate in band.**
  - Reinhardt's free-running x_rms (165 ± 5 pm) matches thermal (161 pm) to 0.8σ.
  - Jain's measured reheating rates scatter on both sides of recoil theory (20.9 vs 15.5; 29.4 vs 38.2; 19.4 vs 21.5; 38.0 vs 35.8 ×10³ s⁻¹). The authors attribute the scatter to gas and the paraxial approximation.
  - Tebbenjohanns' and Magrini's rates are backaction-dominated by design.
- **Out-of-band note:** Dania 2024 (1.28 kHz) reports total heating 3.3×10⁴ phonon/s, against a gas estimate of 2.6×10³. The authors attribute it to electric-field noise. Not converted here.
- **Coincidence check:** Gisler (1.406 MHz) and Seis (1.486 MHz) agree within ×1.6. That is expected (same material, similar thickness and dilution-refrigerator heating level), not a coincidence of physics.
- **Zheng review-vs-paper discrepancy (Stage 33):** resolved in favour of the paper (§54.2).

## 54.7 Caveats
1. All bounds are on narrow (±0.1 %) lines at each experiment's mode frequency. A release line between them is unconstrained, which is why the decade coverage barely moves.
2. The beam formula (§54.1.2) is new. It is validated by its exact plate limit and by numerical integration, not against a published beam/CSL number. The sphere formula *is* validated against a published bound (ratio 1.04).
3. Nanoladder: the 140 mK / 190 zN / Q values are read as belonging to the 25.22 kHz diamond device in Table I. Internally, F_th from T, Q and m gives 169–176 zN, consistent with 190 (+42/−33). The larger of the two heating estimates is used.
4. Tier B rows (Norte, Gisler) assume the mode is at bath temperature. They are not bounds on measured noise.
5. Levitated-sphere bounds assume a rigid homogeneous sphere and use −1σ size or quoted mass where it matters. Γ/λ ∝ m·f(R/r_C), so smaller masses give conservative results.

## 54.8 Verdict
- The 8.2 kHz – 4.6 MHz window now has eleven more converted lines, including the first in 8.2–121 kHz. Every new tier-A point, however, is at 10⁻⁶–10⁻⁴: far above Adler's band and 2–4 orders above the existing membrane and cantilever points.
- **The physically interesting stretch, 8.2–242 kHz, still has no bound at ≤ 10⁻⁷, and 8.2–25 kHz has nothing converted at all.**
- **Cheapest next step:** a mK cantilever or nanoladder at 10–25 kHz with its mode temperature *measured*. To reach λ ≤ 10⁻⁷ with the nanoladder's Γ/λ ≈ 6.4×10¹⁰ (DERIVED), heating must be ≤ 6.4×10³ s⁻¹, i.e. Q/T ≥ 2×10⁷ K⁻¹ (e.g. Q ≥ 2×10⁶ at 0.1 K). Tao's diamond cantilever (Q/T = 1.6×10⁷ K⁻¹ at 32 kHz) is already close, but its mode temperature was not measured. Alternatively, a mK SiN membrane mode near 20–100 kHz with a measured mode temperature, as in Yuan 2015.
