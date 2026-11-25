# Stage 53 — LISA Pathfinder Δg above 30 mHz → frequency-resolved resonant-release bound λ(f)

**Goal.** Earlier stages used LPF only up to 30 mHz (Helou 0.7–20 mHz ≤ 2.96×10⁻⁸ s⁻¹; Stage 33's interpolated 2018 curve to 30 mHz). Here I take the archive's own Δg time series and convert the measured differential acceleration above 30 mHz into λ(f): a CSL-equivalent rate per nucleon at r_C = 100 nm, as in Stages 31/33. The aim is to push into the empty 0.19–9.3 Hz window (Stage 36).

**Rule.** Every input was fetched this session; the quotes and locations are in `results/stage53_sources.csv`. Labels: **CITED**, **DERIVED**, **ASSUMED**. Nothing was tuned; the binning, the validity cut and the CL were fixed before looking at bounds.

**Files:**
- `scripts/stage53_lpf_highfreq.py` and `scripts/stage53_fig.py`;
- `results/stage53.json`, with segment metadata, archive OIDs and sha256 of every raw file;
- `results/stage53_lpf_lambda_vs_f.csv`, the envelope;
- `results/stage53_lpf_lambda_indicative_above_cut.csv`;
- `results/stage53_lpf_bands_per_segment.csv`;
- `results/stage53_lpf_asd_L4-66.csv`, the measured ASD of the Feb 2017 run;
- `figures/stage53_lpf_lambda_vs_f.png`;
- `data/stage53/lpfsa_deltag_catalogue.json`, the archive's Δg catalogue as fetched, and `data/stage53/LPFSA_Deltag_Release_Notes_doc11850.pdf`.

Raw XML (2.1 GB) is not committed. Re-fetch it by OID; the sha256 of each file is in `stage53.json`.

## 53.1 Where the data live, and how they were reached

- **The COSMOS page is only a portal.** https://www.cosmos.esa.int/web/lisa-pathfinder-archive is empty even when signed in. The archive itself is the **LISA Pathfinder Legacy Archive** at **https://lpf.esac.esa.int/lpfsa/**, a JavaScript (GWT) app. Its `js/config.js` names the backend `/lpfsa-sl/` (CITED, fetched).
- **The backend answers plain HTTPS without login.** `check-user` returns `"USER__USERNAME":"anonymous"` with `ALLOWED_TYPES` including `ANALYSIS_OBJECT` and `DOCUMENT`.
  - Metadata query: `metadata-action?RESOURCE_CLASS=ANALYSIS_OBJECT&...&QUERY=ANALYSIS_OBJECT.ORIGIN='D'`. It returns **337 Δg objects**: Levels L0–L3, 80 segments each, and **17 L4 segments**.
  - Download: `data-action?ANALYSIS_OBJECT.ANALYSIS_OBJECT_OID=<oid>`. It returned the LTPDA XML files with HTTP 200, and every size matches the catalogue's FILESIZE.
  - **The ESA account was not needed.** It was never used, and the password file was not read.
- **Level 4 is what we want.** The archive's *Δg Release Notes* (DOCUMENT_OID 11850, CITED) say L4 = "Δg_noise": Level 3 "cleaned by force glitches … computed exclusively for uninterrupted data-series with the best maximum authority configuration, when no experiments was running". All 17 L4 segments were used: 2.1 GB, Mar 2016 – Jun 2017, 0.86–18.5 d each, fs = 10 Hz (CITED, file header).
- **Not used:** the DOI pages at esdcdoi.esac.esa.int (JS-only, returned HTTP 500/404 to curl) and the HEASARC DRS FITS mirror (DRS telemetry; HEASARC says DRS "did not optimize the delta-g performance").

## 53.2 Method (Stage 33, unchanged) and validation

- **Conversion.** Eq. 13 of Carlesso–Ferialdi–Bassi 2018 in the thick-slab limit gives S_F = 8πħ²λρ²A r_C²/m₀² (two-sided, per cube). The cube is 4.6 cm and 1.928 kg (CITED). The one-sided Δg PSD = 2 × 2 (two independent TMs) × S_F/m². **All measured Δg power is attributed to the effect** (DERIVED); no known noise is subtracted.
- **Formula check (DERIVED).** 5.2 fm s⁻²/√Hz → λ = 2.986×10⁻⁸ vs Helou's 2.96×10⁻⁸ (**ratio 1.009**), the same as Stage 33.
- **Data check (DERIVED vs CITED).**
  - Apr 2016 segment L4-9 (the era of the 2016 PRL): mean ASD over 2–20 mHz = **5.73 fm**, vs the CITED 5.2 ± 0.1. That is ×1.10 in ASD, giving λ = 3.6×10⁻⁸ vs Helou's 2.96×10⁻⁸.
  - Feb 2017 run L4-66 (13.33 d, the Armano 2018 Fig. 1 run, "~13 day"): **1.98 fm**, vs the CITED 1.74 ± 0.05 "above 2 mHz" (×1.14).
  - Our band means include the 1/f rise near 2 mHz and the full segment span, so they come out slightly high, i.e. conservative. The pipeline reproduces the published levels to 10–14 % in ASD.
- **Readout check.** Above 30 mHz, Armano 2018 says Δg is "limited … by the interferometer readout noise of S_x^½ = (34.8 ± 0.3) fm/√Hz … effective Δg ASD of S_x^½(2πf)²" (CITED).
  - In the data, the displacement-equivalent noise ASD/(2πf)² at 0.1–0.3 Hz is 40–52 fm/√Hz for 15 of the 17 segments. L4-9 is 40.2 fm (×1.16 of 34.8).
  - L4-66 and L4-74 (2017) sit at 81–89 fm. That is higher readout noise in those runs, not a signal: the spectrum is f²-shaped. **So above ~50 mHz every bound below is a readout-noise limit, and it rises as f⁴ in λ.**
- **PSD estimate (ASSUMED choices):**
  - Welch with 2000 s Blackman–Harris segments, 50 % overlap, linear detrend (Δf = 0.5 mHz);
  - 1/20-decade bands;
  - 95 % CL from χ² with ≈2K dof per bin (adjacent bins treated as half-correlated);
  - "line" = the worst single 0.5 mHz bin in the band, at 95 % CL.
- **Upper frequency limit (ASSUMED rule, fixed in advance).**
  - The 10 Hz data roll off well below Nyquist (the on-board low-pass is not documented in the fetched handbook).
  - Ground processing has no filter: the history shows only `diff` (5POINT), `split`, `interpmissing` and similar. The 5-point second-derivative stencil is exact to 0.2 % at 1 Hz (DERIVED).
  - A band is kept while ASD/(2πf)² ≥ 1/√2 (−3 dB) of its 0.1–0.3 Hz median.
  - The filter is common to all segments, so the **robust cut is the lowest per-segment −3 dB point: 1.12 Hz**. Per-segment cut points range from 1.12 to 2.24 Hz. Bands above 1.12 Hz are listed only as *indicative*.

## 53.3 Result: λ(f) above 30 mHz (DERIVED, 95 % CL, tightest of 17 segments per band)

| band (Hz) | λ_max band (s⁻¹) | narrow line, 0.5 mHz bin (s⁻¹) | segment | median of 17 segments |
|---|---|---|---|---|
| 0.0317–0.0356 | **1.28×10⁻⁸** (best above 30 mHz) | 1.39×10⁻⁸ | L4-76 (Jun 2017) | 3.1×10⁻⁸ |
| 0.0423–0.0474 | 2.8×10⁻⁸ | 3.6×10⁻⁸ | L4-76 | — |
| 0.089–0.100 | 2.5×10⁻⁷ | 3.0×10⁻⁷ | L4-9 | — |
| 0.100–0.112 | 3.9×10⁻⁷ | 5.1×10⁻⁷ | L4-13 | 5.2×10⁻⁷ |
| 0.178–0.200 | 3.5×10⁻⁶ | 4.6×10⁻⁶ | L4-9 | — |
| 0.200–0.224 (first band fully in window) | 5.6×10⁻⁶ | 7.4×10⁻⁶ | L4-9 | 7.7×10⁻⁶ |
| 0.283–0.317 | 2.2×10⁻⁵ | 2.9×10⁻⁵ | L4-9 | 3.2×10⁻⁵ |
| 0.448–0.502 | 1.3×10⁻⁴ | 1.7×10⁻⁴ | L4-9 | — |
| 0.893–1.002 | 1.4×10⁻³ | 3.2×10⁻³ (1.000 Hz line) | L4-9 | — |
| 1.002–1.125 (robust top) | **1.9×10⁻³** | 2.4×10⁻³ | L4-9 | 2.7×10⁻³ |
| 1.12–2.24 (indicative only) | 3.3×10⁻³ → 0.15 | — | 10 → 1 segments | — |

- **Shape.** From ~50 mHz up, λ_max ∝ f⁴ (readout noise × (2πf)², squared). Each decade in frequency costs about 4 decades in λ.
- **Did we reach 0.19–9.3 Hz? Yes, formally: 0.19 → 1.12 Hz, i.e. 0.77 of the window's 1.69 decades.** These are the first bounds with a stated frequency in that window in this audit. They are weak: 3.5×10⁻⁶ at 0.19 Hz rising to 1.9×10⁻³ at 1.1 Hz. That is ≈1.8 orders (at 0.19 Hz) to ≈3.5 orders (at 1.1 Hz) weaker than the bounds flanking the window: Kovachy atoms ≤ 0.19 Hz at 5.1×10⁻⁸, and Zheng at 9.3 Hz, 6.7×10⁻⁷ (both repo Stage 33 values, *not re-fetched this session*). The 1.12–9.3 Hz part (0.92 dec) is still empty.
- **Below 0.19 Hz.** LPF is tighter than the atom band-average (5.1×10⁻⁸) only below ≈ 55 mHz. Its advantage is that it is frequency-resolved, while the atom bound is an average over its kernel.
- **Consistency with Stage 33.** Stage 33's ASSUMED log-log interpolation gave 1.7×10⁻⁸ at 30 mHz; the measured value is 1.28×10⁻⁸ at 31.6–35.6 mHz (ratio 0.75). Over 2–20 mHz the best segment (L4-74, 1.93 fm) gives 4.1×10⁻⁹, consistent with Stage 33's 3.3×10⁻⁹ from the CITED 1.74 fm.

## 53.4 Caveats (honest)

1. **Pure noise attribution.** The bounds are on the total measured Δg. Above 50 mHz this is known readout noise, so the bound is an instrument limit, not a physics floor. Subtracting the 34.8 fm model is not allowed by the no-tuning rule, and is not reliable, because the readout level varies ×2 between runs.
2. **Filter near the top.** Within the kept bands the on-board low-pass may already attenuate by up to −3 dB in power. So λ near 0.8–1.1 Hz could be **optimistic by up to ×2**. The segment that sets the envelope there (L4-9) has the earliest roll-off. The 17-segment median at 1.0–1.1 Hz is 2.7×10⁻³ and the worst segment is 1.5×10⁻². Above 1.12 Hz nothing is claimed.
3. **Band vs line.** These are band averages (1/20 dec ≈ 12 % wide). A release line narrower than 0.5 mHz is bounded more weakly, by (0.5 mHz)/(line width), as in Stage 33.
4. **Model assumptions.** The fluctuations on the two TMs are taken as independent (37.6 cm ≫ r_C, CITED distance). The thick-slab limit is valid since r_C ≪ 4.6 cm. m₀ = 1 u is ASSUMED.
5. **Envelope = tightest of 17.** This is legitimate: each segment is an independent 95 % upper limit on a time-constant effect, and with ≥ 10³ averages per band the statistical spread is ≈ 1–5 %. But it is not 95 % joint coverage in the strict sense. The 17-segment medians are given as the robust alternative.
6. **Unverified for x12 itself.** The history checked shows the 5POINT `diff` only for the L2 cross-talk regressors. That the same stencil was used for Δẍ is ASSUMED.

## 53.5 Coincidences and flags

- **1.000 Hz line.** In the Feb 2017 run (L4-66) there is a line exactly at 1.0000 Hz, 25.6× the local median PSD in one 0.5 mHz bin (plus window leakage at ±0.5 mHz). An exactly-integer frequency points to the 1 Hz on-board cycle, so it is **instrumental**. It is left in, which makes the 0.89–1.00 Hz line bound 3.2×10⁻³ rather than 1.4×10⁻³. It is not a release candidate.
- **69.5–70 mHz bump.** About 2× the local median in L4-66. It lies inside the 20–100 mHz "sensing crosstalk" region that Level 2 subtracts (Release Notes, CITED), so it is probably a residual. It is also left in.
- **Not coincidences, agreement checks:** the formula reproduces Helou to 1.009; the data reproduce the published levels to ×1.10–1.14; and the measured 30 mHz value agrees with the Stage 33 interpolation to within ×0.75.

## 53.6 Verdict

LPF now gives a measured, frequency-resolved bound from 30 mHz to 1.12 Hz. The best value is **1.3×10⁻⁸ s⁻¹ at 32–36 mHz**, then **2.5–3.9×10⁻⁷ near 0.1 Hz**, **2.2×10⁻⁵ at 0.3 Hz** and **1.9×10⁻³ at 1.0–1.1 Hz**. This covers the lower **0.77 decade of the 0.19–9.3 Hz window, but only weakly**, because LPF's free-fall readout noise rises as f². A resonant release between 0.19 and 1.1 Hz is still allowed up to 10⁻⁵–10⁻³ s⁻¹, and **1.12–9.3 Hz remains empty**. To do better there one needs a sensor whose force noise is flat (not f²-rising) at 1–10 Hz: torsion pendulums such as TOBA, or levitated masses. Those are the Stage 36 data-request targets.
