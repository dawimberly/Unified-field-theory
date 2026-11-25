# Summary for Dan — Stages 30 to 38 (overnight run, Nov 24–25, 2026)

Plain-language recap. Full details, tables and sources are in [HYPERFIELD_RESULTS.md](HYPERFIELD_RESULTS.md); every stage has its own `results/stageNN_sources.csv`, with a quote from each source fetched in-session.

**How to read the numbers.** For "release rate" bounds we use one common yardstick. λ is the extra random-kick rate per nucleon, expressed like the CSL collapse model with a 100 nm smearing length, in s⁻¹. Smaller means a tighter limit. For scale:
- Adler's proposed value is 10⁻⁸ (give or take 100×).
- The original GRW value is 10⁻¹⁶.

---

## The short version
1. **Big G (Stage 30): no.** The ~500 ppm spread between G measurements is not caused by quantum-gravity noise; that noise is ≥ 3.9 orders too small. The weak "bigger test mass → bigger G" trend comes from three old experiments and is tangled up with measurement method.
2. **Your do-si-do idea (Stage 31) survives in two places:**
   - **Angle-dependent release timing** is the strongest survivor. The H₂ data allow an extra term up to ~489 zeptoseconds (conservative; 425 zs with the matched fit), and show a 1.3σ hint that is not significant (Stage 52, which replaced the earlier 485 zs and 1.8σ).
   - **A resonant ("lets go at one frequency") release** is the second.
3. **Can the existing H₂ data settle it (Stage 32)? Probably not.** Re-binning gives about 1.4σ. New data, or the authors' event file, are needed.
4. **Hunting for the resonant release frequency (Stages 33–37).** We now have limits over ~10.8 of the 25.7 decades from 20 µHz to 4 MeV.
   - **The biggest hole (THz → soft X-ray) is mostly untestable with ordinary matter** in this model (Stage 35). This is physics, not a missing experiment.
   - **The one standing anomaly, the 8174 Hz cantilever excess**, can still be a do-si-do line, but only a fairly narrow one (10–470 Hz wide). Part of it looks environmental.
5. **Other do-si-do elements (Stage 38).** The angle element also survives on a second system (helium), with no anomaly. The test there is ~150× weaker than H₂.

---

## Stage by stage

### Stage 30 — Big G vs test-mass size
- 16 CODATA G values compared against test-mass volume.
- The trend is +69 ppm per decade of volume (p = 0.035, about 0.18 after look-elsewhere). It disappears without three old dumbbell experiments (p = 0.38).
- Quantum-gravity (CQ) noise at its maximum allowed size shifts G by ≤ 0.036 ppm, ≥ 3.9 orders short of 500 ppm.
- **Honest negative.** The real pattern is method: dynamic vs compensation vs suspended setups differ at p = 0.001. That is already known.

### Stage 31 — The do-si-do model, all 7 elements
- **Angle (probability laws, lab orientation):** tightly constrained (10⁻² down to 10⁻¹⁸).
- **Angle × release time:** our re-analysis of the Grundmann 2020 H₂ "birth delay" data gives a slope of **441 ± 35.8 (stat) ± 32 (fit model) zs** against the light-travel 246.8 (Stage 52, see `STAGE52_grundmann_504_audit.md`). The residual over the paper's nondipole model is **+171 ± 130 zs (1.3σ)**, so |d₁| < 425 zs (matched fit), or < 489 zs as a conservative bound. **Open.** *(Corrected 2026-11-25. The earlier 504.5 ± 36 zs, +235 ± 128 zs (1.8σ) and < 485 zs came from a Gaussian+constant per-row fit at |cosα| ≤ 0.4. A pure-Gaussian fit reproduces the paper's own points and gives 441.)*
- **Gate threshold:** ≤ 24.86 µeV, with no amplitude gate.
- **Amplitude:** < 2×10⁻²⁰ eV.
- **Openness spread:** intrinsic < 10⁻³⁰.
- **Hysteresis:** |η| < 9×10⁻⁵ (CMB FIRAS).
- **Resonant release:** at that point bounded in only 3 bands.

### Stage 32 — Could re-binning the existing H₂ data decide it?
- A joint angle × bond-length split gives σ(d₁) ≈ 169 zs. The expected significance is ~1.4σ (25 % chance of reaching 2σ).
- It is worth doing because it is immune to the big systematics, but it is not decisive.
- The analysis script is ready (`scripts/stage32_analyze_joint.py`).

### Stage 33 — First frequency map of resonant release
- 21 sources were converted to one yardstick.
- Best numbers:
  - LISA Pathfinder: 3×10⁻⁹ (mHz);
  - cantilevers: 2×10⁻¹⁰ (3.5 kHz);
  - bulk heating: 3×10⁻¹¹ (5.6–13 GHz);
  - X-rays: 10⁻¹²–10⁻¹⁵ (keV).
- Coverage of the Stage 31 gap went from 0.001 to 1.5 decades.

### Stage 34 — LIGO
- Public LIGO noise converted to force noise over 10 Hz – 5 kHz. It reproduces the published LIGO bound to 1 %.
- The O3 detectors are ~4× tighter: **2.2×10⁻⁶ at 35 Hz**.
- LIGO gets weak above ~1 kHz (λ > 1).
- **Honest negative:** germanium 100 keV – 1 MeV was dropped (no usable count table).

### Stage 35 — The THz-to-soft-X-ray hole (13 GHz – 14.5 keV)
- **Key physics finding.** A neutral atom is smaller than the 100 nm kick scale, so the kicks push its electrons and nucleus together. For photon energies below ~12.4 keV, their radiation **cancels**.
  - Heat (phonons) can only carry the kicks up to ~37 GHz before an exponential wall.
  - So most of this window is **untestable with ordinary lab matter**, including qubits, which need ~100 GHz phonons (suppressed by 10⁻¹⁰⁰).
  - **This corrects a Stage 33 suggestion** (cryogenic calorimeters would not help).
- **New limits:**
  - 13–26 GHz: 4×10⁻¹⁰ → 10⁻⁶ (heating tail).
  - 1–8 keV: ~40–55, from the cosmic X-ray background versus free electrons in intergalactic gas (weak).
  - 60–630 GHz: 3×10⁸ (CMB, useless).
- **The cryogenic "low-energy excess"** (an unexplained signal in dark-matter detectors) is **not** a do-si-do candidate:
  - it cannot be produced in this model;
  - it differs from detector to detector;
  - it fades with time and jumps after warm-ups.

### Stage 36 — MHz and GHz resonators; the 0.2–9 Hz hole
- **New formula** for thin vibrating membranes, checked against our two validated limits.
- **Best new number: 2.9×10⁻⁹ at 1.486 MHz** (Copenhagen soft-clamped membrane, Seis 2022). That is inside Adler's range and ~9 orders better than the trapped ions that were the only MHz data before.
- Also:
  - 1.2×10⁻⁶ at 1.14 MHz (Rossi 2018);
  - 1.6×10⁻⁶ at 5.05 GHz (ETH quantum acoustic resonator, 150× better than the Stage 33 HBAR point).
- **Honest negative:** these are single frequencies, not bands. The 0.19–9.3 Hz window is still empty: no published in-window bound was found.

### Stage 37 — The 8174 Hz cantilever excess
- **2017:** unexplained extra force noise of 1.87 aN²/Hz. As a white noise it reads as λ = 10⁻⁷·⁷.
- **2020:** the same hardware, retuned to 3532.7 Hz with a heavier layered mass, saw nothing (λ < 2×10⁻¹⁰).
- A do-si-do line at 8174 Hz must therefore drop by ≥ 100× within 4.6 kHz. That means a Lorentzian width (half-width) of 10–470 Hz, or a Gaussian σ ≤ 1530 Hz.
- **No other existing experiment can see it.** LIGO is 6.4 orders short. Harmonics are untestable.
- **Leans instrumental:**
  - the excess rose 2.8σ with the cryocooler running;
  - the 2020 data sit 1.9σ below what a non-mass-coupled floor of the same size would give.
- **Verdict: open, not supported.**

### Stage 38 — Other elements with new public data
- **Angle element, P₂ term, on the helium atom** (ETH RABBITT 2016): |d₂| ≲ 15 as. Consistent with no effect, ~150× weaker than H₂.
- **Dropped:**
  - nuclear-decay "spread" tests (no numerical limit fetchable);
  - gate and hysteresis (no new data).

---

## Coverage of the resonant-release frequency axis (after Stage 36)
| range | covered |
|---|---|
| whole span 20 µHz – 3.8 MeV (25.66 decades) | **10.84 decades** (formal, any λ) |
| Stage-31 gap 20 mHz – 19 keV (20.36 decades) | 6.54 decades formal; 3.92 with λ ≤ 1 s⁻¹ |
| 13 GHz – 14.5 keV (8.43 decades) | 2.37 formal; **0.45 with λ ≤ 1** (mostly untestable, Stage 35) |
| 4.6 MHz – 5.6 GHz (3.08 decades) | single lines only (5.05 GHz) |
| 0.19 – 9.3 Hz (1.69 decades) | **0** |

Tightest limits by region:
- LPF 3×10⁻⁹ (mHz);
- cantilever 2×10⁻¹⁰ (3.5 kHz);
- **membrane 2.9×10⁻⁹ (1.5 MHz, new)**;
- bulk heating 3×10⁻¹¹ (5.6–13 GHz);
- X-rays 10⁻¹²–10⁻¹⁵ (keV).

## Honest negatives (things that did not work or were wrong)
- Big G: no new-physics signal (Stage 30).
- Existing H₂ data cannot decide the angle question alone (Stage 32).
- Germanium 100 keV – 1 MeV: no usable public tables (Stage 34).
- **Stage 33's calorimeter idea was wrong:** neutral matter does not radiate in this model below ~12 keV (Stage 35).
- THz → soft X-ray is structurally untestable in the lab for this model (Stage 35).
- The cryogenic low-energy excess is not this effect (Stage 35).
- 0.19–9.3 Hz: no usable published data (Stage 36).
- The 8174 Hz excess shows signs of an environmental contribution (Stage 37).
- Nuclear decay "spread" limits: not quantified in any fetchable source (Stage 38).

## Open questions
1. Does the 1.3σ H₂ angle residual (Stage 52; 1.8σ with the older Stage 31 fit) grow or vanish with bond-length-split data? (Needs the Frankfurt event file.)
2. Is the 8174 Hz floor mass-coupled? Re-measure near 8.17 kHz with a different mass geometry.
3. Does the do-si-do picture *predict* a release frequency f₀? Without one, the resonant element can only be cornered, not tested. Candidates so far: 8174 Hz (Stage 37 window), or none.
4. Would a non-mass-proportional coupling, or a sub-nm smearing length, change Stage 35? Both are outside the current parametrization; sub-nm is excluded by others.

## Possible data requests (none sent; drafts on request)
| who | what | why |
|---|---|---|
| Grundmann / Dörner (Frankfurt) | event-level H₂ data with angle and bond length (KER) | Stage 31/32 test. *Already requested earlier; the reply-check routine is watching for an answer.* |
| Vinante (Southampton/Trento) | raw 2017 (8100–8240 Hz) and 2020 (3515–3550 Hz) spectra; any later runs near 8.17 kHz | Stage 37: is the floor flat, and is it mass-coupled? |
| Schliesser group (Copenhagen) | membrane occupation/heating vs temperature for other modes and devices (Seis 2022, Rossi 2018) | turn the 1.5 MHz point into a MHz band (Stage 36) |
| Chu group (ETH) | populations and T₁ for all HBAR modes, not just 404 (12.5 MHz comb around 5 GHz) | a line scan across GHz (Stage 36) |
| Majorana / LNGS Ge groups | binned 100–1000 keV spectrum plus detection-efficiency simulation | the 100 keV – 1 MeV hole (Stage 34) |
| Keller group (ETH) | tabulated angle-resolved He delays with errors, both hemispheres | a d₁ and a tighter d₂ (Stage 38) |
| Oosterkamp/Fuchs (Leiden); Ando group (TOBA, Tokyo) | force/torque noise of low-frequency (sub-10 Hz) modes, with geometry | the empty 0.19–9.3 Hz window (Stage 36) |
| — (public, no request needed) | ESA LISA Pathfinder archive data above 30 mHz | possibly fills 0.03–1 Hz; the published text does not give the range |

## Commits (branch `hyperfield-stages`)
| stage | commit |
|---|---|
| 30 | a669cc6, 1cd6153 (verification) |
| 31 | d94a1f2 |
| 32 | b83dd49 |
| 33 | 599e4f7 |
| 34 | c3f83a1 |
| 35 | dab73df |
| 36 | 2f6a701 |
| 37 | c948f0c |
| 38 | c6c4001 |

Figures:
- `figures/stage33_frequency_coverage.png`
- `figures/stage34_frequency_coverage.png`
- `figures/stage35_frequency_coverage.png`
- `figures/stage36_frequency_coverage.png`
- `figures/stage37_8174Hz.png`
