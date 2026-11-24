# Hyperfields, Branes and a Locked Universe: Ten Numerical Tests of Dan Wimberly's Picture

*Write-up of the Hyperfield project, Stages 1–10. Draft dated 24 Nov 2026. Detailed source: `HYPERFIELD_RESULTS.md`, raw numbers in `results/*.json`, scripts in `scripts/`, figures in `figures/`. Not pushed anywhere.*

## Abstract
Dan's picture: the world is a 13-dimensional "crystalline hyperspherical field". Atomic probability clouds are 3D shadows of higher-dimensional spherical fields (hyperfields), and particles are field states "locked" onto our 3D space.

We turned the picture into ten concrete calculations, with no fitted parameters and every input labelled DERIVED, ASSUMED or CITED.

**What held up:**
- **Hydrogen as a 4D sphere:** hydrogen's orbitals are exactly 4D spherical harmonics. That is Fock's 1935 result, reproduced to machine precision.
- **Fields in different dimensions:** "not every field lives in every dimension" is exactly the peer-reviewed brane-world idea.
- **The lock:** a topological kink locks exactly one chiral fermion onto a 3-brane.
- **The diagonal slide:** it selects the E8 lattice (or the 24-cell) without tuning.
- **Nucleon mass:** nucleons really are about 91% locked field energy.

**What failed:**
- Open extra dimensions ruin hydrogen by 7–34%.
- The E8 lattice torus does not give the Standard Model.
- Standard-Model induced gravity needs a cutoff at the Planck scale.
- The 6 TeV lock transition is invisible to LISA.

**Lock-specific results.** These come from the lock model, but none is a confirmed new discovery:
- a hysteresis floor on the transition speed β/H;
- selection of E8 by the slide-plus-churn rule;
- a strong-lock gravitational-wave signal that BBO or DECIGO could falsify.

## 1. The original idea (Dan)
1. **13D crystalline hyperspherical field.** Spacetime has extra dimensions, organised like a crystal (a lattice) and like a sphere.
2. **Probability clouds as hyperfields.** An electron's orbital in an atom is a 3D view of a symmetric higher-dimensional field.
3. **Diagonal slide.** Push a cube "top-back-left to bottom-front-right through the centre". In maths: glue the half body-diagonal (½,…,½) into the cubic lattice.
4. **Multi-viscosity churn.** The hypersphere churns like a ball of fluids with different viscosities, one per quantum field.
5. **Particles as locked field states.** Matter is field energy locked onto our 3D brane, with hysteresis (a first-order "lock").

## 2. Method and rules
- **No tuning.** No parameter is adjusted to match data. Where a number has to be chosen, it is labelled ASSUMED and scanned.
- **Labels.** DERIVED means computed here. ASSUMED means an input. CITED means taken from literature (verified on the web where marked).
- **Coincidences flagged.** Numerical matches with no mechanism are called out, not used.
- **Honest negatives.** A failed test is reported as a result.
- **Engineering.** Short checkpointed runs, everything reproducible with `.venv/bin/python scripts/<stage>.py`.

## 3. Stage results

### Stage 1: Hydrogen as a 4D sphere
- **Question:** are hydrogen clouds shadows of a 4D sphere?
- **Result (DERIVED):**
  - Fock's map turns momentum space into the 3-sphere S³.
  - The kernel eigenvalues are 2π²/n to 1e-15.
  - Energies are −13.605693, −3.401423, −1.511744, −0.850356 eV.
  - The degeneracy is n².
  - Overlap with the real orbitals is 1 − F ≤ 3e-16 for 1s…4f.
- **Verdict: survives, but no new number.** This is Fock (1935) and the known hidden SO(4) symmetry. The cloud is a stereographic projection in momentum space, not a flat slice.
- Figure: `figures/stage1_slices_vs_actual.png`

### Stage 2: More open dimensions (D = 3…13, 25, 26)
- **Question:** do hydrogen clouds come from a D-dimensional sphere with D > 3?
- **Result (DERIVED):**
  - The scale-free ratio (1S–2S)/(1S–3S) is wrong by −7.1% (D=4) up to −33.9% (D=25). It is measured to about 1e-12.
  - Shell-2 state counts become D+1.
  - The 1s size is 9 a₀ at D=13.
  - With Gauss-law forces there are no bound atoms at all.
  - Photon-visible compact dimensions must have R ≲ 4e-17 m.
- **Verdict: ruled out for open dimensions.** Extra dimensions must be compact, or invisible to the photon. Closest published relative: Bars' two-time physics (13D, signature 11+2).
- Figure: `figures/stage2_spectra_vs_D.png`

### Stage 3: Cost scoping
- **Result:**
  - Helium in any D is cheap: three internal coordinates, about 1 s on a laptop.
  - A full Standard-Model lattice in 13+1D needs 222 PB at L = 8.
  - The gauge coupling has mass dimension −5 there, so there is no continuum limit even in principle.
- **Verdict: scoping only.** A direct 13D field simulation is impossible and ill-defined without a UV completion.

### Stage 4: Brane worlds
- **Question:** what if only some fields live in the extra dimensions?
- **Result:**
  - Hydrogen cannot test gravity-only extra dimensions: the shift is 5.5e-27.
  - **Survivors:**
    - ADD with n ≥ 3 at M_D ≳ 6–9 TeV
    - RS1 with k·r_c ≲ 10.5
    - RS2 with ℓ < 44 µm
    - normal-branch DGP with r_c ≳ 3.5/H₀
  - **Excluded:** ADD n = 1, ADD n = 2 as a hierarchy fix, and self-accelerating DGP.
- **Verdict: survives.** Dan's "not every field in every dimension" is exactly ADD/RS/DGP. No new number.
- Figure: `figures/stage4_add_R_vs_n.png`

### Stage 5: The lock (gravity in the bulk, Standard Model on a kink)
- **Question:** can particles be "locked" field states on a brane?
- **Result (DERIVED):**
  - **The lock:** a kink with winding Q = 1 gives exactly one chiral zero mode, ψ₀ ∝ cosh^{−s}(μy), and the unlocked states are gapped at gv.
  - **Gravity:** it deviates from Newton only below 4e-12 m (n=4) to 5e-16 m (n=9).
  - **Short-range force strength:** the Yukawa strength is quantized, α = (4/3) × (first-shell count), e.g. 24 for Z⁹, 13.33 for S⁹ and 8/3 for E8⊕S¹.
  - **Brane thickness:** below 5.5e-21 m.
  - **Hysteresis:** a φ⁶ wall potential gives a window −1 < a < 1/2.
- **Verdict: survives (consistent, partly derived).** The lattice choice is observable only in principle. The "imperfect lock" needs a gapless bulk.
- Figure: `figures/stage5_geff.png`

### Stage 6: Diagonal slide and churn
- **Question:** does the slide or churn pick a dimension and lattice without tuning?
- **Result (DERIVED):**
  - **The slide gives:**
    - the 24-cell/D4 at d = 4, because |h|² = 1 matches the nearest neighbours;
    - E8 at d = 8 (even unimodular, kissing 240).
    - It fails at d = 9 and d = 13: D₉⁺ and D₁₃⁺ do not close, and Z⁹ ∪ (Z⁹+h) has kissing 18.
  - **Energy minimization:** E8 emerges robustly from minimizing the energy of 240 points on S⁷. The d = 9 analog is not special.
  - **Crystallization:** it is first-order with r_s/r_c = 9/8.
  - **Churn bound:** brane velocity v < 3.7e-8 c (GW170817) makes churn unobservable.
  - **Multi-viscosity layering:** it depends on assumed definitions. Not robust.
- **Verdict: survives, with a twist.** The rule picks E8 (8 extra dimensions), not 9. "13D" would need 3+1+8+1, which is numerology.
- Figures: `figures/stage6_slide_lattices.png`, `figures/stage6_thomson.png`

### Stage 7: Protons, neutrons and "Einstein-leaning" systems
- **Result:**
  - **Locked field energy:** about 91% of the nucleon mass is field energy, about 1% quark masses (lattice QCD).
  - **Proton pressure:**
    - repulsive core of about 2.5e35 Pa at 0.1 fm;
    - sign change at 0.58 fm (lattice) / 0.6 fm (experiment);
    - the gluon layer is the larger one.
  - **Extra-dimension gravity at 1 fm:** E8's α = 320 at 1 fm is 1e-26 to 1e-39 of the strong force.
  - **Equivalence principle:** MICROSCOPE forces the lock to gravitate universally to 5e-12.
  - **Neutron stars:** twin stars need chosen parameters, and on SLy they conflict with 2 M☉ pulsars.
- **Verdict:** the locked-field language fits real QCD, but that is standard physics. Pressure does not make gravity repel. No new number.
- Figures: `figures/stage7_proton_pressure.png`, `figures/stage7_ns_MR.png`

### Stage 8: Sakharov induced gravity
- **Question:** can Standard-Model fields induce 1/G with the right sign and size?
- **Result:**
  - The one-loop N_eff = 1 − 24ξ, so the sign is right only for Higgs ξ < 1/24.
  - The size needs a cutoff Λ ≈ 3–6 M_Pl.
  - The finite (FFZ) version needs about 62+ heavy bosons at 4–75 M̄ and a vacuum-energy tuning of 1e-120.
- **Verdict: ruled out as a bridge, no prediction.** Induced gravity re-expresses 1/G and predicts nothing.
- Figure: `figures/stage8_induced_G.png`

### Stage 9: Does the E8 torus give the Standard Model?
- **Result (DERIVED):**
  - **Symmetry:** a flat E8 torus gives U(1)⁸ and a non-chiral spectrum.
  - **E8 gauge symmetry** needs strings at the self-dual radius.
  - **Couplings:** bulk gauge couplings would be about 1e-36.
  - **KK shells:** degeneracies 240 : 2160 : 6720, spaced 1.9–3.6 GeV, which is unobservable.
  - **Branchings** were checked: 248 → E6×SU(3), SU(5)×SU(5) and SO(10)×SU(4).
- **Verdict: ruled out.** It adds nothing beyond heterotic E8×E8 or Lisi (Distler–Garibaldi no-go).
- Figure: `figures/stage9_e8_torus.png`

### Stage 10: Gravitational waves from the lock transition
- **Question:** can LISA see the first-order lock transition at T* ≈ 6 TeV?
- **Result:**
  - **The model's numbers:** with natural O(1) inputs the φ⁶ lock gives α = 1e-4–0.12 and β/H = 400–9,800. LISA's SNR is then ~0.
  - **Hysteresis floor:** β/H ≥ S_c(4a₁−3), which is about 130 at a₁ = 1.
  - **Stage-5 headline:** 0.035 Hz, which means β/H ≈ 35, with α = 0.1 gives LISA SNR 3.3, below detection.
  - **BBO/DECIGO:** they detect the strong lock (v ≈ 2 T_c, a₁ = 1) at SNR 170–380 / 14–35.
  - **Existing limits:** LIGO-Virgo O3 and NANOGrav constrain nothing here.
- **Verdict: no LISA detection predicted.** BBO or DECIGO could falsify the strong-lock half; the weak half cannot be tested by anything.
- Figure: `figures/stage10_gw.png`

## 4. Summary table
| Stage | Question | Key number | Verdict |
|---|---|---|---|
| 1 | Hydrogen as 4D sphere | 1 − F ≤ 3e-16; E_n = −13.6057/n² eV | Survives (Fock 1935, no new number) |
| 2 | Open D > 3 spheres | ratio off by 7–34% | Ruled out |
| 3 | Compute cost | 222 PB at L=8; non-renormalizable | Scoping |
| 4 | Brane worlds | ADD n ≥ 3, M_D ≳ 6–9 TeV | Survives (known physics) |
| 5 | Kink lock | 1 chiral zero mode; α = (4/3)N₁ | Survives, partly derived |
| 6 | Diagonal slide / churn | selects E8 (240) or D4 (24); not 9 | Survives; picks 8, not 9 |
| 7 | Nucleons, neutron stars | 91% field energy; node at 0.6 fm | No new number |
| 8 | Induced gravity | ξ < 1/24; Λ ≈ 3–6 M_Pl | Ruled out as bridge |
| 9 | E8 torus → SM? | U(1)⁸, non-chiral, g ≈ 1e-36 | Ruled out |
| 10 | GW from lock | LISA SNR ~0; BBO up to ~380 | No LISA signal; BBO-falsifiable |

## 5. What is genuinely Dan's vs already known
**Genuinely Dan's.** These are model-specific statements, not discoveries yet:
1. **Hysteresis floor on the transition speed.** In the φ⁶ lock, β/H ≥ S_c(4a₁ − 3), about 130 at a₁ = 1. The lock prefers fast transitions, which is why LISA misses it.
2. **E8 selection by the diagonal slide plus churn.** A no-tuning rule picks E8 (or D4) and rejects 9 and 13.
3. **A BBO/DECIGO-falsifiable strong lock.** If v ≳ 1.5 T_c at T* ≈ 6 TeV, BBO or DECIGO should see SNR of order 10–400.

**Already known.** The pieces Dan's picture rests on:
- **Hydrogen as a 4D sphere:** Fock 1935, Bargmann 1936, Bars' two-time physics.
- **Brane worlds:** ADD, RS, DGP.
- **Fermion localization on a kink:** Jackiw–Rebbi and Rubakov–Shaposhnikov.
- **E8 optimality on S⁷:** Cohn–Kumar.
- **Proton mass and pressure:** lattice QCD; Burkert et al. 2018.
- **Induced gravity:** Sakharov, FFZ.
- **E8 gauge symmetry:** heterotic strings, Lisi.
- **Gravitational-wave templates:** Caprini et al. 2016.

## 6. Coincidences flagged (not evidence)
- **1/24:** the Higgs-ξ sign flip at 1/24 (Stage 8) equals the repo's FFZ target 1/24. They are unrelated quantities.
- **6 TeV:** this is the ASSUMED ADD M_D choice. The gravitational-wave peak landing in the BBO/DECIGO sweet spot follows from that choice.
- **The "3" in E6×SU(3):** it is the holonomy of a 6D compact space, not the generation count.
- **13 = 3+1+8+1:** writable, not derived.
- **S⁹ degeneracies 10, 54, 210 = the Stage-2 D=9 shell counts:** a mathematical identity.
- **D4 and Bars' SO(4,2):** the same group appearing twice.
- **ADD n=2 radius vs dark-energy length:** a well-known coincidence.

## 7. Unverified items
- **Detector noise curves:** the LISA (Robson–Cornish–Liu), galactic-foreground and DECIGO/BBO (Yagi–Seto) coefficients were coded from memory. The noise minima match the published curves. The ET/CE curves are crude.
- **Lattice values:** kissing numbers and densities of the best known lattices for d = 9–13.
- **Citations recalled from memory, not re-checked:** Cohn–Kumar, the QGP η/s value, Larsen–Wilczek / Kabat / Visser, Rubakov–Shaposhnikov.
- **Neutron-star inputs:** SLy parameters and the NICER constraint values.
- **Collider and Higgs numbers:** the n = 8 collider bound and the LHC Higgs-ξ limit.
- **Stage-10 approximations:** thin-wall nucleation, and the wall-energy efficiency κ_v not corrected for v_w < 1.

## 8. Next steps
1. **BBO/DECIGO:** turn the strong-lock region into a forecast with a full bounce calculation for 1/2 < a₁ < 3/4, and pin the ET/CE curves.
2. **ChatGPT import:** bring Dan's earlier ChatGPT physics conversations into `docs/chatgpt/` with `chatgpt_import/convert_chatgpt_export.py`, to check nothing in the original idea was missed.
3. **Optional:** a predictive version of induced gravity, meaning a symmetry that fixes the heavy spectrum.
