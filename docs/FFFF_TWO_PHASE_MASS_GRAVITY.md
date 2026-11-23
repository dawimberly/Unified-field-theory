# Two-phase mass/gravity (Cassini-safe)

Runnable: `bots/gpt/research_batches/2026-11-21-ffz/review/ffz_two_phase.py`  
Test: `.../test_ffz_two_phase.py` (hooked in `scripts/run_checks.py`).  
Stacked on three-phase (#28).  
**No tuning.** DERIVED vs ASSUMED labeled.

Dan’s two-phase version: light = E/B wave (massless; no gravity phase); matter carries a locked mass/gravity pair. Nonzero \(T^\mu_\mu\) switches the gravity phase on.

## (1) Formalization

| Tool | What it gives | Mass↔gravity lock? |
| --- | --- | --- |
| de Broglie / Compton clock \(e^{-imc^2 t/\hbar}\) | Phase rotation at \(\omega_C=mc^2/\hbar\) | **No** (single U(1) clock) |
| Dirac mass \(m(\bar\psi_L\psi_R+\mathrm{h.c.})\) | Genuine L↔R locked pair; zitterbewegung \(\sim\hbar/(2mc)\) | **No** (chiral lock ≠ gravity) |
| Higgs vev / condensate phase | Masses from EWSB; \(\xi_H R\|H\|^2\) optional | \(\xi_H\) **assumed**; not a gravity wave phase |
| Dilaton / trace coupling | Sources \(\propto T^\mu_\mu\); off for classical EM | **Derived** for that channel only |
| EP \(m_i=m_g\) | Geometric in GR (universal free fall) | **Not** derived as two-phase locking (**assumed** reading) |

**Honest:** Dirac L–R and Compton clocks are real two-component / phase structures for *mass*. They do **not** derive a locked mass↔gravity wave pair. EP equality is metric geometry, not a phase-lock condition.

## (2) Avoid the #28 clash

| Sector | Coupling | Photons |
| --- | --- | --- |
| Metric \(g_{\mu\nu}\) | Universal (full \(T^{\mu\nu}\)) | Light bending & Shapiro **OK**; Cassini \(\lvert\gamma-1\rvert<2.3\times10^{-5}\) intact |
| Extra dilaton \(\varphi\) | \(\alpha_T\,\varphi\,T^\mu_\mu\) | Off for classical EM (**derived**); on for massive matter |

\(\alpha_T\) is a **free** fifth-force strength (**assumed**).

| Bound | Value | Role |
| --- | ---: | --- |
| MICROSCOPE \(\lvert\eta\rvert\) | \(\lesssim1\times10^{-15}\) | Composition-dependent dilaton |
| Eöt-Wash \(\lvert\eta\rvert\) | \(\lesssim10^{-13}\) | Torsion-balance EP |
| Cassini → \(\omega_{\mathrm{BD}}\) | \(\gtrsim4.3\times10^{4}\) | Universal scalar admixture \(\lesssim1/(2\omega)\sim10^{-5}\) |

None of these closes the FFZ \(h^2 z/K\) gap (**derived**: different channel).

## (3) Induced-gravity \(a_1\) sign

| Count | \(a_1\) | Sign OK? | Label |
| --- | ---: | --- | --- |
| UV SM (12 massless vectors) | \(-43/12\) | **no** | **derived** (#27) |
| IR: W±,Z as Proca (\(+\tfrac12\)), γ+gluons massless | \(-1/12\) | **no** | arithmetic **derived**; IR for Sakharov **assumed** wrong venue |
| All 12 Proca (fantasy, not SM) | \(+125/12\) | yes | **not SM** |
| #28 shared-Kabat \(\mathbb{Z}_3\) | \(+53/12\) | yes | **assumed** (unchanged by two-phase) |

Dirac L–R locking does **not** remove photon/gluon Kabat contacts (Donnelly–Wall: massless edge modes). Proca removes *its own* contact (**derived**) — already used in the IR row — and still leaves \(a_1=-1/12\).

**No derived two-phase boundary flip** without tuning or dropping the SM massless gauge sector.

## (4) Two-phase vs three-phase

| | Three-phase (#28) | Two-phase (this) |
| --- | --- | --- |
| Claim | E,B,G at 120° of one wave | Light=E/B; matter=mass↔gravity via \(T^\mu_\mu\) |
| Derived from GEM/KK/duality? | **No** | Partial (Dirac, Compton, trace) |
| Photons / Cassini | Clash if gravity phase off for light | **Avoided** if metric universal |
| \(a_1\) flip | Shared-Kabat **assumed** | IR Proca still \(-1/12\); no derived flip |
| EP | — | Geometric; not phase-lock |
| Ledger gap | Untouched | Untouched |

**Improvement:** Cassini-safe by construction. **Still negative** for a derived mass–gravity phase lock and for flipping SM \(a_1\).

## Plain statement

Two-phase is the better framing for photons (universal metric + mass-only dilaton), but mass↔gravity as a locked phase pair is not derived, EP is not a phase lock, and induced \(a_1\) stays wrong-sign (\(-43/12\) UV, \(-1/12\) IR Proca) without tuning.
