# Three-phase idea (E, B, gravity) + Majorana \(a_1\) reconciliation

Runnable: `analysis/ffz/ffz_three_phase.py`  
Test: `.../test_ffz_three_phase.py` (hooked in `scripts/run_checks.py`).  
Stacked on SM-roster (#27).  
**No tuning.** DERIVED vs ASSUMED labeled.

## (0) Reconcile #22 vs #27: Majorana gaugino \(a_1\)

| Observable | Maxwell+ghosts | Majorana | Net \(N=1\) U(1) | Where used |
| --- | ---: | ---: | ---: | --- |
| Geometric entropy / EE (Solodukhin) | \(-2/3\) | \(+D/6=+1/3\) | \(\boldsymbol{-1/3}\) | **#22** |
| Sakharov induced \(1/G\) (Adler–Visser) | \(-2/3\) | \(+1/12\) (Weyl weight) | \(\boldsymbol{-7/12}\) | **#27** |

**Derived:** both Maxwell numbers agree; the fermion weights differ because EE counts every on-shell dof as \(+1/6\) in the area law, while Sakharov weights Weyl fermions as \(+1/12\) in \(a_1=C_{\mathrm{AV}}/6\) units (\(N_{\mathrm{Dirac}}=N_{\mathrm{Weyl}}/2\)). Ratio EE/Sakharov for Majorana \(=4\).

**For induced gravity (\(\kappa=0\) Sakharov reading): use \(+1/12\) and net \(\boldsymbol{-7/12}\).**  
#22’s \(+1/3\) / \(-1/3\) are correct for EE/Kabat bookkeeping but **wrong if read as \(1/G\)**. Qualitative #22 conclusions (Kabat not canceled by gaugino; multiplet wrong-sign; \(\Lambda^2\neq h^2z/K\)) still hold with \(-7/12\).

## (1) Three-phase formalization

Dan’s idea: electricity, magnetism, gravity as three phases of one wave; gravity “unrealized” for massless photons, on for massive matter.

| Framework | Supports E/B/G as 120° phases? | Notes |
| --- | --- | --- |
| GEM (Mashhoon) | **No** | \(E_g,B_g\) are linearized metric potentials, not a third Maxwell phase |
| Kaluza–Klein | **No** | \(g_{MN}\to g_{\mu\nu}+A_\mu+\varphi\); dilaton can miss photons (\(T^\mu_\mu=0\)) but metric does not |
| EM duality / Weyl | **No** | Duality mixes \(E\leftrightarrow B\) only; no gravito root of unity |
| SU(3)/\(\mathbb{Z}_3\) | **No** (analogy only) | Cube roots sum to 0 (**derived** algebra); map to (E,B,G) **assumed**, not from Einstein–Maxwell |

### Trace channel (“unrealized for photons”)

Maxwell \(T_{\mu\nu}\) is **traceless** in \(d=4\) (**derived**): \(T^\mu_\mu(\mathrm{EM})=0\).  
A dilaton / scalar-tensor coupling \(\propto T^\mu_\mu\) therefore **vanishes for pure EM** and turns on for nonrelativistic mass (\(T^\mu_\mu\sim-\rho c^2\)).

### Clash with experiment / GR

Photons **do** gravitate via the full \(T^{\mu\nu}\) (light bending \(4GM/(c^2 b)\), Shapiro delay). Equivalence principle: null geodesics feel the metric sourced by all stress-energy. Identifying the “gravity phase” with *only* the trace channel and dropping metric coupling for photons is ruled out by light deflection and Cassini.

**Honest:** no derived three-phase wave; the only clean “off for photons” knob is trace coupling, which is **not** GR gravity.

## (2) Induced-gravity sign test (\(\mathbb{Z}_3\) / boundary)

SM baseline (#27): \(a_1=-43/12\) (wrong sign).

| Option | \(a_1\) | Sign OK? | Label |
| --- | ---: | --- | --- |
| Shared Kabat among \(\mathbb{Z}_3\) triples (\(a_{1,\mathrm{vec}}=0\)) | \(+53/12\) | yes | **assumed** (not Donnelly–Wall: each U(1) has its own edge modes) |
| Triple fermion copies | still needs unjustified multiplicity | — | **assumed** = tuning |
| Free boundary \(\kappa_b\) | continuous | — | **tuning** |

**Derived without extra assumption:** does **not** flip the sign. Shared-Kabat would, but that is a free change of contact multiplicity.

## (3) Testable predictions and bounds

| Prediction (if gravity phase = photon–metric decoupling) | Bound | Status |
| --- | --- | --- |
| \(\gamma\to0\) (no / reduced light bending & Shapiro) | Cassini \(\lvert\gamma-1\rvert<2.3\times10^{-5}\) (Bertotti et al. 2003) | **ruled out** |
| Small admixture \(\varepsilon\) in photon–metric coupling | \(\lvert\varepsilon\rvert\lesssim2.3\times10^{-5}\) | room only at that level |
| Gravitational birefringence from \(\mathbb{Z}_3\) phases | — | **not derived** from the idea; astro bounds exist in literature (not re-derived) |
| Massive vs photon EP violation | LLR/EP \(\eta\lesssim10^{-13}\) (published; **assumed** here) | any differential coupling tightly constrained |

## Plain statement

For Sakharov \(1/G\), the \(N=1\) multiplet \(a_1\) is \(\boldsymbol{-7/12}\) (#27 correct; #22’s \(-1/3\) is EE). The three-phase E/B/G picture is not derived from GEM, KK, or duality; trace-only “unrealized photons” conflicts with light bending and Cassini. \(\mathbb{Z}_3\)/boundary contact does not flip SM \(a_1\) without tuning.
