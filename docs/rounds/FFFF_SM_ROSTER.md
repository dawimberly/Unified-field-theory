# Full Standard Model roster in induced gravity (Sakharov, \(\kappa=0\))

Runnable: `analysis/ffz/ffz_sm_roster.py`  
Test: `.../test_ffz_sm_roster.py` (hooked in `scripts/run_checks.py`).  
Stacked on particle-plasticity (#26).  
**No tuning.** DERIVED vs ASSUMED labeled.

Standing ledger: \(C\approx0.051317\), \(c_Y=1/24\), GAP \(\approx+0.009651\pm0.0024\) (units \(h^2 z/K\)).

## (1) One-loop \(a_1\) per species

| Species | \(a_1/\int R\) (d=4) | Label |
| --- | ---: | --- |
| Real scalar | \(1/6-\xi\) | **derived** (Bunch–Parker / HK) |
| Weyl fermion | \(+1/12\) | **derived** (Adler–Visser weight: \(N_{\mathrm{Dirac}}=N_{\mathrm{Weyl}}/2\) in units of \(1-6\xi\)) |
| Dirac | \(+1/6\) | **derived** |
| Massless vector + ghosts + Kabat | \(-2/3\) | **derived** (Kabat 1995; Larsen–Wilczek; Solodukhin; Donnelly–Wall edge modes) |
| Massive vector (Proca) | \(+1/2\) | **derived** (\(3\times1/6\), no Kabat); UV SM uses massless count |

Higgs \(\xi_H\): **assumed** (free vs conformal \(1/6\)); dependence shown below.

### SM degrees of freedom (UV, \(\Lambda\gg v\))

| Sector | Count |
| --- | ---: |
| Higgs real scalars | \(4\) |
| Weyl fermions | \(45\) (\(+3\,\nu_R\) optional \(\to48\)) |
| Massless gauge vectors | \(12\) (\(8+3+1\)) |

### Known wrong-sign / conformal issue

SM alone does **not** induce a healthy Newton constant: the \(-4 N_1\) vector weight dominates. Literature: Visser, *Sakharov’s induced gravity: a modern perspective* (gr-qc/0204062); Frolov–Fursaev–Zelnikov; Adler; Zee. Reproduced numerically below.

## (2) Net \(\Lambda^2\) coefficient (sign)

Adler–Visser bracket \(C_{\mathrm{AV}}=6\,a_1^{\mathrm{net}}=N_0(1-6\xi_H)+N_{\mathrm{Dirac}}-4 N_1\).

| Roster | \(a_1^{\mathrm{net}}\) | \(C_{\mathrm{AV}}\) | Sign OK for \(G_{\mathrm{ind}}\)? |
| --- | ---: | ---: | --- |
| SM, \(\xi_H=0\) | \(-43/12\approx-3.583\) | \(-43/2\) | **no** |
| SM, \(\xi_H=1/6\) | \(-17/4=-4.25\) | \(-51/2\) | **no** |
| SM+\(\nu_R\), \(\xi_H=0\) | \(-10/3\approx-3.333\) | \(-20\) | **no** |
| SM+\(\nu_R\), \(\xi_H=1/6\) | \(-4\) | \(-24\) | **no** |
| MSSM, \(\xi=1/6\) | \(-35/12\approx-2.917\) | \(-35/2\) | **no** |

MSSM detail (**derived**): each \(N=1\) vector multiplet contributes \(-2/3+1/12=-7/12\) (Kabat not canceled by gaugino — same qualitative lesson as #22, Sakharov weights); chiral Weyl at \(\xi=1/6\) give \(+1/12\) each (scalars drop out). Net still wrong-sign.

\(\xi_H\) dependence (SM, no \(\nu_R\)): \(a_1=4(1/6-\xi_H)+45/12-8\). Sign flip would need \(\xi_H<-43/48\approx-0.896\) — **unphysical**, not a fix.

**vs WZ:** the WZ Newton test is an *interacting* three-loop \(h^2 z/K\) matching, not free-field \(a_1\). Free MSSM/WZ \(a_1\) being wrong-sign is a different statement from the ledger gap.

## (3) Map onto FFZ normalization

| Question | Answer | Label |
| --- | --- | --- |
| Does the free SM roster change \(0.051317\) vs \(1/24\)? | **No** | **derived** |
| Why? | Free \(a_1\) is a \(\Lambda^2\) renormalization of \(1/G\); the ledger simple pole is \(\propto h^2 z/(K\varepsilon)\) from the Yukawa sector | **derived** (same channel split as gauge-contact #22) |
| When do SM fields enter the ledger channel? | Only through *interacting* diagrams with their Yukawas / gauge couplings | **derived** (scope) |
| Top-Yukawa-dominated proxy | Single-coupling WZ keeps \(c_Y=1/24\) and GAP \(\approx+0.00965\); full multi-Yukawa+gauge re-ledger **not computed** | proxy **assumed** (scope) |

## (4) Gap closure?

| Claim | Result |
| --- | --- |
| Derived \(\delta C\) from free roster | \(0\) — does **not** close GAP |
| Assumed \(\Lambda^2\to z\) @ \(h=1\): \(\delta c=-a_1\) | ratios to GAP \(\sim O(10^2\)–\(10^3)\); **none** equals GAP |
| Coincidences | None identified; wrong-sign SM is literature, not a numerical accident that fixes \(+0.00965\) |

## Plain statement

The full SM (and SM+\(\nu_R\), MSSM) free-field heat-kernel roster reproduces the known wrong-sign induced-\(G\) result and lives in the \(\Lambda^2\) channel. It does **not** move the FFZ three-loop Newton coefficient \(0.051317\) vs \(1/24\), and it does **not** close the \(+0.009651\) gap without tuning.
