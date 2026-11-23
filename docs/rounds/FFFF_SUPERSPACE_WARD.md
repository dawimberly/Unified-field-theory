# Superspace / supergraph rung (SUSY Ward audit)

Runnable: `analysis/ffz/ffz_superspace_ward.py`.  
Ladder entry: `docs/PROVEN_METHODS_LADDER.md` (rung 5).  
**No tuning to \(1/24\).** Full curved-superspace three-loop supergraphs are **not** attempted.

## Scope

| Done (feasible) | Not done |
| --- | --- |
| Non-renormalization ⇒ what is forced for induced \(R\) | Full \(N=1\) curved superspace / SUGRA-background 3-loop supergraphs |
| Wick inventory + boson↔fermion partner map | D-algebra reduction of all 10 classes to supergraphs |
| Two-loop cubic↔Yukawa geometric Ward | |
| Projected higher-pole Ward at \((a_X,a_Y)=(0,-1)\) | |
| Missing / mis-weighted partner audit vs ledger | |

## (A) Does SUSY force the three-loop simple pole \(=1/24\)?

| Claim | Forced by SUSY? |
| --- | --- |
| Soft Y coefficient \(c_Y=1/24\) | **YES** — holomorphic non-renormalization + Martin \(\beta_\Lambda=3/2\) + \(s_Y=-a_Y R/6\) |
| Y-linear higher poles \(1/36,-1/18,+1/24\) | **YES** — soft-source pole equation |
| Jack–Osborn superspace \(\delta=0\) | **YES** for additive *source* beta; does **not** set baseline\(=c_Y\) |
| Baseline (improvement-independent) induced-\(R\) simple \(=1/24\) | **NO** — Kahler / metric matching residue; flat vacuum RG alone does not determine it |

Supergraph D-algebra organizes the flat vacuum sum; a curved metric / Weyl probe is a soft or background spurion that breaks rigid SUSY, so D-algebra does **not** force baseline \(=c_Y\).

## (B) Component Ward checks

### Inventory (10/10)

| Multiset | Multiplicities | Role |
| --- | --- | --- |
| QQ | 4, 16 | pure scalar |
| QCc | 4, 16 | pure scalar |
| QFf | 8 | fermion (QCc tadpole twin) |
| CCcc | 8, 32 | pure scalar |
| CcFf | 8 | fermion (CCcc layout + fermion pair) |
| FFff | 8, 32 | fermion (QQ twin; A=2-cycle, B=4-cycle) |

H tetrahedron **absent** (charge bipartite). Martin E,G,J,K,L present.

### Partner map

| Sector | Boson | Fermion |
| --- | --- | --- |
| tadpole×sunset | QCc | QFf |
| cubic vacuum | CCcc | CcFf |
| quartic / 4-fermion | QQ | FFff(A+B) |
| two-loop seed | cubic sunset | Yukawa sunset |

### Wards that pass

- Two-loop \(zJ\) geometric cancelation (leftover \(z|y|^2 I'/36\)).
- Projected triple: \(1/144+1/48=1/36\).
- Projected double residual \(=0\) (overlap-fixed B).

### Simple-pole aggregates \(L=0\) (bookkeeping only)

| Aggregate | Value |
| --- | ---: |
| QCc+QFf | \(-0.11894\) |
| CCcc+CcFf | \(-0.09424\) |
| QQ+FFff | \(+0.17947\) |
| Full sum (+CT) | \(+0.05141\) |
| Target \(1/24\) | \(+0.041667\) |
| **GAP** | \(\boldsymbol{+0.00975}\) |

Partner aggregates are **not** required to vanish.

## (C) Missing / mis-weighted partners

| Check | Result |
| --- | --- |
| Orphan fermion class | none |
| Multiplicity patterns | match packet / bipartite SUSY structure |
| FFff(A) ×4 spin normalization | present in ledger |
| Gap vs missing-class scale | \(0.0097\) \(\lesssim\) half a typical Yukawa simple (\(\sim0.02\)); **not** a fully missing mult-8 topology |

**FLAGS: none.**

## Plain statement

The superspace/supergraph feasible subset confirms: SUSY **requires** \(c_Y=1/24\) and the higher-pole/flat Wards, but **does not force** the baseline simple pole to \(1/24\). The component ledger respects the partner inventory; no missing/mis-weighted superpartner is identified. The gap **stands at \(+0.00975\)**. This rung does not close it.

Next ladder rung: background-field method with \(R^*\) operation.
