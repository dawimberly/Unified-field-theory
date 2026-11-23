# Surface / boundary / contact terms vs the Newton gap (new-idea rung)

Runnable: `analysis/ffz/ffz_surface_term.py`.  
Dan’s original idea: gravity from mode phase/frequency/amplitude at field-interaction surfaces.  
Stacked on self-consistent loop (`FFFF_SELF_CONSISTENT_LOOP.md`, #19).  
**No tuning to \(1/24\).** Derived vs assumed labeled.

Standing target: GAP \(=+0.00965052\pm0.0024\) (R* ledger vs \(1/24\), units \(h^2 z/K\)).

## (1) Sakharov spectral sum — bulk vs boundary

| Sector | Form | Role at 3 loops, \(\xi=1/6\) | Label |
| --- | --- | --- | --- |
| Bulk modes | \(\sum_n\alpha_n/(2\omega_n)\) ↔ zeta / HK \(a_1\) + interactions | \(a_1^{\mathrm{bulk}}=(1/6-\xi)R\to0\); interacting 3-loop = existing Wick ledger | **derived** |
| Geometric boundary (GHY / extrinsic \(K\)) | Vassilevich \(a_1^\partial\propto\int_{\partial M}\sqrt{h}\,K\) (BC-dependent) | Newton probe has **no** \(\partial M\); Stokes on compact support ⇒ **0** to \(\int\sqrt{g}R\) | formula **assumed** (Vassilevich); vanishing **derived** |
| \(\Box R\) / total derivatives | \(\int\sqrt{g}\,\Box R=0\) (closed / falloff) | Local \(c_\Box\) free (scheme); **0** to integrated EH residue | **derived** |
| Vertex \(\delta\)-contacts | Coincident Wick points / Kähler–metric matching | Forests already in CT/R*; leftover finite \(R\) contact = **scheme** (same class as `FFFF_GAP_PHYSICAL.md`) | forests **derived**; finite coeff **assumed**/free |

## (2) Finite simple-pole pieces the bulk ledger might have dropped

| Candidate | Value (units \(h^2 z/K\)) | Label |
| --- | ---: | --- |
| Bulk \(a_1\) at \(\xi=1/6\) | \(0\) | **derived** |
| Integrated GHY / extrinsic \(K\) | \(0\) | **derived** |
| Integrated \(\Box R\) | \(0\) | **derived** |
| Extra beyond soft dictionary \(s_Y=-a_Y R/6\) | \(0\) | **derived** |
| Vertex \(\delta\) contact (derived coeff) | \(0\) (none fixed) | **derived** (absence); absorbing coeff **scheme** |
| **Total derived surface/contact** | \(\boldsymbol{0}\) | |

## (3) vs GAP \(+0.00965\pm0.0024\)

| Quantity | Value |
| --- | ---: |
| Predicted (derived surface sum) | \(0\) |
| **Ratio predicted / GAP** | \(\boldsymbol{0}\) |
| Free parameter? | Yes — any absorbing finite pure-metric / vertex \(R\) contact is a free scheme coefficient. Derived GHY/\(\Box R\)/\(a_1\) pieces are parameter-free and vanish. |

Numerical proximities (**not derived** as surface terms; within \(\pm1\sigma\) of GAP):

| Expression | Value | ratio to GAP |
| --- | ---: | ---: |
| \(\ln2/72\) | \(0.009627\) | \(0.998\) |
| \(1/104\) | \(0.009615\) | \(0.996\) |
| \(7/720\) | \(0.009722\) | \(1.007\) |

These are recorded only to avoid mistaking coincidence for a surface derivation.

## (4) Optional: Jacobson / entanglement-entropy cross-check

| Check | Result | Label |
| --- | --- | --- |
| Free scalar EE area quad \(\propto(1/6-\xi)\) at \(\xi=1/6\) | \(0\) (matches \(a_1\to0\)) | **derived** |
| Weyl fermion EE | one-loop HK baseline already in packet; not a 3-loop gap source | **derived** (scope) |
| EE \(\leftrightarrow1/G\) finite contact ambiguity | exists for nonminimal scalars; **does not predict** \(+0.00965\) | ambiguity **assumed** (literature); non-prediction **derived** |
| Interacting 3-loop WZ EE | not attempted (same hardness as bulk ledger) | scope |

## Higher poles

Triple/double cancelation **unaffected** (GHY/\(\Box R\) integrate to 0; finite local \(R\) contact is simple/finite matching only).

## Plain statement

A Sakharov / surface / contact analysis of the Newton residue finds **no derived** extra simple pole: geometric boundary and \(\Box R\) drop out of the integrated EH match, \(\xi=1/6\) kills bulk \(a_1\), and vertex contacts that could absorb the gap are the same free scheme contacts already identified in the physical-gap audit. Entanglement-entropy free-field structure is consistent but does not fix GAP. This new-idea rung does **not** explain \(+0.00965\) without tuning a scheme contact to the gap.
