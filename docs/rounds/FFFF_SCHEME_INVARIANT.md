# Scheme-independent test of the Newton gap (new-idea rung)

Runnable: `analysis/ffz/ffz_scheme_invariant.py`.  
Stacked on gauge/Kabat (`FFFF_GAUGE_CONTACT.md`, #22).  
**No tuning to \(1/24\).** Derived vs assumed labeled.

Ledger: simple total \(=+0.051317\pm0.0024\), \(c_Y=1/24\), GAP \(=+0.00965052\pm0.0024\) (units \(h^2 z/K\)).

## (1) Allowed finite CTs and observables

### Allowed finite counterterms (WZ+Yukawa, \(\xi=1/6\), SUSY)

| Sector | Operator | Allowed? | Absorbs GAP? | Label |
| --- | --- | --- | --- | --- |
| Matter | \(Z_z\) / wave-function / Yukawa CT | yes | no | **derived** (already in ledger) |
| Matter | \(\xi R\|\Phi\|^2\) with \(\xi\neq1/6\) | **no** | no | **derived** (SUSY fixes \(\xi\)) |
| Matter | non-holomorphic curvature–Kähler \(R\|\Phi\|^2\) | **no** | no | **derived** (SUSY/holomorphy) |
| Matter | FI / D-term | **no** | no | **derived** (no gauge in pure WZ) |
| Vacuum | \(\int\sqrt{g}\,R\) (finite \(\kappa\)) | **yes** | **yes** | op **derived** (curved-space renorm); finite \(\kappa\) **scheme** |
| Vacuum | \(R^2\), \(W^2\), Euler | yes | no | **derived** (not linear-\(R\) Newton) |
| Vacuum | \(\Box R\) | yes | no | **derived** (total derivative) |

Packet: finite shift \(\kappa\) of the pure-metric coefficient moves its beta by \(4\kappa h^2\sum z/K\).

### Observables

| ID | Observable | Value | SI? | Label |
| --- | --- | ---: | --- | --- |
| O1 | Higher poles at \((0,-1)\) | triple \(=0\), double \(=0\) | **yes** | **derived** |
| O2 | Soft \(c_Y\) | \(1/24\) | **yes** | **derived** (Martin \(\beta_\Lambda\)) |
| O3 | \(D\equiv c_0-c_Y\) at \(\kappa=0\), \(L=0\) | \(+0.00965052\pm0.0024\) | **no** | value **derived**; SI status **derived** (κ shifts \(D\)) |
| O4 | \(D+\kappa\) with \(\kappa=-\mathrm{GAP}\) | \(0\) | no | **derived** (removal, not a prediction) |
| O5 | \(D/c_Y\) at \(\kappa=0\) | \(0.2316\) | **no** | same scheme dependence as \(D\) |
| O6 | \(R\)-simple / 3-loop \(a\)-anomaly ratio | — | — | **infeasible** (3-loop anomaly not in packet) |

## (2) Does GAP survive as SI nonzero?

| Question | Answer |
| --- | --- |
| Survives as scheme-independent nonzero? | **No** |
| Absorbed by? | finite vacuum \(\int\sqrt{g}\,R\) with \(\kappa=-\mathrm{GAP}=-0.00965052\) |
| Matter CT of existing SUSY WZ absorbs it? | **No** |
| Vacuum EH absorbs it? | **Yes** |

## (3) Missing-operator diagnosis

| Claim | Result | Label |
| --- | --- | --- |
| Missing matter operator forced? | **No** | **derived** |
| Minimal absorber | finite vacuum \(\int\sqrt{g}\,R\) (\(\kappa=-\mathrm{GAP}\)) | **derived** |
| Forced vs optional | operator **forced** by curved-space renormalizability; finite \(\kappa\) **optional** (scheme) | **derived** |

Candidates ruled out as “genuinely missing”:
- \(R\|\Phi\|^2\) beyond \(\xi=1/6\) — optional SUSY-breaking, not required
- curvature–Kähler / superspace curvature — not forced (SUSY kinetic already conformal)
- FI/D-term — out of pure WZ
- \(R^2/W^2\) — do not absorb the linear-\(R\) Newton residue

## Higher poles

Triple/double remain SI zeros; finite \(\kappa\) does not reopen them.

## Plain statement

The \(\sim+0.01\) Newton gap is a **scheme-dependent matching residual**: it does **not** survive as a scheme-independent nonzero quantity once the allowed finite vacuum Einstein–Hilbert counterterm is admitted. No missing matter Lagrangian operator is forced. Choosing \(\kappa=-\mathrm{GAP}\) is a scheme statement for \(G_{\mathrm{ind}}\), not evidence of new physics beyond the vacuum CT sector already required by curved-space renormalizability.
