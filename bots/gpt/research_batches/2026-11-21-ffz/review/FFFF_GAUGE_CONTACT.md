# Gauge / Kabat / Donnelly–Wall contact vs the Newton gap (new-idea rung)

Runnable: `ffz_gauge_contact.py`.  
Dan’s surface idea via gauge fields: add an \(N=1\) vector multiplet to WZ.  
Stacked on ξ-running (`FFFF_XI_RUNNING.md`, #21).  
**No tuning to \(1/24\).** Derived vs assumed labeled.

Standing GAP \(=+0.00965052\pm0.0024\) (units \(h^2 z/K\)).

## (1) Gauge-sector \(a_1\) and contact prescription

| Piece | \(a_1/\int R\) (d=4) | Label |
| --- | ---: | --- |
| Maxwell + ghosts | \(D_s/6-c_{\mathrm{Kabat}}=2/6-1=\boldsymbol{-2/3}\) | **derived** (Solodukhin/Kabat) |
| Kabat contact \(c\) | \(\boldsymbol{1}\) | **derived** (HK) |
| Donnelly–Wall reading | edge modes ≡ Kabat contact (one negative surface scalar) | coeff **derived**; interpretation **assumed** |
| Majorana gaugino | \(D/6=1/3\), contact \(=0\) | **derived** |
| Net \(N=1\) U(1) multiplet | \(\boldsymbol{-1/3}\) | **derived** |
| SU(\(N\)) | \(-(N^2-1)/3\) | **derived** |

**Prescription used:** Donnelly–Wall edge modes for the statistical meaning of Kabat’s term; coefficient from heat kernel.

## (2) SUSY cancelation? \(g^2\) vs \(h\)

| Question | Answer | Label |
| --- | --- | --- |
| Does the gaugino cancel the Kabat contact? | **No** (contact is gauge/edge-mode only) | **derived** |
| Net multiplet remainder | \(-1/3\) per U(1) | **derived** |
| Is that remainder in the \(h^2 z/K\) ledger channel? | **No** — free \(a_1\) renormalizes \(1/G\) with a \(\Lambda^2\) piece, mass-independent at leading order | **derived** |
| \(N=2\) relation \(\tilde g=g^2/K=h\) | Fixes \(g\) vs \(h\), but **does not** remove the \(\Lambda^2\leftrightarrow z\) choice | relation **assumed**; limitation **derived** |

Assumed conversion \(\Lambda^2\to z\) at \(h=1\): \(\delta c=+1/3\) in ledger units (ratio to GAP \(\sim35\)) — **not a prediction**.

## (3) vs GAP; \(\ln2\)

| Quantity | Value |
| --- | ---: |
| Predicted (derived MS / \(h^2 z/K\) channel) | \(0\) |
| **Ratio predicted / GAP** | \(\boldsymbol{0}\) |
| Assumed \(\Lambda^2\to z\) @ \(h=1\) | \(+0.333\) (ratio \(\sim34.5\); not derived) |
| \(\ln2\) from continuum U(1) edge modes? | **No** |
| \(\ln2/72\) reproduced? | **No** |

Continuum Maxwell / Donnelly–Wall edge modes give a continuous flux entropy, not \(\ln2\). A \(\mathbb{Z}_2\) gauge theory would produce \(\ln2\), but that is not the \(N=1\) vector multiplet. \(\ln2/72\) is **not** derived here.

## (4) Higher poles

Free one-loop vector multiplet: no triple/double. Pure-WZ ledger unchanged (gauge was out of model). Full interacting gauge+WZ re-ledger **not attempted**.

## Plain statement

The Kabat contact is real and SUSY does **not** cancel it against the gaugino (net \(a_1=-1/3\) per U(1)), but that effect lives in the \(\Lambda^2\) renormalization of \(1/G\), not in the WZ three-loop \(h^2 z/K\) simple pole. Derived ratio to GAP \(=0\). \(\ln2/72\) is not reproduced. This gauge-contact rung does **not** explain the gap.
