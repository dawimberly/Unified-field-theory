# Is the ~+0.01 simple-pole gap physical?

Runnable: `ffz_gap_physical_audit.py`. **No tuning to \(1/24\).**  
Prior ledger: `FFFF_SIMPLE_GAP_AUDIT.md` (gap \(\approx+0.00965\pm0.0024\)).

Dan order: **(3) physical meaning → (1) missing piece / tighten → (2) \(a_Y\) scan**.

---

## (3) Physical meaning

### (a) Is \(1/24\) a theorem?

| Claim | Status | Label |
| --- | --- | --- |
| Soft Y coefficient \(c_Y=1/24\) | **Yes — RG-fixed.** Martin \(\beta_\Lambda=(3/2)h^2 z s_Y/K\), engineering weight pole \(=-\beta/(2\cdot3)\), \(s_Y=-a_Y R/6\). Lower loops reproduce \(1/6\) and \(1/24\). | **derived** |
| Y-linear higher poles \(1/36\), \(-1/18\) | **Yes — fixed** by the soft-source pole equation independently of the three-loop \(\beta\) coefficient. | **derived** |
| Baseline (improvement-independent) simple \(=1/24\) | **No — not a t Hooft / RG theorem.** Packet: flat vacuum RG alone does not determine the baseline; the equality is the projected Newton **matching test**. | **matching conjecture / test** |

So a finite residual simple pole in the **baseline** is **not forbidden** by RG consistency. If the baseline graphs+CTs are complete and correct, the gap can be physical (scheme-dependent). If something is still wrong in the baseline, the gap is an error — but RG does not force that conclusion by itself.

### (b) Interpretation as a \(G_{\mathrm{ind}}\) shift

Ledger total \(\approx+0.05132\) vs target \(1/24\approx+0.04167\):

\[
\frac{\mathrm{gap}}{1/24}\approx +23.2\%.
\]

In \(\Gamma_{\mathrm{ind}}\supset -R/(16\pi G_{\mathrm{ind}})\), this is a **\(\sim+23\%\) shift of the three-loop matching piece** of \(1/G_{\mathrm{ind}}\) relative to the Y-cancelation residue (not a 23% shift of the full \(G_{\mathrm{ind}}\), which also includes lower loops and thresholds).

**Scheme dependence:** yes. A finite pure-metric \(R\) contact redefinition moves the associated beta together with the matching condition (projection checkpoint: shift \(\kappa h\sum z/K\) changes the beta by \(4\kappa h^2\sum z/K\) from \(\beta_h=3h^2\), \(\beta_z=hz\)). That can absorb the entire simple residual. Higher poles are **not** absorbable this way — and they already cancel after the FFff(B) overlap fix.

**Verdict (3):** gap is **allowed** as a scheme-dependent three-loop matching residual; it is **not** an automatic proof of a graph error.

---

## (1) Missing piece / topology / weak entries

### Topology coverage

Re-ran the Wick inventory (packet multiplicities). All ten oriented classes map onto the ledger:

| Multiset | Classes / mult | Ledger entry |
| --- | --- | --- |
| QQ | 2 / 4,16 | QQ exact \(+1/16\) |
| QCc | 2 / 4,16 | QCc doubled + tadpole pole \(-1/36\) |
| QFf | 1 / 8 | QFf derived |
| CCcc | 2 / 8,32 | CCcc exact \(-5/72\) |
| CcFf | 1 / 8 | CcFf QMC |
| FFff | 2 / 8,32 | FFff(A)×4 + FFff(B) overlap-only |
| CT | — | \(Z_z\) + two-loop matter CT (**derived**) |

**Omitted classes with size estimates:**

| Class | Status | Size vs gap |
| --- | --- | --- |
| H tetrahedron | absent (charge bipartite) | 0 |
| gauge / gaugino | out of model | 0 |
| extra soft holomorphic masses | packet sets to 0; Y already projected | 0 |
| FFff(B) non-overlap finite beyond field \(c_0+\Delta c_0(W/U)\) | **assumed** zero (overlap-only) | soft-linear excluded (\(\sim+0.046\) if wrongly kept); empirical \(2\times\Delta c_0\) shifts simple by \(\sim-0.009\) (**not derived**) |
| pure-metric finite contact | scheme, not a Wick graph | can absorb whole gap |

**No omitted Wick topology of size \(\sim0.01\).** The only \(\sim0.01\)-sized handle still labeled **assumed** is B’s non-overlap finite content (empirical \(2\times\Delta c_0\)), which is not a new graph and is not derived.

### Tightened weak entries

| Piece | New value | Error | Label |
| --- | --- | --- | --- |
| QCc doubled simple | \(-0.016497423297\) | \(\sim10^{-11}\) | **derived** (faces+tail) |
| CcFf simple | \(-0.024786\) | \(\pm1.0\times10^{-4}\) | **reproduced** QMC (power 20, 8 seeds; \(2\times\mathrm{SEM}\)) |
| FFff(B) field \(c_0\) | \(191.72\pm1.36\) | \(\Rightarrow\) simple \(\pm0.00059\) | **reproduced** numerical |
| FFff(B) \(\Delta c_0=64/3\) | exact | — | **derived** (overlap \(W/U\)) |
| FFff(B) beyond that | 0 in ledger | — | **assumed** (overlap-only) |

None of these moves the gap by \(\sim0.01\).

---

## (2) \(a_Y\) / \(a_X\) scan

Needed after (3)+(1): gap still unexplained by a missing Wick graph, and not RG-forbidden.

Higher poles: \(\mathrm{residual}_{\mathrm{higher}}\propto(a_Y+1)\) once baseline higher poles cancel Y at the projected point (X additive higher poles vanish in the packet through this order). So **triple+double cancel only at** \((a_X,a_Y)=(0,-1)\).

At that unique point the simple residual stays \(\approx+0.00965\neq0\).

| Locus | Result |
| --- | --- |
| Higher poles = 0 | only \((0,-1)\) |
| Higher + simple = 0 | **empty** |
| Simple = 0 ignoring higher | \(a_Y\approx-1.232\) (reopens higher poles by factor \(\approx-0.232\)) |

**No simultaneous solution.**

---

## Bottom line for Dan

1. **\(1/24\) as \(c_Y\):** derived (RG). **Baseline \(=1/24\):** matching test, not a pole-equation theorem — residual simple pole **allowed** and **scheme-absorbable** (\(\sim+23\%\) of the three-loop \(1/G_{\mathrm{ind}}\) piece).
2. **Missing Wick graph of size \(\sim0.01\):** none found; ledger covers the inventory. Remaining assumed handle: B non-overlap finite (empirical \(2\times\Delta c_0\)), not derived.
3. **\((a_X,a_Y)\) scan:** no point satisfies higher poles and simple together.

Next ideas (not done here): either accept a finite pure-metric contact / scheme statement for \(G_{\mathrm{ind}}\), or derive (not fit) any extra B finite content that preserves pl\(=0\).
