# FFff(B) non-overlap finite \(\Delta c_0\) — first principles

Runnable: `ffz_ffff_B_nonoverlap_finite.py`. **Not tuned to \(1/24\).**

Prior status: field \(c_0\) + overlap \(\Delta c_0(W/U)=64/3\) on the ledger; anything beyond that was **assumed** zero. Empirical \(2\times(64/3)\) would close the \(\sim+0.01\) gap but was not derived.

## Question

Does the coupling-weight vertex copy contribute a **second** \(\Delta c_0=64/3\) (or any other value) through **non-overlap** sectors (power/log, single-soft finite), while remaining consistent with the rigor result that only soft-soft ll enters the overlap addback for the slope?

## Derived answer

\[
\boxed{\Delta c_0^{\mathrm{non\text{-}overlap}} = 0 \quad \text{(exact)}}
\]

**Label: derived** from sector Taylor / inclusion-exclusion (same regularization that fixes \(\Delta a_1=-16\)).

### Forest argument

On a double-log sector the overlap addback that restores the double-subtracted corner is

\[
\varepsilon\int_0^1\mathrm{d}x\,\mathrm{d}y\,x^{\varepsilon u-1}y^{\varepsilon v-1}G(0,0)
=\frac{G(0,0)}{\varepsilon\,u\,v}.
\]

It depends only on the soft-soft leading value. Soft-linear is removed by single-axis subtractions and never enters that addback (`FFFF_B_OVERLAP_RIGOR.md`).

The same addback, with \(W^{-3\varepsilon}U^{\varepsilon}=1+\varepsilon(-3\ln W+\ln U)+\cdots\), supplies the **one** finite piece

\[
\Delta c_0^{\mathrm{overlap}}(W/U)=\frac{64}{3}.
\]

By construction that integral has:

| Channel | In forest addback? |
| --- | --- |
| soft-soft leading \(\to\Delta a_1=-16\) | yes |
| soft-soft \(W/U\to 64/3\) | yes |
| soft-linear (ll) | **no** |
| single-soft-only | **no** |
| power/log | **no** (\(\Delta\mathrm{pl}=0\)) |

So the non-overlap finite of the forest copy is **exactly zero** — not a second \(64/3\).

### Full-operator diagnostic (not the forest term)

If one instead added a **global** extra vertex/coupling operator (already excluded by pl\(=0\)):

| Piece | Value | \(=\;64/3\)? |
| --- | ---: | --- |
| ll soft-linear | \(-105.303856\) | no (\(\Delta\approx-126.6\)) |
| pl \(W/U\) finite | \(-48.538020\) | no (\(\Delta\approx-69.9\)) |
| pl residue | \(-48\) | would spoil double pole |

None of these is \(64/3\). Empirical \(2\times(64/3)\) has **no** home in non-overlap sectors under the regularization that fixes the slope.

## Updated simple pole vs \(1/24\)

B uses \(c_0=\mathrm{field}+64/3+0\) with field \(c_0\approx191.72\pm1.36\):

| Quantity | Value |
| --- | ---: |
| B simple | \(-0.03259\pm0.00059\) |
| Projected total | \(+0.05145\pm0.0021\) |
| Target \(1/24\) | \(+0.041667\) |
| **GAP** | \(\boldsymbol{+0.0098\pm0.0021}\) |

**Gap does not close.**

## Status labels

- **Derived:** forest non-overlap \(\Delta c_0=0\); overlap \(W/U=64/3\); soft-linear/pl excluded from the correction.
- **Reproduced:** numerical soft-linear \(\approx-105.304\), pl \(W/U\approx-48.54\) (diagnostic only).
- **Not derived:** empirical second \(64/3\).
