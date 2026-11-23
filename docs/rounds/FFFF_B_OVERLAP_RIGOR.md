# Overlap rigor: soft-soft only, from the regularization

Runnable: `ffz_ffff_B_overlap_rigor.py`.

This note derives, from sector Taylor inclusion-exclusion (not from the \(1/24\) target), why the extra vertex-weight copy that supplies \(\Delta a_1=-16\) and \(\Delta c_0=64/3\) does **not** include soft-linear \(\approx-105.304\). It also checks the coupling-weight residue against \(\lambda_F\).

## Sector inclusion-exclusion

On a double-log sector the two soft parameters \((x,y)\) are subtracted by Taylor operators at each axis. For a reduced density \(G(x,y)\),

\[
R=G(x,y)-G(x,0)-G(0,y)+G(0,0).
\]

Through bilinear order \(R\) is \(O(xy)\): soft-linear terms are removed by the single-axis subtractions. The **overlap addback** that restores the double-subtracted corner is

\[
\varepsilon\int_0^1\!\mathrm{d}x\,\mathrm{d}y\,x^{\varepsilon u-1}y^{\varepsilon v-1}G(0,0)
=\frac{G(0,0)}{\varepsilon\,u\,v},
\]

which depends only on the soft-soft leading value \(G(0,0)\). Soft-linear never enters that addback.

A single log axis with an \(O(\varepsilon)\) numerator gives \(\varepsilon/(\varepsilon u)=1/u\) from the leading piece and a finite soft-linear piece with **no** \(1/\varepsilon\) pole. So soft-linear is not an overlap pole channel.

## What the missing overlap copy contributes

If the field-weight-only assembly under-counts the overlap addback by one \(O(\varepsilon)\) soft-soft copy with leading density \(\rho\) (ll residue \(\int\rho=-16\)), the correction is exactly that addback:

| Piece | Value | In overlap addback? |
| --- | --- | --- |
| \(\Delta a_1\) (soft-soft leading) | \(-16\) | yes |
| \(\Delta c_0\) from \(W^{-3\varepsilon}+U^{\varepsilon}\) on that density | \(64/3\) | yes (same leading \(\rho\)) |
| Soft-linear \(\approx-105.304\) | — | **no** |
| Power/log shift | — | **no** (overlap is ll) |

Global doubling of the vertex operator would also move power/log (vertex/coupling pl residue exactly \(-48\)) and soft-linear; that overshoots the double-pole budget (the known \(-48\) pl shift from global \(\lambda_F\) doubling). The forest structure selects soft-soft leading only.

## Coupling weight vs \(\lambda_F\)

Packet: Yukawa Weyl weight \(\varepsilon\), field weight \(\lambda_F=(4-d)/2\). At \(d=4-2\varepsilon\), \(\varepsilon=\lambda_F\).

A local Weyl variation \(\delta y=-\varepsilon\sigma(x_v)y\) with curvature probe \(\sigma\sim x_v^2\) inserts the diagonal position moment \(T_{vv}\) at each Yukawa vertex. The field-weight vertex piece is \(-\lambda_F\sum_i T_{ii}\). Therefore, as operators at this order,

\[
A_{\text{coupling}}=A_{\text{vertex}}.
\]

The runnable check asserts that identity and finds both ll residues equal to \(-16\) (pl residues equal to each other as well). **Equality is an operator identity**, not a fit to the Newton target.

## Status labels

- **Derived from regularization:** overlap addback = soft-soft leading only; soft-linear excluded; \(\Delta a_1=-16\); \(\Delta c_0=64/3\); \(A_{\text{coupling}}\equiv A_{\text{vertex}}\) at this order.
- **Still assumed / scheme-level:** that the missing assembly piece is exactly one such overlap copy of this operator (rather than a different forest operator with the same soft-soft leading but different non-overlap footprint). Global alternatives that change pl are excluded by the double-pole budget, which supports but does not uniquely prove the identification.
