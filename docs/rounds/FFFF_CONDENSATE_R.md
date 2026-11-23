# Condensate / phase-locking as a finite \(R\) contact (new-idea rung)

Runnable: `analysis/ffz/ffz_condensate_r.py`.  
Ladder: `docs/PROVEN_METHODS_LADDER.md` (**new idea**, not an established peer-reviewed method).  
**No tuning to \(1/24\).** Derived vs assumed labeled below.

## (1) Background vev — SUSY?

| Item | Status | Label |
| --- | --- | --- |
| \(W=mXY+yX^3/6\) ⇒ \(F_X=mY+(y/2)X^2\), \(F_Y=mX\) | SUSY vacuum \(X=Y=0\) unique | **derived** |
| SUSY-preserving flat direction with \(\langle\phi\rangle\neq0\) | **None** | **derived** |
| Background \(\langle X\rangle=v\) (phase-locked condensate) | Requires soft and/or F-term breaking | **assumed** |

## (2) Induced shift of the \(R\) coefficient

### Tree level (\(\xi=1/6\))

Jordan-frame contact \(\xi R|X|^2\) with \(v^2\equiv|\langle X\rangle|^2\):

\[
\delta_{\mathrm{tree}}=\frac{\xi v^2}{h^2 z/K}=\frac{K v^2}{6\,h^2 z}
\qquad\text{(units \(h^2 z/K\))}.
\]

| Property | Result | Label |
| --- | --- | --- |
| Form of \(\delta_{\mathrm{tree}}\) | \(K v^2/(6 h^2 z)\) | **derived** |
| Pole structure | Finite matching / simple residue only; **no** \(1/\varepsilon^3,1/\varepsilon^2\) | **derived** |

### One / two loops (masses \(m^2\sim |y|^2 v^2=h K v^2\))

| Piece | Result | Label |
| --- | --- | --- |
| \(\delta_{1\mathrm{loop}}\) estimate (\(c=1\)) | \(K v^2/(12 h z)\) | structure **derived**; \(c=1\) **assumed** |
| Tree vs 1-loop | \(\sim 2/h\) → tree dominates at weak \(h\) | **derived** |
| Full 2-loop re-ledger with split masses | Not re-run; UV higher poles IR-mass blind | no re-ledger **assumed**; UV blindness **derived** |

## (3) Is \(v\) fixed? Gap-absorbing value

| Condition | Fixes \(v\)? |
| --- | --- |
| \(a_Y=-1\) projection | **No** (source weight, not a vev) |
| Yukawa RG fixed point | **No** (fixes \(h_*\), not \(v^2/z\)) |
| Coleman–Weinberg for SUSY \(W\) | **No** (minimum at \(X=Y=0\)) |
| Soft scale \(m_{\mathrm{soft}}\) | **External input** (\(v\sim m_{\mathrm{soft}}/\|y\|\)) |
| Gap itself | **Tuning**, not a prediction |

**Verdict:** \(v\) is **free** (or external). No independent model condition predicts the gap-absorbing value.

### Value that absorbs GAP \(=+0.00975\)

\[
\frac{K v^2}{6 h^2 z}=\mathrm{GAP}
\quad\Rightarrow\quad
\frac{v^2}{z}=3.703\times10^{-4}\,h^2.
\]

| \(h\) | \(v^2/z\) | \(\|v\|/\sqrt{z}\) |
| ---: | ---: | ---: |
| \(1\) | \(3.703\times10^{-4}\) | \(1.924\times10^{-2}\) |
| \(0.1\) | \(3.703\times10^{-6}\) | \(1.924\times10^{-3}\) |

Independent prediction of this value? **No.**

## (4) Simple pole vs pole structure

- Tree \(\xi v^2\): **simple / finite matching residue only**.
- Soft mass loops: may shift \(m^2/\varepsilon\) pieces; **triple/double of the \(h^2\) ledger unchanged**.
- Higher-pole cancelation at \((a_X,a_Y)=(0,-1)\) **stands**.

## Plain statement

A condensate can **parametrize** the \(\sim+0.01\) gap as a Jordan-frame \(\xi v^2\) contact (same scheme class as the pure-metric contact in `FFFF_GAP_PHYSICAL.md`), but it does **not explain** it: \(v\) is not fixed by the Newton-test ingredients, and choosing the gap-absorbing \(v\) is tuning. This new-idea rung does not close the gap without tuning.
