# FFff(B) two-log overlap assembly fix

## What was wrong

Independent extraction of the field-weight-only FFff(B) kernel gave projective slope \(224\) and bare double \(-1/96\) at \(L=0\). Against corrected A (\(5/288\)) and the counterterm budget, the projected double residual vs Y was \(-1/144\), with uncanceled \(+L/48\) in the simple pole. The triple pole still cancelled.

## The \(-16\) term

Sector-splitting the numerator \(\varepsilon\)-derivative (`ffz_ffff_B_slope16_hunt.py`) shows that the fermion-vertex Weyl-weight piece on **double-log** sectors (two overlapping logarithmic endpoints) has residue **exactly \(-16\)**. That is the full \(\Delta a_1\) needed for slope \(224\to 208\) and double \(-1/96\to -1/288\).

Global doubling of \(\lambda_F\) also shifts power/log by \(-48\) and overshoots. The fix is therefore an **overlap assembly** correction, not a global rescaling of \(\lambda_F\).

## Fix applied

Runnable: `ffz_ffff_B_overlap_fix.py`.

\[
a_1 \;=\; 224 + (-16) \;=\; 208
\qquad\Longrightarrow\qquad
\text{B double}\;=\;-\frac1{288}\quad(L=0).
\]

With A double \(5/288\), bare subtotal \(5/144\), CT \(-5/48\):

\[
\frac5{144}+\frac5{288}-\frac1{288}-\frac5{48}-\Bigl(-\frac1{18}\Bigr)\;=\;0.
\]

Scale logs in the double still cancel (\(-L/48-L/16+L/12=0\)). The residual-driven \(+L/48\) in the simple pole is gone.

## Physical reading

Packet notes assign the Yukawa *coupling* Weyl weight \(\varepsilon\), equal at this order to the field weight \(\lambda_F=(4-d)/2\). The field-weight-only kernel implements \(-\lambda_F\sum M_{ii}\). On the two-log overlap, a second soft-soft \(O(\varepsilon)\) factor with the same residue \(-16\) was missing from the endpoint assembly. This matches the parallel audit’s suspect (assembly of \(O(\varepsilon)\) terms where two endpoints overlap).

**Derived:** vertex-ll residue \(-16\); arithmetic that \(\Delta a_1=-16\) clears the ledger.  
**Assumed in the fix:** that missing factor is exactly one extra copy of that soft-soft piece (coupling-weight / overlap forest), not a different operator with a different power/log footprint. Global kernel operators tried (`d C_{ii} N`, etc.) did not give \(\mathrm{ll}=-16\) with \(\mathrm{pl}=0\); hence the correction stays in the assembly step.

## What is still open

- Independent derivation of the coupling-weight operator that reproduces \(\mathrm{ll}=-16\) without a hand assembly step.
- Projected simple pole vs \(1/24\) — see `FFFF_B_SIMPLE_LEDGER.md` (vertex-ll \(W/U\) finite supplies exact \(\Delta c_0=64/3\); gap not closed).
- Finite-\(\varepsilon\) absolute normalization vs ChatGPT’s \(\sim 1.00333\).
