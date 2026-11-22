# Soft-linear audit of the extra vertex-ll copy

Runnable: `ffz_ffff_B_soft_linear_audit.py`.

## What was computed

On FFff(B) double-log sectors, for one extra copy of the fermion-vertex Weyl-weight density:

| Channel | Value | Notes |
| --- | --- | --- |
| Leading soft-soft residue | \(-16\) | Supplies \(\Delta a_1\); exact |
| \(W^{-3\varepsilon}+U^{\varepsilon}\) finite | \(64/3\) | Exact; already in the simple ledger |
| Soft-linear with \(U^{-5}\) product rule | \(\approx -105.304\) | Quadrature-stable (n=12…40); U-rule shift exactly \(+30\) vs raw |
| Soft-linear raw (no U rule) | \(\approx -135.304\) | |

The soft-linear measure factor is the standard limit
\(\varepsilon/((\varepsilon u+1)\varepsilon v)\to 1/v\) of a single soft Taylor term.

## Why soft-linear is excluded from the overlap fix

The overlap assembly that fixed the double pole is a **two-log** correction: both soft axes logarithmic, where two \(O(\varepsilon)\) endpoint factors overlap (`FFFF_B_OVERLAP_FIX.md`). Soft-linear is a **single**-soft channel (one log endpoint × regular). Global \(\lambda_F\) doubling is already ruled out because it shifts power/log by \(-48\). By the same region logic, soft-linear is outside the assembly fix and must not be appended to \(c_0\).

If it is wrongly included,

\[
c_0 \approx 192+\tfrac{64}{3}-105.304 \approx 108
\quad\Rightarrow\quad
\text{B simple}\approx +0.013,\ \text{gap}\approx +0.055,
\]

which is much worse than the overlap-only ledger and is not a credible finish-line.

## Overlap-only status

With \(a_0=-48\), \(a_1=208\), \(c_0=192+64/3\):

- B simple \(\approx -0.03271\)
- B needed \(\approx -0.04173\) (using A simple \(4\times 0.03724\))
- Gap \(\approx +0.00903\)

Closing that gap by A alone requires A's projective finite coefficient \(c_0^A\approx 12.33\) instead of the reported \(\sim 14.93\) (\(\Delta\approx -2.6\)). That is inside the qualitative uncertainty stated in `FFFF_A_REGULATED.md` (“finite coefficient is more sensitive”), but it is not yet a controlled re-extraction.

## Next math step

Re-extract FFff(A) projective \(c_0\) / simple constant with tighter \(\varepsilon\) fits and cutoff control. Do not fold soft-linear into B’s overlap \(c_0\) unless an independent derivation shows the assembly really doubles the single-soft Taylor channel.
