# FFff(A) simple / projective \(c_0\) re-extraction

Runnable: `ffz_ffff_A_simple_reextract.py` (no tuning to \(1/24\)).

## Method

Regulated integrand as in `FFFF_A_REGULATED.md` (analytic long tails, short-face Taylor). Fit

\[
y(\varepsilon)=\varepsilon\Bigl(I-\frac2{\varepsilon^2}\Bigr)
=-\frac{10}3+c_0\,\varepsilon+b\,\varepsilon^2+\cdots
\]

with intercept constrained to \(-10/3\) (supported by free intercepts \(\approx-3.333\)). Systematics: cuts \(\{0.5,1,2\}\), seeds \(\{12,41,77\}\), \(\varepsilon\in\{0.01,0.0075,0.005,0.00375,0.0025,0.001875\}\), Sobol power 17.

## Result (independent of Newton target)

| Quantity | Value |
| --- | --- |
| Projective \(c_0\) | \(15.11\pm 0.58\) |
| Cut×seed range | \([14.54,\ 15.70]\) |
| Three-smallest-\(\varepsilon\) check | \(15.20\pm 0.85\) |
| Packet prior | \(\sim 14.93\) |
| A simple \(L=0\) pre-\(\times4\) | \(0.03740\pm 0.00050\) |
| A simple \(L=0\) post-\(\times4\) | \(0.1496\pm 0.0020\) |

Error is the max of sample std and half-range over cut×seed.

## Relation to the \(1/24\) gap (informational only)

Overlap-only B (\(c_0^B=192+64/3\)) plus this A leaves projected simple gap \(\approx+0.0097\) vs \(1/24\). Nulling that gap would require projective \(c_0^A\approx 12.33\), which lies **outside** the controlled interval \(15.11\pm 0.58\). The re-extraction does **not** move A toward 12.33; the simple-pole test remains open for other reasons (missing graphs, further B finite pieces, or a different assembly operator).
