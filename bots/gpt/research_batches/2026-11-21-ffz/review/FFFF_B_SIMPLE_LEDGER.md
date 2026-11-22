# Projected Newton simple-pole ledger

Target at \(a_X=0\), \(a_Y=-1\): improvement-independent three-loop simple pole \(=1/24\) in units \(h^2 z/K\), at \(L=0\). Higher poles already cancel after the FFff(B) overlap fix (`FFFF_B_OVERLAP_FIX.md`).

Runnable: `ffz_ffff_B_simple_ledger.py`.

## Known \(L=0\) simple pieces

| Piece | Value | Status |
| --- | --- | --- |
| Scalar aggregate (QQ+QCc+CCcc) | \(-0.05121964552\) | Packet (`ffz_qcc_full_regulated_check.md`) |
| QFf | \(-0.07466151326053\) | Exact finite (`ffz_qff_curvature_pole.md`) |
| CcFf | \(\approx -0.02480\) | Review QMC (checkpoint \(-0.0247\)) |
| FFff(A) | \(\approx 0.148962\) | \(4\times\) regulated inference \(0.03724\) (spin fix) |
| CT constant | \(\approx 0.085120\) | \((C_0+6)/18+(1-\zeta(2))/24\) |
| **Subtotal without B** | \(\approx 0.083401\) | |
| **B simple needed** | \(\approx -0.041734\) | \(1/24-\)subtotal |

\(C_0=-15/2+2\sqrt3\,\mathrm{Cl}_2(\pi/3)\). Scale logs \(L,L^2\) cancel separately against the FFff higher-pole budget and CT (`COUNTERTERM_BUDGET.md`).

## FFff(B) conversion

Projective \(I=a_0/\varepsilon^2+a_1/\varepsilon+c_0\) with overlap-fixed \(a_0=-48\), \(a_1=208\):

\[
\text{B simple}_{L=0}
= \Bigl(-\frac{163}{1728}+\frac{\pi^2}{64}\Bigr)-\frac{c_0}{2304}.
\]

Field-weight regulated fits give \(c_0\approx 192\) (`FFFF_B_REGULATED.md`). At that value, B simple \(\approx -0.02345\).

## Extra vertex-ll finite: \(\Delta c_0=64/3\)

The same soft-soft vertex density that supplies \(\Delta a_1=-16\) contributes, through \(W^{-3\varepsilon}\) and \(U^{\varepsilon}\) on double-log sectors,

\[
\Delta c_0 = \int_{\mathrm{ll}}\rho\,(-3\ln W+\ln U_{\mathrm{red}}) = \frac{64}{3}
\]

exactly (quadrature recovers \(21.333\ldots\) to machine precision; residue still \(-16\)).

| \(c_0\) scenario | B simple | Gap vs needed |
| --- | --- | --- |
| \(192\) (field-weight only) | \(\approx -0.02345\) | \(\approx +0.01828\) |
| \(192+64/3\) (derived W/U) | \(\approx -0.03271\) | \(\approx +0.00903\) |
| \(192+2\cdot(64/3)\) (empirical; **not derived**) | \(\approx -0.04197\) | \(\approx -0.00023\) |

## Soft-linear (audited, excluded)

See `FFFF_B_SOFT_LINEAR_AUDIT.md`. Soft-linear of the extra vertex copy is quadrature-stable at \(\approx-105.304\) (U-rule shift \(+30\)), but it is a **single**-soft channel. The overlap fix only corrects two-log soft-soft, so soft-linear is **not** added to \(c_0\). Including it opens a \(\sim+0.055\) gap.

## What is still open

1. Scheme-level identification that the missing piece is exactly one overlap copy of this operator (see `FFFF_B_OVERLAP_RIGOR.md`; global doubling excluded by pl \(=-48\)).
2. Projected simple versus \(1/24\): FFff(A) re-extraction gives \(c_0^A=15.11\pm 0.58\) (`FFFF_A_SIMPLE_REEXTRACT.md`), **not** \(\approx 12.33\). Gap \(\approx+0.01\) is not closed by A.
3. Any remaining finite pieces or missing graphs beyond the present ledger.

Homogeneity still predicts \(+L/96+3L^2/32\) in bare B simple from the fixed double/triple.
