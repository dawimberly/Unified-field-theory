# Fixed-mass density-arrangement scan

**Conversation date range** $2026-10$ (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records fixed-mass density-arrangement scan. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Completed 45 cases; 0 failed solves. Same continuous mass, radius, ambient density and scalar coupling in each comparison. Only core/mantle density arrangement changes. Model: $n=1$ chameleon scalar, mean body density contrast 1000.

| mu | Core radius / radius | Core density / mean | Exterior charge change |
|---|---:|---:|---:|
| **DERIVED** 0.1 | 0.3 | 0.25 | +0.3484% |
| **DERIVED** 0.1 | 0.3 | 0.5 | +0.2338% |
| **DERIVED** 0.1 | 0.3 | 1.0 | +0.0000% |
| **DERIVED** 0.1 | 0.3 | 1.5 | -0.2404% |
| **DERIVED** 0.1 | 0.3 | 1.9 | -0.4377% |
| **DERIVED** 0.1 | 0.6 | 0.25 | +2.9051% |
| **DERIVED** 0.1 | 0.6 | 0.5 | +2.0486% |
| **DERIVED** 0.1 | 0.6 | 1.0 | +0.0000% |
| **DERIVED** 0.1 | 0.6 | 1.5 | -2.6302% |
| **DERIVED** 0.1 | 0.6 | 1.9 | -5.2718% |
| **DERIVED** 0.1 | 0.8 | 0.25 | +7.2304% |
| **DERIVED** 0.1 | 0.8 | 0.5 | +5.4336% |
| **DERIVED** 0.1 | 0.8 | 1.0 | +0.0000% |
| **DERIVED** 0.1 | 0.8 | 1.5 | -7.9870% |
| **DERIVED** 0.1 | 0.8 | 1.9 | -15.7119% |
| **DERIVED** 0.2 | 0.3 | 0.25 | +0.1149% |
| **DERIVED** 0.2 | 0.3 | 0.5 | +0.0771% |
| **DERIVED** 0.2 | 0.3 | 1.0 | +0.0000% |
| **DERIVED** 0.2 | 0.3 | 1.5 | -0.0790% |
| **DERIVED** 0.2 | 0.3 | 1.9 | -0.1437% |
| **DERIVED** 0.2 | 0.6 | 0.25 | +0.9794% |
| **DERIVED** 0.2 | 0.6 | 0.5 | +0.6889% |
| **DERIVED** 0.2 | 0.6 | 1.0 | +0.0000% |
| **DERIVED** 0.2 | 0.6 | 1.5 | -0.8862% |
| **DERIVED** 0.2 | 0.6 | 1.9 | -1.8049% |
| **DERIVED** 0.2 | 0.8 | 0.25 | +2.6014% |
| **DERIVED** 0.2 | 0.8 | 0.5 | +2.0067% |
| **DERIVED** 0.2 | 0.8 | 1.0 | +0.0000% |
| **DERIVED** 0.2 | 0.8 | 1.5 | -5.3628% |
| **DERIVED** 0.2 | 0.8 | 1.9 | -17.4865% |
| **DERIVED** 0.3 | 0.3 | 0.25 | +0.0732% |
| **DERIVED** 0.3 | 0.3 | 0.5 | +0.0491% |
| **DERIVED** 0.3 | 0.3 | 1.0 | -0.0000% |
| **DERIVED** 0.3 | 0.3 | 1.5 | -0.0503% |
| **DERIVED** 0.3 | 0.3 | 1.9 | -0.0915% |
| **DERIVED** 0.3 | 0.6 | 0.25 | +0.6259% |
| **DERIVED** 0.3 | 0.6 | 0.5 | +0.4398% |
| **DERIVED** 0.3 | 0.6 | 1.0 | +0.0000% |
| **DERIVED** 0.3 | 0.6 | 1.5 | -0.5629% |
| **DERIVED** 0.3 | 0.6 | 1.9 | -1.1432% |
| **DERIVED** 0.3 | 0.8 | 0.25 | +1.6726% |
| **DERIVED** 0.3 | 0.8 | 0.5 | +1.2875% |
| **DERIVED** 0.3 | 0.8 | 1.0 | +0.0000% |
| **DERIVED** 0.3 | 0.8 | 1.5 | -3.5840% |
| **DERIVED** 0.3 | 0.8 | 1.9 | -17.6088% |

## Controls and scope

**DERIVED** Core density ratio 1 is the uniform control. Positive changes mean a larger exterior scalar charge than the uniform sphere. The mantle density is set to $1000(1-f^3\mathrm{ratio})/(1-f^3)$, enforcing equal continuous mass. All sampled densities are positive.

**DERIVED** Every successful solve checks adaptive collocation convergence, positive scalar field, equal continuous mass and integrated exterior flux balance. The two most extreme responses were rerun at tenfold tighter tolerance and doubled outer radius; each charge changes by less than 0.02%. Charge is extracted five ambient Compton lengths beyond the surface.

**DERIVED** This calculation concerns an added scalar force. It does not test intrinsic spin, piezoelectricity, chemistry at fixed density, a varying Einstein gravitational constant, or quantum-gravity unification. The parameters were sampled, not optimized over all profiles. The two extrema are extrema of this finite scan only.

## Next step

**ASSUMED** Use density-matched materials as a null control: this density-only model predicts identical scalar profiles for identical density profiles regardless of material name. A material/spin claim needs a specified additional coupling and independent laboratory constraints before assigning numerical predictions.

## Reproduction

**ASSUMED** Run python C:\workspace\physics\scripts\screening_density_scan.py (NumPy and SciPy required). Full results include solver diagnostics; refinements.json contains extreme-case reruns.
