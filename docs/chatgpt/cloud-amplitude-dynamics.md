# Conservative amplitude-cloud evolution

**Conversation date range** 2026-10 (exact conversation dates UNVERIFIED)

**DERIVED** This checkpoint records conservative amplitude-cloud evolution. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**UNVERIFIED** Result: the tested amplitude-cloud pairs disperse; no persistent localized bound pair is demonstrated during this interval. About 89% of initial energy is near the original centers, falling to about 3% by the end. The main pair energy RMS radius grows from 3.69 to 10.14. The initial amplitudes 0.15, 0.30 and 0.45 give the same qualitative outcome, as do the refined time step and grid. Quieting locally is explained by outward spreading of field energy, not established binding or gravitational radiation suppression.

**DERIVED** This is actual Hamiltonian time evolution of an exact restricted sector of the existing Candidate 0 action. It is not gradient descent and contains no damping, bath, gravity, preferred relative phase, or forced settling.

**DERIVED** With n and theta constant, the equation is rho_tt = laplacian(rho) - lambda rho(rho^2-v^2). The tested trajectories remain rho>0. The initial field is rho=v-A sum exp(-|x-x_i|^2/(2 sigma^2)), with zero initial velocity. These are prepared amplitude depressions, not stationary particle solutions or Hopfions. Topological charge is zero.

**DERIVED** v=1, lambda=0.5, sigma=1.5, center separation=6, box side=48, duration=16; all values are dimensionless free inputs. The massive amplitude frequency about vacuum is 1. Six-neighbor spatial differences and velocity Verlet evolve a periodic lattice Hamiltonian. Boundary energy is monitored; the box is large enough to avoid significant returning radiation during the tested interval.

| Run | Initial energy | Final localized fraction | Max relative energy drift |
|---|---:|---:|---:|
| **DERIVED** pair_main | 2.5201435 | 0.029524 | 0.00138 |
| **DERIVED** single_control | 1.2644154 | 0.009798 | 0.00141 |
| **DERIVED** pair_small | 0.66251074 | 0.029307 | 0.00144 |
| **DERIVED** pair_large | 5.3933511 | 0.030042 | 0.00134 |
| **DERIVED** pair_half_dt | 2.5201435 | 0.029557 | 0.000347 |
| **DERIVED** pair_fine | 2.5326183 | 0.029541 | 0.000618 |

**DERIVED** Initial pair energy minus two isolated controls: -0.0086873241. This comparison concerns overlapping prepared configurations, not optimized binding energy.

**DERIVED** Temporal localization difference: 3.26857e-05; spatial difference: 1.70008e-05.

**DERIVED** Localization is the energy fraction inside the union of radius-3 balls around the original centers. This is a diagnostic, not a definition of binding: a moving bound object could leave those regions. RMS radius and field slices provide additional dispersion diagnostics. The single control tests dispersion unrelated to interaction.

**DERIVED** The energy-derived quadrupole is recorded for inspection only. No gravitational-wave power is reported: these relativistic field waves need not satisfy the slow-motion assumptions of the elementary quadrupole formula, and the finite periodic field box is not an isolated asymptotically flat gravitational solution. Scalar field radiation is not gravitational radiation.

**DERIVED** A finite-duration dispersing run cannot prove that no bound state, oscillon, or topological composite exists for other initial data. These calculations do not evolve the n field and therefore do not answer the two-Hopfion binding question. They do test whether these amplitude-only initial clouds form persistent localized structures without external damping.

**DERIVED** Numerical consistency checks: 21/21 passed.

**DERIVED** Next full test requires a constrained time-evolution solver for n, including the Skyrme term's field-dependent kinetic matrix, resolved individual Hopfions, topology diagnostics and paired isolated controls. Adding a hand-written drag term would not solve that problem.

**DERIVED** Local action source: experiments/hydrogen_cloud_workbench/CANDIDATE_ACTION.md. No new physical force or empirical constant is inferred.
