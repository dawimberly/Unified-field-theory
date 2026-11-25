# First real-time knotted-cloud tests

**Conversation date range** $2026-10$ (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records first real-time knotted-cloud tests. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** This is the first constrained time-evolution implementation of the orientation field n together with the amplitude rho in Candidate 0. The phase theta is constant, which is a consistent sector. No drag, bath, external force, preferred-phase force, or gravitational dynamics is added. The calculations are classical field dynamics, not quantum particles or a derivation of gravity.

## What changed from the amplitude test

**DERIVED** The knotted field and its full Skyrme kinetic energy are now active. In lattice notation $b_i=D_i n \times n, F_{0i}=\dot n\cdot b_i$. The orientation kinetic matrix is

**DERIVED** $M=\rho ^2 I+\kappa  \sum _i b_i b_i^T.$

**DERIVED** The acceleration includes the time derivative of M and enforces $n \cdot n_{ddot}=-$|n_dot|$^{2}$. Spatial forces are obtained by varying the discrete Lagrangian, including explicit occurrences of n in F. This avoids substituting a simpler wave equation for the stated action.

## Initial data and numerical method

**DERIVED** Compact smooth SU(2) quaternion lifts create a single Hopf seed and a product pair. Their intended continuum Hopf numbers are 1 and 2. Each lift has radial profile $f(r)=\pi [1-10s^{3}+15s^{4}-6s^{5}]$ for $s=r/R<1$, and zero beyond R. The target field n is the lifted rotation of a constant unit vector. The second lift is optionally conjugated by an internal z rotation, with angles 0 and pi tested separately.

**DERIVED** rho starts at 1 and all velocities start at zero. These are prepared seeds, not relaxed single-soliton ground states. Pair overlap energy is not a binding energy relative to minimized isolated objects. The chosen relative orientations are initial conditions; the model has not selected them dynamically.

**DERIVED** Free dimensionless parameters are $v=1, \lambda =64, \kappa =1$, compact support $radius R=5$, center separation 7, periodic box side 20, and $N=64$ per axis $(dx=0.3125)$. $\Lambda =64$ is a stiff-amplitude test choice to keep rho away from its polar boundary; it is not inferred from physics. It differs from the earlier $\lambda =1$ static example and $\lambda =0.5$ scalar-only experiment. No physical mass, length or time calibration is supplied.

**DERIVED** Spatial derivatives use forward lattice differences with their exact adjoints in the force. Projected RK4 evolves the nonseparable kinetic system, with $dt=0.02$. Normalizing n and projecting velocity back onto its tangent plane correct integrator errors; no energy-minimization step is used during time evolution.

## Completed trajectories

| Run | Duration | Initial energy | Final energy | Max relative drift | Initial approximate Hopf Q | Final approximate Hopf Q |
|---|---:|---:|---:|---:|---:|---:|
| **DERIVED** single64 | 1.0 | 468.943045 | 468.942205 | $1.79\times10^{-06}$ | 0.858259 | 0.886056 |
| **DERIVED** pair64_phase0 | 0.5 | 956.991733 | 956.990927 | $8.43\times10^{-07}$ | 1.707411 | 1.738358 |
| **DERIVED** pair64_phasepi | 0.5 | 918.864515 | 918.863777 | $8.03\times10^{-07}$ | 1.718954 | 1.750851 |

**DERIVED** These are short transient runs. They do not establish long-time binding, stable isolated particles, or reduced gravitational radiation. The single trajectory is a control for seed reorganization. Pair trajectories can also reorganize because the initial seeds are not stationary.

**DERIVED** Final kinetic-energy fractions (all start at zero): single64: 29.46%; pair64_phase0: 12.12%; pair64_phasepi: 12.56%. The structures are dynamically reorganizing, not demonstrated to have settled into quieter states.

## Numerical validation and topology limitations

**DERIVED** The semi-discrete instantaneous energy derivative, the unit-vector acceleration constraint, and an independent finite-difference directional derivative of static energy pass before evolution. Every accepted trajectory checks total energy drift below 0.2%, unit-vector error below $1\times10^{-12}$, and tangent-velocity error below $1\times10^{-11}$. These are implementation checks, not tests of physical truth.

**DERIVED** An independent coarse-grid integrator refinement gives a successive-state error ratio of 13.440 when dt is halved twice, consistent with approaching fourth-order time integration. That coarse grid is not used to validate topological charge.

**DERIVED** Hopf charge is estimated by centered-difference F and a Fourier Coulomb-gauge inverse curl. It is not an exactly conserved lattice integer. On initial seeds, refining $N=40$ to $N=96$ improves the single estimate from about 0.677 to 0.934 and the phase-zero pair from 1.336 to 1.864. $N=64$ still underestimates the expected values by about 14%. The longitudinal part of the estimated magnetic-like two-form and nearest-neighbor dot products are recorded as resolution diagnostics. This two-form is not identified with electromagnetism.

**DERIVED** Consequently small changes in estimated charge are not proof of topological unwinding or protection. Spatially converged time-evolving topology and relaxed single-soliton controls remain required before trusting pair binding. Initial-grid refinement alone is not dynamical grid convergence.

## Gravity status

**DERIVED** The total field energy's trace-free second moment is recorded as a quadrupole diagnostic. No gravitational-wave power is calculated. Neither scalar/orientation waves nor changes in the approximate Hopf number are automatically gravitons. This action has no dynamical metric; a weak-field Einstein radiation calculation would import gravity as a reference and require a justified isolated-source regime and full stress-energy accounting.

**DERIVED** The immediate next numerical gate is a spatially resolved isolated knot that remains stable over several intrinsic time scales, followed by pair trajectories with identical parameters and separated controls. A persistent lower-energy composite and outgoing field energy would then be meaningful evidence for binding within this candidate. Quantum-gravity unification would still require the missing gravitational and gauge sectors.
