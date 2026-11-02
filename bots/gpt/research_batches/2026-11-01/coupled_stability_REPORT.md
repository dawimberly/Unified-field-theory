# Coupled field stability test
1 November 2026. Candidate 0, amplitude plus constrained orientation, constant phase.

## Outcome
The requested numerical test was run. It has not established a stable, localized continuum field solution. The coarse relaxation approaches the uniform vacuum, and grid-sensitive loss of winding prevents treating the results as physical instability or measurement collapse. Do not interpret energy minimization as the particle physically disappearing.

## What was tested
The existing Candidate-0 periodic cubic-lattice static energy was minimized with both rho and all tangent directions of n free. Unit n was enforced by normalization with its exact chain-rule gradient. Parameters: v=1, lambda=64, kappa=1, box side 20, compact Hopf seed radius 5; grids N=32 and N=48. These are dimensionless choices inherited from earlier work, not measured particle parameters. The intended continuum seed has charge one, but the coarse initial estimates already show substantial under-resolution.

The rho lower bound was 0.05 to avoid the singular polar region. Whether it became active is reported below. Periodic boundaries differ from an isolated infinite-space solution; this is a computational gate, not a continuum theorem. No topology penalty or winding constraint was inserted to force success. The optimizer was L-BFGS-B, capped at 600 iterations per grid. Relaxation is an energy-search method, not Hamiltonian time evolution.

| Grid | Initial energy | Final energy | Initial Hopf estimate | Final Hopf estimate | Maximum physical static-force residual | Iterations |
|---|---:|---:|---:|---:|---:|---:|
| 32^3 | 412.651118 | 2.16178947e-09 | 0.544566 | 1.12601e-25 | 7.27575e-05 | 600 |
| 48^3 | 452.720281 | 49.4029469 | 0.762382 | 0.00418006 | 10.5248 | 600 |

N=32: optimizer success=False; message=STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT; minimum rho=0.99999995; floor active=False; minimum final neighboring orientation dot product=1; elapsed=122.1 seconds.

N=48: optimizer success=False; message=STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT; minimum rho=0.678601442; floor active=False; minimum final neighboring orientation dot product=-0.949338695; elapsed=406.6 seconds.

The force residual is the maximum of the amplitude derivative and projected orientation derivative, divided by lattice volume h^3. Optimizer convergence in a normalized coordinate chart is not sufficient by itself; the physical residual must also be small. Neither unfinished optimization nor small energies should be promoted to a stationary-knot result.

## Independent checks completed
- Static energy was cross-checked against the existing real-time solver at the initial seed on each grid.
- Directional energy derivatives agreed with the analytic gradient. Individual errors are in the JSON files.
- The full amplitude-orientation constrained Hessian-vector product was derived, including the sphere-constraint curvature term -(n dot grad_n E) delta n.
- Hessian-vector finite-difference relative errors: amplitude 1.11e-10, orientation 1.92e-10; bilinear symmetry error 1.61e-16.
- All 648 amplitude/orientation vacuum eigenvalues on an independent 6^3 lattice matched the Fourier prediction within 4.83e-13. Negative eigenvalues: 0; orientation constant-mode zeros: 2.
- At constant theta the quadratic phase contribution is integral rho^2 |grad delta theta|^2 / 2, so it is nonnegative, with no quadratic cross terms with rho or n. Its omission from the static minimizer does not hide a negative phase-energy direction for positive rho.
- Five perturbation directions, at two finite-difference sizes, were probed at each endpoint: amplitude, orientation, mixed, translation-like, and internal rotation. These are directional curvature diagnostics, not a complete knot eigenvalue spectrum. At a nonstationary endpoint they are not normal-mode stability evidence.

## Exact vacuum result
At rho=1 and constant n and theta, the continuum amplitude branch has omega^2=k^2+128. The two orientation branches and phase branch have omega^2=k^2. Thus the homogeneous vacuum is linearly stable in this model. Zero wave-number orientation and phase modes express the continuous internal degeneracy; they are not growing instabilities. A stable vacuum is not a localized particle.

## Why the winding diagnostic matters
For a smooth continuum field n, an allowed smooth evolution or homotopy cannot simply erase the protected Hopf class under the required boundary assumptions. The finite lattice supplies vectors only at nodes and does not guarantee a smooth topology-preserving interpolant. Winding can slip through unresolved cells. The rho field does not remove this protection in the written model because n is an independent field, defined even when rho vanishes.

Consequently, disappearance of the numerical Hopf estimate does not establish physical unwinding. Its absolute value near the vacuum is tiny; the relative longitudinal-B diagnostic is not meaningful when B is essentially zero. Agreement that energy decreases is weaker than spatial convergence of the field or its topology.

## Decision
The current lattice calculation fails the prerequisite for a claimed full stability spectrum of a nontrivial object: a resolved, stationary, winding-preserving solution. The Hessian implementation is ready and checked, but calculating and interpreting a complete knot spectrum at an unresolved or nonstationary endpoint would give false precision.

This rejects the claim that our existing prepared seed and simple size minimum already demonstrate stable matter. It does not reject every continuum solution of Candidate 0, and it supplies no evidence of gravity-driven collapse. A dissipative optimization trajectory is not a physical decay prediction.

The next justified step is a topology-controlled discretization or an adaptive/refined stationary solve with a validated interpolant and residual, followed by the generalized constrained eigenproblem H delta q=omega^2 M delta q. M is the positive kinetic matrix in the existing action when rho>0. Under those conditions a genuine negative constrained Hessian direction would imply a linear instability; positive sampled curvatures cannot prove the absence of all such directions.

## Reproduction and files
Run coupled_stability_test_2026_11_01.py --N 32 --iterations 600, then repeat with --N 48. Run coupled_hessian_check_2026_11_01.py and coupled_vacuum_spectrum_2026_11_01.py, then this report generator. NumPy, SciPy, Matplotlib and threadpoolctl are required. The solver import accepts the copied dependency in the archive directory via Python's usual script-directory search; it also includes the original workspace scripts directory for the current environment.

Original dependency: C:/workspace/physics/scripts/knotted_cloud_dynamics.py. Original model: C:/workspace/physics/experiments/hydrogen_cloud_workbench/CANDIDATE_ACTION.md. No original research sources were edited. Saved endpoint NPZ files are candidate arrays, not certified stationary solutions. No new empirical data, gravity solver or quantum dynamics was used.