# Smooth, topology-preserving field follow-up
1 November 2026. Candidate 0, infinite flat three-space, v=kappa=1 and lambda=64.

## Main result
The earlier grid's artificial unwinding was avoided by switching to globally smooth, topology-preserving trial fields. The radius/amplitude-only optimum is demonstrably not stationary once smooth orientation deformations are admitted: an explicit small change lowers its energy by approximately 0.20330045, consistently across three integration grids.

The enlarged family reaches energy 280.5673060785, a 9.7838% reduction from the radius/radial-amplitude control. This is a better candidate configuration in the same continuum action, not a certified stationary solution of the full field equations or a new physical particle. The omitted-mode and integration checks below determine how far the result can be trusted.

## What changed, and what did not
The previous computations supplied vectors at lattice nodes, allowing the interpolating field to become unresolved and lose its winding. This follow-up defines the field analytically at every point. Numerical quadrature integrates its energy; it does not define the topology. Spacetime is infinite flat R3 for this test, rather than a finite periodic box. The action and its couplings are unchanged, but the geometry, seed family and integration method differ, so absolute energies must not be compared directly to the old periodic compact-seed runs.

There is no new force, topology penalty, noise term, imposed collapse rule, gravitational dynamics or physical scale calibration. The constant phase is retained. Its quadratic perturbation energy is nonnegative for positive amplitude and has no mixed quadratic terms at a constant-phase background.

## Construction
Start with the degree-one rational Hopf texture n_H(y), y=x/R. Let a(y) be a finite combination of smooth polynomial-times-Gaussian basis functions. Apply its Cayley rotation:

    n(y)=[(1-a.a)n_H+2a(a.n_H)+2(a cross n_H)]/(1+a.a).

For real finite coefficients the denominator is at least one and n has unit length everywhere. The basis functions vanish at infinity. Replacing a by t*a continuously, 0<=t<=1, gives a smooth homotopy to the original Hopf map, preserving its unit Hopf class (with the inherited orientation convention). This does not force arbitrary dynamics to retain topology; it respects the smooth topological sector already present in the written independent-n model. It does restrict the trial family and cannot represent every field in that sector.

Amplitude is rho=exp[-sum c_b H_b(y)], hence strictly positive and tending to one at infinity. The first five H_b are radial Gaussians; the extended solve adds polynomial-Gaussian terms, permitting nonspherical amplitude. Spatial derivatives of both fields are analytic. There is no node-wise renormalization or singular interpolant.

After changing variables to y, the static energy is

    E=R*A+(1/R)*B+64*R^3*C,
    A=integral [|grad_y rho|^2+rho^2|grad_y n|^2]/2,
    B=integral sum_(i<j) [n.(partial_i n cross partial_j n)]^2/2,
    C=integral (rho^2-1)^2/4.

R is dimensionless. Once a is nonzero, it is a reference coordinate scale rather than a directly measured object radius. The angular field can change the object's shape and effective extent independently.

Radial integration uses r=tan(t) with Gaussian quadrature in 0<t<pi/2. Polar cosine uses Gaussian quadrature and azimuth uses a periodic trapezoidal rule. Integration-grid refinement is separate from variational-basis refinement.

## Completed minimizations

| Family | Parameters | Dimensionless energy |
|---|---:|---:|
| Radius and radial amplitude | 6 | 310.9944410337 |
| Linear orientation modes | 18 | 280.7897240288 |
| Quadratic orientation modes | 36 | 280.7745686302 |
| Cubic orientation and quadratic amplitude modes | 76 | 280.5673060785 |

The extended endpoint has maximum parameter gradient 3.14858e-05; active parameter bounds: False. Its final optimizer termination is CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH. Minimum sampled amplitude is 0.83613734. Parameter bounds limit the search; no conclusion extends beyond the tested family or their range. No model parameter was fitted to an observed particle mass.

## Explicit countercheck of the radius-only claim
At its own optimum the six-parameter radius/amplitude control has a small retained gradient, but a large omitted orientation gradient. The deformation used here is a_z(y)=epsilon*y_z*exp(-|y|^2/4), with epsilon=0.01 and other orientation coefficients zero. It is smooth, localized and topology preserving. Numerical results:

| Integration grid | Energy change |
|---|---:|
| [28, 12, 24] | -0.203300710453 |
| [40, 18, 36] | -0.203300456100 |
| [56, 24, 48] | -0.203300446974 |

This is a numerically resolved first-variation counterexample, not an interval-certified proof. It is stronger evidence against calling the radius-only configuration an equilibrium than the earlier lattice unwinding, because no topology change is involved.

## Integration and omitted-mode checks

| Integration grid | Energy | Maximum parameter gradient |
|---|---:|---:|
| [40, 18, 36] | 280.5672877894 | 0.000432055 |
| [56, 24, 48] | 280.5673060877 | 0.00108452 |
| [72, 30, 60] | 280.5673061703 | 0.00108286 |
| [72, 30, 60] | 280.5673061612 | 3.14306e-05 |

These rows include pre-refinement results if a refined reoptimization was required. The final label is extended_order3_refined; only rows identified as after-refinement are evaluations of that new endpoint when applicable. Comparisons across different endpoints must not be mistaken for quadrature errors.

Embedding the final parameters into an even larger family gives maximum NEW amplitude-mode gradient 9.23337e-08 and maximum NEW orientation-mode gradient 10.2393. These directions were absent from the optimization. Non-negligible new-mode gradients mean the configuration is still not stationary in the larger family, regardless of a positive Hessian in the smaller one. The larger-family evaluation uses its own reported integration grid and is not itself a further optimized solution.

## Restricted Hessian
The Hessian was calculated by central differences of the analytic parameter gradient, with step 0.0001. Its relative antisymmetry is 2.04e-08. Its eigenvalues range from 0.12029688 to 7764.1028, with 0 below -1e-5 and 0 of absolute size below 1e-5.

These are curvatures in the finite chosen parameter space, not frequencies of the full continuum field. A positive restricted Hessian does not prove dynamical stability. A negative eigenvalue at a nonstationary endpoint is not by itself a normal-mode instability of a stationary solution. We have not solved the complete generalized field eigenproblem with its kinetic metric. A second full Hessian-difference step was not used. The smallest-eigenvalue direction was independently checked using symmetric energy differences at steps 0.001 and 0.0005; details are in smallest_mode_direct_check.json. These curvatures are not continuum oscillation frequencies.

## Implementation controls
The undeformed frozen-amplitude Hopf energy matches 16*pi^2*(R+1/R), with maximum reference integral relative error 3.55e-15. Unit-vector error is 1.89e-15, and the spatial tangent condition error is 1.78e-15. Analytic parameter gradients were compared against independent complex-step gradients; maximum check error is 6.37e-12. A random directional finite-difference check also passed.

A separately written quaternion formula and centered spatial differences checked 20 sample locations in the deformed field. The largest orientation-derivative discrepancy was 1.19e-09; the largest amplitude-derivative discrepancy was 2.01e-10. These are implementation controls, not empirical tests.

## Scientific conclusion
We can now lower the energy without losing topology numerically. The old radius-only equilibrium claim fails a specific smooth deformation test. The expanded configurations are better candidates, but the explicit omitted-mode gradients and quadrature checks remain the appropriate gate for further work. The mathematically justified next step is targeted basis enlargement or a full stationary solve initialized from the smooth candidate, followed by a full constrained stability calculation after stationarity is established.

This work addresses the stable-field-structure question only. It does not derive quantum measurement, gravity, Einstein's equations, particle spin/statistics or an observed mass. The current action still has no dynamical metric. Making a better classical knot is not, by itself, progress on the origin of gravity. No user choice is required to interpret that limitation.

## Reproduction
Run smooth_hopf_variational_2026_11_01.py, smooth_hopf_derivative_check_2026_11_01.py, smooth_hopf_radial_countercheck_2026_11_01.py and smooth_hopf_extended_2026_11_01.py, then smooth_hopf_omitted_check_2026_11_01.py, smooth_hopf_smallest_mode_check_2026_11_01.py and this packaging script. Dependencies: Python, NumPy, SciPy, Matplotlib, threadpoolctl. The source now uses the checked analytic parameter gradient; original coarse runs initially used complex-step gradients, so iteration histories can differ on reproduction while endpoint comparisons and checks remain the relevant criteria. Reproduction overwrites only this follow-up's output directory.

## Direct omitted-mode perturbation check

The strongest newly added direction was tested with three small step sizes on three integration grids. These perturbations preserve topology by the same smooth construction.

| Grid | Step | Energy change |
|---|---:|---:|
| [40, 18, 36] | -0.01 | -0.0712392776394 |
| [40, 18, 36] | -0.001 | -0.00992683686439 |
| [40, 18, 36] | -0.0001 | -0.00102071001925 |
| [56, 24, 48] | -0.01 | -0.0712503069456 |
| [56, 24, 48] | -0.001 | -0.00992790552579 |
| [56, 24, 48] | -0.0001 | -0.00102081650778 |
| [72, 30, 60] | -0.01 | -0.071250305814 |
| [72, 30, 60] | -0.001 | -0.00992790438721 |
| [72, 30, 60] | -0.0001 | -0.00102081638352 |