"""Package the completed static stability gate, including numerical limitations."""
from pathlib import Path
import json,zipfile,shutil,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent
out=root/'coupled_stability_2026_11_01'
results=[json.loads((out/f'N{N}_results.json').read_text()) for N in [32,48]]
hc=json.loads((out/'hessian_implementation_checks.json').read_text())
vc=json.loads((out/'vacuum_spectrum_check.json').read_text())
rows=[]
for r in results:
 rows.append(f"| {r['N']}^3 | {r['initial_energy']:.9g} | {r['final_energy']:.9g} | {r['initial_topology']['hopf_estimate']:.6g} | {r['final_topology']['hopf_estimate']:.6g} | {r['stationarity_max_force']:.6g} | {r['iterations']} |")
fig,ax=plt.subplots(figsize=(8,4.5))
for r in results:
 hist=[dict(iteration=0,energy=r['initial_energy'])]+r['history']
 ax.semilogy([x['iteration'] for x in hist],[max(x['energy'],1e-16) for x in hist],label=f"{r['N']} cubed grid")
ax.set(xlabel='Energy-minimization iteration (not physical time)',ylabel='Dimensionless lattice energy',title='Coupled field relaxation: resolution gate')
ax.legend();fig.tight_layout();fig.savefig(out/'relaxation_energy.png',dpi=160);plt.close(fig)
report='''# Coupled field stability test
1 November 2026. Candidate 0, amplitude plus constrained orientation, constant phase.

## Outcome
The requested numerical test was run. It has not established a stable, localized continuum field solution. The coarse relaxation approaches the uniform vacuum, and grid-sensitive loss of winding prevents treating the results as physical instability or measurement collapse. Do not interpret energy minimization as the particle physically disappearing.

## What was tested
The existing Candidate-0 periodic cubic-lattice static energy was minimized with both rho and all tangent directions of n free. Unit n was enforced by normalization with its exact chain-rule gradient. Parameters: v=1, lambda=64, kappa=1, box side 20, compact Hopf seed radius 5; grids N=32 and N=48. These are dimensionless choices inherited from earlier work, not measured particle parameters. The intended continuum seed has charge one, but the coarse initial estimates already show substantial under-resolution.

The rho lower bound was 0.05 to avoid the singular polar region. Whether it became active is reported below. Periodic boundaries differ from an isolated infinite-space solution; this is a computational gate, not a continuum theorem. No topology penalty or winding constraint was inserted to force success. The optimizer was L-BFGS-B, capped at 600 iterations per grid. Relaxation is an energy-search method, not Hamiltonian time evolution.

| Grid | Initial energy | Final energy | Initial Hopf estimate | Final Hopf estimate | Maximum physical static-force residual | Iterations |
|---|---:|---:|---:|---:|---:|---:|
'''+ '\n'.join(rows)+'\n\n'
for r in results:
 report+=f"N={r['N']}: optimizer success={r['optimizer_success']}; message={r['message']}; minimum rho={r['minimum_rho']:.9g}; floor active={r['floor_active']}; minimum final neighboring orientation dot product={r['final_topology']['min_neighbor_dot']:.9g}; elapsed={r['elapsed_seconds']:.1f} seconds.\n\n"
report+='''The force residual is the maximum of the amplitude derivative and projected orientation derivative, divided by lattice volume h^3. Optimizer convergence in a normalized coordinate chart is not sufficient by itself; the physical residual must also be small. Neither unfinished optimization nor small energies should be promoted to a stationary-knot result.

## Independent checks completed
- Static energy was cross-checked against the existing real-time solver at the initial seed on each grid.
- Directional energy derivatives agreed with the analytic gradient. Individual errors are in the JSON files.
- The full amplitude-orientation constrained Hessian-vector product was derived, including the sphere-constraint curvature term -(n dot grad_n E) delta n.
'''
report+=f"- Hessian-vector finite-difference relative errors: amplitude {hc['amplitude_hvp_relative_error']:.3g}, orientation {hc['orientation_hvp_relative_error']:.3g}; bilinear symmetry error {hc['hessian_symmetry_relative_error']:.3g}.\n"
report+=f"- All {vc['modes']} amplitude/orientation vacuum eigenvalues on an independent 6^3 lattice matched the Fourier prediction within {vc['max_eigenvalue_error']:.3g}. Negative eigenvalues: {vc['negative_modes']}; orientation constant-mode zeros: {vc['zero_modes']}.\n"
report+='''- At constant theta the quadratic phase contribution is integral rho^2 |grad delta theta|^2 / 2, so it is nonnegative, with no quadratic cross terms with rho or n. Its omission from the static minimizer does not hide a negative phase-energy direction for positive rho.
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
'''
(out/'REPORT.md').write_text(report)
shutil.copyfile(Path(r'C:\workspace\physics\scripts\knotted_cloud_dynamics.py'),out/'knotted_cloud_dynamics.py')
shutil.copyfile(Path(r'C:\workspace\physics\experiments\hydrogen_cloud_workbench\CANDIDATE_ACTION.md'),out/'CANDIDATE_ACTION.md')
scripts=['coupled_stability_test_2026_11_01.py','coupled_hessian_check_2026_11_01.py','coupled_vacuum_spectrum_2026_11_01.py','package_coupled_stability_2026_11_01.py']
for s in scripts:shutil.copyfile(root/s,out/s)
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file() and p.name!='manifest.json'}
(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
with zipfile.ZipFile(root/'coupled_stability_2026_11_01.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in out.iterdir():
  if p.is_file():z.write(p,p.name)
print(report)
