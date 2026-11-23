"""Archive the topology-preserving follow-up and its scope limits."""
from pathlib import Path
import json,zipfile,shutil,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent;out=root/'smooth_hopf_2026_11_01'
load=lambda name:json.loads((out/name).read_text())
base=load('radial_control.json');s1=load('order1.json');s2=load('order2.json');ext=load('extended_summary.json');counter=load('radial_nonstationarity_countercheck.json');ic=load('implementation_checks.json');dc=load('independent_spatial_derivative_checks.json');omitted=load('extended_omitted_mode_countercheck.json')
final=load(ext['final_label']+'.json')
rows=[('Radius and radial amplitude',6,base['energy']),('Linear orientation modes',18,s1['energy']),('Quadratic orientation modes',36,s2['energy']),('Cubic orientation and quadratic amplitude modes',ext['parameter_count'],ext['energy'])]
lines=['| Family | Parameters | Dimensionless energy |','|---|---:|---:|']+[f'| {name} | {n} | {E:.10f} |' for name,n,E in rows]
ref=['| Integration grid | Energy | Maximum parameter gradient |','|---|---:|---:|']
for rr in ext['quadrature_refinement']:ref.append(f"| {rr['grid']} | {rr['energy']:.10f} | {rr['max_gradient']:.6g} |")
min_ev=min(ext['restricted_hessian_eigenvalues']);max_ev=max(ext['restricted_hessian_eigenvalues'])
reduction=100*(base['energy']-ext['energy'])/base['energy']
report=f'''# Smooth, topology-preserving field follow-up
1 November 2026. Candidate 0, infinite flat three-space, v=kappa=1 and lambda=64.

## Main result
The earlier grid's artificial unwinding was avoided by switching to globally smooth, topology-preserving trial fields. The radius/amplitude-only optimum is demonstrably not stationary once smooth orientation deformations are admitted: an explicit small change lowers its energy by approximately 0.20330045, consistently across three integration grids.

The enlarged family reaches energy {ext['energy']:.10f}, a {reduction:.4f}% reduction from the radius/radial-amplitude control. This is a better candidate configuration in the same continuum action, not a certified stationary solution of the full field equations or a new physical particle. The omitted-mode and integration checks below determine how far the result can be trusted.

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

'''+ '\n'.join(lines)+f'''

The extended endpoint has maximum parameter gradient {ext['max_parameter_gradient']:.6g}; active parameter bounds: {ext['bound_active']}. Its final optimizer termination is {final['message']}. Minimum sampled amplitude is {final['min_rho']:.8g}. Parameter bounds limit the search; no conclusion extends beyond the tested family or their range. No model parameter was fitted to an observed particle mass.

## Explicit countercheck of the radius-only claim
At its own optimum the six-parameter radius/amplitude control has a small retained gradient, but a large omitted orientation gradient. The deformation used here is a_z(y)=epsilon*y_z*exp(-|y|^2/4), with epsilon=0.01 and other orientation coefficients zero. It is smooth, localized and topology preserving. Numerical results:

| Integration grid | Energy change |
|---|---:|
'''
for rr in counter['rows']:report+=f"| {rr['grid']} | {rr['energy_change']:.12f} |\n"
report+='''
This is a numerically resolved first-variation counterexample, not an interval-certified proof. It is stronger evidence against calling the radius-only configuration an equilibrium than the earlier lattice unwinding, because no topology change is involved.

## Integration and omitted-mode checks

'''+ '\n'.join(ref)+f'''

These rows include pre-refinement results if a refined reoptimization was required. The final label is {ext['final_label']}; only rows identified as after-refinement are evaluations of that new endpoint when applicable. Comparisons across different endpoints must not be mistaken for quadrature errors.

Embedding the final parameters into an even larger family gives maximum NEW amplitude-mode gradient {ext['new_amplitude_gradient']:.6g} and maximum NEW orientation-mode gradient {ext['new_orientation_gradient']:.6g}. These directions were absent from the optimization. Non-negligible new-mode gradients mean the configuration is still not stationary in the larger family, regardless of a positive Hessian in the smaller one. The larger-family evaluation uses its own reported integration grid and is not itself a further optimized solution.

## Restricted Hessian
The Hessian was calculated by central differences of the analytic parameter gradient, with step {ext['hessian_difference_step']}. Its relative antisymmetry is {ext['hessian_symmetry_relative_error']:.3g}. Its eigenvalues range from {min_ev:.8g} to {max_ev:.8g}, with {ext['negative_restricted_eigenvalues_below_minus_1e_5']} below -1e-5 and {ext['zero_like_restricted_eigenvalues']} of absolute size below 1e-5.

These are curvatures in the finite chosen parameter space, not frequencies of the full continuum field. A positive restricted Hessian does not prove dynamical stability. A negative eigenvalue at a nonstationary endpoint is not by itself a normal-mode instability of a stationary solution. We have not solved the complete generalized field eigenproblem with its kinetic metric. A second full Hessian-difference step was not used. The smallest-eigenvalue direction was independently checked using symmetric energy differences at steps 0.001 and 0.0005; details are in smallest_mode_direct_check.json. These curvatures are not continuum oscillation frequencies.

## Implementation controls
The undeformed frozen-amplitude Hopf energy matches 16*pi^2*(R+1/R), with maximum reference integral relative error {ic['frozen_hopf_integral_relative_error']:.3g}. Unit-vector error is {ic['unit_vector_max_error']:.3g}, and the spatial tangent condition error is {ic['tangent_derivative_max_error']:.3g}. Analytic parameter gradients were compared against independent complex-step gradients; maximum check error is {ic['analytic_vs_complex_step_max_error']:.3g}. A random directional finite-difference check also passed.

A separately written quaternion formula and centered spatial differences checked 20 sample locations in the deformed field. The largest orientation-derivative discrepancy was {max(x['orientation_max_error'] for x in dc):.3g}; the largest amplitude-derivative discrepancy was {max(x['amplitude_max_error'] for x in dc):.3g}. These are implementation controls, not empirical tests.

## Scientific conclusion
We can now lower the energy without losing topology numerically. The old radius-only equilibrium claim fails a specific smooth deformation test. The expanded configurations are better candidates, but the explicit omitted-mode gradients and quadrature checks remain the appropriate gate for further work. The mathematically justified next step is targeted basis enlargement or a full stationary solve initialized from the smooth candidate, followed by a full constrained stability calculation after stationarity is established.

This work addresses the stable-field-structure question only. It does not derive quantum measurement, gravity, Einstein's equations, particle spin/statistics or an observed mass. The current action still has no dynamical metric. Making a better classical knot is not, by itself, progress on the origin of gravity. No user choice is required to interpret that limitation.

## Reproduction
Run smooth_hopf_variational_2026_11_01.py, smooth_hopf_derivative_check_2026_11_01.py, smooth_hopf_radial_countercheck_2026_11_01.py and smooth_hopf_extended_2026_11_01.py, then smooth_hopf_omitted_check_2026_11_01.py, smooth_hopf_smallest_mode_check_2026_11_01.py and this packaging script. Dependencies: Python, NumPy, SciPy, Matplotlib, threadpoolctl. The source now uses the checked analytic parameter gradient; original coarse runs initially used complex-step gradients, so iteration histories can differ on reproduction while endpoint comparisons and checks remain the relevant criteria. Reproduction overwrites only this follow-up's output directory.
'''
report+='\n## Direct omitted-mode perturbation check\n\nThe strongest newly added direction was tested with three small step sizes on three integration grids. These perturbations preserve topology by the same smooth construction.\n\n| Grid | Step | Energy change |\n|---|---:|---:|\n'
for rr in omitted['rows']:report+=f"| {rr['grid']} | {rr['step']:.6g} | {rr['energy_change']:.12g} |\n"
report_path=root/'smooth_hopf_followup_2026_11_01.md';report_path.write_text(report)
fig,ax=plt.subplots(figsize=(8,4.5));ax.bar(range(len(rows)),[r[2] for r in rows],color=['#8b949e','#388bfd','#388bfd','#2ea043']);ax.set_xticks(range(len(rows)),['Radius + amplitude','18 parameters','36 parameters',f"{ext['parameter_count']} parameters"]);ax.set(ylabel='Dimensionless continuum energy',title='Energy decreases within a topology-preserving family');fig.tight_layout();fig.savefig(out/'family_energy.png',dpi=150);plt.close(fig)
scripts=['smooth_hopf_variational_2026_11_01.py','smooth_hopf_derivative_check_2026_11_01.py','smooth_hopf_radial_countercheck_2026_11_01.py','smooth_hopf_extended_2026_11_01.py','smooth_hopf_omitted_check_2026_11_01.py','smooth_hopf_smallest_mode_check_2026_11_01.py','package_smooth_hopf_2026_11_01.py']
for name in scripts:shutil.copyfile(root/name,out/name)
shutil.copyfile(report_path,out/report_path.name)
shutil.copyfile(Path(r'C:\workspace\physics\experiments\hydrogen_cloud_workbench\CANDIDATE_ACTION.md'),out/'CANDIDATE_ACTION.md')
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file() and p.name!='manifest.json'};(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
with zipfile.ZipFile(root/'smooth_hopf_followup_2026_11_01.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in out.iterdir():
  if p.is_file():z.write(p,p.name)
print(json.dumps(dict(report=str(report_path),energy=ext['energy'],reduction_percent=reduction,extra_amplitude_gradient=ext['new_amplitude_gradient'],extra_orientation_gradient=ext['new_orientation_gradient'],hessian_min=min_ev,negative_modes=ext['negative_restricted_eigenvalues_below_minus_1e_5'])))
