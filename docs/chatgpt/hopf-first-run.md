# First Hopf stability experiment

**Conversation date range** 2026-10-25 to 2026-10-30 (archive window)

**ASSUMED** This checkpoint records first Hopf stability experiment. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**UNVERIFIED** Run: 2026-10-25. Exploratory numerical research, not evidence for a new physical theory.

## Model

**DERIVED** Unit-radius compact spatial S3, with a unit vector n in S2 at every mesh vertex. Dimensionless coefficients $a=b=1$. The seed is the analytic charge-one Hopf map after a stereographic dilation of 2. Relaxation is projected gradient descent with an energy-decreasing line search, capped at 400 steps. It is not physical time evolution.

**CITED** Energy: $E=\int[|Dn|^2/2+(|Dn|^4-\operatorname{tr}((Dn\,Dn^T)^2))/4]\,dV$. The control omits the quartic term. This is a reference Faddeev-Skyrme-type energy, not the original Possibility Cloud action. — eq./page UNVERIFIED.

**DERIVED** The tetrahedra use chordal piecewise-linear finite elements; vertex vectors have unit norm, but the interpolated field need not. This can permit artificial unwinding. Normalized spherical-face fluxes and a solved discrete potential give an approximate Hopf charge; it is not an exact protected lattice invariant.

**CITED** Reference: Faddeev and Niemi, Knots and Particles (1996), https://arxiv.org/abs/hep-th/9610193. No published result is claimed to have been replicated quantitatively. — eq./page UNVERIFIED.

## Results

| Nodes | Tetrahedra | Run | Initial energy | Final energy | Initial Q estimate | Final Q estimate | Final gradient RMS | Energy concentration |
|---:|---:|---|---:|---:|---:|---:|---:|---:|
| **DERIVED** 120 | 600 | hopf | 145.56 | $6.78688\times10^{-15}$ | 0.983813 | -0.000000 | $6.97\times10^{-09}$ | 0.529 |
| **DERIVED** 120 | 600 | perturbed_hopf | 145.874 | $1.34415\times10^{-14}$ | 0.984137 | -0.000000 | $6.68\times10^{-09}$ | 1.071 |
| **DERIVED** 120 | 600 | no_quartic_control | 50.6896 | $1.45625\times10^{-14}$ | 0.983813 | 0.000000 | $7.6\times10^{-09}$ | 1.101 |
| **DERIVED** 840 | 4800 | hopf | 227.698 | 216.846 | 0.999059 | 0.999704 | 0.00288 | 1.028 |
| **DERIVED** 840 | 4800 | perturbed_hopf | 230.735 | 216.847 | 0.999085 | 0.999703 | 0.003 | 1.028 |
| **DERIVED** 840 | 4800 | no_quartic_control | 64.641 | $6.3776\times10^{-13}$ | 0.999059 | 0.000000 | $8.31\times10^{-09}$ | 1.06 |
| **DERIVED** 6480 | 38400 | hopf | 256.93 | 236.955 | 0.999947 | 0.999974 | 0.00258 | 1.293 |

**DERIVED** Energy concentration is $V\,\int (e^{2})/(\int e)^{2}$: 1 means uniform energy density; larger values mean concentration. It is unreliable once energy is essentially zero.

## Interpretation and limits

- **DERIVED** The base 120-node mesh loses the winding and relaxes to a constant field, even with the quartic term. The base lattice is inadequate for this seed and discretization.
- **DERIVED** Compare refined runs above before judging persistence. Charge near one at the final iteration is finite-run evidence only; a nonzero gradient means relaxation is unfinished.
- **DERIVED** A nonuniform map on a compact sphere is not automatically an isolated localized particle. Near-uniform energy density would mean that the structure has spread across the sphere.
- **DERIVED** Keeping charge on two meshes is not proof of continuum convergence, a local energy minimum, or nonlinear stability. A Hessian analysis, longer runs, more perturbations, improved topology preservation, and further refinement are outstanding.
- **DERIVED** Perturbations use one deterministic random seed (20261025), so this is not a statistical robustness study.
- **DERIVED** Uniform-field energy and analytic energy-gradient checks passed on every mesh; the largest finite-difference discrepancy is listed below. Saved energy traces were checked for monotonic decrease.
- **DERIVED** No emergent gravity, two-body force, inverse-square law, quantum behavior, physical mass, or experimental agreement was tested.

## Mesh and derivative verification

| Nodes | Chordal volume | Derivative error |
|---:|---:|---:|
| **DERIVED** 120 | 16.69252677 | $1.22\times10^{-08}$ |
| **DERIVED** 840 | 18.85946501 | $3.96\times10^{-09}$ |
| **DERIVED** 6480 | 19.51116208 | $4.36\times10^{-08}$ |

**DERIVED** The unit S3 continuum volume is $2\,\pi ^{2} = 19.73920880$. Chordal volume should approach it with refinement.

## Reproduce

**DERIVED** Requires Python, NumPy, SciPy, Matplotlib, and threadpoolctl. From the physics directory:

```powershell
python scripts/hopf_stability_experiment.py
python scripts/summarize_hopf_experiment.py
```

**DERIVED** The simulation overwrites only its own experiments/hopf_first_run outputs. NPZ files preserve initial/final fields and mesh; JSON preserves energy traces and topology diagnostics.
