# Matter organization and the Newtonian pull

Date: 2026-11-19. Baseline: `fbfa641` on `main`. This is a calculation in the existing Candidate-0 / phase-lock setting. It does not modify the action, select a measure, or claim a new force.

## Plain result

Stable arrangements of the same constituents differ in the Newtonian pull only through ordinary stress-energy. On the lattices below, every pair that is stable under the repository Hamiltonian is either distinguished by its energy and multipoles, or is a winding/anti-winding pair with the same energy density, the same multipoles, and the same even pair correlator.

The written action does not define a further organization-dependent coupling. One toy coupling, with a single free coefficient α, is written below and labeled as an assumption. It can assign different active masses to two same-stress phase patterns, but the pattern that would show the split is not a stable critical point. Once the full active mass is matched, that toy leaves no residual pull. The Qin et al. (2026) force benchmark limits the assumed coefficient only in the reading where α changes the active mass and is left out of the matched inertia.

No Earth or Moon interior model is used.

## What the repository already fixes

Three existing lines were read before any coupling was added.

Candidate 0 is a classical field theory for an amplitude ρ, a unit orientation n, and a phase θ. Its static energy on flat three-space, in the normalization fixed by the November 1 notes, is

```
E = ∫ d³x [
    (1/2) |∇ρ|²
  + (1/2) ρ² |∇n|²
  + (1/2) ρ² |∇θ|²
  + (κ/2) Σ_{i<j} (n · (∂_i n × ∂_j n))²
  + (λ/4) (ρ² − v²)²
]
```

The benchmark values used throughout this repository are v = 1, λ = 64, κ = 1. They are dimensionless. The November 1 stability note gives the vacuum spectrum at ρ = v: the amplitude branch has ω² = k² + 128, and the two orientation branches and the phase branch have ω² = k². The zeros at k = 0 are the global rotation and the global phase. The same notes, and the October 30 cross-reference, state that this action is written on an imposed flat metric. It has no dynamical metric. A coherence score is not an energy, and a mass defect is not a change in G.

The S³ language in the project denotes more than one object. The internal target of the Hopf field is an S³, projected to the orientation n on S². A round spatial S³ appears only as a compact domain in a published minimizer used as a mathematical check. Four-dimensional spacetime is a third object. This calculation uses ordinary flat space. The microscopic variables are θ and n on lattice sites. A global shift of θ, and a global rotation of every n, leave the energy and the pair correlations unchanged. The observables used here are therefore

```
C_{ij} = cos(θ_i − θ_j),    S_{ij} = n_i · n_j.
```

The Possibility Cloud note says, as a qualitative hypothesis, that a tighter resonant lock amplifies gravity beyond the mass defect. That sentence is not a term in the action. The October 30 cross-reference already records the obstruction: released binding energy lowers the mass only if it leaves the system, and the defect itself does not change G.

The composite-gravity line in this repository is a regulated auxiliary-metric determinant for free scalars. Its stress tensor is the metric variation of that matter action, including the contact terms catalogued in the October 30 stress and source reports. Those contacts are an ambiguity in the composite operator, not a phase-lock amplifier. Stochastic composite gravity is recorded as a comparison branch and is not adopted here. Nothing in those calculations supplies a Newton constant that depends on a lock correlator once the stress-energy is fixed.

So the ordinary gravitational theory for this matter model, if gravity is Einstein gravity, is the Newtonian field of the stress-energy belonging to E. The sections below compute that field for matched arrangements, then introduce one toy because the action does not contain anything else.

## Lattice transcription

At ρ = v the potential is minimized. On a cubic lattice of spacing a the two-derivative terms become

```
E = J h,
J = v² a,
h_θ = Σ_{bonds} (1/2) wrap(θ_i − θ_j)²,
h_n = Σ_{bonds} (1 − n_i · n_j).
```

`wrap` returns the principal value in (−π, π]. This is the standard nearest-neighbor transcription of ∫ (1/2) ρ² (|∇θ|² + |∇n|²) d³x. The overall factor J is not fixed by the dimensionless continuum parameters, because a is not fixed. Rest energy ε₀ per site is likewise an external normalization. Every fractional monopole below is proportional to J/ε₀. The worked example sets J/ε₀ = 1 so the bookkeeping is visible. That choice is not a material binding fraction.

The Skyrme density built from n · (Δn × Δn) is identically zero when every spin is collinear. It does not change the energy difference between the collinear arrangements. The two-derivative term is what decides their stability.

Two lattices are used, because they are the smallest ones on which the two questions separate.

- A 4-site ring is the smallest lattice with a smooth winding sector and with two phase patterns that share every nearest-neighbor |Δθ|.
- The 8 vertices and 12 edges of a cube are the smallest three-dimensional lattice. The graph is the same as the periodic 2×2×2 lattice. Sites are embedded as a geometric cube in R³ so that bond energy sits on ordinary edge midpoints.

## Stability

Differentiate U(ρ) = λ(ρ² − v²)²/4. Then U''(ρ) = λ(3ρ² − v²). At the repository values,

```
U''(v) = 128,    U''(0) = −64.
```

The arrangements are compared at ρ = v, where the amplitude curvature is positive and the continuum gap is ω² = k² + 128. The negative curvature at ρ = 0 is the existing low-amplitude diagnostic. It is not a property of these arrangements.

On the ring, h_θ = (1/2) θᵀ L θ with L the cycle Laplacian. The Hessian is L. Its eigenvalues are 0, 2, 2, 4. The zero is the global phase. There is no negative eigenvalue, so the quadratic energy is locally stable. The smooth critical points are the uniform windings. For four sites the steps are 0 and ±π/2.

| Arrangement | h_θ | C₂ | winding | max \|force\| | stationary | stable |
|---|---:|---:|---:|---:|---|---|
| locked, all θ = 0 | 0 | +1 | 0 | 0 | yes | yes |
| circulating, step +π/2 | π²/2 = 4.93480220054 | −1 | +1 | 0 | yes | yes |
| anticirculating, step −π/2 | π²/2 | −1 | −1 | 0 | yes | yes |
| folded, steps +π/2, −π/2, −π/2, +π/2 | π²/2 | 0 | 0 | π | no | no |

C₂ is the mean of cos(θ_i − θ_{i+2}). A global phase shift leaves h_θ and C₂ unchanged, which the script checks. Opposite windings match in h_θ and in C₂. The folded pattern matches the circulating nearest-neighbor energy bond by bond and differs in C₂. It is not a critical point: the largest |∂h/∂θ| equals π.

On the cube the orientation energy of the aligned state is 0. The transverse Hessian in each of the two tangent directions is the hypercube Laplacian, with eigenvalues 0, 2, 2, 2, 4, 4, 4, 6. A finite-difference Hessian agrees to 3.3×10⁻⁸. The zero in each transverse direction is a global rotation. The phase Hessian on the same graph is that same Laplacian and contributes one further zero. The kernel has dimension 3, which is the lattice form of the continuum count (two rotations and one phase). With ρ = v the kinetic term of θ and of transverse n is positive, so the aligned cube is linearly stable.

The Néel arrangement puts n = ±ẑ on the two sublattices. Every bond has n_i · n_j = −1, so h_n = 24. The torque on each spin vanishes because every neighbor is antiparallel, so the point is stationary. A common lab-frame tilt of size ε changes the energy by −24 ε² to leading order (numerical curvature −24.0000005647). The direction is unstable. A single flipped spin has h_n = 6 and is not a second vacuum.

Global rotations leave both the aligned energy and the Néel energy unchanged to numerical zero in the script.

## Ordinary baseline

Constituent rest energy ε₀ sits on the sites. Gradient energy sits on the bond midpoints. With c kept explicit,

```
M = (N ε₀ + E) / c²,    E = J h.
```

In the worked example c = 1, a = 1, ε₀ = 1, and J/ε₀ = 1. Masses are then pure numbers. The Newtonian potential is −Σ m/|x − r|, checked against the multipole formula

```
Φ = −M/R − (n · Q · n) / (2 R³)
```

about the center of mass, through the quadrupole. A two-mass control recovers Q_zz = 4 m d² and agrees with the direct potential at R = 100 to a relative error below 10⁻⁷. Direct and quadrupole potentials for the cube arrangements agree at R = 80 to a relative error of about 10⁻⁹.

| Arrangement | h | M | fractional monopole vs aligned | \|Q\| about the CM |
|---|---:|---:|---:|---:|
| aligned | 0 | 8 | 0 | 0 |
| Néel | 24 | 32 | 3 | 0 |
| single flip | 6 | 14 | 0.75 | 1.5746719775 |

The Néel-versus-aligned difference is entirely the monopole

```
ΔM = ΔE / c² = 24 J / c²,
(M_Néel − M_aligned) / M_aligned = 3 (J/ε₀).
```

Cubic symmetry puts the center of mass at the geometric center and sets the quadrupole to zero. At any distance the fractional potential difference equals the fractional mass difference. With a realistic J/ε₀ the same identity holds and the fraction is smaller by that factor. The example value 3 is the choice J = ε₀, not a predicted mass excess.

The single flip has the same constituents and a different distribution. Its center of mass in the geometric frame is (1/7, 1/7, 1/7), so the dipole about the cube center is nonzero. The quadrupole eigenvalues about the center of mass are −0.642857142857, −0.642857142857, and 1.285714285714. Spreading the same extra mass 6 uniformly over the eight sites matches the total mass and gives zero quadrupole.

Along the body diagonal, the fractional potential difference between the localized flip and that smeared, mass-matched distribution is

| R/a | common origin (dipole still present) | centers of mass matched |
|---:|---:|---:|
| 10 | 0.0258266875908 | 0.000422061256407 |
| 40 | 0.00625314958209 | 2.82218584773×10⁻⁵ |
| 80 | 0.00310972719659 | 7.11701331824×10⁻⁶ |

The ratio of the R = 40 value to the R = 10 value is 0.242 in the common-origin column, the 1/R dipole scaling. After the centers of mass coincide the same ratio is 0.0669, against 1/16 = 0.0625 for a pure quadrupole; the remainder is the next multipole. Matching the monopole removes the 1/R piece. Matching the center of mass removes the dipole. Matching the quadrupole leaves a difference that falls two powers of distance faster than the monopole. For the Néel pair there is no multipole left to match: the distributions agree once ΔM = ΔE/c² is included.

Opposite windings on the ring have the same |Δθ| on the same bonds, so they have the same bond masses and the same multipoles. Their Newtonian mass channel agrees. They differ by the sign of the phase current. That sign is ordinary momentum in T_{0i} and would source the ordinary gravitomagnetic field. This comparison is about the Newtonian mass channel. The dimensionless model does not fix an internal speed, so no gravitomagnetic force is quoted.

## Why the action does not add a further piece

Take two stable states. If their energies differ, the far-field monopole differs by ΔE/c², and any mismatch in where that energy sits appears in the multipoles just computed. If their energies and bond-wise |Δθ| agree, as for opposite windings, the even correlator C₂ agrees as well: on the ring a stable state is a uniform step Δ, the energy fixes Δ², and C₂ = cos(2Δ). A function of the nearest-neighbor energy already determines C₂ on the stable set.

The folded pattern is the counterexample that makes the correlator independent of the bond energy: every bond still has |Δθ| = π/2, so every local function of a single bond agrees, while C₂ is −1 for the circulation and 0 for the fold. The identity that sees the reversal is sin(Δ₁) sin(Δ₂), which is +1 along the circulation and −1 on the reversing steps. That pattern is not stationary. Holding it would take an external constraint. The constraint does work, and that work belongs in the energy budget of whatever apparatus holds the fold.

In the long-wave limit the independence disappears. For a slow step δ,

```
C₂ ≈ cos(2δ) = 1 − 2 δ² + O(δ⁴) = 1 − 4 × (δ²/2) + O(δ⁴),
```

and δ²/2 is the bond energy. The script’s check at δ = 0.01 gives C₂ + 4×(bond energy) − 1 = 6.7×10⁻⁹. A continuum operator built from (∇θ)² therefore repeats the gradient energy. The sign-sensitive split exists at finite spacing, at gradients of order one radian per lattice step, which is the cutoff of this lattice rather than a second continuum source.

Coupling C₂ minimally, as another term in the matter Lagrangian, would add its own contribution to the stress tensor. That contribution is ordinary mass-energy of the new term. It is removed by matching the full mass, including the new term. It is not a residual on top of that match.

## One labeled toy

The missing coupling, if one wants a pull that survives matching of the Candidate-0 energy, is an active-mass term that is not the metric variation of that energy. The single-coefficient toy used here is

```
ρ_active(x) = ρ_ordinary(x) + α C₂(x) / c².
```

α has units of energy per volume. C₂ is dimensionless, so α/c² is a mass density and α/(ρ c²) is a pure number. The ordinary Hamiltonian, including rest energy and J h, is left unchanged. The toy is an assumption. It is not derived from Candidate 0, from the S³ Hopf projection, or from the auxiliary-metric determinant.

Checks on this coupling:

- Double counting against the nearest-neighbor energy: on the circulating/folded pair the nearest-neighbor energies agree and C₂ differs by 1, so α is not a rescaling of J for that pair. In the long-wave expansion above, C₂ repeats the gradient energy. Both statements are true, on different classes of configurations. The stable configurations are in the second class.
- Units close as above.
- Conserved stress tensor: the lattice C₂ is computed from coordinate differences on a fixed flat background. It does not depend on a dynamical metric, so metric variation of −α C₂ does not produce the mass density written in the toy. Promoting the separations to proper distances would make some stress tensor, and that tensor is not calculated here and is not equal to a hand-added α C₂/c². Embedding the toy in one conserved, minimally coupled stress tensor fails. The failure is recorded in `flags` in the results file.
- Equilibrium: a minimally coupled potential large enough to matter would also contribute pressure. For a wafer in mechanical equilibrium the virial identity returns the active mass to the total energy. The unequilibrated combination ρ + 3p/c² is not used as a prediction for these wafers.
- If the full active mass, including α, is what the torsion balance and the laboratory scale both measure, the paired difference vanishes.

The toy therefore predicts an extra pull only when a C₂ difference is maintained at fixed ordinary inertial mass and α is allowed to change the active mass alone. No stable pair computed here supplies that C₂ difference.

## Scaling at fixed density

Hold the lattice spacing and the rest-mass density fixed, and grow the sample. Let β = α/(ρ c²) be fixed. For a source swap against an unchanged test mass, with the same geometric kernel as the Newtonian force,

```
δF/F = β × (mean correlator of the source).
```

A short-range correlator whose mean differs by a constant ΔC between two preparations then gives δF/F = β ΔC at every N. A sample-spanning lock with Φ = |N⁻¹ Σ exp(iθ)|² = 1 gives the same N-independent mean. For independent uniform phases the exact moments are

```
E[Φ] = 1/N,    Var(Φ) = (N − 1)/N³,    rms(Φ) = √(N − 1) / N^(3/2).
```

The mean shift of a random sample-spanning pattern falls as 1/N. The rms falls as 1/N, not as 1/√N. A Monte Carlo of 20,000 draws on the cube (seed 0) gives Φ mean 0.12423 against 1/8 = 0.125, and rms 0.116093 against 0.116927.

The mean of the nearest-neighbor cosine over independent uniform phases is 0. Its rms on the cube is a fluctuation of the spatial average. The Monte Carlo rms is 0.204772. An estimate that treats the 12 bonds as independent, √((1/2)/12) = 0.204124, is close on this small graph; the script reports both and uses the mean, which is zero, as the shift. The 1/√N_bonds factor is the fluctuation under that independence assumption. It is not applied to the mean shift of a prepared arrangement.

The static Hamiltonian does not assign a correlation time. If, as an extra assumption, the correlator were a stationary process with exponential autocorrelation time τ, an average over duration T ≫ τ would multiply the variance by 2τ/T. For the purely illustrative values τ = 1 s and T = 100 h that factor is 5.56×10⁻⁶. It is not applied to any mean in the table.

With β = 1 the fractional mean shifts are

| N | short-range mean, ΔC = 1 | locked, sample-spanning | random sample-spanning mean | random sample-spanning rms |
|---:|---:|---:|---:|---:|
| 8 | 1 | 1 | 0.125 | 0.116927 |
| 64 | 1 | 1 | 0.015625 | 0.0155024 |
| 512 | 1 | 1 | 0.00195312 | 0.00195122 |
| 4096 | 1 | 1 | 0.000244141 | 0.000244111 |

A short-range mean and a locked sample-spanning mean stay constant as the sample grows at fixed density. The random sample-spanning mean shrinks as 1/N. Nothing in the table puts a 1/√N factor on the mean.

## Qin et al. (2026) in experimental units

The benchmark is Qin et al., Nature Communications (2026), DOI 10.1038/s41467-026-78187-4: two silicon wafers, each about 7.6 mg, face-to-face gap 0.654 mm, reported force 471 ± 23 aN (1σ), and a 5.3% uncertainty on the inferred G. The paper’s calculated 16ω_d torque at that gap is 11.17 aN·m and the force arm is 23.694 mm. Their ratio is 471.427 aN, which is the published 471 aN figure. That number is the modulated torque divided by the arm, including 16 source wafers and the spokes. It is not the static attraction of one face-on pair. The spokes are 0.53 aN·m, a fraction 0.04745 of the calculated torque.

Wafer dimensions used here are the paper’s 2.598 × 2.597 × 0.487 mm³. With mass 7.6 mg the density is 2312.99 kg/m³. The static pair force was integrated on cell-centered grids (8×8×4, 12×12×6, 16×16×8). The last interval changes the close-gap force by 0.13%.

| gap | point-mass force | extended force | extended/point | 50 μm skin / extended | skin mass fraction |
|---:|---:|---:|---:|---:|---:|
| 0.654 mm | 2961.16 aN | 1070.7 aN | 0.361583 | 0.12209 | 0.102669 |
| 1.554 mm | 925.437 aN | 561.186 aN | 0.606402 | 0.117804 | 0.102669 |
| 50 mm | 1.51242 aN | 1.5105 aN | 0.99873 | 0.10356 | 0.102669 |

At 50 mm the extended body agrees with the point mass to 0.13%. At the experimental gap the finite size changes the absolute force by a factor 0.36, so geometry matters for the absolute Newtonian force.

It cancels in the fractional toy signal when the extra active density is proportional to the ordinary density throughout the wafer. Multiplying the source mass by 1.02 multiplies the extended force by 1.02 to the convergence of the grid. Then

```
δF/F = α ΔC₂ / (ρ c²)
```

at every gap, and the extended/point factor drops out. The distance ratio of the extended force between 1.554 mm and 0.654 mm is 0.524129, against 0.312525 for point masses. A bulk fractional shift does not change that ratio.

The distance dependence changes if the same coefficient is supported only in a skin. A 50 μm slab on the near face, at the wafer density, holds mass fraction 0.102669. Its share of the force is 0.12209 at 0.654 mm and 0.117804 at 1.554 mm. Replacing the force fraction by the mass fraction would be wrong by about 19% at the close gap. The skin-force grid change at that gap is 0.17%, smaller than the gap-to-gap change in the fraction. This skin is an alternate support of the same α, used only to apply the geometry correction the distance dependence requires. The headline translation uses bulk support, because C₂ of a uniform preparation fills the sample.

Under bulk support, a source swap, and |ΔC₂| = 1, the published 1σ force figure gives

```
|α| / (ρ c²) ≤ 23/471 = 0.0488323,
|α| ≤ 1.01513×10¹⁹ J/m³.
```

If both wafers carry the same ΔC₂, the fractional force is about twice as large and the coefficient limit is 0.0244161. The 5.3% G uncertainty is the absolute inference. A paired swap cancels the common G calibration; the force figure 23 aN is the resolution the pair still has to beat. The spoke fraction is not folded into the limit. If only part of the attractor is reorganized, the torque fraction is smaller than |α|/(ρ c²) by the share of the torque that part produces, and the paper’s budget isolates the spokes (4.7%) without isolating one wafer.

For scale only: a cohesive energy of 4.63 eV per silicon atom over 28.085 u c² is a fractional rest energy 1.76981×10⁻¹⁰. A source-swap force difference of that size on the 471 aN channel would be 8.3×10⁻⁸ aN. The Qin benchmark does not approach chemical binding. The cohesive energy is an external comparison. It is not an input to the lattice and it is not an interior model of a planet.

What the paired comparison must resolve, for bulk support and a source swap, is

```
|δF| = 471 aN × |α ΔC₂| / (ρ c²).
```

Resolving 23 aN limits |α|/(ρ c²) at the 0.049 level when |ΔC₂| = 1. That limit is a limit on the assumed coefficient in the active-mass-only reading. In the reading where the stress tensor is conserved and minimally coupled, matching the full active mass removes the difference, and the same number is not a bound on a new force.

The assumptions that control the estimate are: the value of ΔC₂ that can actually be prepared and held; whether α is excluded from the matched inertial mass; bulk versus skin support; and which part of the attractor carries the reorganization. The first of these is empty for the stable states of the Hamiltonian computed here.

## Worked numerical example

Repository parameters v = 1, λ = 64, κ = 1. Illustrative ratio J/ε₀ = 1, with c = a = ε₀ = 1.

1. Locked ring: h = 0, C₂ = 1, stable. Circulating ring: h = π²/2, C₂ = −1, winding +1, stable. The mass difference is (π²/2) J/c². Anticirculating: same mass, same C₂, winding −1. Folded: same nearest-neighbor h, C₂ = 0, force π, not a solution.
2. Aligned cube: h_n = 0, stable, quadrupole 0. Néel cube: h_n = 24, stationary, curvature −24, quadrupole 0. The entire far-field difference is ΔM = 24 J/c².
3. Single flip: h_n = 6. Against the same extra mass smeared uniformly, the fractional potential along the body diagonal is 2.58×10⁻² at R = 10a before the centers of mass are matched, and 4.22×10⁻⁴ after.
4. Toy, only if a held ΔC₂ is inserted by hand at fixed ordinary mass: δF/F = [α/(ρ c²)] ΔC₂ for bulk support. Qin’s 23 aN on 471 aN limits |α|/(ρ c²) by 0.0488323 per unit |ΔC₂| in a source swap. No stable pair in this calculation has a ΔC₂ left after the ordinary mass is matched.

## Reproduction

From the repository root, with Python 3.10+ and the NumPy/SciPy pins in `requirements.txt`:

```sh
python3 scripts/matter_organization_gravity.py
```

The script writes `experiments/matter_organization_gravity_2026_11_19/results.json` and exits nonzero if an internal check fails. This run used Python 3.12.3, NumPy 2.5.5, and SciPy 1.18.1. G and c are `scipy.constants` values, recorded in the JSON.

### Script output

```
Matter organization and Newtonian pull
python 3.12.3, numpy 2.5.5, scipy 1.18.1
amplitude U''(v) = 128; U''(0) = -64
ring: locked E,C2 = 0, 1; circulating E,C2,w = 4.93480220054, -1, 1; anticirculating E,C2,w = 4.93480220054, -1, -1; folded E,C2,max|F| = 4.93480220054, 0, 3.14159265359
ring Hessian eigenvalues [-8.586881206085195e-16, 2.0, 2.0000000000000004, 3.9999999999999996]; negative 0, zero 1
cube orientation energies aligned, Neel, flip = 0, 24, 6
Neel tilt curvature / eps^2 = -24.0000005647 (expected -24); vacuum Hessian max abs error 3.331e-08
baseline masses aligned, Neel, flip = 8, 32, 14
baseline quadrupoles aligned, Neel, flip = 0.000e+00, 0.000e+00, 1.5746719775
ΔM(Neel-aligned) = 24; fractional monopole 3
flip vs smeared at R=10.0: common-origin 0.0258266875908, CM-matched 0.000422061256407
flip vs smeared at R=40.0: common-origin 0.00625314958209, CM-matched 2.82218584773e-05
flip vs smeared at R=80.0: common-origin 0.00310972719659, CM-matched 7.11701331824e-06
random-phase Monte Carlo: Phi mean 0.12423 (exact 0.125), rms 0.116093 (exact 0.116927); short-range mean -0.0013879, rms 0.204772, iid-bond rms estimate 0.204124
Qin density 2312.99 kg/m^3; torque/lever 471.427 aN; spoke fraction 0.0474485
gap 0.654 mm: point 2961.16 aN, extended 1070.7 aN, extended/point 0.361583, grid change 1.332e-03, skin/extended 0.12209 (mass fraction 0.102669)
gap 1.554 mm: point 925.437 aN, extended 561.186 aN, extended/point 0.606402, grid change 8.244e-04, skin/extended 0.117804 (mass fraction 0.102669)
gap 50 mm: point 1.51242 aN, extended 1.5105 aN, extended/point 0.99873, grid change 3.445e-06, skin/extended 0.10356 (mass fraction 0.102669)
distance ratios 1.554/0.654 mm: extended 0.524129, point 0.312525, skin fraction near 0.12209, far 0.117804
uniform source rescaling 1.02 -> force ratio 1.02
toy bound |alpha|/(rho c^2) at |ΔC2|=1, source swap, 1σ: 0.0488323; |alpha| 1.01513e+19 J/m^3; both-wafers reading 0.0244161
external Si cohesive fraction 1.76981e-10; force resolution for that fraction 8.3358e-08 aN
N=8: short-range mean 1, locked spanning 1, random spanning mean 0.125, random spanning rms 0.116927
N=64: short-range mean 1, locked spanning 1, random spanning mean 0.015625, random spanning rms 0.0155024
N=512: short-range mean 1, locked spanning 1, random spanning mean 0.00195312, random spanning rms 0.00195122
N=4096: short-range mean 1, locked spanning 1, random spanning mean 0.000244141, random spanning rms 0.000244111
assumed fluctuation average: variance factor 5.55556e-06 (not applied to the mean)
On the lattices computed here, arrangements that are stable under the repository Hamiltonian differ gravitationally only through ordinary stress-energy. Opposite windings match in energy, multipoles, and the even pair correlator. The Néel arrangement is stationary but unstable, and its far field differs from the vacuum by ΔM = ΔE/c² with vanishing extra multipoles. The one toy coupling can split a same-stress pair, but that pair is not stable, and the coupling is not part of the action. Matched full active mass leaves no residual.
wrote /workspace/experiments/matter_organization_gravity_2026_11_19/results.json
```
