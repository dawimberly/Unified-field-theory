# Dimension, field concentration, and gravity: mechanism audit
1 November 2026. Research checkpoint for Dan Wimberly.

## Finding
Three spatial dimensions are special in several precise senses, but the fluid and gravity thresholds are not yet one mechanism. The strongest project-specific result is that positive quadratic and quartic gradient energies balance under spatial dilation only for 2<d<4, selecting d=3 among integer dimensions. This permits a finite-size minimum along a dilation family; it does not prove a stable particle, gravity, objective collapse or dimensional selection of the universe.

Completed work: inspected primary fluid/gravity literature and the existing Candidate 0 action, derived the dimension counts and scaling conditions below, and executed independent small algebraic controls. No nonlinear field simulation or formal proof-checker run was performed. The prior Navier-Stokes paper's theorem and introductory mechanism were inspected in the previous lookup; its 166-page proof was not independently verified.

## 1. Distinguish spatial and spacetime dimensions
Write d for spatial dimensions and D=d+1 for spacetime dimensions. An internal three-component orientation n is not three extra spatial dimensions. A compactified spatial R3 is not the same object as the field's internal S2 or spacetime.

## 2. What changes in fluid flow
For incompressible constant-density flow in three spatial dimensions,

    partial_t u+(u.grad)u = -grad p + nu Laplacian u + f,
    div u=0.

Pressure is divided by the constant density. Taking the curl yields

    D_t omega=(omega.grad)u+nu Laplacian omega+curl f.

For strictly planar u=(u_x(x,y),u_y(x,y),0), omega=(0,0,omega_z), so (omega.grad)u=omega_z partial_z u=0. In three-dimensional flow it need not vanish. Under periodic or sufficiently decaying conditions,

    (1/2) d/dt integral |omega|^2
       = integral omega_i S_ij omega_j
         -nu integral |grad omega|^2 + integral omega.curl f,
    S_ij=(partial_i u_j+partial_j u_i)/2.

Stretching can amplify rotation, but its sign depends on alignment with strain. It is not necessarily growth and does not alone prove blowup. Smooth externally forced blowup also requires keeping curl f and all forcing assumptions; the force cannot be discarded.

A local control is u=(-a*x/2-Omega*y, Omega*x-a*y/2, a*z). Its divergence vanishes, omega_z=2 Omega, stretching_z=2a Omega, and omega.S.omega=4a Omega^2. Positive a amplifies the local rotation; negative a attenuates it. These are kinematic evaluations at a specified instant, not a finite-energy Navier-Stokes blowup construction: the affine field has infinite total energy on R3.

The 3D cross product is special notation. In higher dimensions vorticity is an antisymmetric two-form and deformation can still amplify it. Do not infer that nonlinear concentration exists only in d=3. Conversely, singularity formation in some other 2D or 1D equation is not excluded by smoothness of 2D incompressible Navier-Stokes.

## 3. What changes in Einstein gravity
A massless gravitational plane wave in D-dimensional flat spacetime has a symmetric traceless polarization tensor on the D-2 transverse directions. Hence

    N=(D-2)(D-1)/2-1=D(D-3)/2, D>=3.

For D=3,4,5,6 this gives 0,2,5,9 polarizations. Ordinary 2+1 Einstein gravity has global/geometric effects but no local vacuum gravitational waves. In D=3 the Riemann tensor is determined algebraically by Ricci curvature, with no independent Weyl tensor. D=4 admits independent vacuum tidal curvature and radiation. Higher dimensions also admit gravity.

This counting is not a vortex-stretching equation. Linear gravitational waves exist without a nonlinear instability. It also concerns ordinary Einstein gravity; modifying the action can change degrees of freedom in low dimensions.

## 4. An established bridge, with an important dimension mismatch
Bredberg, Keeler, Lysov and Strominger construct a near-horizon/hydrodynamic correspondence relating an incompressible fluid in p spatial plus one time dimension to a vacuum Einstein geometry in p+2 spacetime dimensions. The construction uses boundary, regularity and expansion assumptions; it is not an unrestricted equivalence of all fluid and gravitational solutions.

Thus a planar p=2 fluid corresponds to D=4 gravity, while p=3 corresponds to D=5. The presence of 4D gravitational radiation does not require vortex stretching in its boundary fluid. This defeats the simplest proposed identification of the two minimum-dimension thresholds, though not every possible deeper connection.

The paper explicitly considers forcing as well. Applying a particular forced singular fluid solution would require checking its boundary forcing dictionary and whether the hydrodynamic expansion remains controlled near the singular time. The blowup drives gradients large, exactly where a derivative expansion may fail. No naked spacetime singularity or quantum outcome follows just by substituting the fluid solution.

## 5. A direct calculation in our existing field model
Candidate 0 has positive static terms with two spatial derivatives, four spatial derivatives, and a nonnegative potential. Under x -> x/s, retaining the same profile shape and field amplitudes,

    E_d(s)=A s^(d-2)+B s^(d-4)+C s^d,
    A=E_2>0, B=E_4>0, C=E_0>=0.

The coefficients refer to a fixed admissible finite-energy reference configuration. This compares scaling in each dimension, not the same topological sector across all dimensions: pi_d(S2) changes with d. Coupling dimensions can change too. The exponents and sign test remain valid within a specified action in each d.

Ignoring the potential temporarily, stationary size requires

    s_*^2=(4-d)B/[(d-2)A].

A positive finite solution requires 2<d<4. Only integer d=3 qualifies. At d=3, E(s)=A s+B/s has s_*=sqrt(B/A) and positive second derivative 2B/s^3. The quartic term resists shrinking and the quadratic term resists spreading.

Restoring C gives in d=3

    E'(s)=A-B/s^2+3C s^2,
    E''(s)=2B/s^3+6C s>0.

There is exactly one positive minimum along this scaling family. For C>0,

    s_*^2=2B/[A+sqrt(A^2+12BC)].

With dimensionless A=B=C=1, s_*=0.6589829633. This value is an illustration, not a physical length prediction. Four parameter controls gave derivative residuals below 9e-16 and positive dilation Hessians.

In d=2, E=A+B/s^2+C s^2: a positive potential term can also stabilize size. Therefore d=3 is special for the TWO-gradient competition, not a universal uniqueness theorem. At d=4, E=A s^2+B+C s^4 decreases toward its infimum as s shrinks; for d>4 all terms decrease on shrinking. Neither result proves actual finite-time dynamical collapse. Higher derivatives, charge, boundaries, time dependence or other interactions change the test.

This generalizes the project's earlier 3D Derrick audit. It is an established scaling method applied to the candidate, not a newly discovered theorem.

## 6. Does the model contain fluid stretching or gravity?
Its smooth phase gradient has zero curl locally: curl grad(theta)=0. Identifying phase gradient with a nonrelativistic velocity therefore does not produce a smooth vortical fluid away from defects. The relativistic conserved phase current rho^2 partial^mu theta is not automatically an incompressible Navier-Stokes velocity; normalization and a physical dictionary matter.

The orientation texture does provide F_ij=n.(partial_i n cross partial_j n). In three dimensions B_i=(1/2)epsilon_ijk F_jk has vanishing divergence for smooth fields, by the closedness of the pulled-back S2 area form. But divergence-free does not imply its ideal-fluid transport law. To obtain vortex stretching one would need a derived velocity v and transport equation partial_t B=curl(v cross B), or justified corrections. Candidate 0's field equations have not been shown to supply that identification. Its Skyrme kinetic matrix also differs from fluid kinetic energy.

The written vacuum still has scalar excitations and no dynamical metric. The new size argument does not create a propagating spin-2 mode or remove the earlier light-bending failure. Since the user's current program allows changing quantum dynamics while retaining GR, an alternative is to couple the matter model consistently to supplied Einstein dynamics. That would import gravity rather than derive its origin, and is a separate stated route.

## 7. Interpretation and next mathematical gate
The evidence supports: dimensionality changes available deformation, radiation and stabilization mechanisms. It does not support gravity existing only in our dimension or a Navier-Stokes singularity causing quantum collapse.

Keep the agreed interpretation-neutral objective: identify dynamics and observables before deciding whether any event is collapse. The immediate gate in Candidate 0 is a stationary or explicitly time-dependent coupled solution with controlled topology, followed by its full constrained stability operator. A dilation minimum alone leaves deformation, amplitude, phase and translation modes unchecked. Treating a prepared knot as a stationary solution would invalidate a normal-mode stability claim.

If investigating the fluid connection specifically, first derive or disprove the texture transport dictionary. Without it, running fluid blowup profiles through the texture code would test a new imported equation instead of our action. There is currently no derived common mechanism connecting the three dimensional facts.

## Sources
- Local model: C:/workspace/physics/experiments/hydrogen_cloud_workbench/CANDIDATE_ACTION.md
- Local prior audit: C:/workspace/physics/experiments/cloud_amplitude_dynamics/STATIC_SECTOR_AUDIT.md
- Local coupled runs: C:/workspace/physics/experiments/knotted_cloud_dynamics/REPORT.md
- Bredberg et al., From Navier-Stokes To Einstein, JHEP07(2012)146: https://arxiv.org/abs/1101.2451 ; PDF inspected, especially the setup and expansion structure.
- Carlip, Quantum Gravity in 2+1 Dimensions, Living Reviews 8:1 (2005): https://arxiv.org/abs/gr-qc/0409039 ; dimensional result checked against the accessible review record and indexed text.
- Vorticity derivation, UCL teaching notes: https://uclnatsci.github.io/Electromagnetism-Fluids-and-Waves/fluids/vorticity.html
- Recent forced Navier-Stokes construction: https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf

Reproduction: python dimension_mechanism_checks_2026_11_01.py. Uses only Python standard library. JSON records dimension counts, eight local strain/rotation controls and four restricted size minima. These checks establish arithmetic consistency, not empirical confirmation.