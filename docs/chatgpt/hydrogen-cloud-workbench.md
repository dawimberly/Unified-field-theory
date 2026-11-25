# Hydrogen benchmark and a concrete cloud-action audit

**Conversation date range** 2026-10-25 to 2026-10-30 (archive window)

**ASSUMED** This checkpoint records hydrogen benchmark and a concrete cloud-action audit. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Prepared for Dan Wimberly, 2026-10-25. This report develops and tests one explicit candidate based on the supplied amplitude, phase and topology ideas. It records failed requirements as well as successful calculations. It is not a unified theory or an independent prediction of particle properties.

## Main findings

**CITED** The standard reduced-mass, nonrelativistic Coulomb calculation gives hydrogen binding 13.5982872641 eV. The NIST ionization reference is 13.598434599702 eV. The difference, 0.0001473356 eV, is not claimed as a numerical failure: relativistic, radiative, recoil and finite-size refinements are outside the leading Schrodinger model. — eq./page UNVERIFIED.

**DERIVED** For the selected cloud action, relaxing the amplitude of a fixed Hopf texture and its radius yields a restricted-family minimum at $R=1.1004012$, energy 251.23527833 in the explicitly dimensionless $v=\kappa =\lambda =1$ example. This is below the frozen-amplitude family's minimum 315.82734083. It demonstrates amplitude-energy reorganization in that chosen action, not a particle mass or atomic binding energy.

**DERIVED** The leading amplitude-mediated interaction is scalar. A localized-source approximation gives a screened attractive potential, not a derived electromagnetic interaction. Fitting its ground energy does not generally reproduce the hydrogen excited spectrum. Its charge-sign behavior is also different from electrostatics.

**DERIVED** Most decisively, the specified vacuum contains amplitude, phase and internal-orientation scalar waves. The texture two-form is not an independent photon field, and no dynamical metric is present. This candidate therefore does not currently meet the electromagnetic or gravitational requirements for unification.

## Completion of the requested steps

| Requested work | Completed calculation | What remains unestablished |
|---|---|---|
| **DERIVED** Define cloud amplitude and phase | Explicit fields, symmetries, vacuum, topology and one action in CANDIDATE_ACTION.md | These definitions are hypotheses, not measured microscopic structure |
| **DERIVED** Derive interaction versus separation | Integrate out the linear amplitude field; verify finite-source cross energy | The full nonlinear, quantum two-Hopf interaction and electron/proton assignments are not derived |
| **DERIVED** Check for stable binding | Solve the effective radial quantum Hamiltonian; negative levels lie below its zero-at-infinity threshold | A bound level in this effective potential is not a demonstrated physical atom or stable full-field composite |
| **DERIVED** Compare binding and spectrum to hydrogen | Converged Coulomb reference and screened-potential tests, including an explicit one-level fit | No independently predicted alpha, particle mass, or hydrogen spectrum from the cloud |
| **DERIVED** Keep gravity in the same framework | Derive the matter stress-energy and audit the vacuum degrees of freedom | No emergent Einstein equation or spin-2 gravitational mode |

## 1. Standard hydrogen reference

**CITED** Inputs are the measured CODATA 2022 electron and proton masses and fine-structure constant. Quantum mechanics, the Coulomb law, and $\hbar$ are imported. The reduced mass is $\mu =me\,mp/(me+mp)$. Define $a_\mu =\hbar /(\mu *c\,\alpha ), E_\mu =\mu *c^{2}\,\alpha ^{2}$. Then — eq./page UNVERIFIED.

**DERIVED** $H/E_\mu = -\nabla^{2}/2 - 1/r$, with r measured in a_mu.

**DERIVED** Numerically $a_\mu =5.29465409463\times10^{-11} m$ and $E_\mu =27.196574528273 eV$. Analytically $E_n=-E_\mu /(2n^{2})$. The radial solver uses $u=rR, u(0)=u(120)=0$, a centered finite-difference kinetic term, and radial spacings 0.02,0.01,0.005. Richardson extrapolation removes the leading second-order spacing error. All quoted energies use one separated-particle threshold $E=0$.

| Level | Extrapolated numerical energy / E_mu | Analytic energy / E_mu | Difference / E_mu |
|---|---:|---:|---:|
| **DERIVED** 1s | -0.49999999984 | -0.50000000000 | $1.63\times10^{-10}$ |
| **DERIVED** 2s | -0.12500000000 | -0.12500000000 | $1.93\times10^{-12}$ |
| **DERIVED** 3s | -0.05555555556 | -0.05555555556 | $1.17\times10^{-14}$ |
| **DERIVED** 2p | -0.12499999999 | -0.12500000000 | $8.99\times10^{-12}$ |
| **DERIVED** 3p | -0.05555555555 | -0.05555555556 | $1.37\times10^{-12}$ |
| **DERIVED** 3d | -0.05555555555 | -0.05555555556 | $4.46\times10^{-12}$ |

**DERIVED** The numerical 1s radial wavefunction agrees with $u(r)=2\,r\,\exp\left(-r\right)$ with squared-overlap fidelity 0.999999999987. Its mean radius is 1.500006250 a_mu and its virial residual $2T+V$ is $-6.25\times10^{-06} E_\mu$.

**DERIVED** The leading nonrelativistic 2p-to-1s wavelength is 121.56844562 nm. Actual Lyman-alpha lines have resolved fine structure; agreement at this approximate level is a reference check, not precision spectroscopy or a new prediction.

## 2. One candidate action, not a collection of independent potentials

**DERIVED** The complete definitions, field equations, static energy, stress-energy tensor and limitations are in CANDIDATE_ACTION.md. The fields are amplitude rho, periodic phase theta, and an independent internal unit vector n with target S2. The background is locally Lorentzian; that spacetime structure is assumed, not emergent. Finite-energy n configurations compactify spatial R3 to S3 for Hopf classification; this is not an assertion of a globally spherical physical universe.

**DERIVED** In natural units, the matter Lagrangian is

**DERIVED** $L = -(\partial \rho )^{2}/2 - \rho ^{2}[(\partial \theta )^{2}+(\partial n)^{2}]/2 - \lambda (\rho ^{2}-v^{2})^{2}/4 - \kappa F_{munu} F^munu/4$,
**DERIVED** $F_{munu} = n \cdot (\partial _\mu n \times \partial _\nu n)$, $n \cdot n=1$.

**DERIVED** rho is a field amplitude, not automatically a normalized probability density. Quantum probabilities require a state of the fields. Standard field quantization is an assumption; a nonperturbative quantum soliton spectrum has not been calculated. The free inputs v,lambda,kappa are not fixed by the notes.

## 3. Actual amplitude relaxation in a topological ansatz

**DERIVED** For the explicit degree-one stereographic Hopf profile n_R and constant theta, frozen $\rho =v$ gives $E(R)=16\,\pi ^{2}\,(v^{2}\,R+\kappa /R)$. Its finite-radius minimum within that family is $R=\sqrt{\kappa }/v$. Numerical derivatives independently checked the local quadratic/quartic invariants of n_R, and quadrature checked their exact integrals.

**ASSUMED** Next $\rho =v\,f(r/R)$ was solved as a radial boundary-value problem with $f\prime(0)=0$ and $f(\infty )=1$, keeping n_R fixed. The energy includes amplitude gradients, amplitude potential, texture gradients and the four-derivative texture energy. The tail of the quadratic texture energy is integrated analytically with $f=1$ outside the numerical domain. Outer-domain extension checks the resulting approximation.

| Radius R | Frozen-amplitude energy | Relaxed-amplitude energy | Central $\rho /v$ |
|---:|---:|---:|---:|
| **DERIVED** 0.25 | 671.13309927 | 646.95051700 | 0.00393833 |
| **DERIVED** 0.5 | 394.78417604 | 353.54615629 | 0.00641511 |
| **DERIVED** 1.0 | 315.82734083 | 252.74552259 | 0.01373586 |
| **DERIVED** 2.0 | 394.78417604 | 315.94363518 | 0.04528394 |
| **DERIVED** 4.0 | 671.13309927 | 603.70988711 | 0.25232151 |
| **DERIVED** 8.0 | 1283.04857214 | 1244.62710353 | 0.74775374 |

**DERIVED** Optimizing the radius gives $R=1.10040118, E=251.235278334$, central $\rho /v=0.015714978$. Increasing the outer radial domain from 70 to 100 changes E by $9.455\times10^{-07}$. Nearby radii 5% above and below have higher energy within this restricted family.

**DERIVED** This is a constrained variational result. The n texture is not fully relaxed; a Hessian or general stability analysis was not performed. The finite-radius optimum and lower amplitude energy do not establish an electron, proton, spin, charge, or binding between two textures. Converting its dimensionless energy to a measured mass would require an independently fixed physical scale.

## 4. Interaction energy derived from the same action

**DERIVED** Expanding $\rho =v+\sigma$ gives amplitude mass $m_\sigma =\sqrt{2\,\lambda }\,v$. A weak, prescribed texture acts as a source J approximately equal to v*(|grad n|$^{2}+$|grad theta|$^{2}$). Integrating out the static quadratic amplitude field gives

**DERIVED** $V_\sigma = -\int [J1(x) \exp\left(-m_\sigma *|x-y|\right)/(4\,\pi *|x-y|) J2(y) d3x d3y]$.

**DERIVED** For point-source idealizations this becomes $V(d)=-q1\,q2\,\exp\left(-m_\sigma *d\right)/(4\,\pi *d)$. Numerical Fourier integration for Gaussian sources agrees with their analytic finite-width expression. This verifies the sign and separation dependence of the stated quadratic approximation.

**DERIVED** The full nonlinear texture interaction also contains texture deformation, phase/orientation fields and overlap effects. The Hopf ansatz has algebraic tails, so it is not an exactly compact source and its full asymptotic response must not be asserted to be a pure Yukawa field. The point-source spectral tests are therefore a diagnostic of an amplitude-mediated attraction channel, not the completed two-soliton calculation.

**DERIVED** Same-sign scalar sources attract, and opposite-sign scalar sources repel. If one simply labels their signs as proton/electron electric charges, this reverses electrostatic behavior. Assigning both sources the same sign can produce attraction, but then does not explain why like electric charges repel. This problem survives taking the scalar mediator mass to zero.

## 5. Stable effective bound states and the hydrogen spectral test

**DERIVED** To test the shape of the amplitude-channel potential, use the standard relative-motion Schrodinger equation with measured electron/proton reduced mass and potential

**DERIVED** $V/E_\mu = -C \exp\left(-s\,r\right)/r$,

**ASSUMED** where C is its strength divided by the measured Coulomb strength and $s=m_\sigma *a_\mu$ in natural units. These units import alpha; they do not derive it. First $C=1$ was checked. Then for each s, ONLY C was fitted so the numerical 1s energy equals -0.5 E_mu, allowing the excited spectrum to test the potential. This fit is deliberately labeled calibration. It targets the nonrelativistic reference, not the more precise NIST energy; omitted QED physics must not be absorbed into a supposedly new force. — eq./page UNVERIFIED.

| Screening s (rows include FITTED calibration) | 1s with $C=1 / E_\mu$ | Fitted C for $1s=-0.5$ | Calibrated 2s / E_mu | Calibrated 2p / E_mu | $2s-2p$ separation (eV) |
|---:|---:|---:|---:|---:|---:|
| **FITTED** 0.0 | -0.49998750 | 1.00001250 | -0.12500234 | -0.12500339 | $2.8331095\times10^{-05}$ |
| **FITTED** 0.01 | -0.49006201 | 1.00998849 | -0.11770225 | -0.11765525 | -0.0012781796 |
| **FITTED** 0.05 | -0.45180396 | 1.04944595 | -0.09199424 | -0.09095538 | -0.028253453 |
| **FITTED** 0.1 | -0.40704564 | 1.09792247 | -0.06616073 | -0.06264427 | -0.095635789 |
| **FITTED** 0.2 | -0.32679643 | 1.19275932 | -0.02998532 | -0.02007637 | -0.26948945 |
| **FITTED** 0.5 | -0.14810680 | 1.46712459 | 0.00028299 | 0.00070025 | -0.011348055 |
| **FITTED** 1.0 | -0.01028204 | 1.90844776 | 0.00034812 | 0.00070101 | -0.0095974274 |

**DERIVED** The nonrelativistic Coulomb reference has both 2s and 2p at -0.125 E_mu. Screening changes their values and splits that degeneracy even after calibrating 1s. The tiny $s=0$ calibration change in C compensates finite radial discretization, not a physical coupling correction. Near-zero positive box eigenvalues are continuum-in-a-box levels, not bound states. Negative levels close to threshold need larger-box localization tests before interpreting them as physical shallow states; no phase-count claim is based on them.

**DERIVED** A massless scalar limit can reproduce the mathematical Coulomb radial spectrum if its strength is supplied, but that does not fix its charge-sign structure, make it a photon, determine alpha, or supply gravity. The numerical screened series is a sensitivity diagnostic; it is not a precision experimental bound on mediator mass.

## 6. Binding-energy mass decrease and energy accounting

**CITED** Using the NIST reference binding energy and the standard conversion gives $\Delta M=2.42414115546\times10^{-35} kg$ per ground-state hydrogen atom, a fraction $1.44851721\times10^{-08}$ of the separated proton-plus-electron mass. This is a standard mass-energy calculation, not a cloud-model prediction. — eq./page UNVERIFIED.

**DERIVED** Formation must export the released energy or leave it in other degrees of freedom. Photon energy, recoil and any other channels share the energy balance. If all outgoing energy stays inside a larger enclosing system, it remains in that system's total energy. No graviton branch, new gravitational charge, universal discrete gravity step, or distinction called gravity $A/B$ is calculated here.

## 7. Electromagnetic and gravitational audit

**DERIVED** Around $\rho=v$, $\theta=\mathrm{constant}$, $n=\mathrm{constant}$, the perturbative waves are one massive amplitude scalar and three massless scalars. n has an internal target index; it is not a spacetime vector photon. The composite F begins quadratically in small orientation fluctuations, and its squared term begins at fourth order. Numerical amplitude scaling confirms the power 4. Thus the action has no independent quadratic Maxwell vector wave around this vacuum.

**DERIVED** There is also no dynamical metric or perturbative massless spin-2 gravitational wave. Coupling a matter action to a background metric defines stress-energy, but does not derive Einstein's equations. Adding an Einstein-Hilbert term would import the missing gravitational dynamics. A different phase or more sophisticated model could have additional collective modes; the present audit does not rule out every emergent-gravity idea.

**DERIVED** A more direct gravity check follows from the phase equation. Its time and spatial kinetic coefficients are both multiplied by $\rho ^{2}$. In $3+1$ dimensions its effective metric is $g_{eff}=(\rho /v)^{2}\,\eta$, conformal to the assumed flat background. Numerical reconstruction verifies the densitized inverse metric. Geometric-optics null rays are unchanged, so simply calling phase waves light would not produce Einstein light bending. This statement concerns the principal wave characteristics; finite-wavelength scattering from amplitude gradients can still occur.

**DERIVED** Adding a gauge field by naively gauging theta has its own problem: nonzero $\rho =v$ would give that gauge field a mass, so long-range electromagnetism would require an unbroken gauge sector. Merely renaming a phase wave or texture two-form is insufficient. An extended gauged-Hopfion model exists in the literature, but its independent gauge field is an additional structure, not something this action already derives.

### Rotating-phase background: a useful change, but still different cones

**DERIVED** The original notes include intrinsic rotation, so the same action was also checked with $\theta =\Omega *t$ and $\rho 0^{2}=v^{2}+\Omega ^{2}/\lambda$. This homogeneous charge-carrying state is not the zero-charge vacuum and is not automatically literal spatial rotation. Coupled amplitude and phase fluctuations give a gapless branch with $c_{phase}^{2}=(\lambda *v^{2}+\Omega ^{2})/(\lambda *v^{2}+3\,\Omega ^{2})$. Integrating out the gapped amplitude at low momentum increases the phase inertia coefficient. This is a propagation effect derived within the candidate, not imposed as a separate wave-speed formula.

| Phase rate $\Omega (v=\lambda =1)$ | Low-momentum phase-wave speed | Internal-orientation wave speed |
|---:|---:|---:|
| **DERIVED** 0.0 | 1.000000000 | 1.000000000 |
| **DERIVED** 0.2 | 0.963624112 | 1.000000000 |
| **DERIVED** 0.5 | 0.845154255 | 1.000000000 |
| **DERIVED** 1.0 | 0.707106781 | 1.000000000 |
| **DERIVED** 2.0 | 0.620173673 | 1.000000000 |

**DERIVED** For example $\Omega =1$ gives phase speed $1/\sqrt{2}$, while orientation waves remain at 1. The rotating background therefore makes analogue propagation more interesting, but it does not establish a single geometry shared by all modes. The dispersion equation and its small-momentum limit were numerically checked. No spatially inhomogeneous rotating solution or gravitational field equation was obtained.

## 8. Fine-structure constant and incidence angle

**DERIVED** No measured alpha was used in the cloud action or the dimensionless texture calculation. It WAS used in the hydrogen reference and to define the diagnostic atomic units. No stationary interaction angle, electromagnetic current normalization or value of alpha is derived. A phase difference, spatial incidence angle, Hopf integer and electromagnetic coupling are distinct quantities until a calculation connects them.

## Verdict and next research decision

**DERIVED** **Retain:** the explicit energy functional, its amplitude-relaxation calculation, its stress-energy accounting, and the hydrogen numerical reference.

**DERIVED** **Reject as a complete hydrogen/unification model in its current form:** Candidate 0 lacks the required derived particle assignments and perturbative gauge/gravitational degrees of freedom. Its isolated amplitude channel is not electromagnetism. Matching one energy by fitting its strength is not progress toward predicting constants.

**DERIVED** **Do not spend the next round tuning lambda or kappa to atomic constants.** The next mathematical requirement is a justified microscopic quantum action or phase that yields an unbroken electromagnetic gauge sector and a universal gravitational geometry/dynamics. The present action does not supply that bridge. Deriving those structures is not completed by writing down the known Maxwell and Einstein terms; such additions should be treated as a reference effective theory, not the claimed emergence.

**ASSUMED** The failure of this particular candidate does not falsify every cloud or topology framework. It identifies why this explicit minimal version cannot yet carry the proposed unification. A full quantum theory has not been solved, and the final microscopic mechanism remains open.

## Reproduce and inspect

**DERIVED** 122 numerical assertions passed. They check implementation and approximations, not experimental confirmation. Checks include Coulomb grid extrapolation, wavefunction/virial checks, calibrated screened-state refinement, scalar Green-function integration, Hopf derivative invariants, amplitude boundary-value convergence and restricted radius stability.

**ASSUMED** Run `python scripts/hydrogen_cloud_workbench.py`. Dependencies: NumPy, SciPy, Matplotlib, threadpoolctl. Both Python scripts must be in the same scripts directory. CANDIDATE_ACTION.md is an authored companion rather than generated by the script. JSON contains all parameters and numerical checks; CSV files preserve radial states and amplitude profiles. Only this workbench's own outputs are overwritten.

## Primary sources and provenance

- **CITED** NIST/CODATA 2022 fundamental constants: https://physics.nist.gov/cuu/Constants/Table/allascii.txt . Measured inputs, not outputs. — eq./page UNVERIFIED.
- **CITED** NIST hydrogen ionization reference: https://physics.nist.gov/cgi-bin/Elements/elInfo.pl?element=1 . — eq./page UNVERIFIED.
- **CITED** NIST hydrogen spectral lines: https://physics.nist.gov/PhysRefData/Handbook/Tables/hydrogentable3.htm . Leading reference only; resolved fine structure is outside the solver. — eq./page UNVERIFIED.
- **CITED** Faddeev and Niemi, Knots and Particles: https://arxiv.org/abs/hep-th/9610193 . Motivation for texture models, not replication. — eq./page UNVERIFIED.
- **CITED** Samoilenka and Shnir, Magnetic Hopfions in the Faddeev-Skyrme-Maxwell model: https://doi.org/10.1103/PhysRevD.97.125014 . Explicitly gauged model for comparison. — eq./page UNVERIFIED.
- User theory notes: docs/POSSIBILITY_CLOUD_THEORY.md and archive/MASTER_PHYSICS_DOCUMENT.md. Candidate 0 is a new proposed formalization, not a claim that all its terms were already derived in those notes. — eq./page UNVERIFIED.
