# Candidate 0 amplitude, phase and Hopf action

**Conversation date range** 2026-10-25 to 2026-11-01

**ASSUMED** Candidate 0 supplies a flat-background field action with amplitude, phase and orientation. **DERIVED** The archive develops stability, symmetry and energy checks in this model. **UNVERIFIED** The relationship to observed particles or gravity remains open.

## Research record

# Candidate 0: amplitude, phase, and Hopf texture

**ASSUMED** Date: 2026-10-25. This is a proposed effective field model, not an established theory of particles or a derivation of quantum gravity. Its purpose is to make the supplied cloud ideas precise enough to test and reject when necessary. No physical parameter is inferred by fitting a numerical coincidence.

## 1. Fields, geometry, and probability

**ASSUMED** Work locally on a Lorentzian background with signature (-,+,+,+). The metric is initially a background, not an emergent or dynamical field. Use natural units $\hbar =c=1$ in this document.

| Symbol | Definition | Status |
|---|---|---|
| **ASSUMED** $\rho (x) \ge 0$ | Real medium amplitude, mass dimension 1 | Proposed field |
| **ASSUMED** theta(x), period 2 pi | Dimensionless collective phase | Proposed field |
| **ASSUMED** $n(x), n \cdot n=1$ | Three-component INTERNAL orientation; target S2 | Proposed order parameter |
| **ASSUMED** g_munu | Lorentzian background metric | Imported spacetime structure |
| **ASSUMED** $v>0$ | Vacuum amplitude | Free dimensional input |
| **ASSUMED** $\lambda >0$ | Amplitude-potential coefficient | Free dimensionless input |
| **ASSUMED** $\kappa >0$ | Four-derivative texture coefficient | Free dimensionless input |

**ASSUMED** rho is not automatically a normalized single-particle probability density. In a field quantization the state is formally a wavefunctional Psi[rho,theta,n], with probabilities assigned to configurations after specifying its Hilbert space and constraints. Quantum mechanics is assumed; it is not obtained from classical random fluctuations.

**ASSUMED** The three fields are independent subject to $n \cdot n=1$. In particular n remains defined if rho vanishes. If instead n were merely a parametrization of a complex field that becomes zero, its topology could cease to be protected; that is a different model and must not be silently substituted.

**ASSUMED** Finite-energy configurations approach $\rho=v$, $\theta=\mathrm{constant}$, $n=n_\infty$ at spatial infinity. Compactifying R3 by its boundary point gives a domain topologically S3 for classifying the n texture. This is NOT an assumption that the physical universe is a finite round S3. The n map can carry an integer Hopf invariant in pi3(S2). The phase alone does not supply that Hopf invariant.

## 2. One candidate matter action

**ASSUMED** Define

**ASSUMED**

$$
F_{\mu\nu}=\mathbf n\cdot(\partial_\mu\mathbf n\times\partial_\nu\mathbf n),\qquad
U(\rho)=\frac{\lambda}{4}(\rho^2-v^2)^2.
$$

**ASSUMED** The action is

**ASSUMED**

$$
S_0=\int d^4x\sqrt{-g}\left[
-\frac12\partial_\mu\rho\,\partial^\mu\rho
-\frac{\rho^2}{2}\left(\partial_\mu\theta\,\partial^\mu\theta+
\partial_\mu\mathbf n\cdot\partial^\mu\mathbf n\right)
-U(\rho)-\frac{\kappa}{4}F_{\mu\nu}F^{\mu\nu}
\right].
$$

**CITED** This is an amplitude/phase extension of the type of texture energy investigated in Faddeev-Skyrme models. It is a chosen candidate, not a result already present in the archive. Its symmetries include global phase shifts, global internal rotations of n, and background spacetime covariance. It has no independent electromagnetic gauge potential and no Einstein-Hilbert term. As an effective field theory it also has no demonstrated ultraviolet completion. — eq./page UNVERIFIED.

**ASSUMED** For static fields in flat space,

**ASSUMED**

$$
E=\int d^3x\left[\frac12|\nabla\rho|^2+
\frac{\rho^2}{2}(|\nabla\theta|^2+|\nabla\mathbf n|^2)
+U(\rho)+\frac{\kappa}{4}F_{ij}F_{ij}\right].
$$

**ASSUMED** All displayed static contributions are nonnegative. Binding would mean E_two-configurations is less than the sum of two isolated energies, not that any absolute contribution must be negative. Boundary conditions, vacuum subtraction and all parameters must remain the same in that comparison.

## 3. Equations and conserved quantities

**ASSUMED** Let $X=\partial (\theta ) \cdot \partial (\theta )+\partial (n) \cdot \partial (n)$. Then

$$
\Box_g\rho-\rho X-\lambda\rho(\rho^2-v^2)=0,
$$

$$
\nabla_\mu(\rho^2\partial^\mu\theta)=0.
$$

**ASSUMED** The phase equation gives a global conserved current $j_\theta ^\mu =\rho ^{2} \partial ^\mu \theta$. It is not by itself Maxwell's equations or an electric charge assignment.

**ASSUMED** The constrained n equation can be expressed as

**ASSUMED**

$$
P_\perp\nabla_\mu\left[\rho^2\partial^\mu\mathbf n+
\kappa F^{\mu\nu}(\partial_\nu\mathbf n\times\mathbf n)\right]=0,
$$

**ASSUMED** where $P_{perp}=I-n n^T$ projects onto variations tangent to S2. This expression follows by varying F under tangent variations of n.

**ASSUMED** The Hilbert stress-energy tensor, with the constraint imposed, is

**ASSUMED**

$$
T_{\mu\nu}=\partial_\mu\rho\partial_\nu\rho+
\rho^2(\partial_\mu\theta\partial_\nu\theta+
\partial_\mu\mathbf n\cdot\partial_\nu\mathbf n)+
\kappa F_{\mu\lambda}F_\nu{}^\lambda+g_{\mu\nu}\mathcal L_0.
$$

**ASSUMED** It includes gradients, internal motion, interaction terms and stress. Counting amplitude or frequency alone would omit part of the energy. Background covariance defines how this matter could be coupled to gravity, but it does not produce a gravitational field equation.

## 4. A concrete texture energy calculation

**ASSUMED** Use a degree-one Hopf ansatz obtained from the unit spinor

**ASSUMED**

$$
z_1=\frac{2R(x+iy)}{r^2+R^2},\qquad
z_2=\frac{2Rz+i(r^2-R^2)}{r^2+R^2},
$$

$$
\mathbf n=(2\operatorname{Re}z_1\bar z_2,2\operatorname{Im}z_1\bar z_2,|z_1|^2-|z_2|^2).
$$

**ASSUMED** Its invariants are

**ASSUMED**

$$
|\nabla\mathbf n|^2=\frac{32R^2}{(r^2+R^2)^2},\qquad
F_{ij}F_{ij}=\frac{512R^4}{(r^2+R^2)^4}.
$$

**ASSUMED** For theta constant and rho frozen to v,

$$
E(R)=16\pi^2(v^2R+\kappa/R).
$$

**ASSUMED** Within this family the finite-radius minimum is $R=\sqrt{\kappa }/v$, with $E=32 \pi ^{2} v \sqrt{\kappa }$. This is an exact ansatz calculation, NOT a proof of stability against arbitrary deformations. Freezing rho is not a solution of its equation in a nontrivial texture; it is a restricted variational comparison.

**ASSUMED** We then allow $\rho (r)=v f(r/R)$, keeping the n profile fixed, and solve

$$
f''+\frac{2}{s}f'-\lambda v^2R^2 f(f^2-1)-\frac{32f}{(1+s^2)^2}=0,
$$

**ASSUMED** with $f'(0)=0$ and $f(\infty )=1$. The numerical calculation relaxes f and the radius R for $v=\lambda =\kappa =1$. It therefore tests an actual amplitude-energy reorganization in this action. It does not yet calculate a two-soliton binding energy or identify a texture with an electron or proton.

## 5. Derived leading amplitude-mediated interaction

**ASSUMED** Expand $\rho =v+\sigma$. The amplitude excitation has squared mass $m_\sigma ^{2}=2 \lambda v^{2}$. At quadratic order in sigma, a prescribed static texture contributes

**ASSUMED**

$$
E_\sigma=\int d^3x\left[\frac12|\nabla\sigma|^2+
\frac12m_\sigma^2\sigma^2+\sigma J(\mathbf x)\right],
\qquad J\simeq v(|\nabla\mathbf n|^2+|\nabla\theta|^2).
$$

**ASSUMED** Keeping only this linear-response channel gives

$$
(-\nabla^2+m_\sigma^2)\sigma=-J.
$$

**ASSUMED** Minimizing the quadratic field energy gives $E_\sigma ,min=-1/2 \int (J G J)$. For two prescribed sources the cross energy is

**ASSUMED**

$$
V_\sigma=-\int d^3x\,d^3y\;J_1(\mathbf x)
\frac{e^{-m_\sigma|\mathbf x-\mathbf y|}}{4\pi|\mathbf x-\mathbf y|}
J_2(\mathbf y).
$$

**ASSUMED** For localized point-source idealizations $J_i=q_i \delta (x-x_i)$,

$$
V_\sigma(d)=-\frac{q_1q_2}{4\pi}\frac{e^{-m_\sigma d}}d.
$$

**ASSUMED** Gaussian-source integration is independently checked numerically. The sources are prescribed for that diagnostic; their normalization is not derived as an electron or proton charge. The scalar amplitude is therefore a possible attraction channel, but it is not electromagnetism. Like-sign scalar sources attract; applying opposite electric-charge signs would give repulsion, the reverse of electrostatics.

**ASSUMED** This is NOT the full interaction between nonlinear textures. It omits deformation, massless orientation/phase channels, overlap effects and quantum corrections. The displayed Hopf ansatz has algebraic tails; it is not a compact source. Its actual far-field response cannot simply be replaced everywhere by a point-source Yukawa expression. The atomic spectral test below is a diagnostic of the localized amplitude-channel approximation, not a solved Hopf pair.

## 6. Quantum bound-state test and what is imported

**ASSUMED** For two heavy effective objects, adopt the standard nonrelativistic relative-motion equation

$$
[-\hbar^2\nabla^2/(2\mu)+V(d)]\psi=E\psi.
$$

**ASSUMED** This assumes the two objects and their inertial masses. The action above has not derived the electron, proton, their mass ratio, fermionic spin, or their electric charges. The hydrogen reference inserts measured masses and alpha. The screened-potential diagnostic uses that reference to set units and, in one series, fits the 1s energy. Neither operation is an independent prediction of alpha or an atom.

## 7. Vacuum audit: why this candidate is not yet unified

**ASSUMED** Expand $n=(\pi 1,\pi 2,\sqrt{1-\pi 1^{2}-\pi 2^{2}}), \theta =\theta 0+\tau , \rho =v+\sigma$. The quadratic Lagrangian has one massive amplitude scalar and three massless scalars (tau and the two internal orientation fluctuations). F starts at order $\pi ^{2}$, so $F_{\mu\nu}F^{\mu\nu}$ starts at fourth order and supplies NO independent quadratic Maxwell wave. This is checked by explicit amplitude scaling as well as the algebraic expansion.

**ASSUMED** There is no perturbative massless spin-1 photon and no dynamical massless spin-2 metric mode about this vacuum. Calling the composite F an electromagnetic field, or calling a medium ripple a graviton, would not establish those missing properties. This conclusion is about this specified vacuum and action; it is not a no-go theorem for every emergent model or every possible phase.

**ASSUMED** There is also a direct propagation check. For phase perturbations on a static amplitude profile, the principal wave operator is partial_mu[(rho/v)^2 eta^munu partial_nu(theta)]. In $3+1$ dimensions it can be represented by the effective metric $g_{eff}=(\rho /v)^{2} \eta$. This is conformal to the original flat metric, so its geometric-optics null trajectories are unchanged. Varying rho multiplies time and spatial kinetic coefficients equally: the local characteristic phase-wave speed remains c. Finite-wavelength scattering from amplitude gradients is not excluded. However, simply identifying these phase rays with light would not reproduce Einstein light bending. This is a derived kinematic failure of that identification, in addition to the missing dynamical Einstein equation.

**ASSUMED** Adding an independent Maxwell field and an Einstein-Hilbert term would add those structures by assumption. A further trap: naively gauging $\theta via \partial (\theta )-q$ A while rho has nonzero vacuum v gives the vector a mass q v (after canonical Maxwell normalization if $Z=1$). It does not automatically supply long-range electromagnetism. An unbroken U(1) sector would have to be justified.

**ASSUMED** For gravity, adding $S_EH=(M_Pl^{2}/2) \int (\sqrt{-g} R)$ produces the standard Einstein equation with this T_munu, but that imports Einstein dynamics. Quantizing matter on that geometry is not a derivation of a common microscopic origin. Gravity A and B remain undefined; no extra gravitational charge or screening law was inserted.

## 8. Criteria for advancing or rejecting the candidate

### Check of the time-evolving phase in the original intuition

**ASSUMED** Before rejecting a static-vacuum mechanism as sufficient, consider the SAME action on a homogeneous charge-carrying background $\theta =\Omega *t, \rho =\rho 0$. The amplitude equation gives $\rho 0^{2}=v^{2}+\Omega ^{2}/\lambda$. This is a fixed-charge background, not the empty zero-charge vacuum. A time-evolving internal phase is not automatically literal geometric rotation of a fluid sphere.

**ASSUMED** With $\sigma =\delta (\rho )$ and $\tau =\rho 0\,\delta (\theta )$, the quadratic amplitude/phase Lagrangian is

**ASSUMED**

$$
\mathcal L_2=\tfrac12[(\dot\sigma)^2-|\nabla\sigma|^2-M^2\sigma^2+
(\dot\tau)^2-|\nabla\tau|^2]+2\Omega\sigma\dot\tau,
\quad M^2=2\lambda\rho_0^2.
$$

**ASSUMED** Its two dispersion branches satisfy

$$
(k^2+M^2-\omega^2)(k^2-\omega^2)-4\Omega^2\omega^2=0.
$$

**ASSUMED** The gapless phase branch has low-momentum speed

**ASSUMED**

$$
c_\theta^2=\frac{M^2}{M^2+4\Omega^2}
=\frac{\lambda v^2+\Omega^2}{\lambda v^2+3\Omega^2}.
$$

**ASSUMED** This is a real derived change in collective propagation within the same action. Integrating out the gapped amplitude at low momentum increases the phase time-kinetic coefficient relative to its spatial coefficient. However, the two internal-orientation waves still have speed $c_n=1$ at quadratic order. For $v=\lambda =\Omega =1, c_\theta =1/\sqrt{2}$, while $c_n=1$. Thus this background does not automatically supply a universal metric for all excitations. It also supplies no independent perturbative photon or spin-2 field. The numerical work verifies the dispersion polynomial and low-momentum limits. A curved inhomogeneous background and its dynamics have not been solved.

**ASSUMED** The amplitude/texture energy calculation can be retained as a consistency example. Candidate 0 cannot currently serve as a complete hydrogen or unification model: its particle assignments and relevant gauge/gravitational sectors are missing. A fitted screened scalar potential is not a substitute. Further atomic or constant matching is not justified until an action supplies the missing asymptotic degrees of freedom, electric-charge structure and gravitational dynamics.

**ASSUMED** Relevant background, not a claim of replication:

- **CITED** Faddeev and Niemi, Knots and Particles / Nature 387,58 (1997): https://arxiv.org/abs/hep-th/9610193 — eq./page UNVERIFIED.
- **CITED** Samoilenka and Shnir, Magnetic Hopfions in the Faddeev-Skyrme-Maxwell model (2018): https://doi.org/10.1103/PhysRevD.97.125014 . This has an explicitly added gauge field; a composite texture two-form alone is not the same construction. — eq./page UNVERIFIED.
- **ASSUMED** Original local notes: docs/POSSIBILITY_CLOUD_THEORY.md; archive/MASTER_PHYSICS_DOCUMENT.md.
