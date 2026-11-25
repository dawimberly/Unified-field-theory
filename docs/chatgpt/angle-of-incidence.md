# Minimum-energy coupling angle and the fine-structure hypothesis

**Conversation date range** 2025-03-30 to 2026-10-25

**DERIVED** This checkpoint records minimum-energy coupling angle and the fine-structure hypothesis. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Earlier neutrino and fine-structure discussions considered fractions near $1/137.5$ and an incidence-dependent phase match. **UNVERIFIED** The prior suggestion of a privileged neutrino production angle had no derived rate, kinematic threshold or validated data.

**ASSUMED** A proposed geometry condition was $\cos\theta_i\,(\lambda/R)\approx n+\phi_{\mathrm{spin}}/(2\pi)$. **DERIVED** For $E(\theta)=-J\cos\theta-D\sin\theta$ and nonzero $(J,D)$, the minimum is at $\theta_* = \operatorname{atan2}(D,J)$ modulo $2\pi$.

**DERIVED** Run: 2026-10-25. These are hypothetical, dimensionless model calculations, not physical experiments or a derivation of electromagnetic coupling. No measured alpha entered the energies, quantum Hamiltonians, or trajectories. It entered only the final comparison block. Model choices were made after discussing the target, so this is target-blind computation, not a statistically blinded discovery protocol.

## Connection to the supplied notes

**DERIVED** POSSIBILITY_CLOUD_THEORY.md proposes incidence-dependent resonance, phase locking, a harmonic penalty, and an equation `cos(theta_i)*lambda/R approximately n + phi_spin/(2*pi)`. It does not specify an interaction energy between two configurations, calculate lambda/R from a solution, define the incidence surface, or fix an electromagnetic field normalization. The following chosen models test what is possible and what remains unconstrained; they are not presented as consequences already derived from those notes.

## 1. A stable preferred angle can emerge from competing interactions

**DERIVED** Choose `E(theta) = -J cos(theta) - D sin(theta)` with J>0. The minimum is `theta*=atan2(D,J)`, with positive curvature sqrt(J^2+D^2). Here theta is relative orientation, not angle between an incoming velocity and a surface normal. Independent bounded numerical minimization checked the analytic result.

| D/J | Preferred orientation (degrees) |
|---:|---:|
| **DERIVED** 0 | 0.000000000 |
| **DERIVED** 0.001 | 0.057295760 |
| **DERIVED** 0.01 | 0.572938698 |
| **DERIVED** 0.1 | 5.710593137 |
| **DERIVED** 1 | 45.000000000 |

**DERIVED** Multiplying J and D by the same scale leaves the preferred angle unchanged. Independently changing either one changes that angle. Energy minimization selects the angle for GIVEN coefficients; it does not select the coefficients. A small chiral interaction naturally produces a small angle when D/J is small, but its numerical value remains free. D=0 selects zero angle.

## 2. A symmetric model can also select nonzero angles

**DERIVED** For `E=-J cos(theta)+K cos(theta)^2`, theta*=0 for K<=J/2. For K>J/2, the degenerate minima are `theta*=+/-acos(J/(2K))`. Small angles are possible near the bifurcation K=J/2. Both signs have the same energy; no unique handedness is selected.

**DERIVED** To arrange theta*=alpha would require K/J=0.500013313134004, or 2K/J-1=2.66262680089e-05. This is an inverse fit. No mechanism in the tested model pins that distance from criticality.

## 3. Capture angle is different from bound orientation

**DERIVED** We simulated a unit-mass incoming particle in an invented two-dimensional anisotropic interaction well:

**DERIVED** `V(r,beta) = r^-12 - 2 r^-6 - 0.2 exp[-(r-1)^2/(2*0.3^2)] (J cos(beta)+D sin(beta))`,

**DERIVED** with J=1,D=0.1. beta is position direction relative to the target's +x axis. The core is repulsive; the attractive well has a preferred bound orientation. beta is NOT the incidence angle. Incidence was separately measured as `acos(-v dot rhat / |v|)` at the first inward crossing of r=2. Choice of r=2 is an explicit measurement convention, not a derived particle surface.

**DERIVED** Initial positions are (-4,b) for 17 equally spaced b from -1.5 to 1.5, and speeds are 0.35 or 0.7 in the +x direction. We compared zero drag with local drag `-0.6 envelope(r) v`, over at most 60 dimensionless time units, ending early at outward r=5. The finite starting distance is not an asymptotic scattering boundary. The drag is an explicit bath approximation for energy export, NOT an isolated quantum process or radiation derived from the theory.

**DERIVED** Energy-bound at the end means negative mechanical energy and radius less than 2; it does not prove infinite-time stability. Escape, temporary residence, and energy-bound states were recorded separately. Independent accounting of drag work verifies the mechanical-energy changes.

| Incoming speed | Damping | Energy-bound / trajectories | Incidence range among energy-bound trajectories (degrees) |
|---:|---:|---:|---|
| **DERIVED** 0.35 | 0.0 | 0 / 17 | none |
| **DERIVED** 0.35 | 0.6 | 17 / 17 | 0.00105 to 37.90160 |
| **DERIVED** 0.7 | 0.0 | 0 / 17 | none |
| **DERIVED** 0.7 | 0.6 | 17 / 17 | 0.00034 to 45.05857 |

**DERIVED** The angle that minimizes the bound potential is not a universal arrival angle. The results depend on speed, impact parameter, orientation and the assumed energy-export mechanism. This is the useful part of the enzyme analogy: encounter geometry, access to a favorable state, and release of energy are separate questions. The model is not a chemical enzyme simulation.

## 4. Quantization does not fix the free angle

**DERIVED** We quantized an angular rotor: `H = kinetic*Lz^2/2 - J cos(theta)-D sin(theta)` with hbar=1 and kinetic=1. Its ground-state angular distribution has circular mean atan2(D,J); the spectrum depends only on sqrt(J^2+D^2). This spectral identity was checked numerically. Angular momentum cutoffs 20 and 40 agree for D/J=0.01.

**DERIVED** At D/J=0.01 the mean is 0.0099996667 radians, but its circular standard deviation is 0.76877043 radians. A preferred mean need not imply every individual interaction has one exact angle. The angular spread is model-dependent.

**DERIVED** This is ordinary quantum mechanics assumed from the start, with an input potential. It has no dynamical spacetime and does not unify quantum theory with gravity.

## 5. The electromagnetic coupling is not determined by the angle

**DERIVED** In natural units hbar=c=1, a possible electromagnetic extension is `L=-Z F_munu F^munu/4 + g j^mu A_mu`, with a conserved current and Z>0. Canonically normalizing the field gives `e_eff=g/sqrt(Z)` and `alpha_eff=g^2/(4*pi*Z)`. These field/current normalizations would have to be fixed by the underlying theory. They are not fixed by J or D in any tested model.

**DERIVED** We held J=1,D=0.01 fixed and varied g over 0.1,0.3,1 and Z over 0.5,1,2. The same preferred angle then permits nine different electromagnetic couplings. These are algebraic identifiability examples, not nine simulated gauge-field theories. The absence of a unique mapping is the obstruction: orientation alone cannot determine alpha. Assigning alpha_eff=theta* by definition would add an assumption, not derive the coupling.

## Reference comparison, performed after the computations

**DERIVED** The NIST/CODATA 2022 reference is alpha=0.0072973525643, standard uncertainty 1.1e-12. If interpreted as a radian angle, that is 0.4181075036 degrees.

**DERIVED** For the first angular energy, theta*=alpha requires D/J=tan(alpha)=0.007297482098362. Inserting that ratio would reproduce the target by construction. None of the predetermined ratios independently selects alpha. The handful of cases is not an exhaustive search over possible theories.

## Result and next missing derivation

**ASSUMED** Lowest-energy angular selection works in both classical and quantum illustrative models. It has not selected the observed fine-structure constant. The tests separate three ingredients that must be connected: incidence-dependent capture, the final bound orientation, and the canonically normalized electromagnetic coupling. A proposed underlying action must determine all three and also recover gravitational dynamics. Those derivations remain open.

**DERIVED** Useful next work is to specify the microscopic field, its symmetries, conserved charge and two-configuration interaction. Extract J,D and the gauge/current normalization from that one model, rather than fit independently chosen coefficients. The archive's resonance condition alone is underdetermined until its wavelength, radius, mode and spin-phase terms are supplied by solutions.

## Verification and reproduction

**DERIVED** 170 numerical checks passed, including minima, force gradients, trajectory energy accounting and quantum spectra. They check implementation, not experimental confirmation. All trajectories, quantum distributions and parameters are in results.json.

**DERIVED** Run `python scripts/coupling_angle_checks.py` with NumPy, SciPy, Matplotlib and threadpoolctl. Only this script's own output directory and zip are overwritten. The toy models are deterministic; no stochastic angle search or integer-formula fitting was performed.

## Sources

- **DERIVED** Local source: docs/POSSIBILITY_CLOUD_THEORY.md, especially morphology resonance and the resonance-condition sketch.
- **DERIVED** NIST/CODATA 2022 recommended constants: https://physics.nist.gov/cuu/pdf/all.pdf . Used only as the final comparison target.
- **DERIVED** Antisymmetric spin interactions are a physical motivation for sine/cosine angular terms, not a claim that these illustrative coefficients describe a particular material: Cardias et al., First-principles Dzyaloshinskii-Moriya interaction in a non-collinear framework (2020), https://doi.org/10.1038/s41598-020-77219-3 .
