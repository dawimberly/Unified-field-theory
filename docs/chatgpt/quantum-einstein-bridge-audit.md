# Quantum-to-Einstein bridge: derivation and assumption audit

**Conversation date range** $2026-10$ (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records quantum-to-Einstein bridge: derivation and assumption audit. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

## Decision

**ASSUMED** Prioritize known entanglement/geometry derivations and the low-energy effective-field-theory limit over additional custom knot simulations. This work reconstructs selected calculations and identifies the missing microscopic input. It does not solve quantum gravity, simulate a Planck-scale substrate, or reproduce the complete proofs in the references.

## Existing routes compared

| Route | What is controlled | What remains assumed or missing |
|---|---|---|
| **CITED** Donoghue: Quantum gravity as a low-energy effective field theory [1,2] | Quantum corrections to GR can be calculated systematically at low energies. | GR's low-energy fields and couplings are inputs; microscopic origin and Planck-regime completion are not provided. | — eq./page UNVERIFIED.
| **CITED** Jacobson: Local entanglement equilibrium [3] | Conditional relation between first-order entropy stationarity and the semiclassical Einstein equation. | A smooth local geometry, conformal-field assumptions or an extension conjecture, and a finite universal area-entropy coefficient. | — eq./page UNVERIFIED.
| **CITED** Faulkner et al.: Holographic entanglement [4] | The first law plus an entropy/area dictionary constrains linearized gravitational equations near AdS. | The holographic setting and entropy functional are essential; different functionals can yield different gravitational equations. | — eq./page UNVERIFIED.
| **CITED** Faulkner et al.: Nonlinear holographic extension [5] | Einstein dynamics reconstructed to second perturbative order for the specified class of CFT states. | This is not an unrestricted all-orders derivation for our cosmology. | — eq./page UNVERIFIED.

**ASSUMED** The local equilibrium route is selected here because its central coefficients can be independently checked and its assumptions are unusually explicit. This is a practical audit choice, not a claim that it is the uniquely correct theory.

## Reconstructed calculation

**ASSUMED** Use $\hbar =c=k_B=1$, four spacetime dimensions, and a locally flat reference vacuum. Entropy is dimensionless. Denote the Newton constant by G and the Einstein tensor by G_ab. The following is a leading small-ball, first-order reconstruction, not the full tensor or nonconformal-field proof in [3]. — eq./page UNVERIFIED.

**DERIVED** For a three-dimensional ball with small constant sectional curvature k,

**DERIVED** $S_k(r)=\sin (\sqrt{k} r)/\sqrt{k} for k>0,$

**DERIVED** with sinh for $k<0$ and $S_0(r)=r$. Independently calculate

**DERIVED** $V_k(r)=4 \pi  \int_0^r S_k(s)^2 ds,$

**DERIVED** $A_k(r)=4 \pi  S_k(r)^2.$

**DERIVED** Holding volume equal to the flat ball of radius R gives, to leading order,

**DERIVED** delta A at fixed $V = -(4 \pi /5) k R^{4}$.

**DERIVED** At the center of the locally time-symmetric slice, $G_00=3k$. Hence

**DERIVED** delta A at fixed $V = -(4 \pi R^{4}/15) \delta G_00$.

**DERIVED** For a conformal-field vacuum ball, the first-order modular-entropy change for a uniform small energy-density perturbation delta u uses the weight $(R^{2}-r^{2})/(2R)$. Direct integration gives

**DERIVED** $\delta S_{\mathrm{matter}} = 2 \pi \int _{ball} [(R^{2}-r^{2})/(2R)] \delta u dV$

**DERIVED** $= (8 \pi ^{2} R^{4}/15) \delta u$.

**ASSUMED** Now impose the entropy-equilibrium assumption, with a universal coefficient eta:

**DERIVED** $\eta \delta A + \delta S_{\mathrm{matter}} = 0$.

**DERIVED** The result is

**DERIVED** $\delta G_00 = (2 \pi /\eta ) \delta u = 8 \pi G \delta u$,

**DERIVED** so $G=1/(4 \eta )$. In SI units, eta is an entropy-per-area coefficient measured in units of k_B, and $G=c^{3}/(4 \hbar \eta )$.

**CITED** Jacobson: This coefficient matches the local entanglement-equilibrium framework [3]. Requiring the appropriate relation for every local timelike direction is crucial to recovering a tensor equation. Our constant-curvature numerical examples check the scalar geometry coefficient, the modular integral and the equilibrium coefficient; they do not independently prove that all-frame step. Background curvature/cosmological-constant treatment is also outside these flat-reference examples. In particular Lambda's observed value is not derived. — eq./page UNVERIFIED.

## What was actually verified numerically

1. **DERIVED** For finite-dimensional noncommuting density-matrix perturbations, $\Delta <K_0>-\Delta S$ equals relative entropy. Its remainder scales as $\epsilon ^{1.9979}$, confirming the first-order identity's expected quadratic correction. $K_0=-\log\left(\rho _0\right)$ is held fixed. This checks quantum-information algebra, not a CFT or geometry correspondence.
2. **DERIVED** Constant-curvature ball volumes were integrated, their radii adjusted to fix volume, and their areas differentiated independently. Results match the leading area coefficient for $R=0.5, 1$ and 2.
3. **FITTED** For $\eta =0.25, 1$ and 4, the assumed equilibrium reconstructs G approximately 1, 0.25 and 0.0625. Eta was supplied, so no numerical value of G has been predicted.
4. **DERIVED** Three two-qubit states were compared to test whether one simple ordering measure determines entanglement.
5. **DERIVED** A variable-G cosmological substitution was checked against the independent spatial Einstein equation.

**DERIVED** All 33 numerical checks pass. This count is an implementation consistency result; it is not 33 independent empirical tests of emergent gravity. These calculations are inexpensive and reproducible using the supplied script, NumPy, SciPy and Matplotlib.

## Ordering is not a sufficient variable

**DERIVED** The exact two-qubit examples give:

| State | Alignment $<Z_1 Z_2>$ | One-qubit entropy | Entanglement negativity |
|---|---:|---:|---:|
| **DERIVED** Product $up-up$ | 1 | 0 | 0 |
| **DERIVED** Bell $(|\uparrow\uparrow\rangle+|\downarrow\downarrow\rangle)/\sqrt2$ | 1 | ln(2) | 0.5 |
| **DERIVED** Equal classical mixture of $up-up$ and down-down | 1 | ln(2) | 0 |

**DERIVED** Thus equal spin alignment does not determine subsystem entropy. Equal subsystem entropy does not determine entanglement for mixed global states. No value of gravitational coupling can be inferred from these rows: there is no geometry, area law or gravitational dictionary in this two-qubit test. It is a counterexample to treating several meanings of order as interchangeable.

**DERIVED** The eta in the equilibrium framework is a coefficient of the ultraviolet area term, not the entropy of an arbitrary small subsystem and not a phase-alignment score. Under the framework's assumptions it is common to the class of states considered. Ordinary changes in low-energy matter entropy therefore need not change it. If a new microscopic phase changes eta, its gravitational interpretation needs an independent derivation. Conditional on the displayed coefficient relation, a reduction in eta would increase G, not reduce it. This gives no established sign for an undefined ordering mechanism.

## A concrete obstruction to the naive variable-G hypothesis

**ASSUMED** Suppose the proposal is simply

**DERIVED** $G_{ab} + \Lambda g_{ab} = 8 \pi G_{eff}(x) T_{ab}$,

**DERIVED** with Lambda constant, ordinary separately conserved matter, and no additional terms. Taking the covariant divergence gives

**DERIVED** $0 = 8 \pi T_{ab} \partial ^a G_{eff}$,

**DERIVED** using the geometric Bianchi identity and $\nabla ^a T_{ab}=0$. This is generally not satisfied for an arbitrarily varying coupling.

**DERIVED** For an explicit counterexample, take conserved dust in flat $FLRW, \rho =\rho _0 a^{-3}$, and $G(a)=G_0[1+\epsilon \ln\left(a\right)]$. Enforcing the 00 equation sets $H^{2}=8 \pi G(a)\rho /3$. Differentiation yields

**DERIVED** $2 H_\cdot +3 H^{2} = (8 \pi /3) \rho dG/dln(a)$.

**DERIVED** The spatial dust equation would require its left side to vanish, so the prescription fails for epsilon not equal to zero. The supplied script verifies this residual; a constant-G control gives zero. This is an algebraic consistency check, not an observational bound and not a no-go theorem for all variable-coupling theories.

**CITED** Eling, Guedens and Jacobson: Additional covariant dynamics, stress-energy exchange, or geometric terms are needed. In known thermodynamic extensions, modifying the entropy relation can require non-equilibrium entropy production and modified field equations [6]. An arbitrary G(order) fit with these contributions omitted is not a consistent bridge to Einstein gravity. — eq./page UNVERIFIED.

## The microscopic input still required

**DERIVED** The audit does not derive a quantum Hamiltonian, an emergent causal structure, local Lorentz symmetry, universal matter coupling, the numerical value of G, the particle spectrum, or the measured cosmological constant. It assumes local geometry to investigate its equations. A preselected 600-cell supplies geometry; it does not demonstrate that geometry emerges. A graph of qubits is not itself a gravitational theory.

**ASSUMED** A credible next calculation must select a specific microscopic theory and a geometric observable, then derive its area response at a controlled cutoff or renormalization scale. The same low-energy action must account for the gravitational coupling, any extra fields, and conservation laws. If ordering is proposed to alter the coupling, calculate the coefficient change at fixed energy bookkeeping and check the resulting dynamics, rather than equating a change in subsystem entropy with a change in G.

**CITED** Donoghue: The present result supports using entanglement equilibrium as a conditional benchmark. It does not support the assertion that nuclear ordering weakens gravity. The high-energy completion remains open; low-energy quantum corrections can already be treated within effective field theory [1,2]. — eq./page UNVERIFIED.

## Primary references

**CITED** [1] Donoghue, General relativity as an effective field theory: The leading quantum corrections (1994), https://arxiv.org/abs/gr-qc/9405057 — eq./page UNVERIFIED.

**CITED** [2] Donoghue, The Quantum Theory of General Relativity at Low Energies (1996), https://arxiv.org/abs/gr-qc/9607039 — eq./page UNVERIFIED.

**CITED** [3] Jacobson, Entanglement Equilibrium and the Einstein Equation, published version (2016), https://arxiv.org/abs/1505.04753 ; https://doi.org/10.1103/PhysRevLett.116.201101 — eq./page UNVERIFIED.

**CITED** [4] Faulkner et al., Gravitation from Entanglement in Holographic CFTs (2014; revised 2017), https://arxiv.org/abs/1312.7856 — eq./page UNVERIFIED.

**CITED** [5] Faulkner et al., Nonlinear Gravity from Entanglement in Conformal Field Theories (2017), https://arxiv.org/abs/1705.03026 — eq./page UNVERIFIED.

**CITED** [6] Eling, Guedens and Jacobson, Non-equilibrium Thermodynamics of Spacetime (2006), https://arxiv.org/abs/gr-qc/0602001 — eq./page UNVERIFIED.
