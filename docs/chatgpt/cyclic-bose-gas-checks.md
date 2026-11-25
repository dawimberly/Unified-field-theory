# Cyclic interactions and excited quantum order: finite-system test

**Conversation date range** 2026-10-25 to 2026-10-30 (archive window)

**ASSUMED** This checkpoint records cyclic interactions and excited quantum order: finite-system test. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Run: 2026-10-25. This is a three-boson quantum calculation on a finite ring with a finite momentum cutoff, not an experimental reproduction or a demonstration of a new phase of matter.

## Hamiltonian and protocol

**DERIVED** Units: $\hbar =1, \hbar ^{2}/(2m)=1$, ring length $L=2\,\pi$. Bosonic modes have integer momenta $k=-M$,...,M. We work in fixed $N=3$, total momentum zero and even parity. These symmetries are imposed exactly; finite momentum truncation does not preserve the full continuum integrable structure.

**DERIVED** $H(g) = \sum _k k^{2} a_k^dagger a_k + g/(2L) \sum _(p+q=r+s) a_p^dagger a_q^dagger a_r a_s$.

**DERIVED** Start with the noninteracting condensate. Ramp g from 0 to +G; instantaneously switch +G to -G; ramp from -G back to 0. $G=20$. Each ramp has duration $T=15$ with a sine-squared schedule. Exact unitary propagation under the midpoint Hamiltonian is used for each step, primarily $dt=0.08$. This is a finite-G jump, not the paper's ideal crossing through infinite coupling or its measured Feshbach protocol. No continuum interaction-counterterm fit is used.

**CITED** The contact-gas model is motivated by Bastianello et al., https://arxiv.org/abs/2602.17656 and Zeng et al., https://arxiv.org/abs/2602.17657. Their generalized-hydrodynamic thermodynamic results are NOT assumed to follow from this small calculation. — eq./page UNVERIFIED.

## Momentum-cutoff checks

| Cutoff M | Even-parity dimension | Final energy (ground$=0$) | $n(k=0)$ | Particle fraction in outer four modes |
|---:|---:|---:|---:|---:|
| **DERIVED** 8 | 25 | 8.04985383 | 0.99346497 | 0.000435 |
| **DERIVED** 12 | 49 | 8.21353291 | 0.99291793 | $8.28\times10^{-05}$ |
| **DERIVED** 16 | 81 | 8.35083464 | 0.99268823 | $2.85\times10^{-05}$ |
| **DERIVED** 24 | 169 | 8.58679027 | 0.99465393 | $5.71\times10^{-06}$ |

**DERIVED** At largest cutoff, halving dt gives state fidelity 0.9998144318 and energy difference 0.05279014.

**DERIVED** A small population at the momentum boundary is necessary but not sufficient for convergence of a contact gas: high-momentum virtual scattering also matters. Differences between the cutoff results must be retained as a material limitation, especially at strong coupling.

## Protocol controls

| Preparation | Final energy | Zero-mode population | Ground-state overlap |
|---|---:|---:|---:|
| **DERIVED** one cycle, refined dt | 8.53400013 | 0.99426899 | 0.00295384 |
| **DERIVED** one cycle, $T=7.5$ | 8.15118932 | 0.94071664 | 0.00989998 |
| **DERIVED** one cycle, $T=30.0$ | 8.53673899 | 0.98973594 | 0.00018875 |
| **DERIVED** two cycles | 28.73739278 | 0.82629505 | 0.01146127 |
| **DERIVED** sudden attraction control | 11.04407565 | 1.03301304 | 0.01752556 |

**DERIVED** The direct-attraction control is a sudden 0 to -G change followed by the same attractive return ramp; it lacks the repulsive preparation. The two-cycle run uses the same protocol twice, without resetting the state. A high final energy alone is not evidence of order: momentum distributions and correlations are recorded separately.

**DERIVED** In the refined one-cycle output, probability of the specific free Fock pattern (-2,0,2) is 0.91825614. For two cycles the probability of (-4,0,4) is 0.46296574. These are diagnostics of candidate finite-system spectral-flow states, not fractional occupancies or a Fermi-sea proof.

**DERIVED** A canonical ensemble restricted to the same $N=3$, zero-momentum, even-parity sector was matched to the refined cycle's energy. Its inverse temperature is 0.07989324; the momentum-distribution distance from the cycle is 0.42213699. This distinguishes this prepared distribution from that particular thermal reference, not from every possible ensemble.

## Retention and a constraint-breaking perturbation

**DERIVED** We hold the prepared refined-dt state for $t=0$..100 under four time-independent Hamiltonians. The added three-body operator is $W=\int(\psi^{\dagger 3}\psi^3)/6$ with coefficient 0.2. In this $N=3, P=0$ sector its matrix is rank one. It generically changes the continuum contact-gas integrability, but finite cutoff already perturbs integrability. This is not a controlled continuum integrability-breaking scaling study. W preserves particle number, total momentum and parity.

**DERIVED** Occupation distance is sum_k |$n_k(t)-n_k(0)$|/(2N). It is 0 for identical distributions and at most 1. It is not full-state fidelity, not a thermalization measure, and not proof of a finite lifetime.

| Hold coupling g | Three-body coefficient | Maximum occupation distance | Mean distance in latter half of hold |
|---:|---:|---:|---:|
| **DERIVED** 0.0 | 0.0 | $3.0149199\times10^{-16}$ | $1.5451357\times10^{-16}$ |
| **DERIVED** 0.0 | 0.2 | 0.0056875026 | 0.002482083 |
| **DERIVED** 1.0 | 0.0 | 0.10194868 | 0.051448685 |
| **DERIVED** 1.0 | 0.2 | 0.11143555 | 0.05681887 |

**DERIVED** At $g=0$ without W, every individual momentum occupation is conserved exactly by the Hamiltonian; their persistence is therefore an expected control, not a novel discovery. The full state can still develop relative phases, so conserved occupations do not prove all correlations are stationary.

**DERIVED** The tested three-body perturbation causes only a small redistribution at $g=0$, while restoring two-body interactions changes the occupations more substantially. This is not evidence that every perturbation destroys the state, nor that observed changes measure a protection lifetime.

**DERIVED** At nonzero g and/or W, occupations may evolve even while total energy remains constant. Such changes do not by themselves prove destruction of quantum order. The finite isolated system can show recurrences and does not have a bath that would force it to its ground state.

## Relation to the larger physics proposal

**DERIVED** This tests access to excited states through cyclic quantum interactions and distinguishes conserved observables from energy minimization. It does not establish a fractional Fermi sea, thermodynamic criticality, emergent fermionic statistics, an incidence-angle law, the electromagnetic fine-structure constant, or Einstein gravity. $N=3$ is too small for a thermodynamic phase claim. No statistics parameter is identified with the fine-structure constant.

**DERIVED** A faithful replication would require a converged continuum calculation or the published generalized-hydrodynamic equations, an appropriate thermodynamic regime, correlation/rapidity diagnostics, and the actual interaction protocol. A claim of protected interacting order also needs additional conserved-charge and finite-size checks.

## Verification and reproducibility

**DERIVED** 39 implementation assertions passed. Checks include contact and three-body normalization, Hermiticity, symmetry basis, unitary norm, temporal refinement, hold energy conservation, and the exactly conserved free-mode occupations. These are numerical checks, not experimental confirmations.

**ASSUMED** Run `python scripts/cyclic_bose_gas_checks.py`. Dependencies: NumPy, SciPy, Matplotlib, threadpoolctl. Parameters, runtime versions, time traces and quantum states are included. The script overwrites only its own experiment outputs.

## Additional temporal convergence check

**DERIVED** A further run at $dt=0.02$ has fidelity 0.9999940565 against $dt=0.04$ and final energy difference 0.00078126427. Its candidate-pattern probability is 0.91829105. Hold and canonical-reference calculations above still use $dt=0.04$. This extra check assesses their input-state temporal error rather than silently replacing that state.
