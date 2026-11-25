# Energy difference when two quantum modes couple

**Conversation date range** 2026-10-25 to 2026-10-30 (archive window)

**ASSUMED** This checkpoint records energy difference when two quantum modes couple. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Run: 2026-10-25. This is an exact calculation in a chosen quadratic quantum model, independently checked by finite-basis diagonalization. It is not a calculation of an actual particle binding energy, a topological transition, or a theory of gravity. Gravity A and B remain undefined and were not introduced.

## Assumptions and energy reference

**DERIVED** The two separate modes have equal modal inertia mu and bare angular frequency omega0. The Hamiltonian is

**DERIVED** $H = (p1^{2}+p2^{2})/(2\mu ) + \mu *\omega 0^{2}\,(q1^{2}+q2^{2})/2 + \eta *\mu *\omega 0^{2}\,q1\,q2$.

**DERIVED** Only |eta|$<1$ has a confined quadratic ground state. q represents a displacement-like generalized coordinate, not a normalized probability density. Quantum wavefunction probabilities are computed from this Hamiltonian. mu is an effective mode inertia, not a mass inferred for a physical particle.

**DERIVED** The same energy reference is used for separated and coupled configurations; the finite change in zero-point energy is retained. No extra coupling-dependent constant is added. A microscopic theory would have to fix any such offset: frequencies alone cannot determine an unspecified constant contribution to the interaction energy. There is no separation coordinate, microscopic mechanism for eta, or derived particle rest-energy scale.

## Exact diagonalization into normal modes

**DERIVED** $Q_{plus}=(q1+q2)/\sqrt{2}$ and $Q_{minus}=(q1-q2)/\sqrt{2}$ give

**DERIVED** $\omega _{plus}=\omega 0\,\sqrt{1+\eta }$, $\omega _{minus}=\omega 0\,\sqrt{1-\eta }$.

**DERIVED** For the quantum ground state, $E_{separate}=\hbar *\omega 0$ and

**DERIVED** $E_{coupled} = \hbar *\omega 0*[\sqrt{1+\eta }+\sqrt{1-\eta }]/2$.

**DERIVED** The available relaxation energy in this model is therefore

**DERIVED** $B = \hbar *\omega 0*{1-[\sqrt{1+\eta }+\sqrt{1-\eta }]/2}$.

**DERIVED** For weak coupling, $B/(\hbar *\omega 0) = \eta ^{2}/8 + O(\eta ^{4})$. This energy shift is a consequence of the chosen cross-coordinate interaction, not a universal result of joining oscillators.

| eta | $\omega _{plus}/\omega 0$ | $\omega _{minus}/\omega 0$ | $E_{coupled}/(\hbar \omega 0)$ | $B/(\hbar \omega 0)$ |
|---:|---:|---:|---:|---:|
| **DERIVED** 0.0 | 1.000000000 | 1.000000000 | 1.000000000 | 0.000000000 |
| **DERIVED** 0.1 | 1.048808848 | 0.948683298 | 0.998746073 | 0.001253927 |
| **DERIVED** 0.2 | 1.095445115 | 0.894427191 | 0.994936153 | 0.005063847 |
| **DERIVED** 0.5 | 1.224744871 | 0.707106781 | 0.965925826 | 0.034074174 |
| **DERIVED** 0.8 | 1.341640786 | 0.447213595 | 0.894427191 | 0.105572809 |
| **DERIVED** 0.95 | 1.396424004 | 0.223606798 | 0.810015401 | 0.189984599 |

## Selected example: eta=0.5

**DERIVED** The two frequencies become 1.2247448714 and 0.7071067812 times omega0. The energy changes from 1 to 0.965925826289 in units $\hbar *\omega 0$. Thus

**DERIVED** $B = 0.034074173711 \hbar *\omega 0$.

**DERIVED** That is a 3.40741737% decrease of these two modes' initial ground energy, NOT a 3.40741737% decrease of an entire particle's unspecified rest mass.

**DERIVED** In $\hbar =\mu =\omega 0=1$ units, final kinetic energy is 0.482962913145, bare potential energy is 0.557677535825, and interaction energy is -0.074714622681. All three are counted. The negative interaction contribution cannot be omitted.

## Amplitude and frequency do not both have to decrease

**DERIVED** Ground-state normal-mode RMS amplitude is $\sqrt{\hbar/(2\mu\omega)}$. The lower-frequency mode consequently has a LARGER quantum RMS amplitude. At $\eta =0.5$, normal-mode RMS amplitudes change from 0.707106781 each to 0.638943104 and 0.840896415. The total ground energy still falls.

**DERIVED** The state develops cross-correlation $<q1 q2> = -0.149429245361$ in these units. For positive eta this anticorrelation lowers the interaction energy. This is not a spatial incidence angle or proof of topological charge.

**DERIVED** The classical formula $E=\mu *\omega ^{2}\,A_{peak}^{2}/2$ applies to a coherent classical oscillation with peak amplitude A_peak. The RMS quantum amplitudes quoted here are fluctuation widths, not that classical peak amplitude. A fixed-amplitude argument cannot be substituted for the quantum ground-state calculation.

## What happens to the released energy?

**DERIVED** Immediately after suddenly switching on the coupling, the initial independent ground state still has expected energy $\hbar *\omega 0$, since its cross-correlation is zero. The lower coupled ground energy does not mean the system instantaneously loses energy. It is initially excited relative to the new ground state.

**DERIVED** We explicitly add a zero-temperature bath for both normal modes, with excitation occupations decaying as $\exp\left(-\gamma *t\right), \gamma =0.2$. The code integrates the occupations and outgoing energy separately and verifies $E_{composite}(t)+E_{exported}(t)=E_{initial}$. The bath's carrier is unspecified: no photon, graviton, or other radiation type is inferred. The damping rate is an input, not a decay prediction.

**DERIVED** For an isolated pair with fixed coupling, energy stays in the pair. If the outgoing energy remains trapped in a larger enclosing system, it remains part of that larger system's total energy. A quasistatic externally controlled change of eta could instead extract work; this model does not describe spontaneous approach and binding of two particles.

## Control: joining need not lower the energy

**DERIVED** For the alternative positive-spring interaction $\eta\mu\omega_0^2(q_1-q_2)^2/2$, the frequencies become omega0 and $\omega 0\,\sqrt{1+2\eta }$. At $\eta =0.5$ its ground energy is $1.207106781187 \hbar *\omega 0$, HIGHER than the uncoupled value by $0.207106781187 \hbar *\omega 0$. An additional attractive contribution would be needed for a net bound-state energy decrease. The form and energy reference of the interaction are essential.

## Standard mass-energy conversion: illustrative scale only

**DERIVED** If one ASSUMES $\hbar *\omega 0=1 eV$, the selected example releases 0.034074173711 eV, or $5.45928449425\times10^{-21} J$. Applying the assumed standard relation $\Delta M=B/c^{2}$ gives a mass decrease of $6.07427319854\times10^{-38} kg$.

**DERIVED** Here $c=299792458 m/s$ and $1 eV=1.602176634\times10^{-19} J$ are exact SI conversions. The 1 eV mode scale is illustrative; it is not inferred from the supplied theory. This calculation does not derive the equivalence principle, a gravitational field, Einstein's equations, or an additional gravitational coupling.

## Numerical checks

**DERIVED** 22 assertions passed: energy decomposition, finite-basis eigenstate residuals, ground-energy convergence at the largest basis, zero initial mean switching work, and energy accounting during relaxation. Finite-basis diagonalizations use 12, 20 and 28 oscillator levels per original mode for $\eta =0.2,0.5,0.8$. These are implementation checks, not physical confirmations.

**ASSUMED** Run `python scripts/binding_mode_energy.py` with NumPy, SciPy, Matplotlib and threadpoolctl. All parameters, exact energies, spectra and relaxation traces are in results.json. The script overwrites only its own experiment directory and output zip.

## Conclusion

**ASSUMED** An explicit cross-coordinate coupling can lower the total quantum ground energy through correlated fluctuations and a changed frequency spectrum. It does not require every amplitude and frequency to decrease. This provides a calculable realization of the proposed energy-reduction idea under stated assumptions. It does not fix the coupling or the physical mode scale and is not yet an atomic, topological, or unified-gravity model.
