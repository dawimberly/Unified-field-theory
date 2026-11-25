# Route 6: numerical target and prerequisites

**Conversation date range** $2026-10$ (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records route 6: numerical target and prerequisites. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

## Objective

**ASSUMED** Test the constrained composite theory's connected composite-metric correlation function in an interacting vacuum. The eventual target is a consistent long-distance gravitational response, not a visually attractive field pattern or an arbitrary tensor correlation.

**ASSUMED** No such simulation has been run in this work. The current code performs a small Gaussian quadrature and scalar-loop controls, not Monte Carlo sampling of the composite-gravity theory.

## Define the theory before selecting an algorithm

**ASSUMED** Specify the full independent-quartic determinant action, integration variables and measure, covariant regulator or its controlled lattice replacement, gauge treatment, counterterms, physical observables and vacuum-selection condition. The formal metric-independent continuum action does not automatically supply a normalized, positive Euclidean sampling distribution. Gauge volume, degenerate determinants, flat directions, branches and any sign problem must be assessed explicitly. A naive Metropolis or HMC run cannot be called a simulation of this theory without resolving those issues.

**ASSUMED** The first executable numerical control, if needed, can be an ordinary O(N) scalar theory with a quartic interaction and its well-defined lattice measure. It would test machinery such as estimator variance, projectors and weak-coupling mass insertions. It would not demonstrate a composite graviton and must be labeled accordingly.

## Observables

1. **ASSUMED** Select and renormalize a composite metric operator h_ab constructed from the underlying fields according to the actual theory. Specify the relational/gauge context in which its correlation is meaningful.
2. **ASSUMED** Measure the connected correlation $C_{abcd}(q)=<h_{ab}(q)h_{cd}(-q)>-<h_{ab}(q)><h_{cd}(-q)>$ in a valid translationally invariant reference setup. On curved backgrounds use the appropriate harmonic or position-space construction rather than importing flat Fourier momenta unchanged.
3. **ASSUMED** Project into TT and scalar sectors using momentum conventions consistent with the discretized derivative. The continuum projector used by this audit cannot be blindly reused with an unrelated lattice derivative or Laplacian.
4. **ASSUMED** Also measure background one-point functions, local contact contributions and Ward-identity residuals. Zero momentum is not defined by the present q-dependent projector and must be handled separately.

## Controls and possible failure signals

- **ASSUMED** Reproduce the free/base-model limit and the weak-quartic analytic result once it exists.
- **ASSUMED** Compare genuinely connected correlations, not disconnected vacuum products.
- **ASSUMED** Repeat with larger physical volume and a sequence of spacings at matched parameters. A finite-volume correlation length or an isolated small eigenvalue is not enough.
- **ASSUMED** Look for the expected low-momentum tensor pole only in the FULL appropriate correlator. A matter stress bubble is a building block of an inverse response and need not have $a 1/q^{2}$ pole itself.
- **ASSUMED** Evaluate spin-zero contamination and positivity/unitarity through appropriate observables. The Euclidean TT projector has rank five; this does not count physical massless helicities.
- **ASSUMED** Distinguish numerical symmetry breaking from a physical interaction effect. Quantify autocorrelation and statistical/systematic uncertainty.
- **ASSUMED** Comparing ordered/disordered states requires fixed microscopic parameters with a specified state or thermodynamic control, plus the correct extra medium tensor structures. Changing lambda alone is not a binding transition.

## Gates before a gravity claim

**ASSUMED** An apparent long-range TT signal must survive volume and continuum controls and be accompanied by the relevant conservation identities, acceptable spectral behavior and universal coupling to distinct probe species. Even that would be a result about a model, not proof that nature uses it. The ultraviolet completion and observational connection remain separate requirements.
