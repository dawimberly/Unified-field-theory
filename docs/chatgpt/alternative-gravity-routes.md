# Three alternative gravity routes: concrete calculations

**Conversation date range** $2026-10$ (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records three alternative gravity routes: concrete calculations. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

## Outcome and choice

**DERIVED** Three known-model calculations were run instead of changing the earlier knot model. All 47 implementation checks pass. These are reproductions and mathematical diagnostics, not new experimental evidence or a solution to quantum gravity.

**DERIVED** Induced gravity is the most direct next analytic target for explaining the origin of an Einstein-like action. Dynamical geometry is the more direct target for explaining the existence of geometry itself, but the pilot here is only $1+1$ dimensional. Compactification provides a concrete mechanism for scale-dependent apparent strength, while assuming gravitational geometry and its action from the start.

**DERIVED** The prior entanglement-equilibrium audit remains a fourth comparison baseline: it relates a gravitational coupling to an assumed area-entropy coefficient but does not determine the microscopic system. That result was not rerun or counted as one of these new calculations.

## 1. Compact extra-dimensional propagation

**ASSUMED** Use an unwarped spatial product $R^{3} x S^{1}$, with circle circumference 2 pi R, and a scalar point-source Green function. Let r be separation in $R^{3}$ at equal circle coordinate. An exact periodic image sum gives

**DERIVED** $V(r) = -C \sum _n [r^2+(2 \pi  R n)^2]^-1$

**DERIVED** $=-C\coth(r/(2R))/(2Rr)$.

**DERIVED** Define $G4 M=C/(2R)$ to normalize the long-distance scalar force. Then

**DERIVED** $g(r)=G4 M/r^2 [coth(x)+x csch^2(x)], x=r/(2R).$

**DERIVED** The same bracket is $1+2 \sum _{n\ge 1} [1+n r/R] \exp\left(-n r/R\right)$. Direct image sums, the mode sum and the closed expression agree numerically. Differentiating gives a force exponent approaching 3 at short distance and 2 at long distance.

| $r/R$ | Force / long-distance 3D extrapolation |
|---:|---:|
| **DERIVED** 0.01 | 400.0000 |
| **DERIVED** 0.1 | 40.0000 |
| **DERIVED** 1 | 4.0053 |
| **DERIVED** 3 | 1.4356 |
| **DERIVED** 10 | 1.0010 |
| **DERIVED** 30 | 1.0000 |

**DERIVED** This is a controlled geometric example of stronger short-distance response. It uses MORE available spatial directions at short scales, rather than fewer. It is a scalar propagation calculation, not the complete tensor potential of a compactified Einstein theory: extra polarizations, radion dynamics and matter localization can affect the actual gravitational answer.

**DERIVED** Separately, reducing the zero-mode Einstein-Hilbert action on a fixed flat circle gives $M4^{2}=(2 \pi R) M5^{3}$, or $G4=G5/(2 \pi R)$, in the corresponding action normalization. Holding the bulk coupling fixed, halving R doubles G4. Thus dimensional shrinkage does not give weakened long-distance coupling in this simplest setup. The scalar C normalization and the Einstein action relation are distinct calculations and are not identified as a complete relativistic model.

**CITED** Arkani-Hamed, Dimopoulos and Dvali: Large-extra-dimension models use related volume dilution to address gravitational weakness [1]. The illustrative one-circle computation here is not a fit to their phenomenology. No compactification radius is inferred, and experimental viability is not evaluated. This approach assumes extra geometry and higher-dimensional gravity; it does not explain their microscopic origin. — eq./page UNVERIFIED.

## 2. Quantum-induced curvature

**DERIVED** For one real scalar on an assumed smooth Euclidean background, specify

**DERIVED** $D=-laplacian+m^2+\xi R,$

**DERIVED** $W=+\frac{1}{2} Tr \ln  D = -\frac{1}{2} \int_{\Lambda _UV^-2}^\infty (ds/s) Tr \exp (-sD),$

**DERIVED** up to the usual normalization/subtractions. With positive scalar curvature for a round sphere, the leading heat-kernel terms are

**DERIVED** $\operatorname{Tr}e^{-sD}=(4\pi s)^{-2}\int\sqrt{g}\,e^{-m^2s}[1+s(1/6-\xi)R+\cdots]$.

**DERIVED** Consequently the induced local term is

**DERIVED** $W_R=-(\frac{1}{32} \pi ^2)(\frac{1}{6}-\xi) I \int \sqrt{g} R,$

**DERIVED** $I=\int_{\Lambda _UV^-2}^\infty ds \exp (-m^2s)/s^2.$

**DERIVED** Comparison with the Euclidean Einstein-Hilbert convention $-\int \sqrt{g} R/(16 \pi G)$ yields

**DERIVED** $\delta (1/G)=(\frac{1}{6}-\xi) I/(2 \pi ).$

**CITED** Solodukhin: For a minimally coupled scalar, $\xi =0$, this becomes $I/(12 \pi )$. These conventions and signs are explicit; the numerical interpretation follows the scalar heat-kernel/induced-Newton-constant treatment in [2], rather than interchanging overall signs between effective-action conventions. General induced-gravity context is reviewed in [3]. — eq./page UNVERIFIED.

**ASSUMED** Set $x=m/\Lambda _UV$. Independent quadrature agrees with

**DERIVED** $I/\Lambda _UV^{2}=\exp\left(-x^{2}\right)-x^{2} E1(x^{2})$,

**DERIVED** with limit 1 at $x=0$.

| $m/\Lambda _UV$ | Curvature integral / massless value |
|---:|---:|
| **DERIVED** 0 | 1 |
| **DERIVED** 0.1 | 0.949671 |
| **DERIVED** 0.3 | 0.741244 |
| **DERIVED** 1 | 0.148496 |
| **DERIVED** 2 | 0.00319823 |

**DERIVED** If this were the only gravitational kinetic contribution, increasing the mass gap at fixed cutoff would reduce $1/G$ and therefore INCREASE G. That is not the same as a prediction about ordering: no microscopic ordering process has been modeled, and there may be bare and other loop contributions. Changing the curvature coupling also changes the answer. At $\xi =1/6$ the displayed R contribution vanishes; at $\xi >1/6$ it reverses sign. A lone wrong-sign contribution cannot be accepted as a healthy Einstein kinetic term simply by taking its absolute value.

**DERIVED** The regulator scale, field species, couplings and any bare terms are inputs. Choosing the cutoff to match measured G would be calibration, not a prediction. The same loop generates vacuum-energy and higher-curvature terms that must also be accounted for. The calculation does not solve the cosmological-constant problem or manufacture spacetime from nothing. Its success is narrower: a specified quantum field determinant produces a curvature response of the Einstein-Hilbert form on a given geometry.

## 3. A sum over causal geometries

**DERIVED** Implement the Euclideanized $1+1D$ causal-dynamical-triangulation strip amplitude, with a circular spatial slice of integer length m evolving to length n. For an entrance-marked boundary the one-step kernel is

**DERIVED** $T_mn = binom(m+n-1,m-1) g^(m+n), m,n>=1, 0<g<\frac{1}{2}.$

**DERIVED** The binomial coefficient counts allowed strip triangulations. It is not a force imposed on pre-existing points. A similarity transformation gives the symmetric operator

**DERIVED** $S_mn=\sqrt{mn}/(m+n) binom(m+n,m) g^(m+n).$

**CITED** Ambjorn and Loll: This is the known model in [4]. For a finite maximum loop length N, diagonalize S. For the infinite-length operator, let — eq./page UNVERIFIED.

**DERIVED** $F=2g/[1+\sqrt{1-4g^2}].$

**DERIVED** The leading eigenvector is proportional to $\sqrt{l} F^l$, its eigenvalue is $F^{2}$, and the associated normalized transfer-ensemble stationary distribution is p(l) proportional to $l F^{2l}$. Its mean length is $(1+F^{2})/(1-F^{2})$. The stationary weight can be obtained through the ground-state transform of the positive Euclidean transfer operator; it is not a Lorentzian measurement probability assigned to a single time-evolved wavefunction.

**DERIVED** The calculation matches the known eigenpair away from strong cutoff effects. Near $g=0.4999$, the analytic mean is 50.002500 lattice units:

| Length cutoff N | Computed mean |
|---:|---:|
| **DERIVED** 64 | 22.682714 |
| **DERIVED** 128 | 38.482599 |
| **DERIVED** 256 | 49.294224 |
| **DERIVED** 512 | 50.002282 |

**DERIVED** This exposes a substantial finite-size bias rather than hiding it. Removing triangulation multiplicities changes the length distribution, confirming that counting geometries has a real role in this ensemble.

**CITED** Loll: Topology is fixed to a cylinder; lattice spacing and causal construction are inputs. The experiment does NOT demonstrate dimensional emergence into $3+1$ dimensions, a de Sitter universe in our dimensions, local gravitons, or recovery of four-dimensional Einstein dynamics. Those are separate questions in higher-dimensional CDT [5]. In two spacetime dimensions the Einstein-Hilbert term is topological, so this pilot cannot determine the ordinary four-dimensional Newton coupling. We have reproduced a known low-dimensional quantum-geometric kernel, not the full CDT program. — eq./page UNVERIFIED.

## What to pursue and reject

1. **DERIVED** Retain induced curvature as a specific bridge from quantum fluctuations to a gravitational action. Next, specify a field spectrum and renormalization/matching conditions and account for every generated term. A loop coefficient alone is insufficient.
2. **DERIVED** Retain compactification as a scale-dependent response mechanism, while recognizing that compact-radius stabilization and full tensor dynamics are required. Reject the claim that shrinking a dimension must weaken gravity.
3. **DERIVED** Retain exact low-dimensional geometry sums as a validated computational starting point. Do not advertise them as four-dimensional Einstein recovery. Moving to $3+1D CDT$ requires a different and significantly larger computation.
4. **DERIVED** Do not identify generic spin alignment, reduced variance, subsystem entropy, spectral dimension and gravitational coupling. The previous audit already showed why these are different quantities.

**ASSUMED** No route has been selected as the experimentally established microscopic theory of our universe. The useful progress is a set of reproducible calculations with explicit missing inputs and two counterexamples to the proposed direction of weakening.

## Reproduce

**ASSUMED** Run alternative_gravity_routes.py in a scripts directory. It writes experiments/alternative_gravity_routes relative to the scripts' parent. Packages: NumPy, SciPy, Matplotlib and threadpoolctl. It runs quickly on an ordinary CPU. Numeric outputs are in results.json and the comparison plot in three_routes.png. The script validates identities and known-model limits, not empirical predictions.

## Primary sources

**CITED** [1] Arkani-Hamed, Dimopoulos and Dvali, The Hierarchy Problem and New Dimensions at a Millimeter (1998), https://arxiv.org/abs/hep-ph/9803315 — eq./page UNVERIFIED.

**CITED** [2] Solodukhin, Newton constant, contact terms and entropy (2015), especially sections 2, 3.1 and 3.5, https://arxiv.org/abs/1502.03758 — eq./page UNVERIFIED.

**CITED** [3] Visser, Sakharov's induced gravity: a modern perspective (2002), https://arxiv.org/abs/gr-qc/0204062 — eq./page UNVERIFIED.

**CITED** [4] Ambjorn and Loll, Non-perturbative Lorentzian Quantum Gravity, Causality and Topology Change (1998), especially equations $8-17$, https://arxiv.org/abs/hep-th/9805108 — eq./page UNVERIFIED.

**CITED** [5] Loll, Quantum Gravity from Causal Dynamical Triangulations: A Review (2020), https://arxiv.org/abs/1905.08669 — eq./page UNVERIFIED.
