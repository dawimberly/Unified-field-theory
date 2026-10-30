# Exact branchwise constraint identity: a sharper matching target

October 29, 2026. Formal application of a standard delta-function change-of-variables identity; not a new quantum-gravity theorem.

## Result

There is an exact way to embed a specified reduced matter path integral into a constrained metric integral on a nonsingular branch. It does not require guessing that an unrestricted auxiliary integral is equivalent.

Let G[phi] solve T(g,phi)=0 uniquely on the selected branch, and let J_T=delta T/delta g be invertible there. Begin with a specified regulated reduced measure mu_red and observable O_red:

    Z_red[j] = integral Dphi mu_red[phi]
               exp(i S[G(phi),phi] + i j O_red[phi]).

Insert the normalized identity

    1 = integral_branch Dg delta[g-G(phi)]
      = integral_branch Dg delta[T(g,phi)] |det J_T|.

Therefore

    Z_red[j] = integral Dphi Dg mu_red[phi]
               delta[T(g,phi)] |det J_T|
               exp(i S[g,phi] + i j O_red[phi]).

This is exact as a finite-regulator change-of-variables identity for real constraints and a single simple root. Absolute determinants and real delta functions do not define complex Lorentzian integration cycles; that case needs its own continuation. Multiple roots require explicit branch restriction or a sum with specified weights.

## What the identity eliminates

A plain flat delta[T] integral yields an inverse Jacobian after metric integration. A determinant-weighted delta[T] integral yields the chosen reduced measure. They are different definitions unless the extra factor is demonstrably irrelevant in the stated prescription.

An unrestricted integral over g is a third definition. Its saddle expansion gives a fluctuation determinant rather than this normalized delta identity. Classical equality of the reduced actions does not identify all three measures.

This establishes a precise matched constrained representation of an already specified reduced measure. It does NOT select mu_red, prove its gauge invariance, construct a regulator, or settle which representation the original informal notation integral_{T=0} denotes. The exact identity is a representation theorem, not an ultraviolet completion.

## Physical source matching

For a composite-metric observable, choose O_red[phi]=O[G(phi),phi]. Inside the determinant-weighted constrained representation it may be replaced by O[g,phi], since the delta constraint enforces equality on that branch. Consequently matched correlators agree by construction, including source derivatives.

No corresponding equality with unrestricted auxiliary-metric correlators has been established. Integrating an independent metric with a source can generate source-dependent local terms; discarding these without matching is not justified.

This sharpens the scalar-pole question: compute the same defined matter/composite observable in the reduced theory or the exactly matched constrained representation. A pole in an unmatched independent-metric correlator is insufficient evidence about that observable.

## Gauge and canonical qualifications

The identity can be formulated first at a finite regulator before a common admissible gauge quotient. A complete continuum gauge-fixed construction still needs its gauge Jacobian, boundary conditions and regulator. The displayed determinant is a configuration-space change-of-constraint Jacobian, not automatically the Senjanovic determinant of the full canonical constraint algebra. First- and second-class constraints must be classified separately.

## Published basis and boundaries

Henneaux, Physics Letters B 238 (1990), 299-304, treats auxiliary elimination with gauge fixing and measure contributions: https://doi.org/10.1016/0370-2693(90)91739-X . Publisher abstract verified; full paper not obtained.

Criado and Perez-Victoria, JHEP 03 (2019), 038, Section 2, explains the necessity of Jacobian and source transformations for generating-functional equivalence: https://arxiv.org/pdf/1811.09413 . Full relevant section inspected in the preceding batch.

Carone, arXiv:1911.13166, footnote 4, explicitly distinguishes imposing the constraint from independent auxiliary-metric integration: https://arxiv.org/html/1911.13166v1 . The construction above addresses that distinction by specifying a measure, rather than asserting it disappears.

Status: the formal matching target is now explicit. The unresolved physical inputs are the reduced regulated quantum measure and its Lorentzian/gauge definition. No spectrum transfer, stability proof, or unification claim follows yet.