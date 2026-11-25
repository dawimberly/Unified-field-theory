# Scalar response matching checkpoint — 30 October 2026

**Conversation date range** 2026-10-30

**DERIVED** This checkpoint records Scalar response matching checkpoint — 30 October 2026. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** This checkpoint advances the free, leading-order auxiliary-metric determinant audit. It does not establish equivalence to constrained composite emergent gravity (CEV), a unitary ultraviolet completion, or unification.

## Definitions and comparison

**DERIVED** Set the physical constituent squared mass to one and write q = Q_E^2 in these units, C = (4 pi)^(-2). Pauli–Villars squared masses are (1,L,2L,5L+1); signed weights cancel the zeroth, first and second mass moments. Matching holds the q and q^2 coefficients of each inverse-response channel at their L=5 values by adding explicit real local polynomials. This changes the zero-bare-term prescription. It is a controlled comparison, not a derived microscopic repair.

**DERIVED** The scalar coefficient A is not itself the normalized scalar projector eigenvalue: K0=2A/3. Its matched coefficients are a=-0.0019490158392333555 and b=-0.0003639828687336804.

## Checks completed

1. **DERIVED** Independent analytic continuation around the physical threshold q=-4 reproduced the adjacent-sheet logarithm for four contour radii and three resolutions, error below 1e-12.
2. **DERIVED** Nine real zeros of the matched scalar/tensor kernels were checked with direct Feynman-parameter Simpson integration at 2,000, 4,000 and 8,000 panels and independent sign brackets. Maximum root residual was 6.4e-14. Matching changes poles; an earlier exclusion argument for the unmodified tensor kernel cannot be transferred to these modified kernels.
3. **DERIVED** A stable series for the subtracted heavy-species contribution verified fixed-q convergence as L increases to 100,000. The error scales as 1/L. This convergence is not uniform across regulator thresholds or the full complex plane.
4. **DERIVED** An independent positive-axis asinh formula verified the large-q coefficient C/16 and the asymptotic constant b-19C/240. At q=1e12, the constant error was about 2.1e-14.
5. **DERIVED** The dispersion-integral root and the elementary-formula root agreed within 9e-14. Independent derivative and dispersion quadrature controls passed.

## Exact restricted result

**DERIVED** For the renormalized physical-species limit with finite real a,b, define

    J(q) = -2 + 2 sqrt(1+4/q) asinh(sqrt(q)/2), q>0,
    R(q) = C/4 [(1-q/2)^2 J(q) - q/6 + 11q^2/60],
    A(q) = a q + b q^2 + R(q).

**DERIVED** The three-subtracted dispersion representation is

    R(q)/q^3 = C/32 integral_0^1 [(3-v^2)^2 v^2 / (4+q(1-v^2))] dv.

**DERIVED** It follows either by substituting s=4/(1-v^2) in the physical cut integral, or by differentiating the elementary logarithm and checking the subtractions. This representation concerns the inverse kernel, not a positive spectral representation of its inverse.

**DERIVED** Let F(q)=A(q)/q^2. If a<0,

    F'(q) = -a/q^2 + C/8 integral_0^1 [(3-v^2)^2 v^2 / (4+q(1-v^2))^2] dv > 0.

**DERIVED** Also F tends to minus infinity at q approaching zero, while

    F(q) = (C/16) log(q) + b - 19C/240 + o(1)

**DERIVED** at infinity. Therefore A has exactly one positive real zero, and it is simple, for every finite b when a<0. This improves numerical existence to an analytic existence-and-uniqueness argument for this specified response. No claim of novelty is made.

**DERIVED** With the reference subtractions, q*=11.3066267794444. The identical constituent-source scalar trace factor 2-q* is approximately -9.3066267794444, so that channel does not cancel this zero. Source normalization and the complete physical observable must still be checked before assigning a physical interpretation.

**DERIVED** Changing b by (-.002,-.001,0,.001,.002) moved the zero to approximately (1403.56,116.729,11.3066,2.38058,1.12513). These examples illustrate the theorem; they are not a parameter fit. The corresponding Euclidean momentum is sqrt(q*) times the constituent mass, not q* times the mass.

## Scope and next unresolved issues

**DERIVED** The uniqueness statement excludes finite signed-PV kernels, additional interactions, higher-order local derivatives, nonlocal ultraviolet modifications and a different constraint/measure prescription. A positive-Euclidean zero in an extrapolated inverse response is not by itself an established physical instability: the EFT validity domain, analytic continuation and actual source observable matter. Finite q^2 adjustments alone cannot erase this zero from the specified extrapolation. Repeating finite scans cannot settle the broader microscopic question.

**DERIVED** Next useful work is to establish an error-controlled domain of validity and inspect the complete sourced response, including the special tuning q*=2 where the identical-source trace vanishes. One must not interpret that channel-specific tuning as removal of the inverse-kernel zero.

## Reproduction

**DERIVED** The companion archive contains the experiment directory and this checkpoint. Run the following scripts from experiments with standard-library Python:

    python matched_regulator_test.py
    python matched_decoupling_test.py
    python scalar_asymptotic_check.py
    python scalar_uniqueness_check.py

**DERIVED** Dependencies are included; no third-party scientific packages are required. Existing JSON result files retain prior runs. The two new analytic-control scripts print their full results as JSON.
