# Scalar pole audit: cancellation in a constituent-matter channel

**Conversation date range** 2026-10-29 to 2026-10-30

**DERIVED** This checkpoint records Scalar pole audit: cancellation in a constituent-matter channel. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Research batch, October 29, 2026 (America/Chicago). Provisional derivation requiring independent review.

## Result and scope

**DERIVED** The round-5 auxiliary-metric scalar response has an additional pole. However, its residue vanishes when saturated with the on-shell, minimally coupled constituent-scalar vertex derived below. A generic external trace source is not the same as this matter source. This revises the interpretation of the earlier generic-source residue calculation; it does not establish vacuum stability or a unified theory.

**DERIVED** Assumptions: the saved round-5 quadratic effective action; leading N; dimensionally continued 4<D<6; flat tuned saddle; minimal curvature coupling xi=0; constituent mass M; Euclidean analytic continuation of on-shell amplitudes; the auxiliary-metric formulation. Equivalence to constrained CEV remains open. All formulas below suppress common external-field normalization factors, which cannot affect pole cancellation.

## Exact matrix contraction

**DERIVED** Use z=Q_E^2, a=(D-2)/4, U=M^2-a z, and

    Pi(z) = (4 pi)^(-D/2) Gamma(2-D/2)
            integral_0^1 [M^2+x(1-x)z]^(D/2-2) dx,
    g = (D-2)/2 * (4 pi)^(-D/2) Gamma(1-D/2) M^(D-2).

**DERIVED** Pi(0)=-g/M^2 and g>0. The round-5 scalar quadratic form is S2=A f^2+2 B f s+C s^2, with

    F = Pi U+g,
    A = -U F/4,
    B = -F/4,
    C = -(Pi+1/lambda)/4,
    det K = F (U/lambda-g)/16.

**DERIVED** Its additional zero is U=lambda*g. At that zero the null vector is v=(1,-lambda*g). Arbitrary external sources can couple to it.

**DERIVED** For the physical constituent vertex, let incoming and outgoing Euclidean continued momenta obey k^2=k'^2=-M^2 and q=k'-k. Then k.k'=-M^2-z/2. The minimal stress vertex is

    t_ab = (k_a k'_b+k_b k'_a)/2 - delta_ab (k.k'+M^2)/2,
    q_a t_ab = 0,
    trace(t) = (1-D/2) k.k' - D M^2/2 = -U.

**DERIVED** Varying the Euclidean action with respect to the covariant metric gives -h.t/2, while s enters as +s phi^2/2. Thus the on-shell scalar pair couples to U f+s: its source direction is J proportional to (U,1). Conservation permits replacing the transverse trace mode by the conformal mode in this contraction. The overall sign of J is immaterial; the relative sign is essential.

**DERIVED** At the extra pole J.v=U-lambda*g=0. More explicitly,

    W_matter = (1/4) (U,1) K^(-1) (U,1)^T
             = [C U^2 -2 B U + A]/(4 det K)
             = -U/(Pi U+g).

**DERIVED** The factor U/lambda-g cancels exactly. For lambda=0, sigma is absent and the same result follows from U^2/(4A). At U=0 that expression tends to zero, rather than diverging. For nonzero lambda the finite limit at U=lambda*g is -lambda/(1+lambda Pi), provided that denominator is nonzero. No statement here covers coincident degeneracies at strong coupling.

**DERIVED** This is one scalar exchange channel between identical-mass constituents, not a complete calculation of the full S-matrix. The direct sigma contribution and its interference are included by the 2x2 contraction. Spin-2 exchange, crossed channels, and regular local terms have not been assembled here. The source derivation and its mapping to the constrained theory require review.

## No other positive Euclidean zero in the remaining denominator

**DERIVED** For this minimal model a=(D-2)/4, set y=z/M^2, r=(D-4)/2 in (0,1), and H(y)=integral_0^1 [1+x(1-x)y]^r dx. Then

    F/g = 1-(1-a y) H(y).

**DERIVED** Concavity gives 0<H(y)<=1+r*y/6. If 0<y<1/a,

    (1-a y) H(y) <= (1-a y)(1+r*y/6) < 1,

**DERIVED** because a>r/6. At y>=1/a, (1-a y)<=0, so F/g>=1. Therefore F(z)>0 for every real z>0. Only z=0 remains as a zero on the nonnegative Euclidean axis. This proof applies to this scalar denominator, not to the entire theory's spectrum, and does not address complex poles or timelike cuts.

## Numerical verification

**DERIVED** scalar_matter_audit.py uses Python standard-library Gamma functions and composite Simpson quadrature. It independently evaluates the original round-4 angular-averaged loop numerators, rather than invoking the saved SymPy expressions. It compares both A and B to the exact round-5 formulas at four dimensions D=4.01,4.2,5,5.5 and four momenta z/M^2=0.09,1,2.89,9. It also checks direct matrix inversion against the reduced matter exchange near the pole for lambda=0,0.05,0.5 at fixed M=1. Holding M fixed is a response-kernel check, not a self-consistent bare-mass gap scan. Some tested couplings near D=4 are outside a strict small-coupling expansion; they test algebra only.

**DERIVED** The first 800-panel run had maximum relative discrepancy 5.66e-11 across the A/B comparisons. Final machine-readable output additionally contains two-sided convergence toward the finite exchange limit and 400/800/1600-panel quadrature checks. The full original SymPy scripts were not run: required packages are absent, and network access blocked installation.

## What remains unresolved

1. **DERIVED** Independently verify the on-shell source, relative metric/sigma sign, and external-leg prescription.
2. **DERIVED** Determine whether an uncoupled scalar pole still produces an instability in gauge-invariant vacuum observables or other matter channels. Absence from this channel is not absence from the theory.
3. **DERIVED** Assemble full leading-N matter scattering with consistent normalization and channel counting.
4. **DERIVED** Resolve constrained-measure equivalence before transferring the auxiliary result to CEV.
5. **DERIVED** Test unequal masses, nonminimal curvature coupling, regulators, and 1/N corrections separately.

**DERIVED** The earlier conclusion of a dangerous pole for a generic trace source remains a property of that source calculation. It cannot automatically be identified with a pole in the on-shell constituent amplitude.

## Inputs and primary context

**DERIVED** Local inputs: C:/workspace/physics/external_reviews/grok_quartic_round4/{common.py,r4_scalar.py,r4_poles.py}; grok_quartic_round5/{DERIVATION.md,r5_residue.py}. These are project research notes, not externally validated results.

**UNVERIFIED** Carone, Claringbold and Vaman, *Composite graviton self-interactions in a model of emergent gravity*, arXiv:1710.09367v2, https://arxiv.org/abs/1710.09367 : defines the constrained scalar model and extracts the graviton from scalar scattering. Section II distinguishes the constrained path integral from an unrestricted auxiliary-metric integral. It does not establish the cancellation claimed in this note.

**UNVERIFIED** Carone, *Composite Gravitons from Metric-Independent Quantum Field Theories*, arXiv:1911.13166, https://arxiv.org/abs/1911.13166 : primary background on tuning, clock-free variants, and the need for ultraviolet completion. It is not evidence that the present model is stable.

**DERIVED** No novelty claim is made.
