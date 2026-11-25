# Scalar threshold crossing and continuation

**Conversation date range** 2026-10 (exact conversation dates UNVERIFIED)

**DERIVED** This checkpoint records scalar threshold crossing and continuation. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Fixed free leading-N PV auxiliary-metric model, constituent M=1, no bare curvature terms. Not constrained CEV or a unitary ultraviolet completion.

**DERIVED** The previous presence/absence of subthreshold scalar poles has a concrete mechanism in the tested family: a pole reaches the physical pair threshold s=4 at Lc=14.557206455684199. A 301-point L scan on [1.05,150] found one threshold sign change; bisection refined it. This is not a proof of global uniqueness.

**DERIVED** The physical-species log near s=4-delta obeys

    J=-2+(pi/2)sqrt(delta)+O(delta).

**DERIVED** Since its scalar vertex factor (1+s/2)^2 tends to 9, the full kernel has

    A(s,L)=A(4,L)+k sqrt(delta)+O(delta), k=9 C pi/8,
    C=(4 pi)^(-2).

**DERIVED** Thus, on the side A(4,L)<0, the nearby pole has binding gap delta approximately [A(4,L)/k]^2. Its connected identical-source residue tends to 4.5 sqrt(delta)/k and vanishes at threshold. These are derived threshold expansions, numerically checked against the full log kernel.

**DERIVED** At relative offset (L-Lc)/Lc=-0.0003 the exact gap is 2.193898636e-7; the gap/prediction ratio is 1.000139567 and the residue/prediction ratio is 1.000139555. No finite-difference estimate at the branch point itself is used.

**DERIVED** Continuing across only the physical constituent cut gives the adjacent-sheet below-threshold expression

    A_II(s,L)=A_I(s,L) - (C/4)(1+s/2)^2 2 pi sqrt(4/s-1).

**DERIVED** For tested positive offsets +0.003,+0.01,+0.03,+0.1, the former first-sheet pole is replaced by a real second-sheet pole below threshold. Example: at +0.01, s_II=3.9997566848986335. On the negative-offset side, the tested nearby roots are on the first sheet instead. This supports a bound-state/virtual-state style threshold trajectory of the response, not a new spacelike instability at the crossing.

**DERIVED** The scattering-theory classification is standard: see the primary research article On the nature of near-threshold bound and virtual states, European Physical Journal A (2021), DOI https://doi.org/10.1140/epja/s10050-021-00413-y . Its bound/virtual pole distinction informs terminology only; it is not evidence that this PV gravity model has a physical hadronic-style unitary S-matrix.

**DERIVED** Limitations: double precision; finite parameter/root brackets; no interval-certified uniqueness; no complex-pole search on the second sheet; no coefficient matching across L; signed regulator fields remain unphysical ultraviolet data; no proof of the microscopic Lorentzian contour or constrained quantum equivalence. A response pole should not be called an observed particle. The old DR spacelike scalar zero is a different calculation and has not been resolved by this threshold result.

**DERIVED** Reproduction: scalar_threshold_transition.py imports spec and C from complex_scalar_audit.py; that file imports the Simpson helper from scalar_matter_audit.py. All dependencies are Python standard library. Results include threshold roots, gaps, residues, and sheet assignments. The archive includes these helpers.

**DERIVED** Next: test the local sheet continuation independently against analytic discontinuity conventions, then examine low-energy matching before attributing the transition to a physical change in one theory.
