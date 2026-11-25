# Pauli-Villars scalar-response control

**Conversation date range** $2026-10$ (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records pauli-Villars scalar-response control. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** This calculation concerns the free, minimally coupled auxiliary-metric model at leading N, per unit N, with $\lambda =0$. It does not establish constrained-CEV equivalence or a unitary ultraviolet completion. It compares regulator prescriptions; it is not a matched scheme transformation of a fixed physical theory.

## Result

**DERIVED** For one fixed signed Pauli-Villars determinant spectrum, the positive-Euclidean-momentum root at $Q^{2}=2 M^{2}$ in the single-mass finite-dimensional prescription is absent at that point. No spacelike zero was found on 500 sampled points from $Q^{2}/M^{2}=0.001$ to 19.9. This is a scan, not a proof of absence.

**DERIVED** The real analytic continuation below the lightest pair-production cut instead has a simple timelike root at $s/M^{2}=2.77851951116$, with $A'(Q^2)=0.00376237642604$ at $Q^{2}=-s$. In the specified real metric-source convention $W=J_f^{2}/(4A)$, its residue is $+66.4473650934 J_f^{2}$. Thus it is a candidate massive scalar with the ordinary Euclidean residue sign, not a demonstrated physical positive-norm particle. It lies below $s=4 M^{2}$.

**DERIVED** The low-energy coefficient $c=G_PV/12=0.00129934389282$ is positive ($M_P^{2}/N=2c$ in the inherited convention). The scalar slope is $-G_PV/8$, consistent with the Einstein scalar constraint term. The derivative was checked independently by symmetric differences.

## Specified regulator

**ASSUMED** Use $z_i=m_i^{2}/M^{2}=(1,5,10,26)$ and $c_i=(1,-15/7,5/4,-3/28)$. The moments $\sum c_i z_i^k$ vanish for $k=0,1,2$. These weights are a chosen determinant regulator, not a claim about physical negative-probability species. Every stress vertex uses its species' own mass.

## Four-dimensional finite response

**ASSUMED** Set $q=Q^{2}/M^{2}, t=x(1-x), C=(4 \pi )^{-2}, U_i=z_i-q/2$. All dimensionful values below use $M=1$.

**DERIVED** $A_PV(q) = \frac{C}{4} \sum _i c_i [ U_i^2 \int_0^1 dx \log (z_i+tq) - U_i z_i \log (z_i) ].$

**DERIVED** An arbitrary common logarithm scale drops out through the PV moments. The logarithmic expression follows from taking the PV sum before the $D \to 4$ limit of the supplied finite-D kernel. Terms polynomial in z_i of degree at most two, including common divergence coefficients, cancel. At $q=0$ the expression vanishes with the stationary-background tuning.

**DERIVED** $G_PV = C \sum _i c_i z_i \log  z_i; A'_PV(0)=-\frac{G_PV}{8}.$

**DERIVED** At the old root $q=2, A_PV=-0.00494520034493$, not zero. The individual masses have different U_i factors; a common single-species $U=0$ root cannot simply be factored out of their sum.

## Separate stress-bubble/contact construction

**DERIVED** For $h_{ab}=f \delta _{ab}$, at $D=4$ the minimal scalar trace vertex is $-k \cdot (k+q)-2z_i$. After shifting $k=l-xq$, write $d_i=z_i+tq$ and $b_i=tq-2z_i$. Angular averaging its square produces

**DERIVED** $L^{2} + [-2 b_i+(1-2x)^{2} q/4] L + b_i^{2}$,

**DERIVED** where $L=l^{2}$. The finite PV sums of the three loop moments have logarithmic pieces proportional to $-3 d_i^{2} \log d_i, +2 d_i \log d_i$, and -log d_i, respectively. Consequently

**DERIVED** $bubble = C \sum _i c_i \int dx [-3 d_i^2 + 2(-2b_i+(1-2x)^2 \frac{q}{4})d_i - b_i^2] \log  d_i.$

**DERIVED** The stationary conformal contact is $C2=-\sum _i c_i p_i=C/4 \sum _i c_i z_i^{2} \log z_i$. The direct diagram kernel C2-bubble/4 agrees with the finite-log expression to maximum absolute error $1.03\times10^{-13}$ across 200 spacelike and subthreshold timelike test momenta. This checks the free minimal quadratic kernel without relying solely on the Weyl-rescaling argument. It does not check interacting regulator fields, the constrained measure, or every Lorentzian source contact.

## Numerical verification

**DERIVED** The script uses 640-point Gauss-Legendre quadrature and compares selected near-threshold results with 1280 points (maximum absolute difference $1.57\times10^{-13}$). Quadrature nodes are cached. The finite-log formula is also compared with the signed finite-D sum at $D=4.01,4.001,4.0001$: errors decrease toward the finite-log answer at every specified test point. Double-precision finite-D sums are not used as the root oracle.

**DERIVED** Only the analytic domain $q>-4$ is evaluated. There is no on-cut continuation, second-sheet search, or broad complex pole search. On the real interval s in [0.001,3.999] one sampled sign-change interval yields the timelike root quoted above. There may be other unsampled or complex roots.

## Scope and next calculation

**FITTED** The regulator masses/weights, background subtraction, and no-bare-curvature condition are physical model choices at finite regulator scale. These results are not a claim that a subtraction convention can remove an observable pole in a fixed matched theory. The absence of the old spacelike root is encouraging only for this stated control.

**ASSUMED** Next: define how the quartic auxiliary field couples to the regulator species, derive all entries of the regulated metric-sigma block and contacts from that same action, and compute the physical conserved-source response under Lorentzian continuation. Replacing only one bubble by a PV sum would not be sufficient. Retain the free result here as a regression control. A positive Euclidean scalar residue alone does not prove unitarity, especially with finite signed regulator determinants.

**CITED** Reference: Chris Li and Diana Vaman, A model of composite gravity with Pauli-Villars regulators, arXiv:2512.06132v1, https://arxiv.org/html/2512.06132v1. Its regulator spin-2 calculation supplies context; it is not cited as proof of this scalar result. Round-5 local derivation supplied the finite-D comparison formula; the current script additionally constructs the $D=4$ minimal trace stress bubble and contact. — eq./page UNVERIFIED.
