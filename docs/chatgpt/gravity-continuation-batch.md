# Gravity continuation: source overlap, complex contours, dispersion

**Conversation date range** 2026-10-29 to 2026-10-30

**ASSUMED** This checkpoint records Gravity continuation: source overlap, complex contours, dispersion. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** October 29, 2026, America/Chicago. Conditional theoretical calculations; no claim of unification or novelty.

## Main conclusions

1. **DERIVED** The dimensionally continued auxiliary-metric scalar pole cancels from identical constituent scattering even with the saved quartic/sigma mixing, but not from general sources. Ninety-six unequal-mass/sigma-coupling cases reproduce the analytic residue, and twelve curvature-response cases retain the pole.
2. **DERIVED** For five existing four-dimensional Pauli-Villars (PV) spectra, twenty finite first-sheet upper-half-plane contours give zero winding for the scalar and spin-2 kernels with their known massless factors removed. Three sampling refinements and adaptive phase subdivision agree. This is numerical evidence in explicit bounded regions, not a certified global zero-exclusion theorem.
3. **DERIVED** Twenty three-subtracted scalar inverse-kernel dispersion reconstructions agree with direct evaluation to maximum absolute error $4.613\times10^{-14}$. This checks the full signed PV continuum, not positivity of a physical ultraviolet completion or a propagator dispersion relation.
4. **DERIVED** The next model-defining question is whether independent geometric dynamics may be introduced, or whether the starting theory must remain matter-only. No modification of the starting action was made in this batch.

## Recovery of prior work

**UNVERIFIED** The continuation recovered existing scalar_source_projection_check, metric_tensor_ward_check and metric_spectral_check reports. They already contained the free identical-constituent cancellation, full free-theory Ward checks, and positive physical-only continuum signs. The preceding scalar_matter_audit is an independent reproduction plus a quartic-mixing extension; the free cancellation itself should not be counted as a newly discovered result. Its added analytic Euclidean denominator argument remains a conditional derivation.

## A. General source overlap

**ASSUMED** Use the saved round-5 DR blocks, with $M^{2}=1, a=(D-2)/4, U=1-a z, F=Pi U+g$:

$$
\begin{gathered}
A=-UF/4, B=-F/4, C=-(Pi+1/\lambda )/4, \\
\det K=F(U/\lambda -g)/16, \\
z_{star}=(1-\lambda *g)/a.
\end{gathered}
$$

**DERIVED** Near z_star the inverse matrix has residue $4 vv^T/(a g), v=(1,-\lambda *g)$, so the connected source quadratic functional $W=J^T K^{-1} J/4$ has residue $(J\cdot v)^2/(ag)$.

**DERIVED** For an O(1) minimally coupled probe of mass squared $\mu ^{2}$ and sigma coupling $\kappa , J=(\mu ^{2}-a z,\kappa )$. Its residue is

$$
[\mu ^{2}-1+(1-\kappa )\lambda *g]^{2}/(a g).
$$

**DERIVED** For $\kappa =1$ this is exactly $(\mu ^{2}-1)^{2}/(a g)$: equal masses decouple; unequal masses do not. $\kappa =1$ assumes an extended quartic sharing sigma with the probe. $\kappa =0$ is a spectator without this quartic. These are explicit probe-model assumptions, not an assertion that both are the same microscopic theory. Probes do not modify the leading-N determinant; O(N) added species would require a new calculation.

**DERIVED** The linearized scalar curvature of $h_{ab}=f \delta _{ab}$ is $R1=(D-1) z f$ in these Euclidean conventions. R1 is invariant under linearized diffeomorphisms about flat space: the double-divergence and Laplacian-trace gauge variations cancel. The connected f covariance is $K^{-1}_{ff}/2$, so

$$
Res <R1 R1> = 2 (D-1)^{2} z_{star}^{2}/(a g) > 0.
$$

**DERIVED** All twelve tested DR parameter choices have $z_{star}>0$. Thus the pole remains in a gauge-invariant linearized curvature response of this auxiliary model. Under the earlier assumed mostly-plus Lorentzian continuation, a positive Euclidean pole means $\omega ^{2}=k^{2}-z_{star}$ and exponential growth at sufficiently small k. This conclusion remains conditional on that Lorentzian contour and on the auxiliary theory being the physical model. It does not transfer automatically to constrained CEV.

**DERIVED** Verification: $D=4.01,4.2,5,5.5; \lambda =0,0.05,0.5$; probe $\mu ^{2}=0,0.25,1,2; \kappa =0,1$. Direct two-sided matrix inversions at $z_{star} +/-1\times10^{-7}$ agree with the formula. Some choices near $D=4$ are not in a strict small-quartic regime; this is a fixed-Hartree-mass algebraic response check, not a controlled full bare-theory/gap-saddle scan.

## B. First-sheet complex search

**DERIVED** Retain the existing PV family $m_i^{2}=(1,L,2L,5L+1), c_0=1$, with three signed regulator coefficients chosen to cancel moments $\sum c_i(m_i^{2})^k=0, k=0,1,2$. No bare curvature terms are added. $C0=(4 \pi )^{-2}$.

**DERIVED** For complex $q=Q_E^{2}$, principal analytic continuation of the Feynman-parameter log gives

$$
\begin{gathered}
J(q,m2)=\int _0^{1} \log (m2+q x(1-x)) dx \\
=\log\left(m2\right)-2+b \log ((b+1)/(b-1)), b=\sqrt{1+4m2/q}.
\end{gathered}
$$

**DERIVED** The scalar kernel is

$$
\begin{gathered}
A(q)=C0/4 \sum _i c_i[(m_i^{2}-q/2)^{2} J(q,m_i^{2}) \\
-(m_i^{2}-q/2)m_i^{2} \log\left(m_i^{2}\right)].
\end{gathered}
$$

**DERIVED** The transverse spin-2 coefficient is

$$
\begin{gathered}
K2(q)=C0/8 \sum _i c_i[\int _0^{1} (m_i^{2}+q x(1-x))^{2} \\
\log (m_i^{2}+q x(1-x)) dx \\
-(m_i^{2})^{2} \log\left(m_i^{2}\right)].
\end{gathered}
$$

**DERIVED** The weighted integral is computed analytically by elementary recurrences, with a convergent local series where cancellation would otherwise reduce precision. Both log formulas are independently checked against 1600-panel Simpson quadrature at twelve complex momentum/mass combinations. Maximum unweighted discrepancy is $3.293\times10^{-11}$, weighted discrepancy $3.883\times10^{-11}$. Conjugation of the log branch is checked.

**DERIVED** For each $L=3,5,10,30,100$, use $R=3.92 L$ and rectangles

$$
\begin{gathered}
-R \le \Re q \le R, \epsilon \le \Im q \le R, \\
\epsilon =0.02 and 0.002.
\end{gathered}
$$

**DERIVED** Count winding of $A(q)/q$ and $K2(q)/q$; dividing removes the known massless zero. Initial 400,800,1600 points per edge are adaptively refined using endpoint/midpoint phase increments until each accepted step is below 0.12 radians. All twenty contours have winding zero at all refinements. A positive-control function with one known enclosed off-axis zero has winding one.

**DERIVED** Limitations: no interval arithmetic or certified bounds between samples; no exclusion in $0<\Im q<0.002$, on the real axis, outside the rectangles, or on other sheets. Reflection gives the corresponding lower-half domains only under the verified conjugate analytic prescription. Separate real-axis analyses contain known subthreshold scalar poles, which these contours intentionally do not enclose. Rectangle corners can have modulus above the first regulator threshold; bounds are on real and imaginary coordinates, not a circular low-energy cutoff. The result says nothing about interacting complex poles or constrained-measure corrections.

## C. Dispersion reconstruction

**DERIVED** At $q=-s-i0$,

$$
\begin{gathered}
\rho _A(s)=\Im A(-s-i0)/\pi \\
=-C0/4 \sum _{open} c_i(m_i^{2}+s/2)^{2} \sqrt{1-4m_i^{2}/s}.
\end{gathered}
$$

**DERIVED** The three-subtracted reconstruction is

$$
\begin{gathered}
A(q)=a_1 q+a2 q^{2}-q^{3} \int _4^\infty \rho _A(s)/[s^{3}(s+q)] ds, \\
a_1=-(C0/8) \sum c_i m_i^{2} \log\left(m_i^{2}\right), \\
a2=(C0/16) \sum c_i \log\left(m_i^{2}\right).
\end{gathered}
$$

**ASSUMED** Use $u=4/s$ and split at every pair threshold $u=1/m_i^{2}$. A sine-squared substitution on each interval smooths the square-root endpoints. Simpson $400/800/1600$ panels are compared. At $q=0.1,1,2,10$ for each of five spectra, the largest final discrepancy from direct analytic A is $4.613\times10^{-14}$.

**FITTED** The full regulator-open continuum is included. Its signed weights are mathematical PV data, not positive physical matter. The subtraction coefficients are fixed analytically, rather than fitted to the evaluation points. This validates analyticity/discontinuity consistency of A, not a full dispersion reconstruction of $1/A$ with all poles.

## D. Why a physical modeling choice is approaching

**DERIVED** The DR auxiliary model has a surviving curvature-response concern. The tested PV models have better bounded diagnostics, but signed regulator determinants are not a physical ultraviolet theory. The original constrained model is a third, not yet proven equivalent, quantum definition.

**ASSUMED** One cannot infer a unique high-energy spectrum from the Einstein coefficient alone. Even matching the curvature-squared coefficient leaves higher-curvature data unspecified: adding $\alpha *q^{3}$ to a scalar inverse kernel leaves its q and $q^{2}$ coefficients unchanged, $yet \alpha =-A(q0)/q0^{3}$ puts a zero at any selected positive q0. This is a mathematical illustration of missing information, not a proposed fix or an action modification implemented here. A covariant operator such as R box R contributes at that order in the scalar quadratic response, with correlated tensor structure that must be derived from its action.

**ASSUMED** There remain routine calculations within each fixed model, but they cannot decide which microscopic theory is intended. A proposed stability repair must state whether independent curvature operators/geometric degrees of freedom are permitted. If the goal is strictly matter-only emergence, the next priority should be a specified constrained path integral, regulator, contour, and measure, followed by deriving its actual Hessian. If independent geometric dynamics are permitted, the task changes to a matched effective-theory or specified ultraviolet-completion study. Neither option has been silently selected here.

## Files and reproduction

**DERIVED** All scripts use Python's standard library. Run from the same directory:

```bash
python complex_scalar_audit.py
python scalar_dispersion_audit.py
python scalar_observable_audit.py
```

**DERIVED** complex_scalar_audit and scalar_observable_audit import scalar_matter_audit. scalar_dispersion_audit imports both helpers. JSON results contain the parameter choices and individual comparisons.

**DERIVED** Primary physical context, not evidence for the original calculations in this batch:

- **CITED** Carone, Claringbold, Vaman, arXiv:1710.09367v2, https://arxiv.org/abs/1710.09367 . Constrained scalar model and matter scattering. — eq./page UNVERIFIED.
- **CITED** Carone, arXiv:1911.13166, https://arxiv.org/abs/1911.13166 . Tuning, clock-free variants, ultraviolet-completion caveats. — eq./page UNVERIFIED.
- **CITED** Li and Vaman, arXiv:2512.06132v1, https://arxiv.org/abs/2512.06132 . PV-regulated composite model and low-momentum graviton pole. Its microscopic regulator definition should not be identified automatically with the signed auxiliary determinant used here. — eq./page UNVERIFIED.

**DERIVED** Status: several additional consistency tests completed. No demonstration of global stability, full unitarity, a realistic matter spectrum, or quantum/Einstein unification.
