# Independent check of Grok round 2

**Conversation date range** $2026-10$ (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records independent check of Grok round 2. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

## Verdict

**DERIVED** The requested contact calculation in s2_contacts.py is correct for the stated Euclidean auxiliary-metric action. The stationary TT contact is $p/4$, not zero. The fixed-effective-mass bubble in s4_bubbles.py and its q-squared coefficient also check out. Under a conventional leading-large-N auxiliary-metric effective-action saddle, a conditional mass-replacement result follows.

**DERIVED** This does not yet establish a complete result for the original constrained CEV quantum theory. The source/measure/gauge equivalence and Lorentzian map are explicit assumptions or open items in DERIVATION.md. Those qualifications should remain attached to the claim and in DISAGREEMENT_LOG.md. The ordinary auxiliary-metric response and a fully normalized composite scattering amplitude must not be silently identified.

## What was inspected and executed

**DERIVED** Read DERIVATION.md, DISAGREEMENT_LOG.md, README.md, outputs.txt and all five scripts from external_reviews/grok_quartic_round2. Several mathematical glyphs in the copied Markdown are corrupted; equations below use ASCII or standard symbols reconstructed from the scripts and context. Original files were left unchanged.

**DERIVED** SymPy and mpmath are absent in the available Python runtime. An isolated dependency installation was attempted and failed. Therefore the five original scripts were NOT rerun, and their supplied outputs are not being presented as this reviewer's execution. Instead an independent Fraction/NumPy/SciPy implementation was executed; all 58 checks pass. Source hashes, results and that implementation are included.

## Contact derivation: confirmed

**DERIVED** Write $a=D/2-1$ and let p denote the isotropic kinetic expectation, not the dimension-dependent a. For fixed P_ab with $\langle P_{ab}\rangle=p\delta_{ab}$ and normalized TT h,

**DERIVED** $\sqrt{g}=1-t^2 tr(h^2)/4+O(t^3),$
**DERIVED** $g_inverse=I-th+t^2 h^2+O(t^3).$

**DERIVED** The quadratic coefficient of $\sqrt{g}[g_{\mathrm{inverse}}:P/2+\langle V\rangle]$ is

**DERIVED** $C_TT=[\frac{p}{2}-\frac{Dp}{8}-<V>/4]tr(h^2).$

**DERIVED** The trace tadpole condition gives $<V\ge -a p$, hence $C_TT=(p/4)tr(h^{2})$. A covariant insertion -c g_ab inside P acts as a potential shift $-Dc/2$; the same substitution leaves the contact $p/4$. This is exactly the algebra tested in s2. Independent exact rational checks and direct matrix finite differences reproduce it for $D=3,4,5,6$.

**DERIVED** Scope of the script: its random H is traceless and its P expectation isotropic. $L1=0$ therefore checks a vanishing TT first variation, not the full off-shell stress identity for arbitrary anisotropic P and traceful H. The separate trace variation does correctly supply the stated background equation. This coverage limit does not invalidate the contact result.

## Bubble and normalization: confirmed within the specified Gaussian saddle

**DERIVED** For $\chi =m_{eff}^{2}>0$,

**DERIVED** $p_\phi =N \Gamma (-\frac{D}{2}) \chi ^(\frac{D}{2})/[2(4\pi )^(\frac{D}{2})],$

**DERIVED** $B2(Q)=N \Gamma (-\frac{D}{2})/[2(4\pi )^(\frac{D}{2})] \int_0^1 [\chi +x(1-x)Q^2]^(\frac{D}{2}) dx.$

**DERIVED** The scalar off-diagonal Wick contraction has its factor 2N. It gives $B2(0)=p_\phi$. With $\gamma 2=(p-B2)/4$ and $p=p_\phi$, its small-Q expansion is

**DERIVED** $\gamma 2(Q)=M_P^(D-2) Q^\frac{2}{8}+O(Q^4),$
**DERIVED** $M_P^{D-2}=N \Gamma (1-D/2) \chi ^{D/2-1}/[6(4\pi )^{D/2}]$.

**DERIVED** Independent numerical Feynman-parameter integration checks this coefficient at $D=3.5,4.02,4.2,5.5$ and $\chi =0.1,1,10$. Positivity of this coefficient near $D=4$ requires approaching from above in this finite-dimensional prescription; it is not true for every $D>4$ across all Gamma-function sign intervals.

**DERIVED** With the definition $\Gamma _{quadratic}=\gamma 2\,t^{2}$ for a normalized shear, the Hessian eigenvalue is 2gamma2. The corresponding normalized metric-mode propagator is $1/(2\gamma 2)$, not literally $1/\gamma 2$. Two conventional metric-source couplings -h.$T/2$ then give Euclidean exchange magnitude $1/(8\gamma 2)$, which tends to $1/[M_P^{D-2}Q^{2}]$. This removes a factor-of-two ambiguity but does not determine the Lorentzian sign or prove equivalence to the full constrained composite amplitude.

## Why the conditional mass-replacement argument works

**DERIVED** For a conventional auxiliary-metric leading-N functional,

**DERIVED** $\Gamma [g,\sigma ]=N/2 \operatorname{Tr} \ln\left(-\Delta _g+m0^{2}+\sigma \right)$
 $+ \int \sqrt{g}[V0-c2-N \sigma ^{2}/(4\lambda )]$
 + specified clock/counterterm terms.

**DERIVED** Eliminating sigma at its stationary point gives the effective metric Hessian as a Schur complement:

**DERIVED** $\Gamma _{hh},eff=\Gamma _{hh}-\Gamma _{hsigma} (\Gamma _{sigmasigma})^{-1} \Gamma _{sigmah}$.

**DERIVED** In a flat isotropic vacuum the mixed tensor is $a \delta _{ab}+b Q_a Q_b$, so its TT projection vanishes. At a regular massive saddle, the TT block is therefore the fixed-sigma Gaussian determinant block plus all local contacts. Enforcing the background condition gives the contact above. This is a functional explanation for the conditional mass substitution, not just a count of a few diagrams. The auxiliary potential has been included.

**DERIVED** Qualifications: the Hubbard-Stratonovich expression is an extremization identity unless its quantum integration contour and normalization are specified. For positive Euclidean quartic coupling, the displayed negative sigma-squared auxiliary term does not itself define a convergent real-sigma probability measure. It must not be mistaken for a physical ghost or a completed Monte Carlo measure.

## What is still missing to identify this with constrained CEV

1. **CITED** Quantum definition/measure. DERIVATION.md replaces the constrained composite construction by $-\ln\int\mathcal D\phi\,\exp(-S_E[g,\phi])$ at external g. This is a legitimate matter generating functional, but the exact map to the original constrained functional integral has not been established. Carone 1911.13166, footnote 4, explicitly distinguishes imposing the stress constraint from merely integrating over an auxiliary metric because a functional determinant can differ. Show that the relevant difference is absent or subleading at the stated order, or state the result as belonging to the auxiliary-metric model only. — eq./page UNVERIFIED.
2. **DERIVED** Euclidean sign map. From the displayed positive Euclidean action and $T=0$, algebraic elimination gives $g_{ab}=-a P_{ab}/V$. At $D=4$ its on-shell density is $-\sqrt{\det P}/V$ on a real positive-metric branch, not automatically the positive determinant density written by analogy with the Lorentzian action. A defined Wick rotation, clock-field contour and potential convention must connect them. This is not a proof that no continuation exists; it is an uncompleted map already acknowledged in the supplied derivation.
3. **CITED** Carone, Claringbold and Vaman: Counterterm/gauge bookkeeping. The original CEV gauge in 1710.09367 Eq. 2.7 has the clock gradient squared $V0/a-c1$. The new text also puts an explicit -c1 delta_ab inside P while writing $X=b x$. Define whether $b^{2}$ is the unshifted coefficient or the already-renormalized clock gradient. Otherwise the same subtraction can be counted twice. The statement $b^{2}=c1$ cannot be identified literally with the original physical clock gradient without this definition; the original tuned gradient vanishes.
4. **DERIVED** Complete source map. The Euclidean normalization magnitude was checked above, but the constrained composite-source identification and Lorentzian continuation remain required for an actual physical scattering residue. Symbolic replacement of m in a free bubble does not by itself supply this map.

## Gap equation: reproduced roots, controlled interpretation

**DERIVED** For $D=4.2, \lambda =0.5, m0^{2}=1$ the two roots are 1.0245266959059707 and 1.657164442808966e16. The relative mass corrections $\lambda G/\chi$ are 0.0239395381 and 1.0 respectively. The remote root is therefore outside a strict small-correction expansion. Its existence as a root of the continued equation is established; interpreting it as a physical vacuum or proving it an artifact of a UV completion is not.

**DERIVED** The safe small-coupling test is |$\lambda G(\chi )/\chi$|$<<1$, with the regulator scale and dimensional coupling defined. The shorthand $\lambda <<8\pi ^{2}$|epsilon| assumes a particular near-four-dimensional scale choice and should not replace that dimensionless test at arbitrary chi.

**DERIVED** A negative quartic coupling gives an unbounded bare potential without additional stabilizing physics. A unique positive gap root in that case does not demonstrate a stable physical vacuum.

## PV toy: uniqueness can be proved, not merely sampled

**DERIVED** For the supplied common-shift offsets $t_i=(0,4,9,25)$, the weights are exactly

**FITTED** $c_i=(1,-\frac{15}{7},\frac{5}{4},-\frac{3}{28}).$

**DERIVED** They satisfy $\sum c_i t_i^k=0$ for $k=0,1,2$. Combining denominators gives the exact identity

**DERIVED** $\sum _i c_i/(u+\chi +t_i)$
**DERIVED** $=900/[(u+\chi )(u+\chi +4)(u+\chi +9)(u+\chi +25)]$.

**DERIVED** For $u\ge 0$ and $\chi >0$ the right-hand side is positive and strictly decreasing in chi. Consequently the regulated tadpole

**DERIVED** $G_PV(\chi )=\int_0^\infty du u\,900/[16\pi ^2 product_i(u+\chi +t_i)]$

**DERIVED** is positive and $G_{\mathrm{PV}}'(\chi)<0$ for every $\chi >0$. Thus $F'=1-\lambda G_{\mathrm{PV}}'>0$ for $\lambda >0$. There is at most one positive symmetric gap solution; for $m0^{2}>0$ the endpoint signs imply exactly one. The rational identity was checked with exact Fraction polynomial arithmetic, and the integral was independently compared with s5's logarithmic formula.

**DERIVED** This proves a statement about this particular common-shift toy. It is not proof for every PV spectrum, a completed interacting PV extension of constrained CEV, or the absence of condensate solutions.

## Condensate limit: correct a nonuniform limit

**DERIVED** The derivative expansion in $Q^{2}/\chi$ assumes $\chi >0$. At $\chi =0, B2(Q)$ does NOT vanish at nonzero Q. Instead,

**DERIVED** $B2(Q)=N \Gamma (-\frac{D}{2})/[2(4\pi )^(\frac{D}{2})]$
 * $(Q^{2})^{D/2} * \operatorname{B} (D/2+1,D/2+1)$.

**DERIVED** This was independently checked by quadrature. The local $Q^{2}$ coefficient inferred from the massive dimensional formula tends to zero, but the massless nonlocal response remains. Scalar radial modes and other regulator scales also require separate treatment. A vanishing massive-formula M_P coefficient in this prescription is not a general proof that gravitational response disappears on every condensate branch.

## Recommended status and next task

- **DERIVED** Contact algebra: independently confirmed.
- **DERIVED** Fixed-effective-mass Gaussian bubble and Euclidean $low-Q$ coefficient: independently confirmed.
- **DERIVED** Mass replacement for the stated isotropic leading-N auxiliary-metric response: conditionally derived and independently checked.
- **DERIVED** Exact same result for the fully defined constrained CEV theory: not yet verified; retain the source/measure/sign/gauge qualifications.
- **DERIVED** Remote dimensional root: reproduced but uncontrolled in the strict expansion.
- **DERIVED** Common-shift PV toy positive-coupling uniqueness: proved for the stated offsets/weights.
- **UNVERIFIED** Spin-zero spectrum, physical unitarity, finite-N effects and observational meaning: unresolved.

**DERIVED** Finish the quantum measure, clock-counterterm and Lorentzian/source normalization map before moving to physical interpretations of gap branches. If the map is established at the stated order, the conditional result can be upgraded. If it is a different auxiliary-metric model, label it that way and assess that model on its own merits.

**CITED** Primary papers consulted: Carone, Claringbold and Vaman, arXiv:1710.09367 (Eqs. 2.7 and 2.10-2.14); Carone, arXiv:1911.13166 (footnote 4); Batz, Erlich and Mrini, arXiv:2011.06541 (background and response/contact distinctions). No novelty claim or external message sent.
