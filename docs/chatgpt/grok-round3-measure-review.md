# Independent review of round-3 measure step

**Conversation date range** $2026-10$ (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records independent review of round-3 measure step. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

## Verdict

**DERIVED** The bare metric dependence is algebraic. On a fixed nonsingular branch with an explicitly chosen ultralocal measure, its local constraint Jacobian and auxiliary-metric fluctuation factors can be treated as scaleless contact loops and dropped in the usual perturbative dimensional-regularization prescription. This is a valid formal argument for that local sector; it is not a numerical consequence of SymPy.

**DERIVED** However, t4_measure.py does NOT establish that this is the complete gauge-fixed constrained measure at the massless tuning. In particular, the supplied clock gauge becomes singular exactly at $g_c=0$. Thus the upgrade should remain a conditional perturbative equivalence proposal, not an independently established result for the fully defined constrained CEV theory. A clock-free formulation or a controlled nonsingular-gauge limiting procedure is needed.

## Files and execution

**DERIVED** Read DERIVATION.md, all four scripts and outputs.txt in external_reviews/grok_quartic_round3. The ASCII equations are now readable. Original files were not modified. The original SymPy scripts could not be rerun because SymPy/mpmath are absent in this runtime; this was already established in the round-2 review. Their supplied outputs are not presented as this reviewer's execution.

**DERIVED** An independent NumPy implementation evaluates the FULL symmetric-metric constraint Jacobian, including off-diagonal components, and checks it against direct finite differences. It also checks the local action Hessian and the clock-gauge Jacobian. All 160 implementation checks pass. These are finite-dimensional local identities, not 160 quantum-equivalence tests. Dimensional regularization's scaleless-integral assignment is not something these numerical checks evaluate.

## What t4 actually checks

**DERIVED** Its first assertion uses $a D=4$ diagonal metric ansatz and verifies that the pointwise Lagrangian Hessian has no symbolic derivatives. This confirms algebraic dependence in that ansatz. Its second assertion checks that one bubble coefficient divided by N is independent of N.

**CITED** It does not calculate the full constraint determinant, its invertibility, the actual constrained surface measure, Faddeev-Popov factors, contour phases, or the transformation of source insertions. The printed claim of equivalence is an interpretation added after the assertions; it is not itself an assertion tested by the script. — eq./page UNVERIFIED.

## Full local metric Jacobian: independent calculation

**ASSUMED** Use Euclidean $S=\int\sqrt{g}\,[g^{-1}:P/2+V]$, $T=P-g[g^{-1}:P/2+V]$. At a stationary isotropic point $g=I, P=\rho I, V=-(D/2-1)\rho$, define $a=D/2-1$ and $n=D(D+1)/2$.

**DERIVED** For an arbitrary symmetric metric variation h at fixed P,V:

**DERIVED** $\delta T=-\rho h+(\rho /2) I tr(h)$.

**DERIVED** The constraint Jacobian J has eigenvalue -rho in the $n-1$ traceless directions and a*rho in the trace direction. Therefore

**DERIVED** $\det J=(-\rho )^{n-1} a\,\rho$.

**DERIVED** The pointwise action Hessian $H=-J/2$ has eigenvalues $\rho /2$ and $-a\,\rho /2$, with

**DERIVED** $\det H=-a\,(\rho /2)^n$.

**DERIVED** Both are nonzero if rho!$=0$ and D!$=2$. This was checked for $D=3,4,5,6$ with both signs of rho, using an orthonormal basis of every symmetric metric component. Trace and shear action Hessians were independently evaluated by finite differences. This strengthens the original diagonal-only check.

**DERIVED** The local factors remain branch/phase dependent for indefinite Hessians. Absolute Jacobians for real constraints and stationary-phase determinants are not automatically the same object.

## The DR statement that is justified

**DERIVED** If an operator really is multiplication by a nonsingular matrix M(x), then, up to a field-independent normalization,

**DERIVED** $\operatorname{Tr} \ln M = \delta ^{D}(0) \int d^Dx tr \ln M(x)$.

**DERIVED** Here $\delta ^{D}(0)=\int d^Dk/(2\pi )^D 1$ is scaleless. The conventional perturbative DR prescription sets it to zero. M(x) may depend on phi and partial phi; that does not put derivatives on the metric delta-function kernel and does not itself invalidate this argument.

**DERIVED** It is log J that vanishes in this prescription, making the regular determinant factor one. The Jacobian determinant is NOT set to zero. A finite, nonzero epsilon does not by itself prevent assigning zero to a scaleless integral: the issue is whether that prescription is consistently part of the definition used for every sector.

**CITED** This is analogous to the familiar removal of Jacobian ghost loops for invertible local perturbative field redefinitions; see Criado and Perez-Victoria, arXiv:1811.09413, section 2. Their equivalence result also keeps track of transformed sources and has invertibility assumptions. It cannot simply be invoked as a theorem for singular gauge fixing or an unspecified constrained surface measure. — eq./page UNVERIFIED.

## Concrete remaining obstruction: the tuned clock gauge

**DERIVED** For the static gauge $F^I=X^I-\sqrt{g_c}x^I$, the infinitesimal diffeomorphism variation on the gauge surface is

**DERIVED** $\delta F^I=+/-\sqrt{g_c} \delta ^I_\mu \xi ^\mu$.

**CITED** The pointwise Faddeev-Popov matrix therefore has determinant proportional to $g_c^{D/2}$. At $g_c=0$ the matrix is zero, not invertible. This is the very tuning used to infer a massless pole. A formula involving ln M or $1/\det M$ on a regular branch cannot be extended to that singular point merely by multiplying it by a DR contact zero. — eq./page UNVERIFIED.

**CITED** Carone 1911.13166, section III.A, explicitly notes that the clock profile cannot implement this gauge fixing when it vanishes. The paper discusses a limiting construction and its unwanted-state issue, and section IV gives a clock-free reformulation. This is not a new no-go argument. It is a known issue that the round-3 equivalence step must accommodate rather than skip. — eq./page UNVERIFIED.

**DERIVED** Possible ways forward:

1. **DERIVED** Work in the clock-free formulation with an explicitly specified regular gauge and constraint measure. Identify common gauge factors and show their treatment is the same in both formulations at the claimed order.
2. **DERIVED** Define the calculation at nonzero g_c, then justify the tuned limit for physical observables, including the scalar sector and sources. Setting DR contact loops to zero before taking a limit is a prescription, not proof that the limiting theory is unitary or gauge equivalent to the exact endpoint.

## Other measure qualifications

- **DERIVED** The notation $\int _{T=0}$ does not, by itself, specify whether the measure is a flat delta-functional, an induced surface measure or a canonical constrained measure. The inverse $\det (\delta T/\delta g)$ follows from a chosen flat-delta prescription with nonsingular roots. State that choice explicitly and match it to the theory being claimed.
- **DERIVED** Multiple stationary branches and contour changes are not regular local determinants. The argument is perturbative around one chosen nonsingular branch; the all-orders claim should be restricted to regular ultralocal perturbative factors on that branch.
- **DERIVED** Before integrating phi, the bare g integral is pointwise. After integrating phi, Gamma_gg contains momentum-dependent matter response. Its determinant is NOT generally ultralocal and cannot all be discarded as delta(0). Its metric loops may be subleading at large N in a regular, consistently gauge-fixed saddle expansion, but that is a different counting argument. $\Gamma =N \gamma$ supports that counting; the t4 N-scaling assertion alone does not verify all its conditions.
- **DERIVED** At fixed regulator and regular large N, some local measure terms may themselves be subleading. Conversely regulator/continuum limits and singular saddles require care. The present checks do not decide the complete PV/cutoff theory either way.
- **DERIVED** A negative dimensionally continued/subtracted expectation of a positive classical quadratic operator permits a formal real background equation. It does not construct a convergent positive real Euclidean sampling measure. The determinant branch, clock/HS contours and Lorentzian continuation remain relevant.

## Status of the other round-3 steps

**DERIVED** The corrected clock dictionary removes the previously identified double subtraction algebraically, given its Euclidean continuation convention. The equations $c1=-p_\phi , g_c=V0/a-c1$ imply $p=p_\phi -g_c=-V0/a$. The stationarity equation then makes c2 the full leading-N potential expectation, and the bubble condition makes $g_c=0$. Those algebraic substitutions are consistent.

**DERIVED** Their interpretation still depends on the gauge and continuation points above. Carrying the overall scattering sign over from the free CEV paper is a convention match, not an independent Lorentzian derivation. A full source-residue conclusion should retain that qualification.

## Recommended conclusion for the disagreement log

**DERIVED** Local ultralocal-Jacobian argument: ACCEPTED on a nonsingular, specified perturbative branch with the usual DR prescription.

**UNVERIFIED** Claim that t4 establishes complete constrained/auxiliary equivalence at the massless clock tuning: NOT ESTABLISHED. The clock-gauge Jacobian is singular there, and the complete constraint/gauge measure has not been supplied.

**DERIVED** The conditional leading-N auxiliary-metric mass-replacement result remains supported. Extending it to constrained CEV should be phrased as a proposal conditional on a completed regular-gauge/measure prescription, not simply conditional on $\delta (0)=0$ alone.

## Single next task

**DERIVED** Specify the clock-free constrained generating functional and its regular gauge fixing, then compare its regulated measure and physical sources with the auxiliary-metric functional at leading N. Separate the demonstrated regular local determinant from any gauge, branch and source factors. That would directly address the missing premise without repeating the already checked massive bubble integral.

**DERIVED** Primary sources:

- **CITED** Carone, Composite Gravitons from Metric-Independent Quantum Field Theories, https://arxiv.org/abs/1911.13166, section III.A, footnote 4, section IV. — eq./page UNVERIFIED.
- **CITED** Criado and Perez-Victoria, Field redefinitions in effective theories at higher orders, https://arxiv.org/abs/1811.09413, section 2. — eq./page UNVERIFIED.

**DERIVED** No quantum equivalence or originality claim. No message was sent to Grok. Original input hashes are recorded in results.json.
