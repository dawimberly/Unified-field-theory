# Response pole checkpoint — 30 October 2026

**Conversation date range** 2026-10-30

**ASSUMED** This checkpoint records Response pole checkpoint — 30 October 2026. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

## Scope

**UNVERIFIED** The calculations concern the free leading-N auxiliary-metric determinant, with real q and $q^{2}$ coefficients fixed to the $L=5$ Pauli–Villars reference before taking the heavy-regulator limit. $q=Q_E^{2}/m_{physical}^{2}$. This matching adds local polynomials and differs from the original zero-bare-curvature prescription. It is not established to be the constrained CEV quantum theory.

**DERIVED** The preceding scalar_matching_checkpoint_2026_10_30.md gives the matching prescription, independent threshold continuation, nine real-root integration controls, fixed-momentum decoupling and scalar asymptotic derivation. The present checkpoint adds analytic classification and exact arithmetic bounds.

## Scalar channel

**DERIVED** Let $C=(4 \pi )^{-2}, w(t)=(3-t^{2})^{2} t^{2}$, and $F0=A/q^{2}$. The limiting response satisfies

$$
F0(z)=a/z+b+(C z/32) \int _0^{1} w(t)/(4+z(1-t^{2})) dt.
$$

**DERIVED** For $a<0$,

$$
\Im F0(z)/\Im z = -a/|z|^{2} + (C/8) \int _0^{1} w(t)/|4+z(1-t^{2})|^{2} dt > 0.
$$

**DERIVED** This is a strict Pick/Herglotz property. It excludes nonreal first-sheet zeros. On each real analytic interval the derivative is positive. The positive interval has exactly one simple zero for every finite b. On (-4,0), a simple zero exists exactly when

$$
b < a/4 + 9C/40.
$$

**DERIVED** Equality places a zero at the branch-point boundary, where the simple analytic-zero argument does not apply. The cut lips at $q<-4$ have a nonzero discontinuity and are not extra analytic real intervals.

**DERIVED** With $a=-0.0019490158392333555$ and $b=-0.0003639828687336804$, the two nonzero first-sheet zeros are approximately

$$
\begin{gathered}
q_{positive} = 11.3066267794444, \\
q_{subthreshold} = -2.60868303692034.
\end{gathered}
$$

**DERIVED** The exact-rational script brackets them respectively in (11.30,11.32) and (-2.61,-2.60). It bounds logarithms by an atanh series and its analytic tail, and bounds the positive dispersion integral by exact cell integrals and endpoint denominators. Floating point is used only to display the resulting bounds.

**UNVERIFIED** The identical-source scalar trace factor $2-q$ is nonzero at the reference positive zero. Tuning b to 0.000816193570175525 moves that zero to $q=2$ and cancels that channel's trace-product pole; it does not eliminate the kernel zero. A conditional curvature numerator proportional to $q^{2}$ remains nonzero. The full constrained physical observable is still unresolved.

## Tensor channel

**DERIVED** For $F2=B/q^{2}, u=-a/3>0$ and $v=b/15$,

$$
F2(z)=u/z+v+(C z/120) \int _0^{1} t^{6}/(4+z(1-t^{2})) dt.
$$

**DERIVED** Unlike the scalar channel, the pole $u/z$ contributes a negative term to the Pick kernel. A positive continuum minus a rank-one term nevertheless gives a useful exclusion result: there can be at most one distinct upper-half-plane zero. If there were two, the divided-difference kernel at these zeros would vanish, forcing a strictly positive-definite continuum Gram matrix to equal a rank-one matrix, which is impossible. This bounds distinct zeros; the local zero count separately addresses multiplicity.

**DERIVED** The full massive formula and independent Feynman-parameter integration locate a conjugate pair at

$$
z = 4.45955827914522 +/- 28.7237236613321 i.
$$

**DERIVED** Upper-half-plane contour counts agree at two resolutions. A finite $L=30$ matched root at $4.14175183501002+31.8836377247863i$ was independently verified by the original signed-PV integral, residual $1.9\times10^{-15}$ or smaller, and a local contour count of one. Stable finite-regulator continuation approaches the limiting pair as L increases to 100,000. This does not contradict the earlier unmodified-kernel Stieltjes exclusion: the matching polynomials change that kernel.

**DERIVED** An exact-rational Rouche certificate subsequently established exactly one simple zero in the disk centered at $4.46+28.72i$ of radius .02. It bounds the complex Simpson remainder by analytic fourth-derivative majorants and bounds the local Taylor remainder by a uniform second derivative. The lower bound on the linear term is $3.4432\times10^{-6}$, exceeding the upper bound $7.7349\times10^{-7}$ on the center error plus nonlinear remainder. All decisive comparisons use rational arithmetic, including coefficient logarithm bounds. The companion implementation also checks its rational complex reciprocal and Simpson remainder bounds against exactly integrable polynomial controls. Combined with the Gram exclusion and real-axis certificates, this gives the complete nonzero first-sheet zero classification of this particular limiting tensor response: the one simple conjugate pair. The physical cut has nonzero boundary discontinuity above threshold and the threshold value is nonzero.

**DERIVED** For $q>0$, $q^2F_2'(q)$ is strictly increasing from -u to infinity, so F2 has one minimum. Numerically it lies at $q=38.8658670925567$, with $F2=2.33471358677353\times10^{-5}$. Exact rational inequalities place the minimum within [20,100] and give a positive lower bound $F2/C>0.00045865$ throughout that interval. Monotonicity extends exclusion to the rest of the positive axis. On $(-4,0)$, $q^2F_2'(q)\le-u+C/90<0$ and $F2(-4+)=-u/4+v-C/600<0$. Exact rational coefficient bounds certify both signs. Thus the reference tensor response has no additional real zeros in these intervals. — eq./page UNVERIFIED.

## Expansion and literature checks

**DERIVED** The origin-centered derivative expansion has radius 4, set by the physical threshold at -4. It converges at $q=1,2,3.5$; at $q=11.3066$ the eighth-order remainder error is about 2, while the sixty-fourth-order error exceeds 1e24. Increasing the Taylor order there does not validate the positive scalar zero. The exact dispersion representation supplies the appropriate continuation. This radius is not, by itself, a physical UV cutoff.

**DERIVED** The large-momentum tensor form is

$$
F2(z) ~ u/z + (C/240) \log z + v - 23C/1800.
$$

**DERIVED** Applying the known Lambert-equation construction to this asymptotic form gives a root at $18.22544+13.89422i$, about 70% away from the full massive root. The published mathematical method is useful, but neglecting massive threshold terms here is quantitatively poor.

**DERIVED** Primary literature used:

- **CITED** Ivanenko et al., Quasi-Herglotz functions and convex optimization: https://pmc.ncbi.nlm.nih.gov/articles/PMC7029951/ — Herglotz positive-measure representations and their signed-measure extensions. The scalar application and tensor Gram argument above are explicit applications to these kernels, not claims made by that paper.
- **CITED** Calmet et al., Physics Letters B 774 (2017) 332–337, https://doi.org/10.1016/j.physletb.2017.09.080 — logarithmic resummed propagator poles, local-coefficient dependence and the Lambert equation. Their model is different; no identification of our zeros with their proposed black-hole precursors follows. — eq./page UNVERIFIED.
- **CITED** Simon, Physical Review D 43 (1991) 3308, https://doi.org/10.1103/PhysRevD.43.3308; Parker and Simon, Physical Review D 47 (1993) 1339, https://doi.org/10.1103/PhysRevD.47.1339 — perturbative treatment and reduction of order in semiclassical gravity. Publisher abstracts were checked. These results motivate checking the approximation and solution prescription before interpreting extra roots; they do not automatically justify removing roots of this induced determinant. — eq./page UNVERIFIED.

## Remaining scientific boundary

**DERIVED** The mathematical response zeros are increasingly constrained. Their physical meaning still requires a specified causal prescription, an approximation valid at those momenta, and equality of the full sourced constrained-CEV generating functional with the auxiliary determinant, including measure, Jacobian and contact terms. A local cancellation in one source channel does not provide that equality. Nor can matching a finite number of derivative coefficients determine an entire microscopic theory.

**DERIVED** No unification, quantum consistency, ultraviolet completion, new particle, black hole or physical instability has been established. The next high-value work is the constrained generating-functional comparison and causal/source analysis, rather than further unconstrained pole scans or divergent Taylor truncations.

## Reproduction

**DERIVED** The companion archive includes this report, the earlier checkpoint, all experiment Python sources, and JSON results. Current scripts use only Python's standard library. Focus on scalar_rational_certificate.py, scalar_pick_check.py, tensor_real_rational_certificate.py, tensor_gram_check.py, tensor_rouche_certificate.py, tensor_subtracted_limit.py and tensor_matched_poles.py. The Rouche script reports success only after its strict rational inequality passes.
