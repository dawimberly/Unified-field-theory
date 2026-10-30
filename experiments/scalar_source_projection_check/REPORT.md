# Scalar source saturation: cancellation in the identical-mass channel

## Main result

The exact finite-D conformal-kernel zero at q=4M^2/(D-2) cancels from the scalar exchange between identical minimally coupled constituent scalars, under the Ward-completed tensor interpretation stated below. Their stress-trace vertex vanishes at the same point. Thus calling this zero a demonstrated observable tachyon in constituent scattering was premature.

The cancellation is channel-specific: an unequal-mass or massless probe has a nonzero trace there. It does not prove that the metric-sector zero is absent or that the vacuum is stable. A new probe here means O(1) species, so it does not modify the leading-N determinant; adding O(N) new matter requires recomputing the kernel.

## Derivation

Use Euclidean q=Q^2 and analytic on-shell p1^2=p2^2=-m_ext^2, Q=p1+p2. E_ab=p1_a p2_b+p2_a p1_b-delta_ab(p1.p2-m_ext^2). Then Q_a E_ab=0 and tr E=2m_ext^2-(D-2)q/2=2U_ext.

For theta_ab=delta_ab-Q_a Q_b/q and P0=theta_ab theta_cd/(D-1), assume the stationary metric Hessian is transverse and write S2=(1/2)h K h. If h_ab=f delta_ab gives S2=A f^2, then K0=2A/(D-1). Coupling h_ab T_ab/2 and integrating the Gaussian gives the connected scalar source functional W0=(tr T)^2/(16A). A cross-exchange coefficient for two separate sources is twice this quadratic coefficient; the sign is unchanged.

In the assumed DR kernel A=-U(PU+g)/4, with U=M^2-(D-2)q/4. For identical constituent sources W0=U^2/(4A)=-U/(PU+g). Since g is nonzero, W0 is regular and zero at U=0. The inverse metric response still has its zero.

## PV source residues

For each previously found subthreshold PV root the identical physical constituent trace is nonzero. The conditional Euclidean W residue is positive, matching the sign of a real source coupled to a canonical positive-kinetic Euclidean scalar. This is stronger than inspecting A′ alone because it includes the source overlap and Gaussian normalization. It remains conditional on the tensor Ward identity and contour/sign conventions of the specified metric model.

| PV L | q at root | Stress trace | Conditional W residue |
|---|---:|---:|---:|
| 3 | -2.09002986 | 4.09002986 | 771.90550895 |
| 5 | -2.77851951 | 4.77851951 | 379.31891583 |
| 10 | -3.74730196 | 5.74730196 | 118.12907468 |

## Checks and limits

159 numerical algebra checks pass, including complex on-shell conservation, source traces, projector saturation, conformal normalization, and cancellation approaching both sides of the DR root. Integer dimensions test tensor algebra; the trace/projector formulas extend algebraically to finite D.

These tests do not calculate the full loop tensor or prove its Ward identities. They also do not establish the constrained CEV measure, Lorentzian integration contour, quartic mixing, 1/N corrections, above-threshold spectral positivity, or stability against arbitrary perturbations. Scalar-channel massless exchange also cannot be diagnosed as a ghost in isolation from the spin-2 Einstein tensor structure.

## Single most decisive next calculation

Compute the full second metric variation, including all seagull/contact terms, and explicitly check Q_a K_abcd=0 at the stationary background using the same regulator. Compare its transverse scalar coefficient with 2A/(D-1). If that map fails, the source verdict above must be revised. Then calibrate Lorentzian exchange against a canonical scalar and study unequal-mass probes.

Reference for the physical scattering setup and universality: Carone, Claringbold and Vaman, Composite graviton self-interactions in a model of emergent gravity, arXiv:1710.09367, sections II and V: https://arxiv.org/html/1710.09367v2 . The source-cancellation result here is our conditional derivation, not attributed to that paper.

Reproduce with python scalar_source_projection_check.py (requires the preceding scalar_full_response_scan/results.json at the documented workspace path, NumPy and SciPy). Input scan results and both scripts are included.