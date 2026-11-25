# Full tensor and contact Ward check

**Conversation date range** $2026-10$ (exact conversation dates UNVERIFIED)

**ASSUMED** This checkpoint records full tensor and contact Ward check. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

## Result

**DERIVED** The full 10 by 10 metric Hessian passes the longitudinal Ward check in 40 momentum/spectrum cases. Maximum absolute longitudinal component: $1.8474\times10^{-13}$. Maximum scalar-projector mismatch with $2A/3$: $1.509\times10^{-13}$. The three subthreshold scalar roots are reproduced directly from the tensor coefficient, with positive connected-source residues.

**DERIVED** This verifies the tensor-projection identity for this four-dimensional PV free determinant. Sixteen additional finite-D tensor tests at $D=4.01,4.2,5,5.9$ check the longitudinal, mixed longitudinal/transverse-trace and transverse-trace bilinears using exact dimensionally integrated loop moments. Their conservation components vanish and the scalar coefficient agrees with the prior full kernel within $1\times10^{-10}$. This supports the channel-specific DR source cancellation within this determinant model. It does not establish constrained CEV equivalence or a Lorentzian UV-complete theory.

## Derivation

**DERIVED** Take Q along axis 3 and $q=Q_E^{2}$, continued to negative q algebraically below the physical pair threshold $q=-4$. For a symmetric polarization h, the linear scalar kinetic vertex is $v_h=k$.h.(k+Q) -(tr h)/2 [k.(k+Q)+z]. Shift $k=l-xQ$ and define $H=h-(tr h)I/2, t=x(1-x), \Delta =z+tq, c_h=-t Q$.h.$Q+(tr h)(tq-z)/2$. Then $v_h=l$.H.$l+(1-2x)l$.H.$Q+c_h$.

**DERIVED** The rotationally averaged v_h v_j is integrated with the four-dimensional PV-summed logarithmic integrals $I_0=-C \log \Delta , I_2=2C \Delta \log \Delta , I4=-3C \Delta ^{2} \log \Delta$. Omitted polynomial pieces vanish after $\sum c_i z_i^k=0, k=0,1,2$. The $l^{4}$ coefficient is $[tr H tr J+2 tr(HJ)]/24$; the $l^{2}$ coefficient is [(1-2x)^2 (H Q).(J Q)+tr H c_j+tr J c_h]/4.

**DERIVED** With $p=-C \sum c_i z_i^{2} \log z_i/4$, the stationary contact Hessian bilinear is $K_{contact}[h,j]=p tr(hj)/2-p tr h tr j/4$. It includes the scalar quadratic metric vertex and the cosmological term tuned to cancel the tadpole. The bubble Hessian is $-\int (v_h v_j)/2$. Contact and bubble are evaluated separately; no transverse projection is imposed on the result.

**DERIVED** All ten real orthonormal symmetric polarizations are included. For Q along axis 3, the four pure gauge directions are h33,h03,h13,h23. Their complete Hessian rows vanish numerically. The bubble alone has nonzero longitudinal rows: omitting contacts fails this test.

**DERIVED** The transverse scalar polarization is $\operatorname{diag} (1,1,1,0)/\sqrt{3}$. Its eigenvalue is $K0=2A/3$, where A is the separately evaluated conformal log kernel. An on-shell identical-constituent source has trace $2-q$; its connected quadratic source residue is (2-q*)^2/[24 K0′(q*)]. This has the sign of a canonically normalized positive-kinetic Euclidean scalar below threshold. Positivity of this isolated residue is not proof of positivity of the whole spectrum.

## Tensor-derived roots

| PV L | q at root | K0 derivative | Connected-source residue |
|---|---:|---:|---:|
| **DERIVED** 3 | -2.090029857 | 0.000902978842 | 771.90551 |
| **DERIVED** 5 | -2.778519511 | 0.00250825095 | 379.318916 |
| **DERIVED** 10 | -3.747301964 | 0.0116509138 | 118.129074 |

## Remaining questions

1. **DERIVED** Does the constrained metric/clock-field measure produce the same Hessian with this regulator? This calculation does not include that measure or independent clock fluctuations.
2. **DERIVED** Does the Lorentzian contour and continuation give the canonical exchange sign when derived from the full microscopic path integral? The Euclidean source sign is now calibrated algebraically, but no contour prescription was established.
3. **DERIVED** What happens on and above the constituent threshold, including complex poles, spectral density and regulator thresholds? No global ghost-free verdict is warranted.
4. **DERIVED** After matching the Einstein and independent curvature-squared coefficients, which regulator differences remain? Zero bare terms across spectra is not a matched physical comparison.

**DERIVED** The next calculation should test the spectral discontinuity above the physical threshold while below the first regulator threshold, with the same full tensor contacts. It can expose negative spectral weight that a single subthreshold pole misses. Constrained-measure equivalence remains a separate requirement.

**CITED** Reference: Carone, Claringbold and Vaman, Composite graviton self-interactions in a model of emergent gravity, arXiv:1710.09367, sections II and V, https://arxiv.org/html/1710.09367v2 . Tensor derivation and numerics above are our calculation, not a claim attributed to that paper. — eq./page UNVERIFIED.

**DERIVED** Reproduce: python metric_tensor_ward_check.py (NumPy and SciPy). Dimensionless constituent $M=1$; results per scalar, overall N suppressed. All input spectra, individual checks, eigenvalues and quadrature refinements are in results.json.
