# Regulated stress bubble, seagull and stationary contact dictionary

Date: 2026-10-30 US Central. Baseline: d21e4d562a7924a053d49475433e414192fb9076.

## Scope

Free minimally coupled Euclidean scalar determinant in D=4, per scalar with signed four-species PV sum. Use the existing covariant regulator prescription, with common ultralocal measure contributions canceled by sum c_i=0; no additional metric-dependent measure term is introduced. This is not an assumed constrained CEV measure. Define q=Q_E^2, C=(4pi)^-2, z_i=m_i^2, and

p=-C sum_i c_i z_i^2 log(z_i)/4.

All polynomials canceled by the three PV moment conditions are understood consistently in tadpoles and bubbles. P[h,j]=tr(hj)-tr(h)tr(j)/2 on symmetric metric variations g=I+h.

## Derivation in the present convention

Start from S_E=integral sqrt(g)[g^{ab} partial_a phi partial_b phi+z phi^2]/2. At flat metric,

T_ab=partial_a phi partial_b phi-delta_ab[(partial phi)^2+z phi^2]/2,
S_h=-h:T/2.

Define Pi[h,j]=<T_h T_j>_connected with the same PV sum; it is a signed regulated correlator, not a positive covariance. For Gamma=-log Z_scalar+lambda integral sqrt(g), fixed-field differentiation gives

Gamma_hj=<S_hj>-Pi[h,j]/4+lambda (sqrt(g))_hj.

This identity follows by differentiating the normalized path integral, not by replacing an operator derivative with a derivative of its expectation value. The bubble factor follows from the two -T/2 insertions.

The regulated coincident contractions are

sum c_i <partial_a phi_i partial_b phi_i>=p delta_ab,
sum c_i <z_i phi_i^2>=-4p.

Hence <T_ab>=p delta_ab. The scalar first derivative is -p tr(h)/2, canceled by lambda=p. For clarity, expansion of the density gives

sqrt(g)=1+tr(h)/2+(tr h)^2/8-tr(h^2)/4+...,
sqrt(g)g^-1=I+tr(h)I/2-h+[(tr h)^2/8-tr(h^2)/4]I-tr(h)h/2+h^2+....

Using the coincident contractions in this *fixed-field* second derivative, rather than differentiating the contractions themselves, yields:

| Contribution | Bilinear coefficient |
|---|---|
| Scalar seagull <S_hj> | p P[h,j] |
| Stationary cosmological Hessian | -p P[h,j]/2 |
| Sum of local contacts | p P[h,j]/2 |
| Scalar connected bubble | -Pi[h,j]/4 |

Thus K=pP/2-Pi/4. This reproduces the repository's normalization: its vertex integral is Pi/2 and its bubble Hessian is minus that integral divided by two. At this fixed prescription the previous quadratic bridge's H=V P^-1/2 agrees with the local matrix for D=4 if V=p, because P^-1=P. This closes the algebraic normalization dictionary inside the specified free determinant. It still does not fix the composite theory's independently regulated coupling or source counterterms.

Normal ordering T at a fixed reference metric leaves its connected two-point function unchanged. A metric-dependent subtraction must, however, be differentiated when computing metric response; it can generate local terms. Therefore a normal-ordering label does not alone justify omitting the seagull or identifying Pi with a metric Hessian.

## Independent new radial checks

With u=ell_E^2, the radial measure is C u du. For barycentric PV weights w_i and c_i=w_i/w_0,

f(u)=sum_i c_i/(u+z_i)=-1/[w_0 product_i(u+z_i)].

This rational representation avoids subtractive cancellation. Define I=integral_0^infinity u^2 f(u)du. Then p=C I/4. The kinetic-plus-mass tadpole cancels pointwise in the signed sum because sum c_i z_i/(u+z_i)=-u f(u).

An independent zero-momentum Wick/angular calculation gives

Pi[h,j](0)=2C{[tr(h)tr(j)/24+tr(hj)/12] J-tr(h)tr(j) I/4},
J=integral_0^infinity u^3[-f'(u)]du=3I.

The last equality follows by integration by parts; boundary terms vanish since f(u) is finite at zero and falls as u^-4. Substitution gives Pi(0)=2pP and K(0)=0 as required. This is a new radial normalization check, not a rerun of the previous finite-q Ward scan.

For L=3,5,10,30,100, independent radial quadrature agrees with the closed logarithmic p to relative error at most 2.23e-16, and the resulting Pi(0) tensor to at most 4.82e-16 in units of p. Actual quadrature tolerances and error estimates are retained; these floating-point errors are empirical agreement, not rigorous interval enclosures. The script and results include every spectrum and separate contact coefficient.

Run `python experiments/stress_contact_dictionary_2026_10_30/check.py`. All new assertions passed. Python/NumPy/SciPy versions are saved. Neither full finite-q stress integrals nor constrained amplitudes were newly calculated here.

## Primary-literature context and limitations

Bzowski, McFadden and Skenderis, JHEP11(2018)153, DOI 10.1007/JHEP11(2018)153, https://inspirehep.net/files/a9e0845c8d7d1acb1809e3504cd3c00e. Inspected the published introduction and Appendix C, especially C.1--C.6: their source-functional conventions retain contact contributions from differentiating operators and metric factors. This supports keeping definitions explicit; their conformal three-point results are not imported as massive PV two-point formulas. The equations and radial checks above are derived here.

The pending bridge to the composite model requires one common regulated source functional, including metric-dependent subtractions and measure. The prior paper's contact cancellation under its conventions is not disproved by a different allocation of our local terms. A physical contour, interactions, 1/N and global scalar spectrum remain open. No health or UV-completion conclusion follows.
