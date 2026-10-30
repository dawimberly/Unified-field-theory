# Scalar finite-index theorem applicability audit
Research batch: 2026-10-29 US Central. Directory date is UTC.
Baseline: cc68498f42a7e7ced47af22f5080adc3eb3dc648.

## Decision
The fixed PV scalar inverse cannot belong to any finite-index generalized Nevanlinna class in the canonical spectral convention. Its already-established negative continuum interval implies arbitrarily many negative squares. Finite-index pole/zero bounds therefore cannot be imported to finish the scalar complex-pole problem.

This eliminates a proposed proof method; it neither proves nor disproves extra complex poles. It is not a new physical instability claim about constrained CEV, whose measure and matched observables remain unresolved.

## Published definition and applicability
Primary source inspected: J. Rovnyak and L. A. Sakhnovich, On the Krein-Langer integral representation of generalized Nevanlinna functions, Electronic Journal of Linear Algebra 11 (2004), 1-15, DOI https://doi.org/10.13001/1081-3810.1117 .
Journal PDF: https://journals.uwyo.edu/index.php/ela/article/download/235/235/235 .
The opening definition uses a reflection-symmetric meromorphic function and the kernel (h(z)-conjugate(h(w)))/(z-conjugate(w)), whose finite Gram matrices must have a uniformly bounded number of negative eigenvalues. The paper's first pages state the associated finite-index framework and representation. We use that definition, not an assumed class membership. The following obstruction is our elementary application to the saved determinant.

## Put the response in the correct sign convention
Let F(q)=1/A(q), with Euclidean q and principal first-sheet cut q<=-4. Set h(z)=F(-z). For an ordinary positive Euclidean propagator F(q)=1/(q+m^2), h(z)=1/(m^2-z) has positive imaginary part in the upper half-plane. Thus the canonical continuum sign is
Im h(s+i0)=pi rho_F(s),
rho_F(s)=Im[1/A(-s-i0)]/pi.
The known massless scalar coefficient and its relation to the Einstein tensor do not enter the obstruction below: all points are strictly inside a regulator-open continuum interval.

## Exact negative points, without repeating a spectrum scan
For z_i=[1,L,2L,5L+1] and s=6L, only the physical and first PV species are above their pair thresholds. The exact first regulator weight is c_1=-5(2L-1)/(4L+1).
Using the existing discontinuity formula,
Im A(-s-i0)=-(C pi/4) S,
S=(1+3L)^2 sqrt(1-2/(3L)) + c_1 (4L)^2/sqrt(3).
Both sides of the following comparison are positive rational numbers:
c_1^2 (4L)^4/3 > (1+3L)^4(1-2/(3L)).
The script checks this strict inequality with exact fractions for L=3,5,10,30,100. It implies S<0, Im A>0, and rho_F=-Im A/(pi |A|^2)<0. Since Im A is strictly nonzero, A cannot vanish at that point.

| L | s | rho_F(s), numerical normalization |
|---|---:|---:|
| 3 | 18 | -0.5319287447588884 |
| 5 | 30 | -0.16126568323687868 |
| 10 | 60 | -0.035956903497867426 |
| 30 | 180 | -0.0037215973270400575 |
| 100 | 600 | -0.0003270231815673326 |

These points are away from every threshold. The analytic boundary value is continuous there, so each lies in an open interval on which Im h(s+i0)<0 and h remains finite. No global interval endpoint calculation is needed for the following proof.

## Why a negative continuum interval gives unbounded negative index
Fix any positive integer n. Choose n distinct real points x_j in one such open interval. Put z_j=x_j+i eta and
K_ij(eta)=[h(z_i)-conjugate(h(z_j))]/[z_i-conjugate(z_j)].
For i=j, eta K_ii=Im h(x_i+i eta), which tends to a strictly negative number as eta decreases to zero.
For i!=j, x_i-x_j stays nonzero, the numerator has a finite limit, and eta K_ij tends to zero.

Therefore eta K tends entrywise (and hence in norm for this fixed n) to a strictly negative diagonal matrix. For sufficiently small positive eta the full n by n matrix is negative definite. Since n was arbitrary, there is no uniform finite bound kappa. This proves exclusion from every finite-index class; testing only n<=32 would not by itself prove it.

The argument needs only the continuous finite boundary values on the chosen interval, not global pole absence. It can therefore be used without assuming the result we were trying to prove.

## Why simple preprocessing does not bypass this obstruction
Subtracting finitely many real isolated-pole terms or real polynomial contact terms from F leaves the continuum imaginary boundary values unchanged away from their finite singularities. The same negative-interval argument survives. A standard squared real-rational factor is nonnegative on the real boundary and cannot reverse the sign throughout such an interval except by degenerating identically; isolated zeros/poles can be avoided.

Likewise adding a real polynomial P(q) to the inverse kernel A, as in the preceding formal low-energy counterterm matching, leaves
Im[1/(A+P)]=-Im A/|A+P|^2
on the cut wherever finite. Such polynomial matching does not cure this scalar regulator-open continuum sign. This does not select a physical UV model or authorize modifying the adopted action.

A finite rational quadrature approximation of a cut can have a finite index that grows with its number of negative-weight nodes. Four regulator species do not imply a finite negative index of the full continuum response. Counting negative continuum directions is also not counting distinct complex poles.

## New execution evidence
Run python experiments/scalar_negative_index_2026_10_30/check.py.
It imports the existing exact scalar evaluation; no older scan or passing test suite is rerun.

- Five exact rational inequalities establish negative sign at s=6L.
- Seventy-five new Hermitian kernel matrices cover five spectra, n=2,4,8,16,32 and eta/L=1e-2,1e-4,1e-6.
- At eta/L=1e-6 every tested matrix has exactly n negative eigenvalues. The largest sampled count is 32 for every spectrum. Hermiticity residual is zero at recorded precision.
- The results include eigenvalue extremes, off-diagonal decay and convergence to the predicted diagonal limit. The analytic arbitrary-n argument supplies the infinite-index conclusion; numerical counts illustrate it.
- Assertions executed successfully on Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0. Code and full results are preserved.

## Next work and limits
Do not spend another batch trying a finite-index theorem on this same full PV scalar response unless the represented function and its hypotheses change explicitly. Global scalar zero control remains open and may need certified slit-plane contour/asymptotic bounds or another analytic representation.

The higher-priority physical comparison is still the imported notes' reduced/constrained quantum measure and matched matter/composite insertions. Inspect the published auxiliary-preserving quantization and field-redefinition results, then write the exact assumptions needed for this model's generating functional. The stochastic construction remains a separate comparison proposal. No consequential physical assumption has been adopted and no user decision is required for this literature/application work.
