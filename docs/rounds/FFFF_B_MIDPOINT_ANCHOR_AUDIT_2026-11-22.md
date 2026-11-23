# FFff(B) midpoint gradient and fixed-source counterterm audit
November 22, 2026. Parent PR head 01513fc8655cada0530ee1da358c2c32fc83d86d. Units h²z/K, L=ln(z/μ²), d=4−2ε. The physical comparison continues to use fixed ξ=1/6 and the original vertex-anchored prescription.

## Midpoint differentiation

Set x0=(u−w)/2, x2=(−u−w)/2, x1=(v+w)/2, x3=(−v+w)/2 at vanishing overall coordinate midpoint. The two spectator differences are r=(u+v)/2+w and t=(u+v)/2−w, with lengths b and d_s. Write their exponential as exp(−r²/(4b)−t²/(4d_s)). For its polynomial spectator tensor P_ij, define q=(r/b+t/d_s)/4 and λ=(1/b+1/d_s)/8. Differentiating both arguments of each spectator and all four vertex factors gives
```
D_ij P = ∂ui∂vj P − q_i∂vj P − q_j∂ui P + (q_i q_j − λ δij)P.
```
In each paired scalar square retain the factor 1/2. Compare at fixed lengths after independently symmetrizing the scalar/fermion lengths within both pairs. Without this symmetrization the gradient of a product of identical scalar kernels distributes the derivative differently from a labeled kinetic propagator.

Five-node tensor Gauss–Hermite integrates the degree-at-most-eight position polynomial in d=2. The test obtains derivatives by extracting polynomial coefficients and analytically differentiating the Gaussian, separately from the review covariance contractions. The normalized flat, mass-independent curvature, and mass-dependent curvature expectations at lengths (1.1,.8,1.3,.7,.5,1.9) are (0.07308345, −0.78296377, 1.11686567) and at (.03,.4,.9,.2,1.1,.7) are (0.47781325, −2.34135333, 2.02165721). Both comparisons agree within 9e−16 after centering the direct covariance and removing its explicit E term. This verifies the differentiated kernel at fixed lengths and integer d=2, with a polynomial differential identity valid componentwise in general d.

The evanescent E split has polynomial ε-derivative (8/3)W N|₄ U(s₀+s₂). Applying the double-projective coefficient functional with the full J valuations gives 42.6666681068, 42.6666666667, 42.6666666667 at 8, 16, 24 Gauss nodes, respectively, supporting 128/3 and the graph double −1/54. E depends only on relative flat covariances and is center invariant.

## The anchoring obstruction survives integration

Simply holding the arithmetic midpoint fixed replaces C by PCP in the insertion moments, where P=I−11ᵀ/4. This operation preserves U and N. I integrated only the difference J_centroid−J_vertex, using its own full finite-d endpoint valuations. Two subtracted extractions (1024 points/sector, seed 23, cutoff 1 and 4096, seed 41, cutoff 0.5) give intercepts 16.01124004 and 15.99518412 and first coefficients 53.67737669 and 54.05936817 in ε²ΔI. Separate endpoint-only Gauss quadrature gives 15.9999702675, 15.99999999999977, 15.99999999999999 at 8, 16, 24 nodes. Thus the raw midpoint replacement shifts the B triple by approximately −1/144. The coefficient 16 is supported numerically. The midpoint-minus-vertex double is approximately −0.0522 to −0.0524 at L=0, using −(ΔI₋₁+25 ΔI₋₂/6)/2304. This spread is a diagnostic, not an error bound. Subtracting E from the raw midpoint therefore gives a gradient double around −0.044. The vertex-normalized +7/864 value remains conditional on restoring the anchoring terms and on the existing vertex-anchored regulated inference.

The raw arithmetic-midpoint prescription fails the triple check and cannot resolve the subleading discrepancy. At fixed arithmetic midpoint c, restore the original vertex anchor by inserting δ(c+(u−w)/2) before differentiating. When moving uᵢ derivatives onto spectators, ∂uᵢδ=∂cᵢδ/2 must be retained; the v derivative leaves the delta invariant for vertex 0. Alternatively derive the corresponding covariant center density with its volume Jacobian before extracting a local R coefficient. Use the triple delta above as the first numerical normalization check. The d=2 derivative test alone does not establish the physical vertex-anchored residue.

## Bare action and counterterms

For W=mXY+yX³/6, normalize the chiral fields canonically before the Weyl redefinition. With only X wavefunction subtraction Z_X, z_B=z/Z_X and h_B=μ^(2ε)h/Z_X³ up to the fixed measure convention. Eliminating the canonical auxiliary fields produces z_B(|X_B|²+|Y_B|²), cubic (m_B* y_B Y_B* X_B²/2 + conjugate), and quartic |y_B|²|X_B|⁴/4. The fermions have off-diagonal mass and an X-only Yukawa vertex. The Weyl weights of cubic, quartic and Yukawa couplings are 1+ε, 2ε, ε. In a renormalized-field presentation the X kinetic and curvature term must be dressed by the same Z_X, under the packet assumption of homogeneous fixed-source subtraction. This leaves the canonical ξ_X=1/6 at a_X=0. The Y scalar keeps its canonical kinetic normalization and its curvature source.

The fixed scalar propagator shift is (a_d−a)/6 times R, a_d=ε/(3−2ε). In σ=−κx², R=4d(d−1)κ, this gives the resolvent term −κ d(4−d)G₂/3 at a=0. For Y, a_{C,Y}=a_Y−a_d since Z_{M,Y}=1. The projected insertion relative to the fixed-ξ baseline remains Δa_{C,Y}=−1 exactly. Dropping −a_d in an absolute conversion instead changes the order-h² baseline by an amount whose magnitude and signs can be read from −a_d times its Y pole polynomial: −1/(108ε²)+1/(81ε) for the translation itself. This cannot be appended as an extra diagram when the fixed scalar propagator already carries the evanescent curvature term. X requires its own Z_M and evanescent conversion. The local-RG argument for homogeneous subtraction is conditional on the packet's covariant kinetic, source and EOM conventions.

The flat mass recursion gives z_B/z=1+h/(2ε)+h²[1/(2ε²)−1/(4ε)]. At a_X=a_Y=0 the heat-kernel coefficient of each scalar is zero and the fixed four-component Dirac trace gives the one-loop baseline proportional to Γ(ε)z^(1−ε)/[6(1−ε)] with the packet sign. The h² insertion after differentiating has coefficient [δ₂−εδ₁²/2]Γ(ε)e^(γ_E ε−εL)/6 with δ₁,δ₂ the mass coefficients above. Its higher poles are 1/12 and −L/12−1/16. The two-loop bare graph is proportional to hz^(1−2ε), giving insertion (2−ε)/ε and higher poles −1/12 and L/6−1/24. The combined triple is zero and the double is L/12−5/48. This checks the mass and coupling insertions. An independent curved local-source projection has yet to test for operators outside the assumed continuation.

The vertex-anchored residual inferred from B is −1/144 in the double and +L/48 in the simple because the bare contribution scales as exp(−3εL). A correction confined to a two-loop or one-loop insertion cannot cancel both with a single double-pole shift and unchanged triples. If both insertion classes were the sole source, their double constants would have to shift by +1/72 and −1/144, respectively, from the two independent equations. These are diagnostics. Next, restore the anchor derivatives at fixed lengths and verify the equivalent projective kernel, then audit the remaining vertex-anchored bare double poles against this paired pole/log check. The physical simple pole versus 1/24 remains open.

## Reproduction

Extract the Python kernel in FFFF_B_PROJECTIVE.md as fffb_midpoint.py and insert `P=np.eye(4)-np.ones((4,4))/4; C=P@C@P` immediately after the covariance construction. In the regulated sector script in FFFF_B_REGULATED.md, form A and B separately for C and PCP (recompute cv and all moments each time), take the differences before constructing J, and keep the same U and N. Use ε=(0,.006,.012,−.006,−.012) and run `10 23 1` and `12 41 .5` with OPENBLAS_NUM_THREADS=1. For endpoint-only quadrature use the ordered-sector functional, applying J|₄ or the E ε-derivative while taking minima from full J and regular-axis powers from the finite-d valuation table. Only the 40 two-log sectors contribute for the midpoint difference. The E check also retains the power/log Taylor coefficients.

Primary references read before tests: [Martin, flat topology Eq. 3.19](https://doi.org/10.1103/PhysRevD.96.096005), [Martin, soft vacuum beta Eq. 4.40](https://doi.org/10.1103/PhysRevD.109.015019), [Jack and Osborn, local sources Eq. 7.5 and Eq. 9.26](https://arxiv.org/abs/1312.0428), and [Bunch and Parker, curved propagator representation](https://doi.org/10.1103/PhysRevD.20.2499).

```python
import itertools
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from fffb_midpoint import kernel

def midpoint(lengths, dim=2):
    # Normalized Gaussian expectation at fixed overall midpoint.
    a,b,c,d,e,f = lengths
    M=np.array([[1,0,0],[.5,.5,1],[0,1,0],[.5,.5,-1]])
    m=M[[0,1,2,3,0,2]]
    lower=np.linalg.cholesky(2*np.linalg.inv(m.T@np.diag(1/lengths)@m))
    nodes,weights=hermegauss(5)
    ids=np.array(list(itertools.product(range(5),repeat=3*dim)))
    wt=np.prod(weights[ids],axis=1)/(2*np.pi)**(3*dim/2)
    coords=np.einsum("ij,njk->nik",lower,nodes[ids].reshape(-1,3,dim))
    u,v,w=coords.transpose(1,0,2)
    g=lambda x,y: np.sum(x*y,axis=-1)
    cp=dim*(4-dim)/3
    qa=(1/a+1/e)/2
    qc=(1/c+1/f)/2
    qs=(1/b+1/d)/8
    # The spectator polynomial is degree at most four in translations.
    # Roots-of-unity coefficient extraction differentiates it without
    # finite-difference truncation. Differentiate the exponential analytically.
    phase=np.exp(2j*np.pi*np.arange(5)/5)
    def spectator(i,j,z1,z2):
        shift=np.eye(dim)[i]*z1+np.eye(dim)[j]*z2
        r=(u+v+shift)/2+w
        t=(u+v+shift)/2-w
        U=u+np.eye(dim)[i]*z1
        V=v+np.eye(dim)[j]*z2
        x0=(U-w)/2; x2=(-U-w)/2
        x1=(V+w)/2; x3=(-V+w)/2
        H=lambda x,y,s: (g(x,x)+g(y,y)+g(x,y)+(dim-1)*s)/3
        T=(r[:,i]*t[:,j]+t[:,i]*r[:,j]-(i==j)*g(r,t))/(4*b*d)
        return np.array([T, -(4-dim)/2*sum(g(x,x) for x in (x0,x1,x2,x3))*T,
                         2*(b*H(x1,x2,b)+d*H(x0,x3,d))*T])
    def hs(r,w,s):
        return (3*g(w,w)/4+g(r,r)/4+dim*s)/3
    pairB=2*sum(s*hs(r,w,s) for r,s in ((u,a),(u,e),(v,c),(v,f)))
    pairA=-cp*(a+e+c+f)
    vals=np.zeros(3)
    r=(u+v)/2+w; t=(u+v)/2-w
    q=(r/b+t/d)/4
    for i,j in itertools.product(range(dim),repeat=2):
        P=np.array([[spectator(i,j,p,k) for k in phase] for p in phase])
        coeff=lambda n,m: np.einsum("a,b,abkn->kn",phase**(-n),phase**(-m),P).real/25
        out=coeff(1,1)-q[:,i]*coeff(0,1)-q[:,j]*coeff(1,0)+(q[:,i]*q[:,j]-qs*(i==j))*coeff(0,0)
        vals += np.array([wt@out[0],wt@(out[1]+pairA*out[0]),wt@(out[2]+pairB*out[0])])
    return vals/4

for lengths in (np.array([1.1,.8,1.3,.7,.5,1.9]),np.array([.03,.4,.9,.2,1.1,.7])):
    cp=2*2/3
    def direct(s):
        n,A,B=kernel(s[None,:],2,0)
        R=np.eye(4)[[0,1,1,0,0,1],1:]-np.eye(4)[[2,2,3,3,2,3],1:]
        U=np.prod(s)*np.linalg.det(R.T@np.diag(1/s)@R)
        return np.array([n[0],A[0]-cp*(s[0]+s[2])*n[0],B[0]])*U
    ds=np.mean([direct(lengths[list(p)]) for p in ((0,1,2,3,4,5),(4,1,2,3,0,5),(0,1,5,3,4,2),(4,1,5,3,0,2))],axis=0)
    result=midpoint(lengths)
    print(lengths,result,ds, np.max(abs(result-ds)))
    np.testing.assert_allclose(result,ds,atol=2e-12,rtol=2e-12)
import sympy as s
e,h,L=s.symbols('e h L')
H=-2*e*h+3*h**2-3*h**3
Zz=1+h/(2*e)+h**2*(1/(2*e**2)-1/(4*e))
trunc=lambda p:s.expand(p).series(h,0,3).removeO().expand()
assert trunc(H*s.diff(Zz,h)+(h-h**2)*Zz)==0
a1=Zz.coeff(h); a2=Zz.coeff(h,2)
one=(a2-e*a1**2/2)*(1/e-L)/6
two=(2-e)/e*(-1/(24*e**2)+(L/12-s.Rational(1,24))/e)
assert s.simplify(one.coeff(e,-3) if one.is_Add else s.expand(one).coeff(e,-3))==s.Rational(1,12)
assert s.simplify(s.expand(one+two).coeff(e,-2)-L/12+s.Rational(5,48))==0
r=-s.Rational(1,96)+s.Rational(1,288)
assert r==-s.Rational(1,144) and -3*r==s.Rational(1,48)
x,y=s.symbols('x y')
assert s.solve([r+x+y,-3*r-2*x-y],(x,y))=={x:s.Rational(1,72),y:-s.Rational(1,144)}
ad=e/(3-2*e)
y3=1/(36*e**3)-1/(18*e**2)+1/(24*e)
assert s.series(-ad*y3,e,0,0).removeO()==-1/(108*e**2)+1/(81*e)
print('Mass recursion, fixed-source counterterms, Y translation, residual and L log passed')

```
