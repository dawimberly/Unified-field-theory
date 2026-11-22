# FFff(B) mixed-pair integration by parts

This is a separate representation of the fixed-d Weyl insertion in FFFF_B_PROJECTIVE.md. It keeps the Yukawa vertex profiles through the differential operations. It does not yet supply an independent evaluation of the gradient part of the pole.

Order the fermion vectors a,b,c,d along edges 02,12,13,03. Distinguish the dimension d_dim from the fourth vector d. Let u=x0−x2 and v=x1−x3 with the two pair midpoints held fixed under each relative derivative. If their separation is w in the orientation m1−m0, the unpaired vectors have separations w+(u+v)/2 and −w+(u+v)/2. Their curved kernels also depend on their own midpoints. Thus derivatives must act on both arguments when expressed in original coordinates.

For a mixed pair with scalar Weyl kernel G and spinor vector b_vec, write F=G²/2. At fixed midpoint the mass perturbation gives δb_vec=−2κz∇Q. The scalar variation is δG=κ(2zQ−c_dim G₂), with c_dim=d_dim(4−d_dim)/3 and G₂ the flat scalar convolution. Hence G b_vec=−∇F+κE, with E=−c_dim G_flat∇G₂. Label the two pairs F₁,F₂ and E₁,E₂. The trace divided by two contracts them with T_ij=b_i d_j+d_i b_j−δ_ij b·d, times all vertex profiles. T includes the curved unpaired vectors and the four Yukawa factors.

The double derivative term integrated at fixed pair midpoints is F₁F₂ ∂_(u_i)∂_(v_j)T_ij after two integrations by parts. The two explicit order-κ terms are κ F₂ ∂_(v_j)(E₁_i T_ij) and κ F₁ ∂_(u_i)(E₂_j T_ij). In these explicit terms all spectator factors are flat. No replacement by flat equations of motion is used. Regulate first in a convergent dimension and then continue. For local-curvature extraction one may hold the overall midpoint of the four vertices fixed; that coordinate is invariant under both u and v derivatives. Fixing vertex 0 before differentiating would instead require accounting for the derivative of the anchoring delta.

## Explicit evanescent pole

Since G₂ has Schwinger measure s ds and G_flat has ds, the explicit E on each paired fermion line is +c_dim s times the original flat mixed-pair vector product, before integration by parts. Thus this part alone is obtained from the ordinary flat trace by multiplying it by c_dim(s₀+s₂), where indices are zero-based fermion-line indices. Apply the existing ordered-sector residue functional to the epsilon derivative of the added polynomial W d_dim(4−d_dim) N U (s₀+s₂)/3. The stripped double-projective residue of that derivative is +128/3 (24-node quadrature). With the radial and trace normalizations this contributes −1/54 to the bare double pole in units h²z/K at L=0. This is an internal split of the unchanged graph. In particular it must not be added to the original regulated graph, which already includes it. Dropping the term would shift the pole by +1/54, a larger shift than the observed 1/144 gap.

The gradient representation should supply the compensating part if the original graph is correct. The current regulated B double would imply +7/864 for that part. Evaluate the midpoint-fixed, differentiated spectator kernel independently before concluding whether the discrepancy comes from B or from the rest of the assembly.

## Runnable algebra checks

Python 3 and SymPy. The Gaussian test checks the IBP algebra in a smooth example, not the graph's Laurent coefficient. The evanescent pole diagnostic is reproducible by replacing the scalar-curvature test polynomial in FFFF_A_ANALYTIC_DOUBLE.md with the positive polynomial specified above and retaining the same full valuations.

```python
import sympy as s
u,v=s.symbols("u v", real=True)
# A one-component smooth test with nonseparable spectator kernel.
# It checks the sign, product derivatives, and both evanescent terms.
F,G,T=1+u*u,1+v+v*v,1+u*v+u*u*v
E=-u*(1+u*u); D=v*(1+v*v)
def der(p,x):
    return s.diff(p,x)-2*x*p
# Gaussian exp(-u^2-v^2) per factor after differentiating the complete product.
def joint(p,x,n):
    return s.diff(p,x)-2*n*x*p
def gauss(p):
    out=0
    for (a,b),coef in s.Poly(s.expand(p),u,v).terms():
        if a%2 or b%2: continue
        out+=coef*s.factorial2(a-1)*s.factorial2(b-1)/6**((a+b)//2)
    return s.expand(out)
direct=gauss(der(F,u)*der(G,v)*T)
reduced=gauss(F*G*joint(joint(T,u,1),v,1))
# F depends on u only and G on v only (including their Gaussian factors).
def moment(poly,n,m):
    out=0
    for (a,b),coef in s.Poly(s.expand(poly),u,v).terms():
        if not(a%2 or b%2):
            out+=coef*s.factorial2(a-1)*s.factorial2(b-1)/(2*n)**(a//2)/(2*m)**(b//2)
    return s.factor(out)
assert moment(der(F,u)*der(G,v)*T,2,2)==moment(F*G*joint(joint(T,u,1),v,1),2,2)
assert moment(-E*der(G,v)*T-der(F,u)*D*T,2,2)==moment(G*joint(E*T,v,1)+F*joint(D*T,u,1),2,2)
# Geometry: for pair midpoints separated by w, the other fermion
# endpoints have differences w+(u+v)/2 and -w+(u+v)/2.
w=s.symbols('w')
r1=w+(u+v)/2
r2=-w+(u+v)/2
assert all(s.diff(p,x)==s.Rational(1,2) for p in (r1,r2) for x in (u,v))
print("Tensor contraction, midpoint derivatives and IBP signs checked")
a,x,y,t,d,z=s.symbols("a x y t d z")
H=(a*x+(1-a)*y)**2+2*d*t*a*(1-a)
# Derivatives on the free Gaussian K_t are included in this expression.
rel=(s.diff(H,x)-s.diff(H,y))/2-(x-y)*H/(2*t)
expected=(x-y)*(x*x+y*y+x*y+(d-1)*t)/3
assert s.factor(s.integrate(-2*t*rel,(a,0,1))-expected)==0
# Flat tensor momentum projection reproduces the established chain numerator.
k,l=s.symbols("k l")
a1,b1,c1=s.symbols("a1 b1 c1") # k²,l²,k·l
p2=a1+b1+2*c1
assert s.expand(2*(a1+c1)*(b1+c1)-p2*c1-(2*a1*b1+(a1+b1)*c1))==0
print("Conformal spinor line and flat Fourier projection passed")

```

- [Bunch and Parker, Feynman propagator in curved spacetime](https://doi.org/10.1103/PhysRevD.20.2499).
- [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), Eq. (3.19).
