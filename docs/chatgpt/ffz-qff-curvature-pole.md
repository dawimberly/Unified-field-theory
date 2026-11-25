# FFZ curved QFf tadpole pole

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**ASSUMED** This checkpoint records the FFZ curved QFf tadpole pole. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**ASSUMED** Continues [the checked QCc scalar sum](ffz-qcc-full-regulated-check.md). This step covers the QFf class only.

**CITED** In the covariant frame the quartic tadpole is constant in R at $\xi =1/6$. It acts as a scalar squared-mass insertion $\delta x=|y|^2T$ in the two-loop Yukawa sunset at fixed fermion mass z, where $T=-z/(K\epsilon )+O(1)$. The mass-shift normalization follows by contracting a pair of the $X^{2}(X^*)^{2}$ slots in $|y|^{2} X^{2}(X^*)^{2}/4$, with four choices. It also reproduces the flat J_SSFFS class of Martin (3.27). — eq./page UNVERIFIED.

**DERIVED** Let G_x be the scalar kernel and $b=-G_z'$ the radial kinetic fermion kernel. Radially integrating the fixed $1/6$ Weyl insertion gives the sunset curvature term $-|y|^2$ $[P/(12d)+(1/12-1/(6d))T_r]$, with $P=\int r^{2}G_x b^{2}$ and $T_r=-\int rG_xG_zb$. This is the original two-loop radial formula with the scalar mass left independent. Its fermion derivative term is $Q=T_r/2$ because $\partial_zb=-rG_z/2$. The integration by parts uses only the radial derivative of the product and holds for unequal x and z. In this paragraph primes on G refer to radial derivatives.

**DERIVED** At $x=z$, the common-mass identities are $P=-zJ/2+dI/3$ and $T_r=-dI/3$, with $I=\int G^{3}$ and $J=\int r^{2}G^{3}$. To differentiate only the scalar, subtract from the common derivative the fermion derivatives. This gives $\partial_xP=\partial_z^{\mathrm{common}}P+(d+2)J/3$ and $\partial_xT_r=\partial_z^{\mathrm{common}}T_r/2-J/4$. With $I\propto z^{d-3}, J\propto z^{d-4}$, the bracket's resulting derivative is

**DERIVED** $D_x=-(d-3)(d-4)I/(72z)+(32-5d)J/(144d).$

**DERIVED** The independent Bessel quadrature below verifies this expression in three convergent dimensions. The QFf curvature term is $-|y|^{4} T D_x$.

**DERIVED** For J, order the three Schwinger lengths as $\rho (1,u,uv)$. The parallel resistance supplies the additional $r^{2}$ moment. The six sectors give

**DERIVED** $K^{2} J = 12d e^{2\gamma _E \epsilon -2\epsilon L} \Gamma (2\epsilon ) \int _{0}^{1} du dv u^\epsilon v (1+v+uv)^{-3+\epsilon }(1+u+uv)^{-2\epsilon }$.

**DERIVED** The integral at $\epsilon =0$ is $1/12$. The derivative of this integral is j₁ in the runnable check, yielding $K^{2}J=2/\epsilon +24j_{1}-1-4L$. The common sunset poles are $K^{2}I/z=-3/(2\epsilon ^{2})+(3L-9/2)/\epsilon$. Substitution into D_x cancels the $1/\epsilon$ terms and the L dependence. Its finite limit is $K^{2}D_x=(24j_{1}-1)/48-1/24+1/18=-0.07466151326053$. Thus the bare QFf result in units $h^{2}z/K$ is

**DERIVED** $-0.07466151326053 / \epsilon$.

**DERIVED** There are no higher poles in this graph. Adding it to the checked scalar sector gives the subtotal $1/(36\epsilon ^{2})+[-0.12588115878-L/12]/\epsilon$. Three Yukawa classes and matter counterterms are still required. No projected Newton conclusion follows from this subtotal.

## Runnable check

**DERIVED** Python 3 with SciPy.

```python
from scipy.integrate import quad,dblquad
from scipy.special import kv,gamma
from math import pi,exp,log
K=16*pi*pi
def check(d):
 nu=d/2-1
 def integrand(r,kind):
  c=r**(-nu)/(2*pi)**(d/2)
  G=c*kv(nu,r)
  b=c*kv(nu+1,r)
  Gx=-r*c*kv(nu-1,r)/2
  weight=2*pi**(d/2)/gamma(d/2)*r**(d-1)
  vals=[G**3,r*r*G**3,Gx*(r*r*b*b/(12*d)-(1/12-1/(6*d))*r*G*b)]
  return weight*vals[kind]
 I,J,D=[quad(integrand,0,60,args=(i,),epsabs=2e-11)[0] for i in range(3)]
 ref=-(d-3)*(d-4)*I/72+(32-5*d)*J/(144*d)
 assert abs(D-ref)<1e-9
 return D
for d in (1.8,2.2,2.6): print(d,check(d))
a=dblquad(lambda t,x:t/(1+t+x*t)**3*(log(x)+log(1+t+x*t)-2*log(1+x+x*t)),0,1,0,1,epsabs=1e-11)[0]
c=(24*a-1)/48-1/24+1/18
assert abs(c+0.07466151326053)<1e-11
print(c)
assert abs(dblquad(lambda t,x:t/(1+t+x*t)**3,0,1,0,1)[0]-1/12)<1e-12
```

**ASSUMED** Next evaluate CcFf, which has the same scalar-cycle layout as one CCcc topology but with a parallel fermion pair. Then compute both FFff classes and assemble the projected higher poles.

- **CITED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), Eq. (3.27).
