# FFZ curved QFf tadpole pole

Continues [the checked QCc scalar sum](ffz_qcc_full_regulated_check.md). This step covers the QFf class only.

In the covariant frame the quartic tadpole is constant in R at ξ=1/6. It acts as a scalar squared-mass insertion δx=|y|²T in the two-loop Yukawa sunset at fixed fermion mass z, where T=−z/(Kε)+O(1). The mass-shift normalization follows by contracting a pair of the X²X*² slots in |y|² X²X*²/4, with four choices. It also reproduces the flat J_SSFFS class of Martin (3.27).

Let G_x be the scalar kernel and b=−G_z' the radial kinetic fermion kernel. Radially integrating the fixed-1/6 Weyl insertion gives the sunset curvature term −|y|² [P/(12d)+(1/12−1/(6d))T_r], with P=∫r²G_x b² and T_r=−∫rG_x G_z b. This is the original two-loop radial formula with the scalar mass left independent. Its fermion derivative term is Q=T_r/2 because ∂z b=−rG_z/2. The integration by parts uses only the radial derivative of the product and holds for unequal x and z. In this paragraph primes on G refer to radial derivatives.

At x=z, the common-mass identities are P=−zJ/2+dI/3 and T_r=−dI/3, with I=∫G³ and J=∫r²G³. To differentiate only the scalar, subtract from the common derivative the fermion derivatives. This gives ∂xP=∂z(common)P+(d+2)J/3 and ∂xT_r=∂z(common)T_r/2−J/4. With I∝z^(d−3), J∝z^(d−4), the bracket's resulting derivative is

D_x=−(d−3)(d−4)I/(72z)+(32−5d)J/(144d).

The independent Bessel quadrature below verifies this expression in three convergent dimensions. The QFf curvature term is −|y|⁴ T D_x.

For J, order the three Schwinger lengths as ρ(1,u,uv). The parallel resistance supplies the additional r² moment. The six sectors give

K² J = 12d e^(2γ_E ε−2ε L) Γ(2ε) ∫₀¹ du dv u^ε v (1+v+uv)^(−3+ε)(1+u+uv)^(−2ε).

The integral at ε=0 is 1/12. The derivative of this integral is j₁ in the runnable check, yielding K²J=2/ε+24j₁−1−4L. The common sunset poles are K²I/z=−3/(2ε²)+(3L−9/2)/ε. Substitution into D_x cancels the 1/ε terms and the L dependence. Its finite limit is K²D_x=(24j₁−1)/48−1/24+1/18=−0.07466151326053. Thus the bare QFf result in units h²z/K is

−0.07466151326053 / ε.

There are no higher poles in this graph. Adding it to the checked scalar sector gives the subtotal 1/(36ε²)+[−0.12588115878−L/12]/ε. Three Yukawa classes and matter counterterms are still required. No projected Newton conclusion follows from this subtotal.

## Runnable check

Python 3 with SciPy.

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

Next evaluate CcFf, which has the same scalar-cycle layout as one CCcc topology but with a parallel fermion pair. Then compute both FFff classes and assemble the projected higher poles.

- [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), Eq. (3.27).
