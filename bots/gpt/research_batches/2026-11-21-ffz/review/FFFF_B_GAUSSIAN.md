# FFff(B) parity-even Gaussian reduction

The one-loop chain has four kinetic vectors. For oriented vertex differences a,b,c,d, the even trace is 2 times [(a·b)(c·d)−(a·c)(b·d)+(a·d)(b·c)]. The overall signed Wick prefactor follows the signed spinor review. The parity-odd trace vanishes for this vacuum response: four differences of four vertices span at most three independent coordinate vectors, and the integrated linear-curvature scalar has no pseudoscalar source. The arbitrary-endpoint kinetic insertion remains parallel to the flat endpoint difference after integrating the split fraction, so it does not add a fourth vector.

For their covariance matrix c_ij per component, the mean trace divided by 2 is d²(c_ab c_cd−c_ac c_bd+c_ad c_bc)+2d c_ac c_bd. To add a Weyl position insertion x_i·x_j, tilt the original vertex Gaussian by exp(λ x_i·x_j). Its first derivative shifts its covariance by C J C, where J has symmetric ij and ji unit entries, or 2 at the ii diagonal entry if i=j. Apply this shift to all covariances in the trace polynomial and add d C_ij times the original polynomial for the normalization of the tilted Gaussian. This avoids dividing by the spin trace where it crosses zero. Use these weighted moments with the scalar and spinor line factors and the vertex factors already established for FFff(A). To evaluate the actual FFff(B) pole, insert the oriented four-cycle differences and the two scalar parallel edges, and resolve the resulting projective endpoints. The check below verifies the Gaussian algebra on a generic covariance. It is not a graph integral or a pole evaluation.

```python
import numpy as np
import sympy as s
rng=np.random.default_rng(42)
L=rng.normal(size=(4,4))*.35+np.eye(4)
C=L@L.T
a,b,c,d=(rng.multivariate_normal(np.zeros(4),C,size=(400000,4)).transpose(2,0,1))
dot=lambda u,v:np.sum(u*v,axis=1)
samples=dot(a,b)*dot(c,d)-dot(a,c)*dot(b,d)+dot(a,d)*dot(b,c)
dimension=4
ab,ac,ad,bc,bd,cd=C[0,1],C[0,2],C[0,3],C[1,2],C[1,3],C[2,3]
expected=dimension**2*(ab*cd-ac*bd+ad*bc)+2*dimension*ac*bd
assert abs(samples.mean()-expected)<5*samples.std()/len(samples)**.5
# Tilt the joint Gaussian by exp(lambda x_i dot x_j); then differentiate at zero.
la=s.symbols("la")
J=np.zeros((4,4));J[0,1]=J[1,0]=1
dc=C@J@C
def trace(c):
 return 16*(c[0,1]*c[2,3]-c[0,2]*c[1,3]+c[0,3]*c[1,2])+8*c[0,2]*c[1,3]
derivative=float(s.diff(trace(C+la*dc),la).subs(la,0))
weighted=samples*dot(a,b)
assert abs(weighted.mean()-(derivative+4*C[0,1]*expected))<5*weighted.std()/len(samples)**.5
print("FFff(B) parity-even trace and weighted coordinate moment passed")

```
