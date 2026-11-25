# FFZ CcFf analytic short-pair subtraction

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**DERIVED** This checkpoint records the FFZ CcFf analytic short-pair subtraction. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Continues [the curved CcFf kernel](ffz-ccff-projective-kernel.md). Fix the sum of the two single-edge lengths to 1. Write p for the scalar-pair total, q for the fermion-pair total, P=pf(1−f), Q=qg(1−g), S=1+P+Q, W=1+p+q. The projective measure is p^(ε−1)q^(ε−2) dp dq df dg dv. The coefficient multiplying it is the smooth numerator H=(qN) S^(−d/2) W^(−k) [A+kB/W], where k=3ε. Thus the short q constant gives a cutoff denominator ε−1 and the q-linear term gives a denominator ε.

**DERIVED** At d=4 only the cubic vertex term in A survives. Anchor at cubic vertex 0. The contracted cycle includes effective edges P,Q and the two single edges. The conditioned variance of cubic vertex 1 is M₁₁=8P(1+Q)/S+4P²Q/[S(1+P)]. The first term follows from the cycle resistance and the second from the covariance with the fermion endpoint difference. Also qN=2(1+P)/S. This yields

**DERIVED** H(p,0)|₄ = −16P/(1+P)³,
**DERIVED** ∫dg ∂qH(p,0)|₄ = −4P(P−2)/(1+P)⁴.

**DERIVED** Both vanish on the scalar corner at d=4. The integrated q-linear coefficient with dp/p is exactly 2, independent of the scalar split after changing p to P. This is one short-face contribution to the projective residue before assembly with the scalar endpoint and evanescent terms.

## Evanescent overlap

**DERIVED** At p=0 anchor the merged cubic vertices. Use the triangle with edges v,1−v,Q. Taking the q expansion of H at finite d before averaging v,g gives

**DERIVED** H₀,c = −d²(d−4)/6,
**DERIVED** H₁,c = −d(d−4)(17d²−72d+4)/72.

**DERIVED** Their leading terms at d=4−2ε are +16ε/3 and −4ε/3. The constant corner is multiplied by 1/[ε(ε−1)] and contributes at finite projective order. The linear corner is multiplied by 1/ε² and contributes −4/(3ε) at projective order under a unit cutoff. The radial Γ(3ε) can promote these to simple and double graph poles respectively. All cutoff terms must be combined with the remaining faces. These pieces alone do not determine the graph's Laurent coefficients.

## Runnable check

**DERIVED** Run with Python 3 and SymPy.

```python
import sympy as s
q,g,v,d=s.symbols('q g v d')
r=g*(1-g); Q=q*r; S=1+Q
# p=0: vertices 0,1 coincide at the root; 03 has length v.
c2=2*(1-v)*(v+Q)/S;c3=2*v*(1-v+Q)/S
c23=(c2+c3-2*Q/S)/2
n=s.Matrix([c2-c23,c23-c3])
C=s.Matrix([[c2,c23],[c23,c3]])
M=d*C+2*n*n.T/(2*Q/S)
M=M.applyfunc(s.cancel)
A=-(4-d)/2*s.trace(M)-d*(4-d)/3
B=s.Rational(2,3)*(v*(M[1,1]+d*v)+(1-v)*(M[0,0]+d*(1-v))
   +q*(M[0,0]+M[1,1]+M[0,1])+(d-1)*q*q*(1-2*r))
k=6-3*d/2
# q times N is d/(2S); measure also carries S^(-d/2) and W^(-k).
raw=d/2*(A+k*B/(1+q))
co0=s.factor(raw.subs(q,0))
co1=s.factor(s.diff(raw,q).subs(q,0)-(k+(1+d/2)*r)*co0)
av=lambda x:s.factor(s.integrate(s.integrate(s.expand(x),(v,0,1)),(g,0,1)))
print('corner constant',av(co0))
print('corner linear',av(co1))
print('4d',av(co0).subs(d,4),av(co1).subs(d,4))
assert av(co0)==-d**2*(d-4)/6
assert av(co1)==-d*(d-4)*(17*d*d-72*d+4)/72
P,Q=s.symbols("P Q",positive=True)
S=1+P+Q
M=8*P*(1+Q)/S+4*P**2*Q/(S*(1+P))
H=-2*(1+P)*M/S**3
assert s.cancel(H.subs(Q,0)+16*P/(1+P)**3)==0
linear=s.factor(s.diff(H,Q).subs(Q,0)/6)
assert s.integrate(s.cancel(linear/P),(P,0,s.oo))==2
```

**DERIVED** Next, integrate the subtracted scalar face and check the long-pair corners with finite d retained. Those terms and the finite projective contribution are needed to obtain the CcFf simple pole and predict cancellation in the full Newton sum.

- **DERIVED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177).
