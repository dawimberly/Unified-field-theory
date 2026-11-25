# FFZ scalar endpoint residues

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**DERIVED** This checkpoint records the FFZ scalar endpoint residues. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Continues [the projective checkpoint](ffz-scalar-projective-checkpoint.md). Work at fixed ξ=1/6 with equal masses and no matter counterterms in the graph contributions below.

1. **DERIVED** Use projective homogeneity to fix the sum of the single-edge lengths to 1. For QCc the single edge has length 1. For CCcc the two single edges have lengths v and 1−v. The two parallel pairs have total lengths p,q and split fractions f,g.
2. **DERIVED** Combining the parallel pairs gives U=pq S, where S=1+pf(1−f)+qg(1−g), and the splitting Jacobian is pq. Set W=1+p+q. With k=N−3d/2, the remaining integrand against p^(ε−1)q^(ε−1) dp dq df dg is S^(−d/2)W^(−k−1)[A+(k+1)B/W], using the homogeneous A,B in the preceding checkpoint evaluated on the unnormalized lengths. CCcc additionally integrates dv. The radial factor is Γ(k+1). For these graphs the UV pair faces are p,q=0. A contracted triangle with both pair totals large does not yield an extra UV divergence in the curvature graph.
3. **DERIVED** Subtract faces by adding and subtracting the value at p=0 with a p<1 cutoff, and likewise for q. Then subtract the overlap once. The integrated corner has weight 1/ε² and each once-subtracted face has weight 1/ε. The mixed remainder is integrable at ε=0. Individual subtractions may use the vertex-fixed graph. Combining the terms uses the same fixed vertex throughout.
4. **DERIVED** At the QCc corner all vertices coalesce. A=−d(4−d)/3 and B=2d/3. Since k+1=3ε, the net numerator is 2dε/3.
5. **DERIVED** At the CCcc corner there are two merged vertices with covariance C=2v(1−v) at the nonroot vertex. A=−d[(6−d)C+(4−d)/3] and B=2d/3. With k+1=1+3ε and ∫C dv=1/3, its integrated numerator is also 2dε/3. Accordingly neither corner has a non-evanescent double projective pole. Both produce +8/(3ε) at projective level.
6. **DERIVED** For QCc on p=0, the conductance graph is two vertices linked by the other parallel pair and the single line. Its covariance is 2qr/(1+qr), r=g(1−g). At ε=0 the correction B drops out. The face after including dq/q gives −8r/(1+qr)^3. Its integral is −4, independent of r away from the measure-zero endpoints.
7. **DERIVED** The total QCc projective pole is therefore (8/3−8)/ε=−16/(3ε). Including Γ(3ε), the Wick weight +1/4 and 1/[4d(d−1)] gives the bare graph double pole −h²z/(108K ε²). It has no triple pole. Its simple pole still needs the finite projective part and must include the evanescent numerator.
8. **DERIVED** For CCcc on p=0 the reduced network is a triangle with edge lengths v,1−v and qr. Retain the two separate Schwinger lengths within the q pair when forming B. Integrating over v gives the following face, including dq/q:
   −8(2q²r²−q²r+3qr−q+r)/[3(q+1)²(qr+1)³].
9. **DERIVED** Integrating this face over q≥0 gives a rational function plus a log. Let r=g(1−g). The result is
   −4r/[3(r−1)]+8(5r−2)/[3(r−1)²]−8(2r²+2r−1)ln r/[3(r−1)³].
   Its g integral is exactly 4/3. Integration by parts with the rational primitive of the log coefficient leaves the integrable rational function 4r(4−r)/[3(1−r)²]. The omitted rational endpoint terms cancel the boundary log. The numerical integral agrees independently.
10. **DERIVED** The CCcc projective residue is 8/3+2(4/3)=16/3. Including Γ(1+3ε), the signed combined Wick factor −5/8 and 1/[4d(d−1)] yields −5h²z/(72K ε), with no double or triple pole. Separately the two CCcc Wick weights −1/8 and −1/2 give simple poles −h²z/(72Kε) and −h²z/(18Kε).

**DERIVED** These results are bare contributions in the baseline curvature sector. The QCc double pole constrains the sum with other graphs and counterterms. Neither Newton constraint can be decided from these graphs individually. The next scalar calculation is the finite once- and twice-subtracted QCc projective remainder, followed by the Yukawa and counterterm sectors.

## Runnable checks

**DERIVED** Python 3 with SymPy and SciPy.

```python
from scipy.integrate import quad
import sympy as s
q,v,r=s.symbols("q v r",positive=True)
S=1+q*r
c1=2*v*(S-v)/S
c2=2*(1-v)*(S-1+v)/S
c12=(c1+c2-2*q*r/S)/2
A=-4*(c1+c2)
B=s.Rational(8,3)*(v*(c1+v)+(1-v)*(c2+1-v)
                  +q*(c1+c2+c12)+q*q*(1-2*r))
face=s.factor(s.integrate(s.expand((A+B/(1+q))/(1+q)/S**2/q),(v,0,1)))
assert s.cancel(face+8*(2*q*q*r*r-q*q*r+3*q*r-q+r)/(3*(1+q)**2*S**3))==0
x=s.symbols("x",positive=True)
rr=x*(1-x)
a=-4*rr/(3*(rr-1))+8*(5*rr-2)/(3*(rr-1)**2)
b=-8*(2*rr**2+2*rr-1)/(3*(rr-1)**3)
primitive=(-16*x**3+24*x**2-16*x+4)/(3*(1-rr)**2)
assert s.cancel(s.diff(primitive,x)-b)==0
regular=s.cancel(a-primitive*s.diff(rr,x)/rr+s.Rational(4,3)*(1/x+1/(1-x)))
assert s.integrate(regular,(x,0,1))==s.Rational(4,3)
assert abs(quad(s.lambdify(x,a+b*s.log(rr),"math"),0,1,epsabs=1e-10)[0]-4/3)<1e-10
assert s.integrate(-8*r/(1+q*r)**3,(q,0,s.oo))==-4
assert (s.Rational(8,3)-8)/3/4/48==-s.Rational(1,108)
assert (s.Rational(8,3)+2*s.Rational(4,3))*(-s.Rational(5,8))/48==-s.Rational(5,72)
print("Scalar face reductions and pole coefficients passed")
```

- **DERIVED** [Heinrich et al., Expansion by regions with pySecDec](https://arxiv.org/abs/2108.10807), projective homogeneity and Cheng–Wu.
- **DERIVED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), scalar vacuum topologies.
