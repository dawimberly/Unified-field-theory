# FFff(A) full soft-corner check

This check supersedes the bilinear short-pair result in the earlier FFff(A) review. Use the connected cycle with adjacent fermion pairs, effective resistances P=qf and Q=rg, and scalar lengths v and 1−v. Here f,g are the products of the two fractions in each pair. The anchor is vertex 0 and the cycle positions are 0,P,P+v,P+v+Q.

For the two endpoint differences Δ₁, Δ₂, let c,e be their variances per coordinate and j their covariance. Gaussian contraction gives E[Δ₁² Δ₂²]=d²ce+2dj². If a_i,b_i are the covariances of vertex i with the two differences, then the additional unnormalized insertion in E[Δ₁²Δ₂² x_i·x_j] is 2d²(e a_i a_j+c b_i b_j)+4dj(a_i b_j+b_i a_j). The remaining piece is d C_ij times E[Δ₁²Δ₂²]. The code divides by the flat Schwinger denominators and multiplies by qr for the projective measure. The signed graph weight remains −1/8 in units |y|⁴, with Wick contractions as in the Martin normalization.

Two spin factors reduce the radial degree by two. With k=6−3d/2, the curvature insertion ρA+zρ²B gives Γ(k−1), and the full projective numerator is S^(−d/2) W^(1−k)[A+(k−1)B/W], including the spin weight in A,B. Alternatively factoring Γ(k) produces the same answer if the A term is divided by k−1. At d=4 the full numerator is the negative of the spin-weighted B term.

Exact d=4 angular averages for the B term are (32/3,0,8,0) for the constant, q, r, qr. Thus the first disputed result is the B mixed derivative −4/27. Its corrected value is zero. The complete finite-d numerator, integrated with unit physical cutoffs and the q^(ε−2)r^(ε−2) measure, gives corner Laurent series 10/(3ε)−169/3. The mixed Taylor term starts at −14ε/3. There is no triple pole from this short-short corner. Multiplication by Γ(3ε−1), the Wick factor, and curvature conversion gives +5/(1728ε²) from this corner under these cutoffs.

The once-soft q-linear and r-linear terms are rational with fifth-order denominators in the effective spectator lengths. The runnable check prints both. In particular the q-linear coefficient contains a term proportional to r² at fixed rg, so the spectator fraction endpoint is nonuniform. The absence of a corner triple pole does not determine the full graph triple pole. The short/long overlaps and FFff(B) still have to be integrated. The combined FFff higher-pole budget remains 1/36 and 1/72−L/12 from the counterterm calculation.

## Runnable check

Python 3 and SymPy.

```python
import sympy as s
q,r,f,g,v,d=s.symbols('q r f g v d')
P,Q=q*f,r*g
S=1+P+Q
t=[0,P,P+v,P+v+Q]
C=s.Matrix(4,4,lambda i,j:2*t[min(i,j)]*(S-t[max(i,j)])/S)
a=C[:,1]; b=C[:,3]-C[:,2]
c=2*P*(S-P)/S; e=2*Q*(S-Q)/S; j=-2*P*Q/S
N=d*(d*(S-P)*(S-Q)+2*P*Q)/(4*S*S)
T=(d*C*N+(d*(a*a.T*e+b*b.T*c)+2*j*(a*b.T+b*a.T))*d/(8*P*Q)).applyfunc(s.cancel)
A=-(4-d)*s.trace(T)/2-d*(4-d)*N/3
B=s.Rational(2,3)*(q*T[1,1]+r*(T[2,2]+T[3,3]+T[2,3])+v*(T[1,1]+T[2,2]+T[1,2])+(1-v)*T[3,3]+N*((d-1)*(q*q*(1-2*f)+r*r*(1-2*g))+d*(v*v+(1-v)**2)))
k=6-3*d/2
W=1+q+r
H=S**(-d/2)*W**(1-k)*(A+(k-1)*B/W)
def at(expr): return s.factor(expr.subs({q:0,r:0},simultaneous=True))
def avg(expr):
 x,y=s.symbols('x y')
 return s.factor(s.integrate(s.integrate(s.integrate(s.expand(expr.subs({f:x*(1-x),g:y*(1-y)})),(x,0,1)),(y,0,1)),(v,0,1)))
print('B', [avg(at(h)) for h in (B.subs(d,4)/S**2,s.diff(B.subs(d,4)/S**2,q),s.diff(B.subs(d,4)/S**2,r),s.diff(B.subs(d,4)/S**2,q,r))])
vals=[avg(at(h)) for h in (H,s.diff(H,q),s.diff(H,r),s.diff(H,q,r))]
print('full',vals)
print('4d',[h.subs(d,4) for h in vals])
assert [h.subs(d,4) for h in vals]==[-s.Rational(32,3),0,-8,0]
ep=s.symbols("eps")
corner=vals[0]/(ep-1)**2+(vals[1]+vals[2])/(ep*(ep-1))+vals[3]/ep**2
print('corner projective Laurent:',s.series(corner.subs(d,4-2*ep),ep,0,1))
H4=H.subs(d,4)
print('qface=',s.factor(s.integrate(s.diff(H4,q).subs(q,0),(v,0,1))))
print('rface=',s.factor(s.integrate(s.diff(H4,r).subs(r,0),(v,0,1))))

```

- [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), Eq. (3.24).
