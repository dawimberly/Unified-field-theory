# FFZ interacting benchmark: two-loop flat-space checkpoint

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**DERIVED** This checkpoint records the FFZ interacting benchmark: two-loop flat-space checkpoint. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

## Scope and conventions

**DERIVED** Each group has two chiral fields with $W=mXY+yX^3/6$, four real scalars and a Dirac pair at squared mass $z=|m|^2$. Set $K=16\pi^2$, $h=|y|^2/K$, $d=4-2\epsilon$, and $L=\ln(z/\mu^2)$. Use the same external-metric dimensional continuation as the preceding checkpoint, including a baseline $\xi_X=\xi_Y=1/6$ in $d$ dimensions. The loop measure is normalized in the $\overline{\mathrm{MS}}$ convention. No additional finite gravitational term is included in the quoted matter thresholds.

**DERIVED** Let $a_X=1-6\xi_X$, $a_Y=1-6\xi_Y$ and $w=2+2a_X+2a_Y$. At matching, the three masses are $z_g=M^2(1,2,4)$ and $w_g=(2,-3,1)$. The previously established equal-Yukawa family is $a_{X,g}=t w_g$ and $a_{Y,g}=(w_g-2)/2-tw_g$. The dimensions and weights here are those of the external-metric model. A full supergravity completion remains an additional physical assumption.

## Independent flat-space calculation of the bare $R$ term

**DERIVED** Work at one vertex at the origin with $g_{\mu\nu}=e^{2\sigma}\delta_{\mu\nu}$, $\sigma(x)=-\kappa x^2$, and $\kappa=R/[4d(d-1)]$. This local Taylor expansion probes the coefficient linear in constant scalar curvature. All propagator corrections below come from flat resolvents with slowly varying sources. The massive derivative expansion uses analytic continuation in $d$.

**DERIVED** After the Weyl field redefinitions, the scalar squared mass is $z e^{2\sigma}$, and its extra curvature source is $(1/6-\xi_d)R$ at this order, where $\xi_d=(d-2)/[4(d-1)]$. The fermion mass is $m e^\sigma$. The cubic, quartic and Yukawa vertices pick up factors of $e^{(1+\epsilon)\sigma}$, $e^{2\epsilon\sigma}$ and $e^{\epsilon\sigma}$, respectively. For momentum $k$, put $G=(k^2+z)^{-1}$, $S=(m-i\gamma\cdot k)G$. The first-order resolvents with one endpoint at the origin are
$$

\delta G/\kappa=-2zG\Delta_kG-\frac{d(4-d)}{3}G^2,\qquad
\delta S/\kappa=-m S\Delta_k S.

$$
**DERIVED** The second term in $\delta G$ is evanescent. The vertex factors and the $d$-dependent normalization are kept before pole subtraction.

**DERIVED** Write $G(r)$ for the position-space scalar propagator and $b(r)=-G'(r)$ for the flat radial vector part of the spinor propagator. Fourier transformation yields
$$

\frac{\delta G}{\kappa}=\frac{4-d}{6}r^2G-\frac{r^3}{3}G',\qquad
\frac{\delta b}{\kappa}=-\frac{d-1}{6}r^2b-\frac{r^3}{3}b'
+\frac{d(d-1)}{3}\partial_z b+\frac{d-1}{3}rG.

$$
**DERIVED** The full scalar and spinor resolvents were cross-checked symbolically in arbitrary $d$. The scalar figure-eight has no geometric $R$ contribution because the conformal scalar coincident variation vanishes. The fermion diagram uses the chiral kinetic contraction. Its flat normalization is fixed by Martin's $f_{FFS}$ formula.

**DERIVED** For $I=\int G^3$ and $J=\int r^2G^3$, the cubic sunset and the fermion–scalar sunset give
$$

F_{R,\mathrm{cubic}}^{(2)}=-\frac{|y|^2zJ}{24d},\qquad
F_{R,\mathrm{Yukawa}}^{(2)}=
|y|^2\left(\frac{zJ}{24d}+\frac{z I'}{36}\right).

$$
**DERIVED** One way to check the factors is to expand the product of propagators and vertices. The cubic product reduces to $(d-1)J/3$ in units of $\kappa$. In the Yukawa product, define $P=\int r^2G(G')^2$, $Q=\int GG'\partial_zG'$ with the derivative on one fermion line, and $T=\int rG^2G'$. The insertion is $P/(12d)+Q/6-T/(6d)$ per $R$ in the integrand. The identities $P=-zJ/2+dI/3$, $Q=-I/2-zI'/6$, and $T=-dI/3$ reduce it to $-zJ/(24d)-zI'/36$. The overall Yukawa Wick factor is $-|y|^2$. Thus the independent flat-source reduction gives the same bare $zI'/36$ as the corrected curved-propagator derivation.

## Group poles and finite Newton thresholds

**DERIVED** Use $I_{\mathrm{bare}}K^2/z=-3/(2\epsilon^2)+(3L-9/2)/\epsilon+O(1)$, $zI'_{\mathrm{bare}}=(1-2\epsilon)I_{\mathrm{bare}}$, and $\delta z=hz/(2\epsilon)$. The baseline remaining local poles, after the one-loop matter insertion, are $hz/K$ times $1/(24\epsilon^2)-1/(24\epsilon)$. The finite single log is $L/12$ in the same units. No $\ln(z)/\epsilon$ term remains.

**DERIVED** An independent soft source for $|X|^2$ or $|Y|^2$ differentiates the flat vacuum integral at fixed cubic vertex and fixed fermion mass. With $x,u$ the squared masses in the $X,Y$ scalar propagators, the renormalized two-loop flat vacuum term in units of $|y|^2/(2K^2)$ is
$$

\mathcal V(x,u)=[J_1(x)-J_1(z)]^2+(2z-x)I_R(z,z,x)-zI_R(x,x,u),
\quad J_1(x)=x[\ln(x/\mu^2)-1].

$$
**DERIVED** At $x=u=z$, $\mathcal V_u=-zI_R'/3$ and $\mathcal V_x=-I_R-zI_R'/3$, where $I_R'$ denotes the common-mass derivative. This gives a direct check of the relative X and Y improvement contributions. The entire group pole after including $\delta a_X=3ha_X/(2\epsilon)$ is
$$

\frac{hz}{K}\frac{w+6a_X}{48}
\left(\frac{1}{\epsilon^2}-\frac{1}{\epsilon}\right).

$$
**DERIVED** It cancels in the equal-Yukawa family because $\sum z w=\sum z a_X=0$. The aggregate finite result is independent of the common single-scale constants. At $\mu=M$,
$$

F_R^{(2)}(t)=\frac{hM^2}{K}\left[
-\frac5{24}(\ln2)^2+\frac{\ln2}{12}
+t\left(\ln2-\frac54(\ln2)^2\right)\right].

$$
**DERIVED** The normalized thresholds for $t=0,1$ are $-0.042332112853, +0.050248800309$. Both assignments have the same one-loop induced $G$. The interacting matter thresholds distinguish them.

**DERIVED** For unequal Yukawas, the pole cancellation at a single scale is the weighted condition $\sum h_gz_g(w_g+6a_{X,g})=0$. Over an open interval, the one-loop solutions have $h_g=h_{g0}/(1-3h_{g0}\ell)$, $a_{X,g}=a_{X,g0}/(1-3h_{g0}\ell)$ and $z_g=z_{g0}(1-3h_{g0}\ell)^{-1/3}$, where $\ell=\ln(\mu/\mu_0)$. Group the terms in $\sum z_gw_g(\ell)$ by distinct nonzero $h_{g0}$. Analytic continuation to each group's branch point shows that the coefficients of the powers $-1/3$ and $-4/3$ must separately vanish for that group. In particular its initial $\sum z_gw_g$ must be zero. No nonempty proper subset of the three weights $M^2(2,-6,4)$ sums to zero. Hence this exact one-loop trajectory requires equal Yukawa magnitudes across the groups, apart from the all-free case. This conclusion is restricted to the specified independent Wess–Zumino groups and mass-independent flow above thresholds.

## Two-loop $R^2$ sector quadratic in the improvements

**DERIVED** Each improvement insertion brings one factor of $R$. Their quadratic local $R^2$ coefficient therefore uses the zero-momentum Hessian of $\mathcal V$, with $(\delta x,\delta u)=-R(a_X,a_Y)/6$. Write
$$

I_R(z,z,z)=z[-3L^\frac{2}{2}+6L+C],\quad
C=-\frac{15}{2}+2\sqrt3\,\mathrm{Cl}_2(\\frac{pi}{3}),\quad
D=-2(C+6)/3.

$$
**DERIVED** The Hessian is
$$

H=\begin{pmatrix}
3L^2-L+D-1 & 2(D-1)/3\\
2(D-1)/3 & L-1/3-2D/3
\end{pmatrix}.

$$
**DERIVED** The improvement-degree-two term is $h(a_X,a_Y)H(a_X,a_Y)^T/(144K)$. It includes the logarithmic and finite terms in this sector. The degree-one and degree-zero terms of the general $R^2$ response are outside this result.

**DERIVED** For an exact derivative check, the bare integral identity
$$

(\partial_1^2-\partial_1\partial_2)I_{\rm bare}
=-\frac{2(d-3)I_{\rm bare}}{3z^2}-\frac{T_2^2}{z},
\qquad T_2=\int (k^2+z)^{-2}

$$
**DERIVED** follows from the $p$-loop integration by parts with $p_\mu$ and the scaling of $I_{\rm bare}$. Subtracting the matching one-loop subdivergence gives $z(\partial_1^2-\partial_1\partial_2) I_R = D-L$. Sector quadrature gave $D=-1.343907238689459$, and symmetric finite differences of Martin's closed-form integral at masses $1,2,4$ reproduced $D-L$.

**DERIVED** In the finite difference $[F_{R^2}(t)+F_{R^2}(-t)]/2-F_{R^2}(0)$, degree-zero and degree-one terms cancel when all assignments share a pure-metric subtraction condition. The two-loop matter-loop result at $t=1,\mu=M$ is $0.26078033947543\,h/K$. To reproduce it, use $(1,-1)H(1,-1)^T=3L^2-D$ and weights $w_g^2=(4,9,1)$. The scale derivative from the two-loop term plus the one-loop running leaves $h a_X^2/(36K)$ per group in the improvement-quadratic sector. It vanishes at $t=0$. This calculation isolates a definite coefficient in a common convention. A finite pure metric $R^2$ term is an independent matching input.

## Reproduction

**DERIVED** An independent Laurent check of the improvement-quadratic Hessian leaves $h a_X^2(1/\epsilon-3/\epsilon^2)/(144K)$ after the one-loop matter counterterms. The bare Hessian poles, in the same normalization, are $3/\epsilon^2-6L/\epsilon$ in the $XX$ entry, zero in $XY$, and $-1/\epsilon$ in $YY$. The counterterm cancels the $L/\epsilon$ terms and removes the $YY$ pole. This agrees with the finite-sector scale derivative above. At $t=0$ the order-$h$ poles in this quadratic sector vanish. The result leaves the degree-zero and degree-one order-$h$ sectors open.

**DERIVED** The leading scalar $R^2$ divergence is already nonzero. The independent $R^2$ part of the one-loop heat-kernel coefficient is $\sum_i a_i^2/72$ before the determinant prefactor. The Dirac contribution has Weyl-squared, Euler and total-derivative terms. In the surviving family the sum of squared improvements is $13-28t+56t^2$, with minimum $19/2$ at $t=1/4$. Even arbitrary splits at fixed FFZ group sums obey the lower bound $(0^2+5^2+1^2)/4=13/2$. This closes the counterterm-free scalar-response branch in the specified weak-coupling spectrum. The trace sum rule probes a conformal-endpoint difference and does not remove this divergence.

**DERIVED** Run this block with Python, SymPy, SciPy and mpmath.

```python
from itertools import combinations
from math import exp, expm1, gamma, log, log1p, pi

import mpmath as mp
import sympy as sp
from scipy.integrate import dblquad, quad
from scipy.special import kv

p, m, d = sp.symbols("p m d")
G = 1 / (p + m*m)

def lap(f, dim):
    return 4*p*sp.diff(f, p, 2) + 2*dim*sp.diff(f, p)

def fourier_slope(f, order, dim):
    return lap((d + order)*f + 2*p*sp.diff(f, p), dim)/3

def check(expr):
    assert sp.factor(expr) == 0, sp.factor(expr)

A, B = m*G, -G
Al, Bl = lap(A, d), lap(B, d+2)
check(-m*(A*Al-p*B*Bl) - (
    (d-1)*Al/6 - fourier_slope(A, 0, d) - d*(d-1)*m*G**2/3
))
check(-m*(A*Bl+B*Al) - (
    (d-1)*Bl/6 - fourier_slope(B, 1, d+2)
    + (d-2)*(d-1)*G**2/3
))
check(-2*m*m*G*lap(G, d)-d*(4-d)*G**2/3 - (
    -(4-d)*lap(G, d)/6 - fourier_slope(G, 0, d)
))

e, L, w, x, y, C = sp.symbols("e L w x y C")
sunset_poles = -sp.Rational(3, 2)/e**2+(3*L-sp.Rational(9, 2))/e
bare = (w*(1-2*e)/72+x/12)*sunset_poles
counterterm = w*(1/e**2-L/e)/24+x*(1/e**2+(1-L)/e)/4
poles = sp.expand(bare+counterterm)
check(poles.coeff(e, -2) - (w+6*x)/48)
check(poles.coeff(e, -1) + (w+6*x)/48)

D = -2*(C+6)/3
H = sp.Matrix([
    [3*L**2-L+D-1, 2*(D-1)/3],
    [2*(D-1)/3, L-sp.Rational(1,3)-2*D/3],
])
a = sp.Matrix([x, y])
v = sp.Matrix([1, -1])
check((v.T*H*v)[0] - (3*L**2-D))
check((6*x*x*L+x*x+y*y)/72 - 2*sp.diff((a.T*H*a)[0]/144, L) - x*x/36)
assert all(sum(c) for n in (1, 2) for c in combinations((2,-6,4),n))

def pole_part(expression):
    series = sp.expand(expression)
    return series.coeff(e,-2)/e**2+series.coeff(e,-1)/e

sunset_first = pole_part((1-2*e)*sunset_poles)
sunset_second = pole_part((1-2*e)*(-2*e)*sunset_poles)
tadpole_second = 1/e-L
derivative_difference = pole_part(-2*sunset_first/3-tadpole_second**2)
derivative_diagonal = sunset_second/9+2*derivative_difference/3
derivative_off_diagonal = derivative_diagonal-derivative_difference
bare_hessian = sp.Matrix([
    [pole_part(2*tadpole_second**2-2*sunset_first/3-sunset_second/3),
     -2*derivative_off_diagonal],
    [-2*derivative_off_diagonal,-derivative_diagonal],
])
curvature_counterterm = -6*x*x*tadpole_second/e+(x*x+y*y)*tadpole_second
check(pole_part((a.T*bare_hessian*a)[0]+curvature_counterterm)
      - x*x*(1/e-3/e**2))

def sector_integrand(t, u):
    scale = 1+u*(1+t)
    V = u*(1+t+u*t)
    return 2*((scale*scale-3*V)/scale/(1+t+u*t)**2-1/(1+t)**2)/u

D_quad = 1-log(2)+dblquad(
    sector_integrand, 0, 1, 0, 1, epsrel=2e-12, epsabs=2e-12
)[0]
mp.mp.dps = 40
C_num = -mp.mpf(15)/2 + 2*mp.sqrt(3)*mp.im(mp.polylog(2,mp.exp(1j*mp.pi/3)))
D_exact = -2*(C_num+6)/3
assert abs(D_quad-D_exact) < 1e-10
ell = mp.log(2)
newton_0 = -5*ell**2/24 + ell/12
newton_1 = newton_0 + ell - 5*ell**2/4
curvature_even = (39*ell**2-14*D_exact)/144
print("Newton thresholds:", float(newton_0), float(newton_1))
print("Sunset derivative constant:", float(D_exact))
print("Even curvature threshold:", float(curvature_even))

def thresholds_at_epsilon(eps):
    euler_gamma = 0.5772156649015329

    def subtracted_sector(v, u):
        w = 1+v
        ratio_log = (1-2*eps)*log1p(u*w)+(eps-2)*log1p(u*v/w)
        return 6*u**(eps-1)*w**(eps-2)*expm1(ratio_log)

    integral = 6*expm1((eps-1)*log(2))/(eps*(eps-1)) + dblquad(
        subtracted_sector, 0, 1, 0, 1, epsabs=1e-11, epsrel=1e-11
    )[0]
    sunset = exp(2*euler_gamma*eps)*gamma(-1+2*eps)*integral
    tadpole_derivative = exp(euler_gamma*eps)*gamma(eps)

    def moment(order):
        return sum(w*expm1(-order*eps*log(z))
                   for z,w in zip((1,2,4), (2,-6,4)))

    double_scale, single_scale = moment(2), moment(1)
    t0 = (1-2*eps)*sunset*double_scale/72 + tadpole_derivative*single_scale/(24*eps)
    difference = sunset*double_scale/12 + tadpole_derivative*single_scale/(4*eps*(1-eps))
    return t0, t0+difference

# Extrapolate the dimensionally continued species sums along a halving sequence.
samples = [thresholds_at_epsilon(0.015/2**j) for j in range(6)]
for power in range(1,len(samples)):
    scale = 2**power
    samples = [tuple((scale*b-a)/(scale-1) for a,b in zip(left,right))
               for left,right in zip(samples,samples[1:])]
assert abs(samples[0][0]-newton_0) < 5e-9
assert abs(samples[0][1]-newton_1) < 5e-9
print("Bare-integral extrapolation:", samples[0])

def radial_integrals(dim, mass_squared):
    mass = mass_squared**0.5
    nu = dim/2-1
    sphere_area = 2*pi**(dim/2)/gamma(dim/2)
    eps = (4-dim)/2

    def integrand(r, diagram):
        scale = (mass/r)**nu/(2*pi)**(dim/2)
        scalar = scale*kv(nu,mass*r)
        vector = mass*scale*kv(nu+1,mass*r)
        vector_prime = -(dim-1)*vector/r-mass_squared*scalar
        scalar_variation = (4-dim)*r*r*scalar/6+r**3*vector/3
        # The Bessel identity gives partial_z vector = -r*scalar/2.
        vector_variation = (-(dim-1)*r*r*vector/6-r**3*vector_prime/3
                            +(dim-1)*(2-dim)*r*scalar/6)
        scalar_sunset = (3*scalar**2*scalar_variation
                         -(1+eps)*r*r*scalar**3)
        fermion_sunset = (scalar_variation*vector**2
                          +2*scalar*vector*vector_variation
                          -eps*r*r*scalar*vector**2)
        return sphere_area*r**(dim-1)*(scalar**3, scalar_sunset, fermion_sunset)[diagram]

    sunset, scalar_part, fermion_part = (
        quad(integrand,0,100/mass,args=(diagram,),epsabs=1e-10,epsrel=1e-10)[0]
        for diagram in range(3)
    )
    direct = (-mass_squared*scalar_part/2-fermion_part)/(4*dim*(dim-1))
    expected = (dim-3)*sunset/36
    return abs(direct-expected)/abs(expected)

max_radial_error = max(radial_integrals(dim,z) for dim in (1.5,2,2.3,2.7)
                       for z in (0.7,1,4))
assert max_radial_error < 1e-8
print("Maximum relative radial check error:", max_radial_error)
```

## Limits and primary references

**DERIVED** The next higher-order Newton diagnostic is the inhomogeneous curvature beta function for $X$. Since the dimensionless interactions are entirely in $X$, a mass-independent subtraction has no corresponding mixing into the curvature coupling of $Y$ in the specified diagonal source sector. With a common Yukawa, a homogeneous $\dot a_X$ preserves both null sums. An additive term $\dot a_X=\gamma a_X+b(h)$ gives $\frac{d}{d\ln\mu}\sum_g w_g=6b(h)$ and $\frac{d}{d\ln\mu}\sum_g z_gw_g=14M^2b(h)$ at the matching point. These identities identify the next calculation in the same scheme. They do not supply $b(h)$ or rule out interaction-corrected matching conditions.

**DERIVED** The independent verification is a flat-resolvent calculation on a quadratic Weyl source. It is a local derivative expansion, so it does not reconstruct the full finite-momentum metric correlator. The benchmark uses the specified $d$-dimensional continuation and matter subtraction. A different finite local gravitational matching changes absolute coefficients. There is no claim here of an interacting UV fixed point, a dynamical supergravity completion, or a derived physical scalar spectrum for constrained CEV.

- **DERIVED** [Frolov, Fursaev and Zelnikov, Statistical Origin of Black Hole Entropy in Induced Gravity](https://arxiv.org/abs/hep-th/9607104)
- **DERIVED** [Martin, Two-loop effective potential for a general renormalizable theory and softly broken supersymmetry](https://arxiv.org/abs/hep-ph/0111209), Eqs. (2.17)–(2.21), (4.12)–(4.15), (7.7)–(7.9).
- **DERIVED** [Mück, Spinor parallel propagator and Green function in maximally symmetric spaces](https://arxiv.org/abs/hep-th/9912059), Eqs. (10), (14) and (15).
- **DERIVED** [Vassilevich, Heat kernel expansion: user's manual](https://arxiv.org/abs/hep-th/0306138), Eq. (4.28).
- **DERIVED** [Prochazka and Zwicky, Flow of the $\Box R$ Weyl anomaly](https://arxiv.org/abs/1703.01239). The endpoint trace sum rule does not fix the absolute pure-metric contact.
