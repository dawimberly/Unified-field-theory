# FFZ finite two-loop even curvature slope

Continues [the fixed-1/6 checkpoint](ffz_fixed_one_sixth_checkpoint.md) and [the original two-loop benchmark](ffz_two_loop_flat_space_checkpoint.md). The prior scheme-qualified obstruction to an additive improvement beta and its two Newton moment consequences carry forward. The higher-loop Newton vacuum simple pole is independent and remains open.

## Finite slope in the specified matter subtraction

For one group use h=|y|²/K, K=16π², z=m², L=ln(z/μ²), and v=tw_g. Let U_X=−vR/6 and U_Y=−U_X. Extract the even family difference in the coefficient of R(q)R(−q). The order-h coefficient of q² is

P₁ = h v²/(216 K z) [L − 0.307317873115].

The three-group sum at μ=M is +0.00359585422706 h t²/(KM²). This is the matter threshold with a common pure-metric subtraction and no independent t-dependent R□R term. The result concerns the Euclidean expansion in q² about zero. It provides a sign and normalization check for subsequent calculations of the nonlocal form factor.

## Derivation

Normalize the two-loop functional by |y|²/2, as in the parent and Martin's scalar and chiral-fermion vacuum diagrams. Let G be the free scalar Green function of mass z and T its coincident limit. Let G_X,G_Y be resolvents with local potentials U_X,U_Y. Up to a source-independent constant the functional reduces to

∫dx (G_X(x,x)−T)² + ∫dxdy [(z−U_X(x))G(x,y)²G_X(x,y)−z G_X(x,y)²G_Y(x,y)].

The chiral fermion numerator is −2∇G·∇G in this normalization. Use 2∇G·∇G=□(G²)−2zG²+2Tδ and integrate □ onto G_X at one endpoint. The resolvent equation produces the endpoint source in the displayed functional. For constant sources the formula reduces to the parent's two-loop potential before subdivergence subtraction. The complete expression keeps the local endpoint insertion at finite momentum.

Expand the resolvents with two opposite-momentum source insertions. In the sunset part the coefficients are −2z per same-line kernel, +z for one distinct-line kernel, and +1 for the kernel with the explicit source at an endpoint. Equal masses allow one representative line or pair for each coefficient. The factors come from 1−3 for same-line insertions and the negative of 1−1−1 for distinct lines. There is no extra factor of two from averaging the Fourier insertions. At q=0 a same-line expansion has Schwinger weight a²/2, a distinct pair has ab, and an endpoint-source expansion has a.

For Schwinger lengths a,b,c, D=ab+ac+bc. The q² slope is minus the effective resistance between the insertion points times the zero-momentum Gaussian. For two points separated by l on line a the resistance is l−l²(b+c)/D. For points distances l and r from their common vertex along lines a,b it is l+r−[(b+c)l²+(a+c)r²+2clr]/D. The latter also follows by inverting the two-loop quadratic form with entries a+c,b+c and off-diagonal c after routing the momentum on the two segments with opposite signs. Integrating these polynomials gives the functions in the runnable check. The endpoint kernel uses the same-line resistance with one point fixed at the endpoint. The figure-eight square gives the derivative of the squared one-loop bubble, −Γ(ε)Γ(1+ε)/(3z) with common scale factors suppressed.

## Endpoint subtraction and cross-check

Scale out total Schwinger length. With dimensional loop measures in the parent convention, the sunset integral has Γ(1+2ε) times the projective kernel [(1+2ε)(2S−D_pair)/scale−E]/scale^(1+2ε), multiplied by the determinant power. S,D_pair,E are the integrated same-line, distinct-line and endpoint resistances. Use the six ordered sectors (1,x,xt). Only sectors with line a largest have a nonzero leading endpoint at ε≠0. Their coefficient is (ε/3)(1+t)^(−2+ε), multiplying x^(ε−1). Its integrated limit is 1/3. This is a finite evanescent term. The remaining sunset integral is ordinary quadrature at ε=0. It plus the endpoint term gives 0.06422737562836.

In units h times the one-loop slope the total bare pole is −1/ε from the figure-eight square. The explicit matter counterterm is +1/ε−1/2 at L=0. The finite −1/2 comes from the derivative of z^(−1−ε) in the mass counterterm. The evanescent sunset and remaining integral give the finite constant 3(0.06422737562836)−1/2. The bare and counterterm factors exp(−2εL) and exp(−εL) give the +L term, independently reproducing the earlier RG check. Numerical dimensional extrapolation below checks the finite subtraction. It shares the diagrammatic representation with the sector result.

## Runnable check

Python 3 with SciPy.

```python
from itertools import permutations
from scipy.integrate import dblquad, quad
from scipy.special import gamma
from math import exp, log
perms=list(permutations(range(3)))
def kernel(a,b,c):
    d=a*b+a*c+b*c
    same=a**3/6-a**4*(b+c)/(12*d)
    diff=a*b*(a+b)/2-(b*a**3*(b+c)/3+a*b**3*(a+c)/3+c*a*a*b*b/2)/d
    end=a*a/2-a**3*(b+c)/(3*d)
    return d,2*same-diff,end
def integral(e):
    def integrand(t,x):
        v=(1,x,x*t)
        result=0
        for i,j,k in perms:
            a,b,c=(v[i],v[j],v[k])
            d,s,p=kernel(a,b,c)
            scale=a+b+c
            result+=(d/x)**(-2+e)*(((1+2*e)*s/scale-p)/scale**(1+2*e)
                )-((e/3)*(1+t)**(-2+e) if i==0 else 0)
        return x**(e-1)*result
    surface=2/3*quad(lambda t:(1+t)**(-2+e),0,1)[0]
    return surface+dblquad(integrand,0,1,0,1,epsabs=3e-10)[0]
i=integral(0)
c=3*i-.5
print(i,c)
for e in [.02,.01,.005,.0025]:
    normalization=exp(2*.5772156649015329*e)
    bare=3*normalization*(gamma(1+2*e)*integral(e)-gamma(e)*gamma(1+e)/3)
    ct=exp(.5772156649015329*e)*gamma(1+e)*(1/e-.5)
    print(e,bare+ct)
print('species', (5*log(2)+(35/4)*c)/216)
assert abs(c+0.3073178731149105)<1e-10
def finite(e,L):
    eg=.5772156649015329
    return 3*exp(2*e*(eg-L))*(gamma(1+2*e)*integral(e)
           -gamma(e)*gamma(1+e)/3) + exp(e*(eg-L))*gamma(1+e)*(1/e-.5)
for L in (0,log(2),log(4)):
    seq=[finite(.02/2**j,L) for j in range(5)]
    for order in range(1,5):
        seq=[(2**order*b-a)/(2**order-1) for a,b in zip(seq,seq[1:])]
    assert abs(seq[0]-(c+L))<2e-7
```

## Next calculation

Use the same momentum-routed representation at finite q² to examine the common-contact-subtracted even response and locate its first spectral threshold. A separate independent check can use direct momentum-space Taylor differentiation and massive vacuum integral reduction. For the Newton branch, compute the order-h² pure-metric R simple pole or its improvement-dependent part and test its mass moments. The present result only addresses curvature degree two and order h.

## Primary papers

- [Martin, Two-loop effective potential for a general renormalizable theory and softly broken supersymmetry](https://arxiv.org/abs/hep-ph/0111209), Eqs. (4.12)–(4.15) and the Wess–Zumino check in Sec. 8.1.
- [Codello and Zanusso, On the non-local heat kernel expansion](https://arxiv.org/abs/1203.2034), Eqs. (4.10)–(4.12) for the independent one-loop normalization.
