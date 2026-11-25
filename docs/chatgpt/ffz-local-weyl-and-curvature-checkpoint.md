# FFZ local Weyl and curvature checkpoint

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**ASSUMED** This checkpoint records the FFZ local Weyl and curvature checkpoint. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Parent: [ffz_two_loop_flat_space_checkpoint.md, version 2](ffz-two-loop-flat-space-checkpoint.md). The Physics folder listing exposed the parent as its most recent FFZ checkpoint. Its threshold and Hessian calculations are carried forward without rerunning their scans.

## Additive improvement beta

**ASSUMED** Use the external-metric theory with $K=16\pi ^{2}, h=$|y|²/K, and $a=1-6\xi$. In mass-independent subtraction the dimensionless source problem reduces to the massless one-chiral X sector. Y has no such vertices. Keep a fixed Yukawa phase. There is no dynamical supergravity assumption.

**DERIVED** Jack and Osborn, Nuclear Physics B 883 (2014), Eqs. (2.47), (2.53), (2.59), express the additive curvature source in terms of the local-coupling source coefficients:
**DERIVED** $\eta =\delta _I B^I-(L_B-\gamma _M)\tau$. The shift by R times a finite function of coupling permits $\tau =0$. For the supersymmetric flat local-source prescription, Eq. (7.5) gives $\delta =0$. Thus the induced additive curvature-source beta vanishes in the associated $\tau =0$ external-metric improvement prescription. This is an all-orders statement at the level of that formal local RG prescription and its supersymmetric flat counterterms. It is conditional on an analytic connection to the benchmark's perturbative matter sector when used to compare schemes.

**DERIVED** Component check in the single-field specialization. Eq. (9.16) directly gives $\delta ^{1}=0$. In Eq. (9.26), take $dy_i=y_i\,d\ln|y|$ and $d\bar y_i=\bar y_i\,d\ln|y|$. With the one-half trace weight of Eq. (5.17), $T=\operatorname{tr}(\bar y y)=h/2$ and $\lambda =h$. The three radial contributions to $\delta ^{2}$ are $-2\lambda T, +Th, +Th$. They cancel. The radial $\beta$ vector therefore produces neither an h² nor an $h^{3} \eta$ in the $\tau =0$ prescription. This component check uses independent diagrammatic coefficients from the same paper as the superspace argument.

**DERIVED** The earlier benchmark holds $\xi =1/6$ at finite d. Its literal pole-only curvature coupling has not been converted to $\tau =0$. Since $\xi _d=(d-2)/[4(d-1)]$ differs from $1/6$ by $-\epsilon /18+O(\epsilon ^{2})$ at $d=4-2\epsilon$, finite source matching requires tracking evanescent terms through subdivergence subtraction. It is not justified to copy the $\tau =0$ cubic coefficient into that scheme.

**DERIVED** To see what survives a finite conversion, write $\beta _h=3h^{2}+u h^{3}+$…, $\Gamma =3h+v h^{2}+$…, and $\beta _a=\Gamma a+b$. Under $a'=a+f(h)$, $b'=b+\beta_h f'-\Gamma f$. For $f=c1 h+c2 h^{2}$, the shift at h² is zero and at h³ it is $(u-v)c1+3c2$. Thus $b_{2}=0$ across these analytic additive prescriptions. More general analytic nonzero multiplicative source normalizations preserve the vanishing. A nonzero cubic coefficient in the original scheme would be removable. Above quadratic order, 3(n−1) multiplies the new coefficient of $h^n$ in f, allowing higher additive terms to be removed recursively as formal series. By contrast, a hypothetical b₂ would force an h ln h term in an invariant offset and could not be removed analytically.

## Consequences for both Newton constraints

**DERIVED** Write $S_j=\sum z_g^j w_g, A_j=\sum z_g^j a_Xg, P_j=\sum z_g^j, j=0,1$. For common h and common mass anomalous rate $\gamma_z$, $S_j'=j\gamma_z S_j+2\Gamma A_j+2bP_j$ and $A_j'=(j\gamma_z+\Gamma)A_j+bP_j$. On the parent's family at matching, the raw slopes are 6b and 14M²b. Both vanish through O(h²), and in the $\tau =0$ coordinates the original algebraic family expressed in those coordinates stays null under homogeneous flow.

**DERIVED** For the same physical sources under $a'=a+f$, the pulled-back constraints are $S_j=S'_j-2fP_j$ and $A_j=A'_j-fP_j$, where primes here label coordinates. Their flows follow the old equations by the chain rule. Conversely, if a is the benchmark coordinate and $a'=a-s(h)$ is homogeneous, the homogeneous null family is specified by $S_j-2sP_j=A_j-sP_j=0$, with $\beta_hs'=\Gamma s+b$. This equation defines a candidate corrected family. Establishing that it has vanishing interacting Newton poles requires matching the local gravitational coefficient and matter sources in the same convention. Setting raw sums to zero separately in two conventions changes the physical curvature vertices. For example, the one-loop complex-scalar threshold from a shift $\Delta a$ is $-\Delta a z(L-1)/(6K)$; reexpressing the same vertex transfers the opposite term into the matching formula. Neither the vanishing of b in some coordinates nor its appearance in others computes the higher-loop pure-metric R counterterm. The finite Newton threshold and its local gravitational matching remain separate.

## Independent curvature-squared form factor

**DERIVED** For the even family difference $[F(t)+F(-t)]/2-F(0)$, both insertions must be the potential shifts $-R \Delta a/6$. In each group the complex scalar determinant yields $-\operatorname{Tr} (GU GU)/2$, with $\Delta a_X=tw_g$ and $\Delta a_Y=-tw_g$. All mixed geometric terms cancel in the even difference at curvature degree two. The finite Euclidean form factor at common subtraction is

**DERIVED** $F_{even}(q^{2})=t^{2}/(36K) \sum w_g^{2} \int _{0}¹ dx \ln [(z_g+x(1-x)q^{2})/\mu ^{2}]$.

**DERIVED** Independent route: the heat-kernel potential form factor in Codello–Zanusso Eqs. (4.10)–(4.12) is $(1/2)\int dx \exp [-s x(1-x)q^{2}]$. With the complex scalar proper-time factor −1 and $\int ds/s$, the momentum-subtracted kernel gives the same logarithm. The local term at $\mu =M$ is $t^{2}\cdot 11 \ln 2/(36K)$, and the derivative at $q^{2}=0$ is $t^{2}\cdot 35/(864KM^{2})$. The coefficient of the determinant pole is $-14t^{2}/(36K \epsilon )$. A common local R² contact changes the absolute constant; the subtracted momentum response is unaffected by an R² contact. Independent higher-derivative gravitational coefficients, if admitted, contribute to analytic momentum powers. The formula is a one-loop result at nonzero momentum. The parent already gives the order-h local even term, $0.26078033947543 h/K$ at $t=1$. Its order-h momentum dependence remains open.

## Runnable verification

**DERIVED** Extract the Python block and run with Python 3, SymPy, and mpmath. This checks the component contraction, finite scheme shift, null-sum covariance, and finite-momentum result against proper time. It uses no paid compute.

```python
import math
import sympy as sp
import mpmath as mp

h, u, v, c1, c2, b, a, w, p, j, gamma_z = sp.symbols(
    "h u v c1 c2 b a w p j gamma_z"
)
trace = h / 2
assert sp.expand(-2*h*trace + h*trace + h*trace) == 0
beta_h = 3*h**2 + u*h**3
Gamma = 3*h + v*h**2
f = c1*h + c2*h**2
shift = sp.expand(beta_h*sp.diff(f, h) - Gamma*f)
assert shift.coeff(h, 2) == 0
assert sp.expand(shift.coeff(h, 3) - ((u-v)*c1+3*c2)) == 0
S_new, A_new, b_new = w+2*f*p, a+f*p, b+shift
assert sp.expand(
    j*gamma_z*S_new+2*Gamma*A_new+2*b_new*p
    - (j*gamma_z*w+2*Gamma*a+2*b*p
       + 2*p*beta_h*sp.diff(f,h)+2*f*j*gamma_z*p)
) == 0
weights = (2,-3,1)
masses = (1,2,4)
assert sum(weights) == sum(z*w for z,w in zip(masses,weights)) == 0
assert sum(masses) == 7
assert sum(sp.Rational(w*w,z) for w,z in zip(weights,masses))/216 == sp.Rational(35,864)
assert math.isclose(sum(w*w*math.log(z) for w,z in zip(weights,masses)), 11*math.log(2))

mp.mp.dps = 25
for q2 in (mp.mpf("0.03"), mp.mpf("0.7"), mp.mpf("6")):
    for z in masses:
        bubble = mp.quad(lambda x: mp.log1p(q2*x*(1-x)/z), [0,1])
        proper_time = mp.quad(lambda s: mp.exp(-z*s)/s *
            mp.quad(lambda x: -mp.expm1(-s*q2*x*(1-x)), [0,1]), [0,1,mp.inf])
        assert abs(bubble-proper_time) < mp.mpf("1e-21")
print("FFZ source, constraints, and proper-time checks passed")
```

## Next calculation

**DERIVED** Convert the $\tau =0$ prescription into the checkpoint's fixed $1/6$ pole-only convention. Determine the finite coefficient of the $q^2\sigma(q)|X|^2$ vertex at orders h and h² with the same evanescent wavefunction and mass-source counterterms. Project away the kinetic equation-of-motion vertex at nonzero external Euclidean momenta. The two finite coefficients determine the benchmark b₃ through the formula above and permit a comparison to the direct curvature simple pole. Then test the interacting Newton counterterm in the matched coordinates. An independent next curvature test is the momentum slope of the order-h even term, obtained from the two-soft-source kernel with opposite nonzero momenta.

## Sources

- **CITED** [Jack and Osborn, Constraints on RG Flow for Four Dimensional Quantum Field Theories](https://arxiv.org/abs/1312.0428), Eqs. (2.47), (2.53), (2.59), (5.17), (7.5), (9.16), (9.26).
- **CITED** [Codello and Zanusso, On the non-local heat kernel expansion](https://arxiv.org/abs/1203.2034), Eqs. (4.10)–(4.12).
- **CITED** [Brown and Collins, Dimensional renormalization of scalar field theory in curved space-time](https://doi.org/10.1016/0003-4916(80)$90232-8$). — eq./page UNVERIFIED.
