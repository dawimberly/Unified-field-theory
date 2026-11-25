# FFZ four-Yukawa curvature kernels and the FFff(B) short-pair corner

**Conversation date range** 2026-11-22

**DERIVED** This checkpoint records the FFZ four-Yukawa curvature kernels and the FFff(B) short-pair corner. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** Continues [the CcFf pole checkpoint](ffz-ccff-regulated-poles.md). Equal masses, fixed ξ=1/6, d=4−2ε. This checkpoint evaluates one FFff(B) corner exactly.

## Gaussian numerator and curvature kernel

**DERIVED** Use cyclic edge directions around the four vertices. For FFff(A) two opposite edges are doubled fermion pairs and the bridges are scalar. For FFff(B) two opposite edges are doubled mixed scalar–fermion pairs and the bridges are fermion. In both cases the Gaussian covariance is twice the inverse anchored conductance matrix. Equivalently in cycle coordinates the per-component covariance of increments is K_ij=2(l_i δ_ij−l_i l_j/S), with effective lengths l and perimeter S. In a parallel pair with lengths s,t use l=st/(s+t).

**DERIVED** Write W for the spin polynomial with the actual Schwinger factors 1/(2s) stripped. Let K now denote the increment covariance within this formula only. The expectation for A is d²K_aa K_cc+2d K_ac², for the two parallel fermion pairs a,c. For B it is d²(K_ab K_cd+K_ad K_bc)+d(2−d)K_ac K_bd. The alternating Weyl trace has an antisymmetric part that vanishes pointwise because the four cyclic coordinate increments sum to zero. Reversing two increments to convert cyclic directions to the F→f convention preserves the four-vector sign. The Wick factors multiplying these polynomials are −|y|⁴/2 for A and +|y|⁴ for B, with the denominators ∏(2s_f) restored.

**DERIVED** For either polynomial expectation T(K), insert a quadratic coordinate moment using the unnormalized rule
**DERIVED** E[W x_i·x_j] = d C_ij T(K) + DT(K)[b_i b_jᵀ+b_j b_iᵀ],
**DERIVED** where C is the coordinate covariance and b_i contains Cov(Δ_e,x_i) for the increments in T. This follows by differentiating the covariance under a Gaussian tilt exp(τ x_i·x_j). The unnormalized rule avoids dividing by a signed trace expectation.

**DERIVED** Define M_ij as the moment divided by T for writing the kernel. The Weyl and curvature-mass weight is
**DERIVED** $A=-(4-d)Σ_i \frac{M_ii}{2}-d(4-d)Σ_scalar \frac{s_e}{3}.$
**DERIVED** The mass-profile weight is
**DERIVED** $B=(\frac{2}{3})Σ_e s_e(M_ii+M_jj+M_ij+d_e s_e),$
**DERIVED** where d_e=d for scalar lines and d−1 for fermion lines. In a gauge with Σu=1, the radial factor is Γ(k)z^(−k)[A+kB] with k=5−3d/2=−1+3ε. The four kinetic factors give degree −2 after Gaussian averaging. Together with the degree-one A insertion this fixes k. For effective lengths fixed by the sum of the two single-edge lengths, write W_total=Σu and use Γ(k)z^(−k)W_total^(−k)[A+kB/W_total]. Include the usual loop measure and divide by 4d(d−1) to convert κ to R.

## FFff(B) double short-pair corner

**DERIVED** Set cyclic effective lengths to (P,v,Q,w), w=1−v. Anchor vertex 0. In the double short limit vertices 0,1 coalesce and vertices 2,3 coalesce at X. Their covariance is C_XX=2vw, and b_X at leading order is (2Pw,2vw,−2Qv,−2vw). Every product of two off-diagonal K entries in T is 4PQvw/S². Applying the tilt gives
**DERIVED** $M_XX=(2d+8)vw-4d/(d+2).$
**DERIVED** At d=4, only the two surviving bridge edges contribute to the mass profile at leading order. They give B=(2/3)(10vw+1/3). A vanishes at d=4 and k=−1.

**DERIVED** For each mixed pair the measure contributes P^(ε−1)[a(1−a)]^(−ε), and analogously for Q. Choose a as the fraction of physical pair length assigned to the fermion, so P/s_f=1−a. At d=4 the spin expectation divided by ∏(2s_f) is 6(1−a)(1−b). The two fraction integrals give 1/4, and the v integral of B is 4/3. Multiplying by the radial residue −1/3, k=−1 and κ/R=1/48 gives **+1/(72 ε³)** in units h²z/K from this local double short corner, with both independent effective-length cutoffs at 1. This is only the leading pole of this localized contribution. Regions with long effective lengths relative to the bridge sum can also contribute to the full bare graph. The finite terms and lower poles of this corner have yet to be combined with those regions.

## Reproduction

**DERIVED** Python 3 standard library. The check obtains the tilt derivative directly from the six off-diagonal entries, compares with the closed conditional moment for multiple rational geometries and dimensions, and checks the integrated leading pole.

```python
from fractions import Fraction as F

for d in (F(3), F(7, 2), F(4), F(9, 2)):
    for j in range(1, 12):
        v, p, q = F(j, 12), F(1, 11), F(1, 13)
        w = 1 - v
        lengths = (p, v, q, w)
        b = (2*p*w, 2*v*w, -2*q*v, -2*v*w)
        # Leading short-pair entries: S=1 at this order.
        def entry(i, j):
            return -2*lengths[i]*lengths[j]
        def derivative(i, j):
            return 2*b[i]*b[j]
        trace = change = F(0)
        for i, j, m, n, coefficient in (
            (0, 1, 2, 3, d*d),
            (0, 3, 1, 2, d*d),
            (0, 2, 1, 3, d*(2-d)),
        ):
            trace += coefficient*entry(i,j)*entry(m,n)
            change += coefficient*(derivative(i,j)*entry(m,n)
                                    + entry(i,j)*derivative(m,n))
        assert trace == 4*p*q*v*w*d*(d+2)
        assert 2*d*v*w + change/trace == (2*d+8)*v*w - 4*d/(d+2)

profile_integral = F(2,3)*(10*F(1,6) + F(1,3))
assert profile_integral == F(4,3)
assert F(-1,3)*(-1)*F(1,48)*6*F(1,4)*profile_integral == F(1,72)
print("FFff(B) covariance tilt and localized triple pole verified")
```

**DERIVED** The next check is to partition the complement using ordered physical-length sectors and isolate the regions with one long bridge or one long edge of a mixed pair. For A, each shrinking fermion pair has P^(ε−2) power counting, so the linear face terms as well as constants enter the subtraction.

- **DERIVED** [Martin, Effective potential at three loops, Eqs. (3.19) and (3.24)](https://doi.org/10.1103/PhysRevD.96.096005).
- **DERIVED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177).
