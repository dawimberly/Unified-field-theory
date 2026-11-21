# FFZ three-loop Newton pole from the Y improvement

Continues [the projected Newton checkpoint](ffz_newton_projection_q4_checkpoint.md) and [the original flat-space checkpoint](ffz_two_loop_flat_space_checkpoint.md). The prior absence of additive matter-source running through h³ in the specified prescription is used here. The new result is a vacuum counterterm and has no implication of additional drift in the two raw moments.

## Y improvement residue

At curvature degree one the coefficient of a_Y can be computed at a flat metric with one soft mass insertion s_Y=−a_Y R/6. All curvature-dependent geometry belongs to a separate term. Set y^{XXX}=y, μ^{XY}=μ^{YX}=m, with conjugates as appropriate, and all soft holomorphic parameters and gauge couplings to zero. Set s_X=0. In Martin's three-loop vacuum beta (4.40), the only linear s_Y terms that survive are two contributions with the soft mass between two μ tensors. In order of appearance the surviving numerical coefficients in units |y|⁴ z s_Y are

- +2 from the μ s μ term in the first tensor group.
- −1/2 from the μ s μ term in the fourth group.
The two contractions immediately before the ζ₃ terms vanish because one μ in each has both indices X.

The ζ₃ contractions vanish because they require a μ with both indices X. Hence β_Λ at order h², linear in s_Y, equals +(3/2) h² z s_Y/K. At zero gauge coupling the flat component check uses the same scalar potential, contractions of the Yukawa chiral kinetic numerators with d-dimensional loop momenta, and fixed spin multiplicities. Applying the result to the benchmark assumes this common flat matter and source subtraction. An alternative continuation with independent evanescent flat source couplings would require matching them first. The calculation does not identify the curved external-metric baseline with a dimensional reduction of the gravitational sector.

The loop simple pole equals minus the vacuum beta divided by twice the loop order. This follows from the engineering weights of the pole polynomial: the vacuum coefficient carries μ^(−2ε), while the order-h^(ℓ−1) polynomial carries a marginal-coupling weight 2(ℓ−1)ε. The independent lower-order normalization checks are β_Λ linear in s_Y = 2z s_Y/K at one loop and −hz s_Y/K at two loops, from (4.37)–(4.38). Their curvature poles are +z a_Y/(6Kε) and −hz a_Y/(24Kε), both agreeing with the original FFZ checkpoint. At three loops the corresponding result is

F_R,pole|Y = +h² z a_Y/(24Kε).

Thus c_Y=1/24 in the original normalization, subject to the flat prescription above. The reverse-index enumeration in the runnable check verifies the two surviving contractions, tests that the cross-mass structure is essential, and reproduces the two-loop Y coefficient. There is no independent three-loop integral evaluation here.

## Higher-pole check for a direct curved calculation

Write the Y-linear vacuum counterterm as zs_Y/K times C(h,ε), with μ^(−2ε) outside. Here s_Y does not run, and β_z/z=h−h²+… follows from the X wavefunction. Solving the pole equation
`[−2ε+(−2εh+3h²+…)∂ₕ+h−h²+…]C=−(2−h+3h²/2+…)`
gives

C = 1/ε + h[1/(4ε²)−1/(4ε)] + h²[1/(6ε³)−1/(3ε²)+1/(4ε)].

The h² higher poles are independent of the new three-loop vacuum beta coefficient. The implied Y-linear loop curvature poles are h²za_Y/K times [1/(36ε³)−1/(18ε²)+1/(24ε)]. The two-loop double pole checks the original FFZ subtraction, and the higher poles provide separate checks for any future direct calculation. In the projected w=0 group, both lower-order Newton pole coefficients vanish along the matter flow. The higher-pole recursion therefore requires the baseline triple and double residues to equal the corresponding Y residues before projection.

## Implication for both Newton constraints

The unweighted condition Σw=0 and mass-weighted condition Σzw=0, together with the X null moments, imply Σza_Y=−7M². The calculated Y part of the aggregate order-h² pole is therefore −7h²M²/(24Kε), independent of t. The question is now whether the improvement-independent coefficient c₀ equals +1/24. The condition tests the fixed-scheme pure-metric residue, even though the matter-source flow preserves both moments to the examined order. A finite pure-metric contact redefinition can move the associated beta coefficient along with its matching condition. A shift in the source weights to cancel a residual pole would change the original curvature assignments.

On the projected point a_X=0,a_Y=−1 there is no need to compute a new soft-source insertion. Compute only the baseline ξ_X=ξ_Y=1/6 curved vacuum pole and subtract 1/24. Equivalently a future calculation at the projected point must reproduce the negative Y term identified here. The baseline calculation requires the genuine metric insertion, including the evanescent curvature and measure terms and all two-loop counterterm insertions. Flat vacuum RG data alone do not determine it.

## Runnable contraction check

Python 3 with the standard library. Each product is encoded by its contracted index labels. The real m,y specialization suffices since the terms under test have only |m|² and |y|⁴.

```python
from fractions import Fraction as F
from itertools import product

def contract(tokens, soft_index, offdiag=True):
    labels = sorted(set("".join(indices for _, indices in tokens)))
    total = 0
    for values in product((0, 1), repeat=len(labels)):
        d = dict(zip(labels, values))
        value = 1
        for kind, indices in tokens:
            loc = tuple(d[i] for i in indices)
            if kind == "y":
                value *= loc == (0, 0, 0)
            elif kind == "m":
                value *= (loc[0] != loc[1]) if offdiag else (loc[0] == loc[1])
            else:
                value *= loc == (soft_index, soft_index)
        total += value
    return total

def terms(soft_index, offdiag=True):
    rows = [
        (2, "ilm kln npq mpq", "ij kr", "jr"),
        (-F(1,2), "ilm knp lmq npq", "ij kr", "jr"),
        (-F(1,2), "ijk ijl mnp mnq", "lq kr", "rp"),
        (-F(1,2), "ijk ijl mnp mnq", "kp lr", "qr"),
        (12, "ijk lpn lmr knq", "jp mq", "ri"),
        (6, "ijk lmn imp knq", "jl qr", "pr"),
        (6, "ijk lmn imp knq", "qp jr", "rl"),
    ]
    return [
        coeff*contract(
            [("y",s) for s in yuk.split()]+[("m",s) for s in mass.split()]+[("s",soft)],
            soft_index, offdiag
        )
        for coeff,yuk,mass,soft in rows
    ]

assert terms(1) == [2,-F(1,2),0,0,0,0,0]
assert sum(terms(1,False)) == 0
two_loop = -contract(
    [("y","ijk"),("y","ijl"),("s","rs"),("m","kr"),("m","sl")], 1
)
assert two_loop == -1
assert -sum(terms(1)) * F(-1,6) / 6 == F(1,24)
assert sum(z*(w-2)/2 for z,w in zip((1,2,4),(2,-3,1))) == -7
assert F(1,4)*4/6 == F(1,6)
assert (-F(1,4)*4-1)/6 == -F(1,3)
assert sum(terms(1))/6 == F(1,4)
print("Y improvement coefficient and lower-loop normalization passed")
```

## Sources

- [Martin, Three-loop effective potential for softly broken supersymmetry](https://doi.org/10.1103/PhysRevD.109.015019), Eqs. (4.37)–(4.40).
- [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), ungauged scalar and fermion vacuum diagram basis.
