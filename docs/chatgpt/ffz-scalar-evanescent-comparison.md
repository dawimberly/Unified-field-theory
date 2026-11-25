# FFZ explicit scalar curvature-mass comparison

**Conversation date range** 2026-11-22

**DERIVED** This checkpoint records the FFZ explicit scalar curvature-mass comparison. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** The [flat-space recursion target](ffz-scalar-evanescent-double-pole-audit.md) for the scalar curvature-mass insertion is −2/27 in the full bare double pole. The four-Yukawa graphs supply 17/216. This independently evaluates all other graph classes through their flat scalar mass derivatives. Couplings and fermion mass are held fixed, and the scalar masses are set equal only after differentiation.

## Bare master residues

**DERIVED** Strip powers of K and retain the leading poles. The tadpole and its negative mass derivative give A(x)→−x/ε and T₂→1/ε. From the published bare-master subtractions, I has leading double pole −(x+y+z)/2, E has triple pole (sum of six pair products)/3, and G has triple pole −w/3−(u+v+y+z)/6 for central mass w. A doubled central denominator in G gives K with leading triple pole 1/3, by a negative mass derivative or divided difference. Likewise the double-propagator sunset has leading double pole 1/2.

## Scalar derivatives of the graph weights

**DERIVED** The table gives coefficients of the bare **flat mass derivative** at the triple pole in units h²z/K. The explicit curvature insertion multiplies each coefficient by ε/18.

| Graph | Coefficient | Derivation |
| --- | ---: | --- |
| **DERIVED** QQ with two tadpoles | −1 | Flat graph is −A²T₂/2 after extracting the quartic couplings. |
| **DERIVED** QQ basketball | −1/2 | Flat graph is −E/8 with four scalar masses shifted together. |
| **DERIVED** QCc with tadpole | −1/2 | Flat graph is z A times one negative mass derivative of I, at fixed cubic coupling. |
| **DERIVED** QCc doubled | −1/4 | Flat graph is z G/4, with all five lines scalar. |
| **DERIVED** Both CCcc classes | 0 | Their bare K triple residues are independent of the masses. |
| **DERIVED** QFf | 0 | Its A(tadpole) multiplies the leading double-pole bracket (x+w)/2−v−y, where x,w are scalar and v,y are fermion masses. Differentiating the tadpole and the bracket gives equal opposite terms at equal masses. |
| **DERIVED** CcFf | −1/2 | Its flat weight is z/2 times the K_SSSSFF combination. Its triple residue is −(x+w)/3−(u+z_s)/6+2(v+y)/3, where x,w,u,z_s are scalar masses and v,y are fermion masses. |

**UNVERIFIED** The sum is **−11/4** in the flat triple pole, or **−11/72** in the explicit scalar curvature-mass double pole. Adding 17/216 from the direct four-Yukawa kernels yields **−2/27**, exactly the independent flat-space recursion target. The test isolates a specified evanescent insertion. The unresolved discrepancy in the full curvature double pole lies in the other terms or their assembly.

## Reproduction

```python
from fractions import Fraction as F

tadpole_qq = -F(1,2)*2
basketball = -F(1,8)*F(6,3)*2
tadpole_qcc = -F(1,2)
doubled_qcc = F(1,4)*(-F(1,3)-4*F(1,6))
cc = 0
qff = F(1,2)*(1-1)
ccff = F(1,2)*(-2*F(1,3)-2*F(1,6))
flat = sum((tadpole_qq,basketball,tadpole_qcc,doubled_qcc,cc,qff,ccff))
assert flat == -F(11,4)
assert flat/18 == -F(11,72)
assert flat/18+F(17,216) == -F(2,27)
print("Explicit scalar curvature-mass comparison passed")
```

- **DERIVED** [Martin and Robertson, Evaluation of the general 3-loop vacuum Feynman integral](https://doi.org/10.1103/PhysRevD.95.016008), Eqs. (2.35)–(2.45).
- **DERIVED** [Martin, Effective potential at three loops](https://doi.org/10.1103/PhysRevD.96.096005), Eqs. (3.17), (3.27).
