# Does \(\xi\) run away from \(1/6\)? (new-idea rung)

Runnable: `ffz_xi_running.py`.  
Stacked on surface/contact (`FFFF_SURFACE_TERM.md`, #20).  
**No tuning to \(1/24\).** Derived vs assumed labeled.

Standing GAP \(=+0.00965052\pm0.0024\) (units \(h^2 z/K\)).

## (1) \(\beta_\xi\) and RG stability

| Argument | Result | Label |
| --- | --- | --- |
| SUSY chiral kinetic \(\int\mathrm{d}^4\theta\,\bar E\,\bar\Phi\Phi\) | Geometric \(\xi\) **not independent**; fixed to \(1/6\) by the multiplet | **derived** (superspace); **assumed** (no SUSY-breaking \(\xi\) deformation) |
| Jack–Osborn local RG: \(\beta_a=\Gamma a+b\), \(a=1-6\xi\) | SUSY \(\tau=0\) prescription (JO 7.5) \(\Rightarrow\delta=0\Rightarrow b=0\) | **derived** (JO + packet); scheme **assumed** (\(\tau=0\)) |
| Component check JO 9.16 / 9.26 | \(\delta^{(1)}=0\); radial \(\delta^{(2)}\): \(-2\lambda T+Th+Th=0\) | **derived** |
| Soft \(\Gamma=h\partial_h(B/h)\) with \(B=6h\gamma\), \(\gamma=h/2-h^2/2+\cdots\) | \(\Gamma=3h-6h^2+\cdots\) | **derived** (packet / JO 5.21) |
| Therefore | \(\beta_\xi=\Gamma(\xi-1/6)\), \(\beta_\xi(1/6)=0\) | **derived** |
| Packet fixed-\(1/6\) additive at \(h^2,h^3\) | both **0** | **derived** (checkpoint runnable) |
| Effect of a new 3-loop \(d_I\) pole (other continuation) | first hits 4d additive beta at \(O(h^4)\), beyond \(O(h^2)\) Newton residue | **derived** (packet scope) |

**Verdict (1):** \(\xi=1/6\) is **RG-stable** in this SUSY WZ setup. It does **not** run away.

## (2) \(\delta\xi\) and \(a_1\propto(1/6-\xi)\)

| Quantity | Value | Label |
| --- | ---: | --- |
| \(\xi_0\) (IC / SUSY) | \(1/6\) | **derived** |
| \(\delta\xi\) | \(0\) | **derived** (\(\beta_\xi(1/6)=0\)) |
| \(a_1\) coeff \((1/6-\xi)\) | \(0\) | **derived** |
| \(\delta c\) in units \(h^2 z/K\) | \(\boldsymbol{0}\) | **derived** |

A non-SUSY independent \(\xi\) deformation would be a **free knob** (tuning), and a pure vacuum \(a_1\) pole lacks the \(z\) weight of the \(h^2 z/K\) test anyway.

## (3) vs GAP; \(\ln2/72\)

| Quantity | Value |
| --- | ---: |
| Predicted \(\delta c\) | \(0\) |
| **Ratio predicted / GAP** | \(\boldsymbol{0}\) |
| \(\ln2/72\) | \(0.009627\) (ratio to GAP \(0.998\)) |
| Does \(\ln2/72\) arise from \(\xi\) running? | **No** |

\(\ln2/72\) is **not** produced by \(\beta_\xi\), JO \(\delta=0\), or \(a_1\propto(1/6-\xi)\). It remains an unexplained numerical proximity.

## (4) Higher poles

Triple/double cancelation **unaffected** (\(\xi\) stays at \(1/6\); soft dictionary and CT unchanged).

## Plain statement

In the SUSY WZ+Yukawa Newton test, \(\xi=1/6\) is protected (superspace + Jack–Osborn \(\beta_\xi=\Gamma(\xi-1/6)\) with additive \(b=0\)). It does not run away, \(a_1\) stays zero, and this rung contributes **ratio 0** vs the \(\sim+0.01\) gap. \(\ln2/72\) is not derived here.
