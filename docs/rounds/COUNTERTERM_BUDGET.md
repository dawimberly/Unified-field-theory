# FFZ baseline counterterm pole budget

In the original fixed-1/6 subtraction, assume the flat supersymmetric matter counterterms and the absence of an additive X curvature counterterm through the required order, as in the packet. This accounting uses units h²z/K and L=ln(z/μ²).

The two-loop bare baseline is proportional to h z^(1−2ε). With δh/h=3h/(2ε) and δz/z=h/(2ε), the combined insertion multiplies that graph by h(2−ε)/ε. Its pole coefficients before the insertion, stripping hz/K, are −1/(24ε²)+(L/12−1/24)/ε. Thus the new counterterm gives −1/(12ε³)+(L/6−1/24)/ε². The existing one-loop mass insertion gives +1/(12ε³)−(L/12+1/16)/ε². Their sum has zero triple pole and double pole (L/12−5/48)/ε².

The required bare FFff(A)+FFff(B) budget is therefore

- Triple pole 1/36.
- Double pole 1/72−L/12.

This uses the independent Y higher-pole target and the bare subtotal 5/144 at L=0. It does not allocate an unfinished corner result to a full graph. With the regulated FFff(A) inference triple 1/576, the remaining bare triple budget for FFff(B) is 5/192.

## Simple-pole counterterm

The equal-mass renormalized sunset has finite constant C=−15/2+2√3 Cl₂(π/3). Reinsert the three one-loop tadpole subdivergences into the finite sunset to get the bare finite term in units z/K²:
−3L²+9L+C−3−3ζ(2)/2. Multiplication by (1−2ε)/36 gives the two-loop curvature graph finite term. The two-loop counterterm contributes twice this finite term minus the simple residue of the original graph. Adding the one-loop insertion from the packet gives the combined simple-pole counterterm

(C+6)/18 + (1−ζ(2))/24 + 7L/48 − L²/8.

As a check, the three-loop bare contribution carries e^(−3εL). The FFff higher-pole budget implies −L/24+L²/8 in its simple pole. Adding the previous bare subtotal's −5L/48 and the counterterm logarithms cancels both powers of L exactly. The full projected simple pole remains open. The FFff graph integrals must supply the remaining bare residues. The higher-pole targets are necessary checks, and the simple-pole cancellation cannot be inferred from them.

## Runnable exact check

```python
from fractions import Fraction as F
two_loop_triple = 2*F(-1,24)
two_loop_double = 2*F(-1,24) - F(-1,24)
assert two_loop_triple + F(1,12) == 0
assert two_loop_double - F(1,16) == -F(5,48)
assert -F(1,18) + F(5,48) - F(5,144) == F(1,72)
assert F(1,36) - F(1,576) == F(5,192)
assert F(1,6)-F(1,12)+F(1,16) == F(7,48)
assert -F(1,6)+F(1,24) == -F(1,8)
```

- [Martin, Two-loop effective potential](https://arxiv.org/abs/hep-ph/0111209), sunset and tadpole subtractions.
