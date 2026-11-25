# Reproduction of the dimension-mechanism checks

**Conversation date range** 2026-11-01

**DERIVED** This checkpoint records reproduction of the dimension-mechanism checks. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

```python
import json
from fractions import Fraction as F
from pathlib import Path
out=Path.cwd()
rows=[]
for d in range(2,7):
    D=d+1
    polar=F((D-2)*(D-1),2)-1
    canonical=F(D*(D-3),2)
    assert polar==canonical
    rows.append(dict(spatial_dimension=d,spacetime_dimension=D,gravitational_polarizations=int(polar),gradient_scaling=d-2,quartic_scaling=d-4,potential_scaling=d))
# u=(-a*x/2-Omega*y,Omega*x-a*y/2,a*z)
examples=[]
for a in [F(-2),F(0),F(1,2),F(3)]:
    for w in [F(1,3),F(2)]:
        divergence=-a/2-a/2+a
        vorticity=2*w
        stretching=a*vorticity
        enstrophy_production=a*vorticity*vorticity
        assert divergence==0
        examples.append(dict(strain=str(a),rotation=str(w),vorticity=str(vorticity),stretching=str(stretching),enstrophy_production=str(enstrophy_production)))
# Three-dimensional E(s)=A*s+B/s+C*s^3.
# Positive A,B, nonnegative C: stationary point is unique in size.
size_cases=[]
for A,B,C in [(1,1,0),(2,3,0),(1,1,1),(2,3,4)]:
    if C==0:
        y=B/A
    else:
        y=2*B/(A+(A*A+12*C*B)**0.5)
    s=y**0.5
    residual=A-B/s**2+3*C*s**2
    hessian=2*B/s**3+6*C*s
    assert abs(residual)<1e-12 and hessian>0
    size_cases.append(dict(A=A,B=B,C=C,stationary_size=s,derivative_residual=residual,size_hessian=hessian))
result=dict(status='Algebraic consistency controls, not empirical tests or PDE simulations',dimension_table=rows,local_affine_flow_examples=examples,three_dimensional_size_controls=size_cases,limits=['Affine examples have infinite total energy on R3 and are local kinematic controls only.','Size minima are restricted dilation minima, not full field stability.','No gravity or collapse dynamics was simulated.'])
(out/'dimension_mechanism_checks_2026_11_01.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
```
