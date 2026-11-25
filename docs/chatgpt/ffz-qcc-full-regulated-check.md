# QCc full regulated check and correction

**Conversation date range** 2026-11-20 to 2026-11-21 (checkpoint sequence)

**DERIVED** This checkpoint records QCc full regulated check and correction. **UNVERIFIED** The results retain the scope and limitations recorded in the source. **UNVERIFIED** Bibliographic details and external claims without a verified paper location remain provisional.

## Research record

**DERIVED** This supersedes the QCc and bare scalar totals in [the finite scalar sum](ffz-qcc-finite-scalar-sum.md).

**DERIVED** The full-integrand test exposed an arithmetic error in the QCc small-pair overlap. At p=q=0 the numerator is −d(4−d)/3 + 3ε(2d/3) = 4dε/3 in d=4−2ε. The earlier checkpoint used 2dε/3. The face script subtracted the correct 4dε/3, which is why the separate face checks passed. The overlap in the assembled projective result was still wrong. The CCcc corner has a different vertex term and is unaffected by this correction.

**DERIVED** The corrected overlap adds 8/(3ε)−4/3 to the previous QCc projective result. The residue becomes 8 and the finite projective coefficient becomes −18.83584915236. After radial integration and the finite-d curvature normalization, the bare doubled QCc result in units h²z/K is

**DERIVED** 1/(72ε²) + [−0.01649742330 − L/24]/ε.

**DERIVED** The updated aggregate of the bare scalar graphs is

**DERIVED** 1/(36ε²) + [−0.05121964552 − L/12]/ε.

**DERIVED** This incorporates the prior numerical QQ simple-pole constant and inherits its numerical status. The baseline Newton test still requires the Yukawa graphs and counterterms.

## Independent full-integrand numerical check

**DERIVED** This check evaluates at nonzero ε before extracting the coefficient. It uses effective pair lengths P,Q instead of physical pair totals as numerical variables and cuts short-face subtractions at P,Q=1, so it does not use the earlier face integrals or finite remainder. The two analytic long-tail models remain the same. Subtract the leading 1/(72ε²), multiply by ε, and extrapolate. At ε=.00125 the two estimates are −.01640531 and −.01642108. Linear extrapolation using the last two ε values yields approximately −.01649. This checks the corrected Laurent coefficient to roughly 10⁻⁵, rather than verifying all displayed digits. SciPy and NumPy are required.

```python
import numpy as np
from scipy.stats import qmc
from scipy.special import beta,gamma
from math import exp
EG=.5772156649015329
def integrand(P,Q,f,g,e):
    d=4-2*e
    p,q=P/f,Q/g
    S,W=1+P+Q,1+p+q
    cx,cy,cxy=2*P*(1+Q)/S,2*Q*(1+P)/S,2*P*Q/S
    A=-d*((6-d)*(cx+cy)/2+(4-d)*W/3)
    B=2*d/3*(p*cx+q*cy+cx+cy+cxy+1+p*p*(1-2*f)+q*q*(1-2*g))
    return S**(-2+e)*(W**(-3*e)*(A+3*e*B/W)-4*d*e/3*(p**(1-3*e)+q**(1-3*e)))
def sample(dim,power,seed):
    return qmc.Sobol(dim,scramble=True,seed=seed).random_base2(power).T
def eval(e,power,seed):
    x,y,a,b=sample(4,power,seed)
    P,Q=x/(1-x),y/(1-y)
    f,g=a*(1-a),b*(1-b)
    d=4-2*e
    c=4*d*e/3
    H=integrand(P,Q,f,g,e)
    hp=integrand(P,0.,f,g,e)
    hq=integrand(0.,Q,f,g,e)
    cutp,cutq=P<1,Q<1
    wp=P**(e-1)/(1-x)**2*f**(-e)
    wq=Q**(e-1)/(1-y)**2*g**(-e)
    mixed=np.mean(wp*wq*(H-cutq*hp-cutp*hq+(cutq&cutp)*c))
    from scipy.special import roots_legendre
    u,w=roots_legendre(900)
    u,w=(u+1)/2,w/2
    p=(u/(1-u))[:,None]
    r=(u*(1-u))[None,:]
    pref=p**(e-1)/(1-u[:,None])**2*r**(-e)
    face=w@(pref*(integrand(p,0.,r,r,e)-(p<1)*c))@w
    f0=beta(1-e,1-e)
    tail=8*d*e/3*beta(1-2*e,1+e)*beta(2*e,2*e)*f0/e
    I=mixed+2*face*f0/e+c*(f0/e)**2+tail
    return e*(I*gamma(3*e)*exp(3*EG*e)/(16*d*(d-1))-1/(72*e**2))
for e in (.01,.005,.0025,.00125):
 print(e,[eval(e,19,seed) for seed in (11,22)],flush=True)
expected=-0.019583843050012695+(-4/9+(8/9)*(7/6))/192
assert abs(expected+0.01649742330)<1e-11
for seed in (11,22):
    assert abs(2*eval(.00125,19,seed)-eval(.0025,19,seed)-expected)<2e-5
```

- **DERIVED** [Heinrich, Sector Decomposition](https://arxiv.org/abs/0803.4177).
