"""Conditional auxiliary-metric determinant audit; not constrained CEV or a physical ghost verdict."""
from pathlib import Path
import json, shutil, zipfile
from functools import lru_cache
import numpy as np
from scipy.special import roots_legendre, gamma
from scipy.optimize import brentq

out=Path(__file__).resolve().parents[1]/'experiments'/'scalar_full_response_scan'
out.mkdir(parents=True,exist_ok=True)
C=(4*np.pi)**-2
@lru_cache(maxsize=None)
def nodes(n):
    x,w=roots_legendre(n); return (x+1)/2,w/2
def dr(q,D,n=320):
    x,w=nodes(n); U=1-(D-2)*q/4
    G=(4*np.pi)**(-D/2)*gamma(1-D/2)
    P=(4*np.pi)**(-D/2)*gamma(2-D/2)*np.dot(w,(1+q*x*(1-x))**(D/2-2))
    return -U*(P*U+(D-2)*G/2)/4
def spectrum(L):
    z=np.array([1.,L,2*L,5*L+1]); mat=np.vstack([z[1:]**k for k in range(3)])
    ci=np.r_[1.,np.linalg.solve(mat,-np.ones(3))]
    assert max(abs(np.dot(ci,z**k))/max(1,np.max(z**k)) for k in range(3))<1e-12
    return z,ci
def pv(q,z,ci,n=320):
    x,w=nodes(n); delta=z[:,None]+q*x*(1-x)
    if np.any(delta<=0): raise ValueError('Outside real subthreshold domain')
    U=z-q/2
    return float(C/4*np.dot(ci,U**2*(np.log(delta)@w)-U*z*np.log(z)))
drrows=[]
for D in [4.01,4.1,4.2,4.5,5.,5.5,5.9]:
    qstar=4/(D-2)
    root=brentq(lambda q:dr(q,D),qstar*.8,qstar*1.2,xtol=1e-12)
    assert abs(root-qstar)<1e-10
    h=1e-5
    derivative=(dr(root+h,D)-dr(root-h,D))/(2*h)
    G=(4*np.pi)**(-D/2)*gamma(1-D/2)
    expected=(D-2)**2*G/32
    assert abs(derivative/expected-1)<1e-7
    drrows.append(dict(D=D,full_spacelike_root=root,factor_root=qstar,slope=derivative,expected_slope=expected))
pvrows=[]
for L in [3.,5.,10.,30.,100.]:
    z,ci=spectrum(L)
    G=C*np.dot(ci,z*np.log(z)); a1=-G/8; a2=C*np.dot(ci,np.log(z))/16
    assert G>0
    approx=-a1/a2
    qs=np.geomspace(.001,4*L*.99,400)
    ss=-np.linspace(.001,3.999,500)
    roots=[]
    for grid in [qs,ss]:
        vals=[pv(q,z,ci) for q in grid]
        for l,r,vl,vr in zip(grid[:-1],grid[1:],vals[:-1],vals[1:]):
            if vl*vr<0:
                root=brentq(lambda q:pv(q,z,ci),l,r,xtol=1e-11)
                h=1e-5
                slope=(pv(root+h,z,ci)-pv(root-h,z,ci))/(2*h)
                err=abs(pv(root,z,ci,640)-pv(root,z,ci,1280))
                roots.append(dict(q2=root,mass_squared=-root,slope=slope,quadrature_difference=err))
    pvrows.append(dict(L=L,z=z.tolist(),weights=ci.tolist(),G=G,einstein_coefficient=G/12,low_q_coefficients=[a1,a2],quadratic_truncation_root=approx,full_roots_in_sampled_domain=roots,spacelike_min=.001,spacelike_max=float(qs[-1]),spacelike_all_negative=all(pv(q,z,ci)<0 for q in qs),quadrature_difference_at_minus3p9=abs(pv(-3.9,z,ci,640)-pv(-3.9,z,ci,1280))))
results=dict(scope='Free minimally coupled auxiliary-metric model, per scalar, M=1, no bare curvature terms, stationary contacts. DR finite D versus four-species PV families. No physical source projection, constrained measure, interactions, or above-cut continuation.',dimreg=drrows,pv=pvrows)
(out/'results.json').write_text(json.dumps(results,indent=2))
lines=['# Full scalar-response and regulator audit','',
'## Outcome','',
'The finite-D minimal auxiliary-metric determinant has an exact spacelike zero at Q^2=4 M^2/(D-2), not merely a zero of a derivative expansion. For five PV spectra, that same spacelike zero is absent on the sampled spacelike intervals. The subthreshold timelike pole is also not regulator-independent. These comparisons change the UV prescription without matching the low-energy curvature coefficients, so they do not establish scheme dependence of a fixed physical theory.',
'','## Equations and assumptions','',
'Set M=1 and q=Q_E^2. In dimensional regularization: U=1-(D-2)q/4; G=(4pi)^(-D/2) Gamma(1-D/2); P(q)=(4pi)^(-D/2) Gamma(2-D/2) integral_0^1 [1+x(1-x)q]^(D/2-2) dx; A(q)=-U(PU+(D-2)G/2)/4.',
'','Thus U=0 is an exact zero and A\u2032(q*)=(D-2)^2 G/32>0 for 4<D<6. This is a conformal inverse-kernel zero, not yet a verified physical scattering pole. A positive slope alone does not establish a positive-norm particle: source overlaps, the full tensor constraints, Lorentzian continuation and sign normalization must be established.',
'','For PV, each species has its own U_i=z_i-q/2. The summed kernel is C/4 sum_i c_i[(z_i-q/2)^2 integral log(z_i+x(1-x)q) dx -(z_i-q/2)z_i log(z_i)], C=1/(4pi)^2. Consequently the shared U factor of the single-mass dimensional-regulator expression is lost. The first three PV moments vanish.',
'','The exact low-q coefficients are a1=-G_PV/8 and a2=C sum c_i log(z_i)/16. The nonzero root of a1*q+a2*q^2 is compared with the full nonlocal kernel. No expansion is used to locate the full roots.',
'','## PV results','',
'| Regulator L | Quadratic root q | Full subthreshold roots q | No sampled spacelike zero |','|---|---:|---|---|']
for r in pvrows:
    lines.append(f"| {r['L']:g} | {r['quadratic_truncation_root']:.6f} | {', '.join(format(x['q2'],'.8f') for x in r['full_roots_in_sampled_domain']) or 'none found'} | {r['spacelike_all_negative']} |")
lines+=['','## What this establishes','',
'Seven finite-D root locations and slopes agree with the analytic factorization. Five PV spectra satisfy moment cancellation and positive induced Einstein coefficient. Each PV spectrum is scanned at 400 spacelike and 500 subthreshold timelike points, with bracketing root refinement. This finite sampling does not prove global root absence. Signed PV weights are not a unitary UV completion.',
'','A pole near the pair-production threshold is sensitive to higher derivative terms; extrapolating an R+R^2 truncation there is not reliable. Conversely, the exact finite-D zero cannot be dismissed solely as a truncation artifact.',
'','## Biggest remaining weakness and next calculation','',
'The conformal kernel has not been converted into a conserved-source, gauge-consistent physical scattering amplitude. Therefore none of these numbers is a final ghost verdict. Next derive the full tensor projector decomposition and the scalar source overlap using the same contact prescription, then analytically continue below threshold and compare the pole residue with a canonical positive-energy scalar exchange. Only after that should interactions or a claimed healthy scalar be added.',
'','For regulator comparisons of the same physical theory, match the Einstein term and independent curvature-squared coefficients; imposing zero bare terms separately in different prescriptions does not perform this matching.',
'','## References and reproduction','',
'Carone, Claringbold and Vaman, Composite graviton self-interactions in a model of emergent gravity, arXiv:1710.09367, Phys. Rev. D 98, 024041 (2018): https://arxiv.org/abs/1710.09367 . Carone, Composite Gravitons from Metric-Independent Quantum Field Theories, arXiv:1911.13166: https://arxiv.org/abs/1911.13166 . These references motivate the program; the conditional kernels and scans above are our calculation, not claims attributed to those papers.',
'','Run python scalar_full_response_scan.py with NumPy and SciPy. Dimensionless masses are in constituent-mass units. No claim of novelty or experimental confirmation.']
(out/'REPORT.md').write_text('\n'.join(lines),encoding='utf8')
shutil.copy2(__file__,out/Path(__file__).name)
with zipfile.ZipFile(out/'scalar_full_response_scan.zip','w',zipfile.ZIP_DEFLATED) as archive:
    for p in out.iterdir():
        if p.is_file() and p.suffix!='.zip': archive.write(p,p.name)
print(json.dumps(results,indent=2))
