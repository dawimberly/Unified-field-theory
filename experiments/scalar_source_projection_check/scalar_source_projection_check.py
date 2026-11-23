"""Source saturation of the previously assumed auxiliary-metric scalar kernel.
The projector-to-conformal map assumes stationary diffeomorphism Ward identities.
This checks that map algebraically, not the full loop tensor Ward identity.
"""
from pathlib import Path
import json, shutil, zipfile
import numpy as np
from scipy.special import roots_legendre, gamma

out=Path(__file__).resolve().parents[1]/'experiments'/'scalar_source_projection_check'
out.mkdir(parents=True,exist_ok=True)
prior=Path(__file__).resolve().parents[1]/'experiments'/'scalar_full_response_scan'/'results.json'
data=json.loads(prior.read_text())
checks=[]
def check(name,value,tol=1e-10):
    assert value<tol,(name,value)
    checks.append(dict(name=name,error=float(value),tolerance=tol))
# Integer-dimensional tests with complex analytically continued on-shell momenta.
for D in [4,5,6]:
    for q in [.3,1.8,-2.7]:
        for m2 in [0.,1.,2.]:
            Q=np.zeros(D,dtype=complex); Q[-1]=np.sqrt(complex(q))
            p1=Q/2; p2=Q/2
            p1[0]=np.sqrt(complex(-m2-q/4)); p2[0]=-p1[0]
            E=np.outer(p1,p2)+np.outer(p2,p1)-np.eye(D)*(p1@p2-m2)
            check('on shell',abs(p1@p1+m2)+abs(p2@p2+m2))
            check('stress vertex conservation',float(np.max(abs(Q@E))))
            trace=2*m2-(D-2)*q/2
            check('stress trace = 2 U_external',abs(np.trace(E)-trace))
            theta=np.eye(D)-np.outer(Q,Q)/q
            P0=np.einsum('ij,kl->ijkl',theta,theta)/(D-1)
            scalar_contract=np.einsum('ij,ijkl,kl',E,P0,E)
            check('scalar projector saturation',abs(scalar_contract-trace**2/(D-1)))
            # K0=2 A/(D-1), S2=(1/2) h K h, h=f delta => S2=A f^2.
            K0=2/(D-1)
            check('conformal-to-transverse normalization',abs(.5*np.einsum('ij,ijkl,kl',np.eye(D),K0*P0,np.eye(D))-1))

x,w=roots_legendre(640); x=(x+1)/2; w=w/2
drrows=[]
for D in [4.01,4.2,5.,5.9]:
    a=(D-2)/4; qstar=1/a
    G=(4*np.pi)**(-D/2)*gamma(1-D/2); g=(D-2)*G/2
    def P(q): return (4*np.pi)**(-D/2)*gamma(2-D/2)*np.dot(w,(1+q*x*(1-x))**(D/2-2))
    # W_scalar=T_trace^2/(16 A); same mass trace=2U, A=-U(PU+g)/4.
    # Therefore W=-U/(PU+g): analytic at U=0 and vanishes there.
    ratios=[]
    for delta in [1e-3,1e-5,1e-7]:
        vals=[]
        for sign in [-1,1]:
            q=qstar+sign*delta; U=1-a*q
            A=-U*(P(q)*U+g)/4
            direct=(2*U)**2/(16*A)
            reduced=-U/(P(q)*U+g)
            check('identical-constituent pole cancellation',abs(direct-reduced))
            vals.append(direct)
        ratios.append(dict(delta=delta,W_left=vals[0],W_right=vals[1]))
    slope=(D-2)**2*G/32
    # Unequal-mass external probe, negligible contribution to leading-N determinant.
    residues={str(m2):(2*m2-(D-2)*qstar/2)**2/(16*slope) for m2 in [0.,1.,2.]}
    drrows.append(dict(D=D,qstar=qstar,near_root=ratios,conditional_W_residue_by_probe_mass_squared=residues))

pvrows=[]
for row in data['pv']:
    for root in row['full_roots_in_sampled_domain']:
        q=root['q2']; derivative=root['slope']; trace=2-q
        residue=trace**2/(16*derivative)
        assert residue>0
        pvrows.append(dict(L=row['L'],q2=q,external_trace=trace,conditional_euclidean_W_residue=residue))
results=dict(scope='Conserved minimally coupled scalar source saturation, conditional on full tensor Ward completion and previous conformal kernel. W=log Z connected source response, not S-matrix normalization.',checks=checks,dimreg=drrows,pv=pvrows)
(out/'results.json').write_text(json.dumps(results,indent=2))
lines=['# Scalar source saturation: cancellation in the identical-mass channel','',
'## Main result','',
'The exact finite-D conformal-kernel zero at q=4M^2/(D-2) cancels from the scalar exchange between identical minimally coupled constituent scalars, under the Ward-completed tensor interpretation stated below. Their stress-trace vertex vanishes at the same point. Thus calling this zero a demonstrated observable tachyon in constituent scattering was premature.',
'','The cancellation is channel-specific: an unequal-mass or massless probe has a nonzero trace there. It does not prove that the metric-sector zero is absent or that the vacuum is stable. A new probe here means O(1) species, so it does not modify the leading-N determinant; adding O(N) new matter requires recomputing the kernel.',
'','## Derivation','',
'Use Euclidean q=Q^2 and analytic on-shell p1^2=p2^2=-m_ext^2, Q=p1+p2. E_ab=p1_a p2_b+p2_a p1_b-delta_ab(p1.p2-m_ext^2). Then Q_a E_ab=0 and tr E=2m_ext^2-(D-2)q/2=2U_ext.',
'','For theta_ab=delta_ab-Q_a Q_b/q and P0=theta_ab theta_cd/(D-1), assume the stationary metric Hessian is transverse and write S2=(1/2)h K h. If h_ab=f delta_ab gives S2=A f^2, then K0=2A/(D-1). Coupling h_ab T_ab/2 and integrating the Gaussian gives the connected scalar source functional W0=(tr T)^2/(16A). A cross-exchange coefficient for two separate sources is twice this quadratic coefficient; the sign is unchanged.',
'','In the assumed DR kernel A=-U(PU+g)/4, with U=M^2-(D-2)q/4. For identical constituent sources W0=U^2/(4A)=-U/(PU+g). Since g is nonzero, W0 is regular and zero at U=0. The inverse metric response still has its zero.',
'','## PV source residues','',
'For each previously found subthreshold PV root the identical physical constituent trace is nonzero. The conditional Euclidean W residue is positive, matching the sign of a real source coupled to a canonical positive-kinetic Euclidean scalar. This is stronger than inspecting A\u2032 alone because it includes the source overlap and Gaussian normalization. It remains conditional on the tensor Ward identity and contour/sign conventions of the specified metric model.',
'','| PV L | q at root | Stress trace | Conditional W residue |','|---|---:|---:|---:|']
for r in pvrows: lines.append(f"| {r['L']:g} | {r['q2']:.8f} | {r['external_trace']:.8f} | {r['conditional_euclidean_W_residue']:.8f} |")
lines+=['','## Checks and limits','',f'{len(checks)} numerical algebra checks pass, including complex on-shell conservation, source traces, projector saturation, conformal normalization, and cancellation approaching both sides of the DR root. Integer dimensions test tensor algebra; the trace/projector formulas extend algebraically to finite D.',
'','These tests do not calculate the full loop tensor or prove its Ward identities. They also do not establish the constrained CEV measure, Lorentzian integration contour, quartic mixing, 1/N corrections, above-threshold spectral positivity, or stability against arbitrary perturbations. Scalar-channel massless exchange also cannot be diagnosed as a ghost in isolation from the spin-2 Einstein tensor structure.',
'','## Single most decisive next calculation','',
'Compute the full second metric variation, including all seagull/contact terms, and explicitly check Q_a K_abcd=0 at the stationary background using the same regulator. Compare its transverse scalar coefficient with 2A/(D-1). If that map fails, the source verdict above must be revised. Then calibrate Lorentzian exchange against a canonical scalar and study unequal-mass probes.',
'','Reference for the physical scattering setup and universality: Carone, Claringbold and Vaman, Composite graviton self-interactions in a model of emergent gravity, arXiv:1710.09367, sections II and V: https://arxiv.org/html/1710.09367v2 . The source-cancellation result here is our conditional derivation, not attributed to that paper.',
'','Reproduce with python scalar_source_projection_check.py (requires the preceding scalar_full_response_scan/results.json at the documented workspace path, NumPy and SciPy). Input scan results and both scripts are included.']
(out/'REPORT.md').write_text('\n'.join(lines),encoding='utf8')
shutil.copy2(prior,out/'input_full_response_results.json')
shutil.copy2(__file__,out/Path(__file__).name)
shutil.copy2(Path(__file__).parent/'scalar_full_response_scan.py',out/'scalar_full_response_scan.py')
with zipfile.ZipFile(out/'scalar_source_projection_check.zip','w',zipfile.ZIP_DEFLATED) as z:
    for p in out.iterdir():
        if p.is_file() and p.suffix!='.zip': z.write(p,p.name)
print(json.dumps(dict(checks_passed=len(checks),pv=pvrows,dr_probe_residues=[dict(D=r['D'],residues=r['conditional_W_residue_by_probe_mass_squared']) for r in drrows]),indent=2))
