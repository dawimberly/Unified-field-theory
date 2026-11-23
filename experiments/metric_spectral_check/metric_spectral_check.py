"""Real-axis cut of the free 4D PV auxiliary-metric response; no UV unitarity claim."""
from pathlib import Path
import json, shutil, zipfile
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import quad

OUT=Path(__file__).resolve().parents[1]/'experiments'/'metric_spectral_check'; OUT.mkdir(parents=True,exist_ok=True)
C=(4*np.pi)**-2
def spec(L):
    z=np.array([1.,L,2*L,5*L+1]); ci=np.r_[1.,np.linalg.solve(np.vstack([z[1:]**k for k in range(3)]),-np.ones(3))]
    return z,ci
def log_integral(s,z):
    if s>4*z:
        b=np.sqrt(1-4*z/s)
        return np.log(z)-2+b*np.log((1+b)/(1-b))-1j*np.pi*b
    if s==4*z: return complex(np.log(z)-2)
    a=np.sqrt(4*z/s-1)
    return complex(np.log(z)-2+2*a*np.arctan(1/a))
def response(s,L):
    z,ci=spec(L); U=z+s/2
    logs=np.array([log_integral(s,m) for m in z])
    A=C/4*np.dot(ci,U*U*logs-U*z*np.log(z))
    beta=np.sqrt(np.maximum(0,1-4*z/s))
    imA=-C*np.pi/4*np.dot(ci,U**2*beta)
    imK2=-C*np.pi*s*s/240*np.dot(ci,beta**5)
    assert abs(A.imag-imA)<1e-9*max(1,abs(imA))
    # rho=(1/pi) Im F(-s-i0) calibrated by F(q)=1/(q+m^2).
    rho_scalar=(1/A).imag/np.pi
    W_rho=(2+s)**2*rho_scalar/16
    return dict(s=s,Re_A=float(A.real),Im_A=float(A.imag),rho_inverse_A=float(rho_scalar),rho_identical_source_W=float(W_rho),Im_spin2_kernel=float(imK2),spin2_inverse_spectral_sign=float(np.sign(-imK2)))

# Explicit full tensor cut: imaginary logarithm -pi only on the open pair interval.
basis=[]
for i in range(4):
    h=np.zeros((4,4));h[i,i]=1;basis.append(h)
for i in range(4):
    for j in range(i+1,4):
        h=np.zeros((4,4));h[i,j]=h[j,i]=1/np.sqrt(2);basis.append(h)
basis=np.array(basis); tr=np.trace(basis,axis1=1,axis2=2); H=basis-tr[:,None,None]*np.eye(4)/2
trH=np.trace(H,axis1=1,axis2=2)
a4=(np.outer(trH,trH)+2*np.einsum('aij,bji->ab',H,H))/24
sv=np.r_[np.ones(3)/np.sqrt(3),np.zeros(7)]
nodes,weights=roots_legendre(64)
def tensor_cut(s,L):
    z,ci=spec(L); im=np.zeros((10,10));q=-s
    for m,c in zip(z,ci):
        if s<=4*m:continue
        b=np.sqrt(1-4*m/s);x=.5+nodes*b/2;w=weights*b/2;t=x*(1-x);delta=m+t*q
        ch=-basis[:,3,3,None]*t*q+tr[:,None]*(t*q-m)/2
        lin=q*np.einsum('ai,bi->ab',H[:,:,3],H[:,:,3])[:,:,None]*(1-2*x)**2/4
        a2=lin+(trH[:,None,None]*ch[None,:,:]+trH[None,:,None]*ch[:,None,:])/4
        poly=-3*a4[:,:,None]*delta**2+2*a2*delta-ch[:,None,:]*ch[None,:,:]
        im+=C*c*np.pi/2*np.einsum('abx,x->ab',poly,w)
    return im

rows=[];checks=[];summaries=[]
for L in [3.,5.,10.,30.,100.]:
    # Stay strictly above physical threshold and below first regulator pair threshold.
    clean=[response(float(s),L) for s in np.geomspace(4.0001,4*L*.999,150)]
    assert min(r['rho_identical_source_W'] for r in clean)>0
    assert all(r['spin2_inverse_spectral_sign']==1 for r in clean)
    high=[response(float(s),L) for s in np.geomspace(4*L*1.001,4*(5*L+1)*10,300)]
    rows.extend(dict(L=L,region='physical_only',**r) for r in clean)
    rows.extend(dict(L=L,region='regulator_open',**r) for r in high)
    for s in [4.2,8.,4*L*1.2]:
        im=tensor_cut(s,L);r=response(s,L)
        ward=float(np.max(abs(im[[3,6,8,9],:])))
        scalar_error=abs(sv@im@sv-2*r['Im_A']/3)
        spin2_error=abs(im[4,4]-r['Im_spin2_kernel'])
        assert max(ward,scalar_error,spin2_error)<1e-8
        checks.append(dict(L=L,s=s,ward=ward,scalar_error=float(scalar_error),spin2_error=float(spin2_error)))
    summaries.append(dict(L=L,physical_interval=[4.0001,4*L*.999],physical_sample_count=len(clean),physical_min_rho_W=min(r['rho_identical_source_W'] for r in clean),above_regulator_negative_scalar_samples=sum(r['rho_identical_source_W']<0 for r in high),above_regulator_negative_spin2_cut_samples=sum(r['spin2_inverse_spectral_sign']<0 for r in high),above_regulator_count=len(high)))
# Independent real-log numerical integration, split at roots, checks the analytic formula.
log_checks=[]
for s,z in [(4.2,1.),(8.,1.),(19.,5.),(21.,5.),(100.,1.)]:
    points=[]
    if s>4*z:
        b=np.sqrt(1-4*z/s); points=[(1-b)/2,(1+b)/2]
    val=quad(lambda x:np.log(abs(z-s*x*(1-x))),0,1,points=points,epsabs=1e-10,limit=300)[0]
    error=abs(val-log_integral(s,z).real)
    assert error<1e-8
    log_checks.append(dict(s=s,z=z,error=error))
result=dict(scope='Real-axis spectral-cut test, free minimally coupled 4D PV determinant, stationary contacts; Euclidean boundary q=-s-i0. Not global spectral positivity, constrained CEV, or interacting unitarity.',summaries=summaries,tensor_cut_checks=checks,log_checks=log_checks,rows=rows)
(OUT/'results.json').write_text(json.dumps(result,indent=2))
lines=['# Continuum spectral check','',
'## Main result','',
'All 750 sampled scalar source spectral weights between the physical pair threshold and first PV pair threshold are positive. The spin-2 inverse-response cut has the positive spectral sign at the same 750 points. Explicit tensor cut integrals reproduce the scalar and spin-2 discontinuities and satisfy the longitudinal Ward check. This closes a real-axis continuum consistency check in the specified free auxiliary-metric model, not a proof of full unitarity.',
'','## Conventions and derivation','',
'Constituent M=1, s is timelike invariant mass squared, q=Q_E^2=-s-i0. A canonical positive Euclidean propagator F(q)=1/(q+m^2) has rho(s)=Im F(-s-i0)/pi positive at its pole. This fixes the spectral sign used here.',
'','The exact Feynman-parameter logarithm has Im integral log(z-s*x*(1-x)-i0) dx = -pi*beta, beta=sqrt(1-4z/s), if s>4z, and zero otherwise. Hence Im A=-C*pi/4 sum c_i(z_i+s/2)^2 beta_i. For the transverse traceless tensor coefficient, Im K2=-C*pi*s^2/240 sum c_i beta_i^5. C=1/(4pi)^2.',
'','Below the first PV threshold only the positive-weight physical scalar contributes. Thus Im A<0 and Im(1/A)>0 analytically throughout this interval; similarly the spin-2 inverse has nonnegative continuum weight wherever finite. This sign conclusion is stronger than the finite scan. The identical-constituent source contributes (2+s)^2/16 times Im(1/A)/pi, which is positive.',
'','The full 10 by 10 tensor cut is independently obtained by integrating the shifted quadratic-vertex polynomial over each open pair interval. Contacts are real and add no discontinuity. Three energies per PV spectrum check the longitudinal rows and the scalar and spin-2 projections. The analytic real logarithm is also checked against independent adaptive integration at five points.',
'','## Results by spectrum','',
'| L | Physical-only s interval | Positive scalar samples | Negative scalar samples above PV threshold | Negative spin-2 cut samples above PV threshold |','|---|---|---:|---:|---:|']
for r in summaries:lines.append(f"| {r['L']:g} | {r['physical_interval'][0]:g} to {r['physical_interval'][1]:g} | 150/150 | {r['above_regulator_negative_scalar_samples']}/300 | {r['above_regulator_negative_spin2_cut_samples']}/300 |")
lines+=['','## What remains open','',
'Positive real-axis continuum weight and isolated-pole residues do not exclude complex poles, poles on other sheets, a bad contour, or failures of a dispersion relation. The scalar massless term cannot be interpreted separately from the Einstein spin-2 tensor structure. Signed PV determinants do not themselves specify a physical unitary UV completion. Above-PV negative weights, if found, concern the regulator-open model; they must not be called low-energy experimental particles.',
'','The next decisive test is a first-sheet complex-pole search and a subtracted dispersion-relation reconstruction below the regulator scale. Agreement there would be stronger than more real-axis parameter scans. The constrained-CEV measure and interactions remain separate gaps.',
'','Reproduce: python metric_spectral_check.py (NumPy and SciPy). No experimental claim or novelty claim. Full spectra, sampled response, integration checks and cut-tensor checks are in results.json.']
(OUT/'REPORT.md').write_text('\n'.join(lines),encoding='utf8')
shutil.copy2(__file__,OUT/Path(__file__).name)
with zipfile.ZipFile(OUT/'metric_spectral_check.zip','w',zipfile.ZIP_DEFLATED) as z:
    for p in OUT.iterdir():
        if p.is_file() and p.suffix!='.zip':z.write(p,p.name)
print(json.dumps(dict(summaries=summaries,max_tensor_ward=max(r['ward'] for r in checks),max_tensor_scalar_error=max(r['scalar_error'] for r in checks),max_tensor_spin2_error=max(r['spin2_error'] for r in checks)),indent=2))
