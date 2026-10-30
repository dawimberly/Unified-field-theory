"""Full 4D PV metric Hessian from scalar bubble plus stationary contacts.
No constrained measure or clock-field action is included.
"""
from pathlib import Path
from functools import lru_cache
import json, shutil, zipfile
import numpy as np
from scipy.special import roots_legendre, gamma
from scipy.optimize import brentq

OUT=Path(__file__).resolve().parents[1]/'experiments'/'metric_tensor_ward_check'
OUT.mkdir(parents=True,exist_ok=True)
C=(4*np.pi)**-2
basis=[]; labels=[]
for i in range(4):
    h=np.zeros((4,4)); h[i,i]=1; basis.append(h); labels.append(f'{i}{i}')
for i in range(4):
    for j in range(i+1,4):
        h=np.zeros((4,4)); h[i,j]=h[j,i]=1/np.sqrt(2); basis.append(h); labels.append(f'{i}{j}')
basis=np.array(basis); tr=np.trace(basis,axis1=1,axis2=2)
H=basis-tr[:,None,None]*np.eye(4)/2
trH=np.trace(H,axis1=1,axis2=2)
angular4=(np.outer(trH,trH)+2*np.einsum('aij,bji->ab',H,H))/24
gauge_indices=[labels.index(s) for s in ['33','03','13','23']]
scalar_vec=np.array([1,1,1,0,0,0,0,0,0,0])/np.sqrt(3)
@lru_cache(maxsize=None)
def nodes(n):
    x,w=roots_legendre(n); return (x+1)/2,w/2
def spec(L):
    z=np.array([1.,L,2*L,5*L+1]); ci=np.r_[1.,np.linalg.solve(np.vstack([z[1:]**k for k in range(3)]),-np.ones(3))]
    return z,ci
def tensor(q,L,n=320,contacts=True):
    z,ci=spec(L); x,w=nodes(n); t=x*(1-x)
    p=-C/4*np.dot(ci,z*z*np.log(z))
    Kcontact=p/2*np.eye(10)-p/4*np.outer(tr,tr)
    loop=np.zeros((10,10))
    for mass2,weight in zip(z,ci):
        delta=mass2+t*q
        ch=-basis[:,3,3,None]*t*q+tr[:,None]*(t*q-mass2)/2
        # l^2 angular terms after shift: (1-2x)^2 (H_a Q).(H_b Q)/4,
        # plus (tr H_a c_b+tr H_b c_a)/4.
        linear=q*np.einsum('ai,bi->ab',H[:,:,3],H[:,:,3])[:,:,None]*(1-2*x)**2/4
        angular2=linear+(trH[:,None,None]*ch[None,:,:]+trH[None,:,None]*ch[:,None,:])/4
        integrand=(-3*angular4[:,:,None]*delta**2+2*angular2*delta-ch[:,None,:]*ch[None,:,:])*np.log(delta)
        loop+=C*weight*np.einsum('abx,x->ab',integrand,w)
    return (Kcontact if contacts else 0)-loop/2
def scalar_A(q,L,n=320):
    z,ci=spec(L); x,w=nodes(n); U=z-q/2
    return C/4*np.dot(ci,U**2*(np.log(z[:,None]+q*x*(1-x))@w)-U*z*np.log(z))

rows=[]; roots=[]
for L in [3.,5.,10.,30.,100.]:
    for q in [0.,.01,.5,2.,10.,-.5,-2.,-3.9]:
        K=tensor(q,L); scale=max(np.max(abs(K)),1e-4)
        ward=np.max(abs(K[gauge_indices,:]))
        scalar=scalar_vec@K@scalar_vec; target=2*scalar_A(q,L)/3
        hI=np.array([1.,1.,1.,1.,0.,0.,0.,0.,0.,0.])
        conformal=hI@K@hI/2
        eig=np.linalg.eigvalsh(K)
        row=dict(L=L,q2=q,ward_absolute=float(ward),ward_scaled=float(ward/scale),scalar_eigenvalue=float(scalar),scalar_expected=float(target),scalar_error=float(abs(scalar-target)),conformal_error=float(abs(conformal-scalar_A(q,L))),eigenvalues=eig.tolist(),bubble_only_longitudinal=float(np.max(abs(tensor(q,L,contacts=False)[gauge_indices,:]))))
        rows.append(row)
    # Direct tensor-based scalar root, independent of log-kernel bracketing.
    qs=-np.linspace(.001,3.999,180)
    vals=[float(scalar_vec@tensor(q,L)@scalar_vec) for q in qs]
    for a,b,fa,fb in zip(qs[:-1],qs[1:],vals[:-1],vals[1:]):
        if fa*fb<0:
            root=brentq(lambda q:scalar_vec@tensor(q,L)@scalar_vec,a,b,xtol=1e-10)
            h=1e-5
            slope=(scalar_vec@tensor(root+h,L)@scalar_vec-scalar_vec@tensor(root-h,L)@scalar_vec)/(2*h)
            residue=(2-root)**2/(24*slope)
            roots.append(dict(L=L,q2=float(root),K0_slope=float(slope),W_residue=float(residue),ward_absolute=float(np.max(abs(tensor(root,L)[gauge_indices,:])))))
maxward=max(r['ward_absolute'] for r in rows)
maxscalar=max(r['scalar_error'] for r in rows)
assert maxward<2e-9, maxward
assert maxscalar<2e-9, maxscalar
assert max(r['conformal_error'] for r in rows)<2e-9
assert all(r['W_residue']>0 for r in roots)
refinements=[]
for L in [3.,100.]:
    for q in [-3.9,2.]:
        err=float(np.max(abs(tensor(q,L,640)-tensor(q,L,1280))))
        assert err<2e-9
        refinements.append(dict(L=L,q2=q,error=err))
drchecks=[]
def dr_bilinear(q,D,th,tj,trhj,hl,jl):
    x,w=nodes(640); t=x*(1-x); delta=1+t*q; pref=(4*np.pi)**(-D/2)
    p=pref*gamma(-D/2)/2
    trH=(1-D/2)*th; trJ=(1-D/2)*tj
    trHJ=trhj+(D/4-1)*th*tj
    c_h=-t*q*hl+th*(t*q-1)/2; c_j=-t*q*jl+tj*(t*q-1)/2
    a4=(trH*trJ+2*trHJ)/(D*(D+2))
    a2=((1-2*x)**2*q*(hl-th/2)*(jl-tj/2)+trH*c_j+trJ*c_h)/D
    I0=pref*gamma(2-D/2)*delta**(D/2-2)
    I2=pref*(D/2)*gamma(1-D/2)*delta**(D/2-1)
    I4=pref*D*(D+2)/4*gamma(-D/2)*delta**(D/2)
    return p*trhj/2-p*th*tj/4-np.dot(w,a4*I4+a2*I2+c_h*c_j*I0)/2
for D in [4.01,4.2,5.,5.9]:
    for q in [.5,1.,2.,-1.]:
        LL=dr_bilinear(q,D,1,1,1,1,1)
        LT=dr_bilinear(q,D,1,D-1,0,1,0)
        TT=dr_bilinear(q,D,D-1,D-1,D-1,0,0)
        x,w=nodes(640); U=1-(D-2)*q/4
        G=(4*np.pi)**(-D/2)*gamma(1-D/2)
        P=(4*np.pi)**(-D/2)*gamma(2-D/2)*np.dot(w,(1+q*x*(1-x))**(D/2-2))
        A=-U*(P*U+(D-2)*G/2)/4
        error=abs(TT-2*A)
        assert max(abs(LL),abs(LT),error)<1e-10,(D,q,LL,LT,error)
        drchecks.append(dict(D=D,q2=q,longitudinal=LL,mixed=LT,transverse_trace_error=error))
results=dict(scope='Free auxiliary-metric determinant, stationary contacts; full 4D PV tensor and finite-D longitudinal/trace bilinears. Not constrained CEV.',basis_labels=labels,rows=rows,roots=roots,max_ward_absolute=maxward,max_scalar_error=maxscalar,quadrature_refinements=refinements,finite_D_checks=drchecks)
(OUT/'results.json').write_text(json.dumps(results,indent=2))
lines=['# Full tensor and contact Ward check','',
'## Result','',
f'The full 10 by 10 metric Hessian passes the longitudinal Ward check in 40 momentum/spectrum cases. Maximum absolute longitudinal component: {maxward:.5g}. Maximum scalar-projector mismatch with 2A/3: {maxscalar:.5g}. The three subthreshold scalar roots are reproduced directly from the tensor coefficient, with positive connected-source residues.',
'','This verifies the tensor-projection identity for this four-dimensional PV free determinant. Sixteen additional finite-D tensor tests at D=4.01,4.2,5,5.9 check the longitudinal, mixed longitudinal/transverse-trace and transverse-trace bilinears using exact dimensionally integrated loop moments. Their conservation components vanish and the scalar coefficient agrees with the prior full kernel within 1e-10. This supports the channel-specific DR source cancellation within this determinant model. It does not establish constrained CEV equivalence or a Lorentzian UV-complete theory.',
'','## Derivation','',
'Take Q along axis 3 and q=Q_E^2, continued to negative q algebraically below the physical pair threshold q=-4. For a symmetric polarization h, the linear scalar kinetic vertex is v_h=k.h.(k+Q) -(tr h)/2 [k.(k+Q)+z]. Shift k=l-xQ and define H=h-(tr h)I/2, t=x(1-x), Delta=z+tq, c_h=-t Q.h.Q+(tr h)(tq-z)/2. Then v_h=l.H.l+(1-2x)l.H.Q+c_h.',
'','The rotationally averaged v_h v_j is integrated with the four-dimensional PV-summed logarithmic integrals I0=-C log Delta, I2=2C Delta log Delta, I4=-3C Delta^2 log Delta. Omitted polynomial pieces vanish after sum c_i z_i^k=0, k=0,1,2. The l^4 coefficient is [tr H tr J+2 tr(HJ)]/24; the l^2 coefficient is [(1-2x)^2 (H Q).(J Q)+tr H c_j+tr J c_h]/4.',
'','With p=-C sum c_i z_i^2 log z_i/4, the stationary contact Hessian bilinear is K_contact[h,j]=p tr(hj)/2-p tr h tr j/4. It includes the scalar quadratic metric vertex and the cosmological term tuned to cancel the tadpole. The bubble Hessian is -integral(v_h v_j)/2. Contact and bubble are evaluated separately; no transverse projection is imposed on the result.',
'','All ten real orthonormal symmetric polarizations are included. For Q along axis 3, the four pure gauge directions are h33,h03,h13,h23. Their complete Hessian rows vanish numerically. The bubble alone has nonzero longitudinal rows: omitting contacts fails this test.',
'','The transverse scalar polarization is diag(1,1,1,0)/sqrt(3). Its eigenvalue is K0=2A/3, where A is the separately evaluated conformal log kernel. An on-shell identical-constituent source has trace 2-q; its connected quadratic source residue is (2-q*)^2/[24 K0\u2032(q*)]. This has the sign of a canonically normalized positive-kinetic Euclidean scalar below threshold. Positivity of this isolated residue is not proof of positivity of the whole spectrum.',
'','## Tensor-derived roots','',
'| PV L | q at root | K0 derivative | Connected-source residue |','|---|---:|---:|---:|']
for r in roots: lines.append(f"| {r['L']:g} | {r['q2']:.9f} | {r['K0_slope']:.9g} | {r['W_residue']:.9g} |")
lines+=['','## Remaining questions','',
'1. Does the constrained metric/clock-field measure produce the same Hessian with this regulator? This calculation does not include that measure or independent clock fluctuations.',
'2. Does the Lorentzian contour and continuation give the canonical exchange sign when derived from the full microscopic path integral? The Euclidean source sign is now calibrated algebraically, but no contour prescription was established.',
'3. What happens on and above the constituent threshold, including complex poles, spectral density and regulator thresholds? No global ghost-free verdict is warranted.',
'4. After matching the Einstein and independent curvature-squared coefficients, which regulator differences remain? Zero bare terms across spectra is not a matched physical comparison.',
'','The next calculation should test the spectral discontinuity above the physical threshold while below the first regulator threshold, with the same full tensor contacts. It can expose negative spectral weight that a single subthreshold pole misses. Constrained-measure equivalence remains a separate requirement.',
'','Reference: Carone, Claringbold and Vaman, Composite graviton self-interactions in a model of emergent gravity, arXiv:1710.09367, sections II and V, https://arxiv.org/html/1710.09367v2 . Tensor derivation and numerics above are our calculation, not a claim attributed to that paper.',
'','Reproduce: python metric_tensor_ward_check.py (NumPy and SciPy). Dimensionless constituent M=1; results per scalar, overall N suppressed. All input spectra, individual checks, eigenvalues and quadrature refinements are in results.json.']
(OUT/'REPORT.md').write_text('\n'.join(lines),encoding='utf8')
shutil.copy2(__file__,OUT/Path(__file__).name)
with zipfile.ZipFile(OUT/'metric_tensor_ward_check.zip','w',zipfile.ZIP_DEFLATED) as z:
    for p in OUT.iterdir():
        if p.is_file() and p.suffix!='.zip': z.write(p,p.name)
print(json.dumps(dict(max_ward_absolute=maxward,max_scalar_error=maxscalar,roots=roots,refinements=refinements),indent=2))
