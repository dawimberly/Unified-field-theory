"""Bounded scalar first-sheet contours, audited against parameter quadrature.
Run: python experiments/complex_scalar_first_sheet_2026_10_30/check.py
Only NumPy/SciPy required. q=Euclidean momentum squared; cut=(-infinity,-4].
"""
from pathlib import Path
import json, platform
import numpy as np
import scipy
from scipy.special import roots_legendre
from scipy.integrate import quad

C=(4*np.pi)**-2
def spectrum(L):
    z=np.array([1.,L,2*L,5*L+1])
    c=np.r_[1.,np.linalg.solve(np.array([z[1:]**k for k in range(3)]),-np.ones(3))]
    return z,c

def logint(q,z):
    # Principal square root and principal atanh, continued from q>0.
    # Compare with the defining principal-log parameter integral below.
    q=np.asarray(q,dtype=complex)
    r=np.sqrt(q[...,None]/(q[...,None]+4*z))
    return np.log(z)-2+2*np.arctanh(r)/r

def A(q,L):
    z,c=spectrum(L); q=np.asarray(q,dtype=complex)
    u=z-q[...,None]/2
    return C/4*np.sum(c*(u*u*logint(q,z)-u*z*np.log(z)),axis=-1)

def direct(q,L,n):
    z,c=spectrum(L); x,w=roots_legendre(n);x=(x+1)/2;w=w/2
    logs=np.log(z[:,None]+q*x*(1-x))@w
    u=z-q/2
    return C/4*np.dot(c,u*u*logs-u*z*np.log(z))

def rectangle(R,h,n):
    vertices=[-R+1j*h,R+1j*h,R+1j*R,-R+1j*R,-R+1j*h]
    return np.concatenate([np.linspace(a,b,n,endpoint=False) for a,b in zip(vertices[:-1],vertices[1:])]+[np.array(vertices[-1:])])

def contour(L,R,h,n):
    q=rectangle(R,h,n)
    for refinement in range(22):
        f=A(q,L)/q; mid=(q[:-1]+q[1:])/2; fm=A(mid,L)/mid
        p1=np.angle(fm/f[:-1]);p2=np.angle(f[1:]/fm)
        bad=(np.abs(p1)>.12)|(np.abs(p2)>.12)|(np.abs(p1-p2)>.03)
        bad|=(abs(f[1:]-f[:-1])>.3*np.minimum(abs(f[1:]),abs(f[:-1])))
        if not np.any(bad):break
        refined=[]
        for j in range(len(mid)):
            refined.append(q[j])
            if bad[j]:refined.append(mid[j])
        refined.append(q[-1]);q=np.array(refined)
    else:raise RuntimeError('Contour refinement did not settle')
    f=A(q,L)/q
    steps=np.angle(f[1:]/f[:-1]); winding=float(steps.sum()/(2*np.pi))
    return dict(L=L,R=R,imaginary_floor=h,points_per_edge=n,winding=winding,
                rounded_count=int(round(winding)),final_points=len(q),refinements=refinement,max_phase_step=float(max(abs(steps))),min_abs_deflated=float(min(abs(f))))

def imA(s,L):
    z,c=spectrum(L); beta=np.sqrt(np.maximum(0,1-4*z/s))
    return -C*np.pi/4*np.dot(c,(z+s/2)**2*beta)

def dispersion(q,L):
    z,c=spectrum(L); a1=-C*np.dot(c,z*np.log(z))/8
    # Two subtractions; map all regulator thresholds and infinity separately.
    edges=list(4*z)+[np.inf]; total=0.; err=0.
    for lo,hi in zip(edges[:-1],edges[1:]):
        val,e=quad(lambda s:imA(s,L)/(s*s*(s+q)),lo,hi,epsabs=1e-11,epsrel=1e-8,limit=300)
        total+=val;err+=e
    return a1*q+q*q*total/np.pi,err*q*q/np.pi

def main():
    audits=[]; contours=[]; dispersions=[]
    for L in [3.,5.,10.,30.,100.]:
        for q in [1+2j,-2+.3j,-8+2j,-50+5j,100+30j]:
            exact=complex(A(q,L)); d1=direct(q,L,512);d2=direct(q,L,1024)
            scale=max(abs(exact),1e-8)
            audits.append(dict(L=L,q=[q.real,q.imag],relative_exact_quad_error=float(abs(exact-d2)/scale),relative_quad_change=float(abs(d2-d1)/scale),conjugation_error=float(abs(A(q.conjugate(),L)-exact.conjugate()))))
        for R in [4*L,20*L,100*L]:
            for h in [.05,.005]:
                for n in [512,1024,2048]:
                    contours.append(contour(L,R,h,n))
        for q in [.5,2.,10.]:
            reconstructed,e=dispersion(q,L); exact=A(q,L).real
            dispersions.append(dict(L=L,q=q,exact=float(exact),reconstructed=float(reconstructed),relative_error=float(abs(exact-reconstructed)/max(abs(exact),1e-10)),integration_error_estimate=e))
    result=dict(scope='Scalar A only, fixed no-bare-curvature PV prescriptions. Numerical bounded evidence, not global root absence or constrained CEV.',versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),parameter_audits=audits,contours=contours,kernel_dispersion=dispersions)
    assert max(x['relative_exact_quad_error'] for x in audits)<1e-8
    assert max(x['conjugation_error'] for x in audits)<1e-10
    assert all(abs(x['winding']-x['rounded_count'])<1e-8 and x['max_phase_step']<.25 for x in contours)
    assert max(x['relative_error'] for x in dispersions)<1e-8
    out=Path(__file__).resolve().parent/'results.json';out.write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(max_quad_error=max(x['relative_exact_quad_error'] for x in audits),counts=sorted(set(x['rounded_count'] for x in contours)),max_phase_step=max(x['max_phase_step'] for x in contours),max_dispersion_error=max(x['relative_error'] for x in dispersions)),indent=2))

if __name__=='__main__':main()
