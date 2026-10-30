"""Once-subtracted inverse scalar dispersion with explicit isolated poles.
Reproduce: python experiments/inverse_scalar_dispersion_2026_10_30/check.py
Uses NumPy/SciPy; previous files are unchanged. No physical health claim.
"""
from pathlib import Path
import json, platform
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import roots_legendre

LD=np.longdouble; CD=np.clongdouble
C=LD(1)/(4*np.arccos(LD(-1)))**2
def spectrum(L):
    z=np.array([1,L,2*L,5*L+1],dtype=LD)
    w=np.array([1/np.prod(z[i]-np.delete(z,i)) for i in range(4)])
    return z,w/w[0]

def A(q,L):
    z,c=spectrum(L);q=np.asarray(q,dtype=CD)
    r=np.sqrt(q[...,None]/(q[...,None]+4*z))
    logs=np.log(z)-2+2*np.arctanh(r)/r;u=z-q[...,None]/2
    return C/4*np.sum(c*(u*u*logs-u*z*np.log(z)),axis=-1)

def boundary(s,L):
    s=LD(s);z,c=spectrum(L);logs=np.empty(4,dtype=CD)
    for i,m in enumerate(z):
        if s>4*m:
            beta=np.sqrt(1-4*m/s)
            logs[i]=np.log(m)-2+beta*(np.log1p(beta)-np.log1p(-beta))-1j*np.arccos(LD(-1))*beta
        elif s==4*m:logs[i]=np.log(m)-2
        else:
            a=np.sqrt(4*m/s-1);logs[i]=np.log(m)-2+2*a*np.arctan(1/a)
    u=z+s/2
    return C/4*np.sum(c*(u*u*logs-u*z*np.log(z)))

def rho(s,L):return float((1/boundary(s,L)).imag/np.arccos(LD(-1)))

def poles(L):
    z,c=spectrum(L);a1=-C*np.dot(c,z*np.log(z))/8
    rows=[dict(q=0.,residue=float(1/a1),kind='massless')]
    grid=np.linspace(-3.999999,-.001,700)
    for a,b in zip(grid[:-1],grid[1:]):
        if A(a,L).real*A(b,L).real<0:
            p=brentq(lambda q:float(A(q,L).real),a,b,xtol=1e-13)
            slopes=[float(A(p+1j*h,L).imag/h) for h in [1e-4,1e-5,1e-6]]
            rows.append(dict(q=p,residue=1/slopes[-1],kind='isolated',derivative_relative_change=abs(slopes[-1]/slopes[-2]-1)))
    return rows

def count_contour(L,left,right,height,n):
    v=[left-1j*height,right-1j*height,right+1j*height,left+1j*height,left-1j*height]
    q=np.concatenate([np.linspace(a,b,n,endpoint=False) for a,b in zip(v[:-1],v[1:])]+[np.array(v[-1:])])
    f=A(q,L)/q;steps=np.angle(f[1:]/f[:-1]);w=float(np.sum(steps)/(2*np.arccos(LD(-1))))
    return dict(L=L,left=left,right=right,height=height,points_per_edge=n,winding=w,count=round(w),max_phase_step=float(max(abs(steps))))

def reconstruct(q,L,R,ps):
    q0=1.;kernel=lambda s:1/(s+q)-1/(s+q0)
    pole=sum(p['residue']*(1/(q-p['q'])-1/(q0-p['q'])) for p in ps)
    z,_=spectrum(L);edges=[float(x) for x in 4*z]+[R];integral=0.;error=0.
    for a,b in zip(edges[:-1],edges[1:]):
        value,e=quad(lambda s:rho(s,L)*kernel(s),a,b,epsabs=2e-10,epsrel=2e-9,limit=300)
        integral+=value;error+=e
    exact=float((1/A(q,L)-1/A(q0,L)).real)
    # Cauchy outer boundary of the slit disk. Constant F(infinity) cancels
    # in this difference; do not silently drop the finite-radius boundary.
    outer=[]
    for n in [128,256,512]:
        t,w=roots_legendre(n);theta=t.astype(LD)*np.arccos(LD(-1))
        zz=R*np.exp(1j*theta)
        value=np.sum(w*(1/A(zz,L))*(zz/(zz-q)-zz/(zz-q0)))/2
        outer.append(complex(value))
    closure=pole+integral+outer[-1].real
    return dict(L=L,q=q,q0=q0,R=R,exact_difference=exact,pole_difference=pole,cut_difference=integral,reconstructed_difference=pole+integral,relative_error=abs(pole+integral-exact)/max(abs(exact),1e-12),quadrature_error_estimate=error,outer_boundary=[dict(nodes=n,real=v.real,imag=v.imag) for n,v in zip([128,256,512],outer)],relative_closure_error=abs(closure-exact)/max(abs(exact),1e-12))

def main():
    pole_rows=[];counts=[];disp=[];asymptotic=[]
    for L in [3,5,10,30,100]:
        ps=poles(L);pole_rows.append(dict(L=L,poles=ps))
        for n in [512,1024,2048]:counts.append(count_contour(L,-3.9999,-.001,.2,n))
        z,c=spectrum(L);limit=3*C/8*np.dot(c,z*z*np.log(z))
        for factor in [100,1000,5000]:
            R=float(factor*L)
            asymptotic.append(dict(L=L,q=R,A=float(A(R,L).real),A_limit=float(limit),relative_difference=float(abs(A(R,L).real/limit-1))))
            for q in [.5,2.,10.]:disp.append(reconstruct(q,L,R,ps))
    result=dict(scope='Fixed free scalar PV auxiliary determinant only. Finite-cut inverse dispersion; no global completeness or physical unitarity claim.',versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,longdouble_mantissa_bits=np.finfo(LD).nmant),pole_rows=pole_rows,subthreshold_counts=counts,dispersion=disp,asymptotic=asymptotic)
    assert [len(p['poles'])-1 for p in pole_rows]==[1,1,1,0,0]
    assert [x['count'] for x in counts]==[1]*9+[0]*6
    assert max(x['max_phase_step'] for x in counts)<.25
    assert max(x['relative_error'] for x in disp if x['R']==5000*x['L'])<1e-4
    assert max(x['relative_closure_error'] for x in disp)<1e-7
    (Path(__file__).resolve().parent/'results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(poles=pole_rows,max_phase_step=max(x['max_phase_step'] for x in counts),max_relative_closure_error=max(x['relative_closure_error'] for x in disp),max_relative_error_by_R={str(k):max(x['relative_error'] for x in disp if x['R']==k*x['L']) for k in [100,1000,5000]}),indent=2))

if __name__=='__main__':main()
