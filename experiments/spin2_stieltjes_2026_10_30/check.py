"""Spin-2 determinant audit and positive-measure representation.
Run from repository root: python experiments/spin2_stieltjes_2026_10_30/check.py
Reference tensor definitions are loaded without running their report generation.
"""
from pathlib import Path
import json, platform
import numpy as np
import scipy
from scipy.special import roots_legendre
from scipy.integrate import quad
LD=np.longdouble;CD=np.clongdouble;PI=np.arccos(LD(-1));C=1/(4*PI)**2

def spectrum(L):
    z=np.array([1,L,2*L,5*L+1],dtype=LD)
    w=np.array([1/np.prod(z[i]-np.delete(z,i)) for i in range(4)])
    return z,w/w[0],w[0]

def kernel(q,L):
    z,c,_=spectrum(L);q=np.asarray(q,dtype=CD)
    r=np.sqrt(q[...,None]/(q[...,None]+4*z))
    I=np.log(z)-2+2*np.arctanh(r)/r;d=z+q[...,None]/4
    return C/15*np.sum(c*(d*d*I-(z*z-z*q[...,None]/8)*np.log(z)),axis=-1)

def direct(q,L,n):
    z,c,_=spectrum(L);x,w=roots_legendre(n);x=(x+1)/2;w=w/2
    delta=z[:,None]+q*x*(1-x)
    return C/8*np.dot(c,(delta*delta*np.log(delta))@w-z*z*np.log(z))

def simplex(q,L,n):
    z,_,w0=spectrum(L);x,w=roots_legendre(n);x=(x+1)/2;w=w/2
    a,b,d=np.meshgrid(x,x,x,indexing='ij');wa,wb,wd=np.meshgrid(w,w,w,indexing='ij')
    Z=((1-a)*z[0]+a*(1-b)*z[1]+a*b*(1-d)*z[2]+a*b*d*z[3]).ravel()
    weight=(wa*wb*wd*a*a*b).ravel();total=0j
    for xx,ww in zip(x,w):
        t=xx*(1-xx);total+=ww*np.sum(weight*t/(Z*(Z+t*q)))
    return -C/(4*w0)*total

def im_kernel(s,L):
    z,c,_=spectrum(L);s=LD(s);beta=np.sqrt(np.maximum(LD(0),1-4*z/s))
    return -C*PI*s*s/240*np.dot(c,beta**5)

def contour(L,R,h,n):
    v=[-R+1j*h,R+1j*h,R+1j*R,-R+1j*R,-R+1j*h]
    q=np.concatenate([np.linspace(a,b,n,endpoint=False) for a,b in zip(v[:-1],v[1:])]+[np.array(v[-1:])])
    for iteration in range(22):
        f=kernel(q,L)/q;m=(q[:-1]+q[1:])/2;fm=kernel(m,L)/m
        p1=np.angle(fm/f[:-1]);p2=np.angle(f[1:]/fm)
        bad=(abs(p1)>.12)|(abs(p2)>.12)|(abs(p1-p2)>.03)
        bad|=abs(f[1:]-f[:-1])>.3*np.minimum(abs(f[1:]),abs(f[:-1]))
        if not np.any(bad):break
        new=[]
        for j in range(len(m)):
            new.append(q[j])
            if bad[j]:new.append(m[j])
        new.append(q[-1]);q=np.array(new)
    else:raise RuntimeError('unresolved contour')
    f=kernel(q,L)/q;phase=np.angle(f[1:]/f[:-1]);winding=float(np.sum(phase)/(2*PI))
    return dict(L=L,R=R,h=h,initial_points_per_edge=n,final_points=len(q),winding=winding,max_phase_step=float(max(abs(phase))))

def dispersion(q,L):
    z,_,_=spectrum(L);edges=[float(i) for i in 4*z]+[np.inf];total=0.;error=0.
    for lo,hi in zip(edges[:-1],edges[1:]):
        value,e=quad(lambda s:float(-im_kernel(s,L)/PI)/(s*(s+q)),lo,hi,epsabs=1e-11,epsrel=1e-9,limit=300)
        total+=value;error+=e
    exact=float(kernel(q,L).real)
    return dict(L=L,q=q,exact=exact,reconstructed=q*total,relative_error=abs(q*total-exact)/abs(exact),quadrature_error=q*error)

def main():
    root=Path(__file__).resolve().parents[2]
    source=(root/'scripts/metric_tensor_ward_check.py').read_text().split('rows=[]; roots=[]')[0]
    source=source.replace('OUT.mkdir(parents=True,exist_ok=True)','')
    ns={'__file__':str(root/'scripts/metric_tensor_ward_check.py')};exec(compile(source,'reference_tensor_definitions','exec'),ns)
    tensor=[];offaxis=[];measure=[];contours=[];dispersions=[];signs=[];low=[]
    for L in [3,5,10,30,100]:
        z,c,_=spectrum(L)
        for q in [-3.9,-2.,.5,2.,10.]:
            full=float(ns['tensor'](q,L,640)[4,4]);k=float(kernel(q,L).real)
            tensor.append(dict(L=L,q=q,full_tensor=full,compact_kernel=k,absolute_error=abs(k-full)))
        for q in [1+2j,-2+.3j,-8+2j,-50+5j]:
            exact=complex(kernel(q,L));d=direct(q,L,1024)
            offaxis.append(dict(L=L,q=[q.real,q.imag],relative_error=float(abs(exact-d)/abs(exact))))
        for q in [1+2j,-2+.3j]:
            values=[complex(simplex(q,L,n)) for n in [12,24,48]];exact=complex(kernel(q,L)/q)
            measure.append(dict(L=L,q=[q.real,q.imag],exact=[exact.real,exact.imag],relative_errors=[abs(v-exact)/abs(exact) for v in values],nodes=[12,24,48]))
        for R in [4*L,100*L]:
            for n in [512,1024,2048]:contours.append(contour(L,R,.005,n))
        for q in [.5,2.,10.]:dispersions.append(dispersion(q,L))
        vals=[float(-im_kernel(s,L)) for s in np.geomspace(4.00001,5000*L,400)]
        signs.append(dict(L=L,samples=400,min_minus_im_kernel=min(vals),negative_samples=sum(v<0 for v in vals)))
        slope=C*np.dot(c,z*np.log(z))/24;quadratic=C*np.dot(c,np.log(z))/240
        estimates=[float(kernel(h,L).real/h) for h in [.01,.005,.0025]]
        low.append(dict(L=L,linear=float(slope),quadratic=float(quadratic),finite_q_slopes=estimates,infinity_limit=float(-C/8*np.dot(c,z*z*np.log(z)))))
    result=dict(scope='Fixed free minimally coupled leading-N auxiliary determinant spin-2 coefficient; proof application is conditional on this model.',versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,longdouble_bits=np.finfo(LD).nmant),tensor_checks=tensor,offaxis_quadrature=offaxis,positive_measure_quadrature=measure,contours=contours,dispersion=dispersions,cut_sign_samples=signs,low_energy=low)
    (Path(__file__).resolve().parent/'results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(max_tensor_absolute=max(x['absolute_error'] for x in tensor),max_offaxis_relative=max(x['relative_error'] for x in offaxis),max_simplex_relative=max(x['relative_errors'][-1] for x in measure),max_dispersion_relative=max(x['relative_error'] for x in dispersions),contour_windings=sorted(set(round(x['winding']) for x in contours)),max_phase_step=max(x['max_phase_step'] for x in contours),negative_cut_samples=sum(x['negative_samples'] for x in signs)),indent=2))
    assert max(x['absolute_error'] for x in tensor)<2e-9
    assert max(x['relative_error'] for x in offaxis)<1e-8
    assert max(x['relative_errors'][-1] for x in measure)<1e-5
    assert max(x['relative_error'] for x in dispersions)<1e-7
    assert all(abs(x['winding'])<1e-8 and x['max_phase_step']<.25 for x in contours)
    assert all(x['negative_samples']==0 for x in signs)

if __name__=='__main__':main()
