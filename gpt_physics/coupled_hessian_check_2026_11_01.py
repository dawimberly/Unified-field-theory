"""Exact constrained Hessian-vector check for the Candidate-0 lattice energy."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh,ArpackNoConvergence
from threadpoolctl import threadpool_limits
threadpool_limits(limits=1)
sys.path.insert(0,str(Path(__file__).resolve().parent))
from coupled_stability_test_2026_11_01 import eg,kd,OUT

def hessian(r,n,h,dr,dn,lam=64.,kap=1.):
    g,b,f=kd.geometry(n,h); dg=[kd.D(dn,i,h) for i in range(3)]
    db=[np.cross(dg[i],n)+np.cross(g[i],dn) for i in range(3)]
    df={(i,j):kd.dot(dn,np.cross(g[i],g[j]))+kd.dot(n,np.cross(dg[i],g[j])+np.cross(g[i],dg[j])) for i,j in f}
    dgr=sum(kd.DT(kd.D(dr,i,h),i,h) for i in range(3))
    dgr+=dr*sum(kd.dot(gi,gi) for gi in g)+2*r*sum(kd.dot(g[i],dg[i]) for i in range(3))+lam*(3*r*r-1)*dr
    dgn=kap*sum(df[i,j][...,None]*np.cross(g[i],g[j])+f[i,j][...,None]*(np.cross(dg[i],g[j])+np.cross(g[i],dg[j])) for i,j in f)
    dld=[-2*(r*dr)[...,None]*g[i]-r[...,None]**2*dg[i] for i in range(3)]
    for i,j in f:
        dld[i]-=kap*(df[i,j][...,None]*b[j]+f[i,j][...,None]*db[j])
        dld[j]+=kap*(df[i,j][...,None]*b[i]+f[i,j][...,None]*db[i])
    dgn-=sum(kd.DT(dld[i],i,h) for i in range(3))
    gn=eg(r,n,h,lam,kap)[2]/h**3
    return dgr,kd.project(dgn,n)-kd.dot(n,gn)[...,None]*dn

def check():
    rng=np.random.default_rng(2011);N=5;h=.7
    r=1+.02*rng.normal(size=(N,)*3);n=rng.normal(size=(N,)*3+(3,));n/=kd.norm(n)[...,None]
    dr=rng.normal(size=r.shape);dn=kd.project(rng.normal(size=n.shape),n)
    norm=np.sqrt((dr*dr).sum()+(dn*dn).sum());dr/=norm;dn/=norm
    hr,hn=hessian(r,n,h,dr,dn)
    ep=1e-5
    def grad_at(e):
        raw=n+e*dn;mag=kd.norm(raw)[...,None];nt=raw/mag
        E,gr,gn=eg(r+e*dr,nt,h)
        return gr/h**3,kd.project(gn,nt)/mag/h**3
    pr,pn=grad_at(ep);mr,mn=grad_at(-ep)
    errr=np.linalg.norm((pr-mr)/(2*ep)-hr)/max(1,np.linalg.norm(hr))
    errn=np.linalg.norm(kd.project((pn-mn)/(2*ep),n)-hn)/max(1,np.linalg.norm(hn))
    assert max(errr,errn)<1e-7
    dr2=rng.normal(size=r.shape);dn2=kd.project(rng.normal(size=n.shape),n)
    hr2,hn2=hessian(r,n,h,dr2,dn2)
    a=np.sum(dr2*hr)+np.sum(dn2*hn);b=np.sum(dr*hr2)+np.sum(dn*hn2)
    sym=abs(a-b)/max(1,abs(a),abs(b));assert sym<1e-12
    result=dict(amplitude_hvp_relative_error=errr,orientation_hvp_relative_error=errn,hessian_symmetry_relative_error=sym)
    (OUT/'hessian_implementation_checks.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))

if __name__=='__main__':check()
