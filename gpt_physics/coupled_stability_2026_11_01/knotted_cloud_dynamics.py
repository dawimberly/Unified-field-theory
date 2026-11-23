"""Constrained real-time Candidate-0 rho,n dynamics, theta constant.

The spatial lattice Lagrangian is varied before time integration. The
Skyrme velocity matrix and its time derivative are retained. No damping.
"""
from pathlib import Path
import argparse,json,time,zipfile
import numpy as np
from scipy.fft import fftn,ifftn,fftfreq
from threadpoolctl import threadpool_limits
threadpool_limits(limits=1)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'experiments'/'knotted_cloud_dynamics'
OUT.mkdir(parents=True,exist_ok=True)
checks=[]

def D(a,i,h): return (np.roll(a,-1,axis=i)-a)/h
def DT(a,i,h): return (np.roll(a,1,axis=i)-a)/h
def dot(a,b): return np.sum(a*b,axis=-1)
def project(a,n): return a-dot(a,n)[...,None]*n
def norm(a): return np.sqrt(dot(a,a))

def geometry(n,h):
    g=[D(n,i,h) for i in range(3)]
    b=[np.cross(gi,n) for gi in g]
    f={(i,j):dot(n,np.cross(g[i],g[j])) for i in range(3) for j in range(i+1,3)}
    return g,b,f

def evaluate(state,h,lam,kappa,need_acc=True):
    rho,p,n,u=state
    g,b,f=geometry(n,h)
    ft=[dot(u,bi) for bi in b]
    rho2=rho*rho
    kinetic=.5*p*p+.5*rho2*dot(u,u)+.5*kappa*sum(fi*fi for fi in ft)
    potential=lam/4*(rho2-1)**2+.5*rho2*sum(dot(gi,gi) for gi in g)
    potential+=.5*sum(D(rho,i,h)**2 for i in range(3))
    potential+=.5*kappa*sum(fi*fi for fi in f.values())
    if not need_acc: return kinetic+potential,kinetic,potential
    # Explicit partial L / partial n, holding n_dot and D_i n fixed.
    ln=kappa*sum(ft[i][...,None]*np.cross(u,g[i]) for i in range(3))
    ln-=kappa*sum(f[i,j][...,None]*np.cross(g[i],g[j]) for i,j in f)
    # partial L / partial(D_i n); exact adjoint of the forward difference.
    ld=[-rho2[...,None]*g[i]+kappa*ft[i][...,None]*np.cross(n,u) for i in range(3)]
    for i,j in f:
        ld[i]-=kappa*f[i,j][...,None]*b[j]
        ld[j]+=kappa*f[i,j][...,None]*b[i]
    rhs=ln+sum(DT(ld[i],i,h) for i in range(3))
    bt=[np.cross(D(u,i,h),n)+np.cross(g[i],u) for i in range(3)]
    mdotu=2*(rho*p)[...,None]*u
    mdotu+=kappa*sum(dot(u,bt[i])[...,None]*b[i]+ft[i][...,None]*bt[i] for i in range(3))
    matrix=rho2[...,None,None]*np.eye(3)
    matrix+=kappa*sum(bi[..., :,None]*bi[...,None,:] for bi in b)
    tangent=project(rhs-mdotu,n)
    accn=np.linalg.solve(matrix,tangent[...,None])[...,0]-dot(u,u)[...,None]*n
    accr=-sum(DT(D(rho,i,h),i,h) for i in range(3))
    accr+=rho*(dot(u,u)-sum(dot(gi,gi) for gi in g))-lam*rho*(rho2-1)
    return (p,accr,u,accn)

def clean(s):
    rho,p,n,u=s
    n=n/norm(n)[...,None]
    return rho,p,n,project(u,n)

def rk4(s,dt,h,lam,kappa):
    def add(a,b,k): return clean(tuple(x+k*y for x,y in zip(a,b)))
    k1=evaluate(s,h,lam,kappa)
    k2=evaluate(add(s,k1,dt/2),h,lam,kappa)
    k3=evaluate(add(s,k2,dt/2),h,lam,kappa)
    k4=evaluate(add(s,k3,dt),h,lam,kappa)
    return clean(tuple(x+dt/6*(a+2*b+2*c+d) for x,a,b,c,d in zip(s,k1,k2,k3,k4)))

def quaternion(x,y,z,R):
    r=np.sqrt(x*x+y*y+z*z)
    s=np.minimum(r/R,1)
    f=np.pi*(1-10*s**3+15*s**4-6*s**5)
    vec=np.stack(np.broadcast_arrays(x,y,z),axis=-1)
    fac=np.sin(f)/np.where(r>0,r,1)
    return np.concatenate([np.cos(f)[...,None],vec*fac[...,None]],axis=-1)

def multiply(q,r):
    s=q[...,0]*r[...,0]-dot(q[...,1:],r[...,1:])
    v=q[...,:1]*r[...,1:]+r[...,:1]*q[...,1:]+np.cross(q[...,1:],r[...,1:])
    return np.concatenate([s[...,None],v],axis=-1)

def seed(N,L,R,separation,pair=False,phase=0):
    a=(np.arange(N)-N//2)*L/N
    x,y,z=np.meshgrid(a,a,a,indexing='ij',sparse=True)
    q=quaternion(x+(separation/2 if pair else 0),y,z,R)
    if pair:
        r=quaternion(x-separation/2,y,z,R)
        # Conjugate the second lift by an internal rotation about z.
        c=np.cos(phase); s=np.sin(phase)
        rx=r[...,1].copy(); ry=r[...,2].copy()
        r[...,1]=c*rx-s*ry; r[...,2]=s*rx+c*ry
        q=multiply(q,r)
    w,qx,qy,qz=[q[...,i] for i in range(4)]
    n=np.stack([2*(qx*qz+w*qy),2*(qy*qz-w*qx),1-2*(qx*qx+qy*qy)],axis=-1)
    n/=norm(n)[...,None]
    return np.ones((N,N,N)),np.zeros((N,N,N)),n,np.zeros_like(n)

def topology(n,h):
    # Independent centered-difference + Fourier Coulomb-gauge estimator.
    # This is approximate; a conserved integer is NOT imposed.
    g=[(np.roll(n,-1,axis=i)-np.roll(n,1,axis=i))/(2*h) for i in range(3)]
    B=np.stack([dot(n,np.cross(g[1],g[2])),dot(n,np.cross(g[2],g[0])),dot(n,np.cross(g[0],g[1]))],axis=-1)
    k=2*np.pi*fftfreq(n.shape[0],h)
    kvec=np.stack(np.meshgrid(k,k,k,indexing='ij'),axis=-1)
    k2=dot(kvec,kvec); denom=np.where(k2>0,k2,1)
    bh=fftn(B,axes=(0,1,2))
    ah=1j*np.cross(kvec,bh)/denom[...,None]
    A=ifftn(ah,axes=(0,1,2)).real
    q=float(np.sum(A*B)*h**3/(16*np.pi**2))
    longitudinal=np.sum(abs(np.sum(kvec*bh,axis=-1))**2/denom)
    rel=float(np.sqrt(longitudinal/max(np.sum(abs(bh)**2),1e-30)))
    nd=min(float(np.min(dot(n,np.roll(n,-1,axis=i)))) for i in range(3))
    return dict(hopf_estimate=q,longitudinal_B_fraction=rel,min_neighbor_dot=nd)

def verify(name,ok,value):
    checks.append(dict(name=name,passed=bool(ok),value=float(value)))
    if not ok: raise AssertionError((name,value))

def selftest():
    rng=np.random.default_rng(15); N=5; h=.8; lam=3.; kap=.7
    rho=1+.02*rng.normal(size=(N,N,N)); p=.03*rng.normal(size=rho.shape)
    n=rng.normal(size=(N,N,N,3)); n/=norm(n)[...,None]
    u=project(.03*rng.normal(size=n.shape),n)
    state=(rho,p,n,u); deriv=evaluate(state,h,lam,kap)
    # dE/dt must vanish for the semi-discrete Euler-Lagrange acceleration.
    eps=1e-6
    plus=clean(tuple(x+eps*y for x,y in zip(state,deriv)))
    minus=clean(tuple(x-eps*y for x,y in zip(state,deriv)))
    rate=(evaluate(plus,h,lam,kap,False)[0].sum()-evaluate(minus,h,lam,kap,False)[0].sum())/(2*eps)
    verify('semi-discrete energy derivative',abs(rate)<2e-6,abs(rate))
    constraint=np.max(abs(dot(n,deriv[3])+dot(u,u)))
    verify('second derivative unit-vector constraint',constraint<1e-11,constraint)
    # Static rho and tangent-n derivatives versus the equation's force.
    static=(rho,0*p,n,0*u); acc=evaluate(static,h,lam,kap)
    dr=rng.normal(size=rho.shape); dn=project(rng.normal(size=n.shape),n)
    plus=clean((rho+eps*dr,0*p,n+eps*dn,0*u))
    minus=clean((rho-eps*dr,0*p,n-eps*dn,0*u))
    fd=(evaluate(plus,h,lam,kap,False)[0].sum()-evaluate(minus,h,lam,kap,False)[0].sum())/(2*eps)
    g,b,f=geometry(n,h)
    ma=rho[...,None]**2*acc[3]+kap*sum(dot(bi,acc[3])[...,None]*bi for bi in b)
    exact=-np.sum(dr*acc[1])-np.sum(dn*ma)
    verify('static variational force directional derivative',abs(fd-exact)<2e-5,abs(fd-exact))
    print(json.dumps({'selftests':checks}),flush=True)

def diagnose(N,L,R,sep,lam,kap):
    for pair in [False,True]:
        s=seed(N,L,R,sep,pair)
        e,k,v=evaluate(s,L/N,lam,kap,False)
        print(json.dumps(dict(N=N,pair=pair,E=float(e.sum()*(L/N)**3),**topology(s[2],L/N))),flush=True)

def run(label,N=40,L=20.,R=5.,sep=7.,lam=64.,kap=1.,dt=.015,end=1.,pair=False,phase=0.):
    h=L/N; s=seed(N,L,R,sep,pair,phase); steps=round(end/dt); dt=end/steps
    records=[]; initialtop=topology(s[2],h); start=time.time()
    a=(np.arange(N)-N//2)*h
    x,y,z=np.meshgrid(a,a,a,indexing='ij',sparse=True)
    rr=x*x+y*y+z*z
    inside=((x+sep/2)**2+y*y+z*z<3**2)|((x-sep/2)**2+y*y+z*z<3**2) if pair else rr<3**2
    edge=(abs(x)>L/2-2)|(abs(y)>L/2-2)|(abs(z)>L/2-2)
    for step in range(steps+1):
        if step%max(1,steps//40)==0 or step==steps:
            e,k,v=evaluate(s,h,lam,kap,False); E=float(e.sum()*h**3)
            qmat=np.array([[float((e*x*x).sum()),float((e*x*y).sum()),float((e*x*z).sum())],
                [float((e*x*y).sum()),float((e*y*y).sum()),float((e*y*z).sum())],
                [float((e*x*z).sum()),float((e*y*z).sum()),float((e*z*z).sum())]])*h**3
            qmat-=np.eye(3)*np.trace(qmat)/3
            records.append(dict(t=step*dt,energy=E,kinetic=float(k.sum()*h**3),
                minimum_rho=float(s[0].min()),unit_error=float(np.max(abs(norm(s[2])-1))),
                tangent_error=float(np.max(abs(dot(s[2],s[3])))),
                localized_fraction=float(e[inside].sum()*h**3/E),edge_fraction=float(e[edge].sum()*h**3/E),
                rms_radius=float(np.sqrt((e*rr).sum()*h**3/E)),quadrupole=qmat.tolist()))
            if s[0].min()<.1: raise RuntimeError('rho approaches polar boundary; trajectory not accepted')
        if step==steps: break
        s=rk4(s,dt,h,lam,kap)
        if (step+1)%10==0:
            print(json.dumps(dict(label=label,step=step+1,steps=steps,seconds=time.time()-start)),flush=True)
    finaltop=topology(s[2],h)
    drift=max(abs(r['energy']/records[0]['energy']-1) for r in records)
    verify(label+' energy drift',drift<.002,drift)
    verify(label+' field constraint',max(r['unit_error'] for r in records)<1e-12,max(r['unit_error'] for r in records))
    verify(label+' tangent velocity',max(r['tangent_error'] for r in records)<1e-11,max(r['tangent_error'] for r in records))
    result=dict(label=label,N=N,L=L,R=R,separation=sep,lambda_=lam,kappa=kap,dt=dt,end=end,pair=pair,phase=phase,
        initial_topology=initialtop,final_topology=finaltop,max_energy_drift=drift,records=records,seconds=time.time()-start)
    (OUT/(label+'.json')).write_text(json.dumps(result,indent=2))
    np.savez_compressed(OUT/(label+'_final.npz'),rho=s[0],n=s[2])
    print(json.dumps({k:result[k] for k in ['label','seconds','max_energy_drift','initial_topology','final_topology']}),flush=True)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--mode',default='check'); parser.add_argument('--N',type=int,default=40)
    parser.add_argument('--label',default='single'); parser.add_argument('--pair',action='store_true'); parser.add_argument('--phase',type=float,default=0)
    parser.add_argument('--dt',type=float,default=.015); parser.add_argument('--end',type=float,default=1.)
    parser.add_argument('--kappa',type=float,default=1.); parser.add_argument('--lambda',dest='lam',type=float,default=64.)
    args=parser.parse_args(); selftest()
    if args.mode=='diagnose': diagnose(args.N,20.,5.,7.,64.,1.)
    if args.mode=='run': run(args.label,N=args.N,pair=args.pair,phase=args.phase,dt=args.dt,end=args.end,kap=args.kappa,lam=args.lam)
