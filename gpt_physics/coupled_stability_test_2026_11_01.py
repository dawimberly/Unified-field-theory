"""Static coupled Candidate-0 stability gate. No gravity or quantum collapse assumed."""
import sys,json,time,argparse
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
threadpool_limits(limits=1)
sys.path.insert(0,r'C:\workspace\physics\scripts')
import knotted_cloud_dynamics as kd
OUT=Path(__file__).resolve().parent/'coupled_stability_2026_11_01'
OUT.mkdir(exist_ok=True)

def eg(r,n,h,lam=64.,kap=1.):
    g,b,f=kd.geometry(n,h)
    gs=sum(kd.dot(gi,gi) for gi in g)
    dens=.5*sum(kd.D(r,i,h)**2 for i in range(3))+.5*r*r*gs+lam/4*(r*r-1)**2+.5*kap*sum(fi*fi for fi in f.values())
    gr=sum(kd.DT(kd.D(r,i,h),i,h) for i in range(3))+r*gs+lam*r*(r*r-1)
    gn=kap*sum(f[i,j][...,None]*np.cross(g[i],g[j]) for i,j in f)
    ld=[-r[...,None]**2*g[i] for i in range(3)]
    for i,j in f:
        ld[i]-=kap*f[i,j][...,None]*b[j]
        ld[j]+=kap*f[i,j][...,None]*b[i]
    gn-=sum(kd.DT(ld[i],i,h) for i in range(3))
    return float(dens.sum()*h**3),gr*h**3,gn*h**3

def run(N,maxiter):
    start=time.time(); L=20.; h=L/N
    r,_,n,_=kd.seed(N,L,5.,7.)
    size=r.size
    def unpack(x):
        rr=x[:size].reshape(r.shape); raw=x[size:].reshape(n.shape)
        mag=np.linalg.norm(raw,axis=-1,keepdims=True)
        return rr,raw/mag,mag
    def objective(x):
        rr,nn,mag=unpack(x); E,gr,gn=eg(rr,nn,h)
        return E,np.r_[gr.ravel(),(kd.project(gn,nn)/mag).ravel()]
    x0=np.r_[r.ravel(),n.ravel()]
    rng=np.random.default_rng(20261101)
    direction=rng.normal(size=x0.size); direction/=np.linalg.norm(direction)
    E0,g0=objective(x0); eps=1e-5
    fd=(objective(x0+eps*direction)[0]-objective(x0-eps*direction)[0])/(2*eps)
    error=abs(fd-g0@direction)
    assert error<2e-6,(fd,g0@direction,error)
    est=kd.evaluate((r,np.zeros_like(r),n,np.zeros_like(n)),h,64.,1.,False)[0].sum()*h**3
    assert abs(E0-est)<1e-10
    top0=kd.topology(n,h)
    history=[]; tick=time.time(); it=0
    def callback(x):
        nonlocal tick,it
        it+=1
        if it%20==0:
            E,g=objective(x); history.append(dict(iteration=it,energy=E,gradient_max=float(abs(g).max())))
        if time.time()-tick>15:
            print(json.dumps(dict(N=N,seconds=round(time.time()-start,1),**history[-1])) if history else str(it),flush=True)
            tick=time.time()
    bounds=[(.05,None)]*size+[(None,None)]*(3*size)
    opt=minimize(objective,x0,jac=True,method='L-BFGS-B',bounds=bounds,callback=callback,options=dict(maxiter=maxiter,maxcor=12,ftol=1e-14,gtol=2e-7,maxls=30))
    rr,nn,_=unpack(opt.x); E,gr,gn=eg(rr,nn,h)
    pgn=kd.project(gn,nn)
    stat=max(float(abs(gr).max()),float(abs(pgn).max()))/h**3
    # Finite directional curvature controls on the constrained manifold.
    probes=[]
    def pert_energy(dr,dn,eps):
        nt=nn+eps*dn; nt/=np.linalg.norm(nt,axis=-1,keepdims=True)
        return eg(rr+eps*dr,nt,h)[0]
    for label in ['amplitude','orientation','mixed','translation_x','internal_rotation_z']:
        dr=np.zeros_like(rr); dn=np.zeros_like(nn)
        if label in ['amplitude','mixed']: dr=rng.normal(size=rr.shape)
        if label in ['orientation','mixed']: dn=kd.project(rng.normal(size=nn.shape),nn)
        if label=='translation_x':
            dr=(np.roll(rr,-1,0)-np.roll(rr,1,0))/(2*h)
            dn=kd.project((np.roll(nn,-1,0)-np.roll(nn,1,0))/(2*h),nn)
        if label=='internal_rotation_z': dn=np.cross(np.array([0.,0.,1.]),nn)
        norm=np.sqrt((np.sum(dr*dr)+np.sum(dn*dn))*h**3)
        if norm<1e-12: continue
        dr/=norm; dn/=norm
        for ep in [1e-3,5e-4]:
            curvature=(pert_energy(dr,dn,ep)+pert_energy(dr,dn,-ep)-2*E)/ep**2
            probes.append(dict(mode=label,epsilon=ep,curvature=curvature))
    result=dict(N=N,L=L,h=h,lambda_=64.,kappa=1.,initial_energy=E0,final_energy=E,iterations=opt.nit,optimizer_success=bool(opt.success),message=str(opt.message),minimum_rho=float(rr.min()),rho_floor=.05,floor_active=bool(rr.min()<.05001),stationarity_max_force=stat,initial_topology=top0,final_topology=kd.topology(nn,h),gradient_directional_error=error,energy_crosscheck_error=abs(E0-est),curvature_probes=probes,history=history,elapsed_seconds=time.time()-start)
    np.savez_compressed(OUT/f'N{N}_stationary_candidate.npz',rho=rr,n=nn)
    (OUT/f'N{N}_results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ['history','curvature_probes']}),flush=True)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--N',type=int,default=32);p.add_argument('--iterations',type=int,default=600)
    a=p.parse_args();run(a.N,a.iterations)

