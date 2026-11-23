"""Topology-preserving continuum variational test of Candidate 0.
Smooth Cayley rotations of a charge-one rational Hopf map; positive amplitude.
No discrete normalization singularities, gravity, or quantum collapse assumed.
"""
from pathlib import Path
import json,time,argparse,itertools
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
threadpool_limits(limits=1)
OUT=Path(__file__).resolve().parent/'smooth_hopf_2026_11_01';OUT.mkdir(exist_ok=True)

def exponents(order):
 return [e for total in range(order+1) for e in itertools.product(range(total+1),repeat=3) if sum(e)==total]

class Model:
 def __init__(self,nr=28,nmu=12,nphi=24,order=1,amp_order=None):
  self.order=order;self.exps=exponents(order)
  t,wt=leggauss(nr);t=(t+1)*np.pi/4;wt*=np.pi/4
  mu,wm=leggauss(nmu);ph=2*np.pi*np.arange(nphi)/nphi
  r=np.tan(t);ct=np.sqrt(1-mu*mu)
  rr,mm,pp=np.meshgrid(r,mu,ph,indexing='ij')
  self.y=np.stack([rr*np.sqrt(1-mm*mm)*np.cos(pp),rr*np.sqrt(1-mm*mm)*np.sin(pp),rr*mm],axis=-1).reshape(-1,3)
  self.w=np.broadcast_to((wt*r*r*(1+r*r))[:,None,None]*wm[None,:,None]*(2*np.pi/nphi),(nr,nmu,nphi)).ravel()
  y=self.y;r2=np.sum(y*y,axis=1);den=1+r2
  z1=2*(y[:,0]+1j*y[:,1])/den;z2=(2*y[:,2]+1j*(r2-1))/den
  dz1=np.zeros((len(y),3),complex);dz1[:,0]=2/den;dz1[:,1]=2j/den;dz1-=z1[:,None]*2*y/den[:,None]
  dz2=2j*y/den[:,None];dz2[:,2]+=2/den;dz2-=z2[:,None]*2*y/den[:,None]
  zz=z1*np.conj(z2)
  self.n=np.stack([2*zz.real,2*zz.imag,abs(z1)**2-abs(z2)**2],axis=-1)
  dz=dz1*np.conj(z2[:,None])+z1[:,None]*np.conj(dz2)
  self.G=np.stack([2*dz.real,2*dz.imag,2*(np.conj(z1[:,None])*dz1-np.conj(z2[:,None])*dz2).real],axis=-1)
  self.widths=np.array([.3,.6,1.2,2.4,4.8])
  self.H=np.exp(-r2[:,None]/self.widths[None,:]**2)
  self.dH=-2*y[:,None,:]*self.H[:,:,None]/self.widths[None,:,None]**2
  gauss=np.exp(-r2/4)
  bas=[];dbas=[]
  for e in self.exps:
   mono=np.prod(y**np.array(e),axis=1)
   dg=np.zeros_like(y)
   for j in range(3):
    if e[j]:
     ee=list(e);ee[j]-=1;dg[:,j]=e[j]*np.prod(y**np.array(ee),axis=1)
   bas.append(gauss*mono);dbas.append(gauss[:,None]*(dg-.5*y*mono[:,None]))
  self.B=np.stack(bas,axis=1);self.dB=np.stack(dbas,axis=1)
  self.amp_exps=[] if amp_order is None else [e for e in self.exps if sum(e)<=amp_order]
  if self.amp_exps:
   ids=[self.exps.index(e) for e in self.amp_exps]
   self.H=np.concatenate([self.H,self.B[:,ids]],axis=1);self.dH=np.concatenate([self.dH,self.dB[:,ids]],axis=1)
  self.na=self.H.shape[1];self.offset=1+self.na
  self.size=self.offset+3*len(self.exps)
  self.grid=[nr,nmu,nphi]
 def fields(self,p):
  c=p[1:self.offset];a=p[self.offset:].reshape(-1,3)
  rho=np.exp(-self.H@c);dr=-rho[:,None]*np.einsum('pbj,b->pj',self.dH,c)
  v=self.B@a;dv=np.einsum('pbj,bk->pjk',self.dB,a)
  vv=np.sum(v*v,axis=1);D=1+vv;vn=np.sum(v*self.n,axis=1)
  M=(1-vv)[:,None]*self.n+2*v*vn[:,None]+2*np.cross(v,self.n)
  n=M/D[:,None]
  vdv=np.einsum('pk,pjk->pj',v,dv)
  dvn=np.einsum('pjk,pk->pj',dv,self.n)+np.einsum('pk,pjk->pj',v,self.G)
  dM=-2*vdv[:,:,None]*self.n[:,None,:]+(1-vv)[:,None,None]*self.G
  dM+=2*dv*vn[:,None,None]+2*v[:,None,:]*dvn[:,:,None]
  dM+=2*np.cross(dv,self.n[:,None,:])+2*np.cross(v[:,None,:],self.G)
  G=dM/D[:,None,None]-M[:,None,:]*(2*vdv/D[:,None]**2)[:,:,None]
  return rho,dr,n,G
 def components(self,p):
  rho,dr,n,G=self.fields(p)
  e2=.5*np.sum(dr*dr,axis=1)+.5*rho*rho*np.sum(G*G,axis=(1,2))
  e4=0*rho
  for i,j in [(0,1),(0,2),(1,2)]:
   F=np.sum(n*np.cross(G[:,i],G[:,j]),axis=1);e4+=.5*F*F
  e0=.25*(rho*rho-1)**2
  return np.array([self.w@e2,self.w@e4,self.w@e0])
 def energy(self,p):
  R=np.exp(p[0]);A,B,C=self.components(p)
  return R*A+B/R+64*R**3*C
 def jac_cs(self,p):
  g=np.zeros(len(p));eps=1e-25
  for j in range(len(p)):
   z=p.astype(complex);z[j]+=1j*eps;g[j]=self.energy(z).imag/eps
  return g


 def jac(self,p):
  R=np.exp(p[0]);rho,dr,n,G=self.fields(p);v=self.B@p[self.offset:].reshape(-1,3);dv=np.einsum('pbj,bk->pjk',self.dB,p[self.offset:].reshape(-1,3));D=1+np.sum(v*v,axis=1)
  G2=np.sum(G*G,axis=(1,2));dr2=np.sum(dr*dr,axis=1)
  gG=R*rho[:,None,None]**2*G;gn=np.zeros_like(n);e4=np.zeros_like(rho)
  for i,j in [(0,1),(0,2),(1,2)]:
   F=np.sum(n*np.cross(G[:,i],G[:,j]),axis=1);e4+=.5*F*F
   gn+=F[:,None]*np.cross(G[:,i],G[:,j])/R
   gG[:,i]+=F[:,None]*np.cross(G[:,j],n)/R
   gG[:,j]+=F[:,None]*np.cross(n,G[:,i])/R
  torque=np.cross(n,gn)+np.sum(np.cross(G,gG),axis=1)
  cc=np.cross(n[:,None,:],gG)
  rotgrad=2/D[:,None]*(torque-np.cross(v,torque))
  tmp=cc-np.cross(v[:,None,:],cc)
  rotgrad+=np.sum(-4*np.einsum('pk,pjk->pj',v,dv)[:,:,None]/D[:,None,None]**2*tmp-2/D[:,None,None]*np.cross(dv,cc),axis=1)
  drotgrad=2/D[:,None,None]*tmp
  ga=self.B.T@(self.w[:,None]*rotgrad)+np.einsum('pbj,pjk,p->bk',self.dB,drotgrad,self.w)
  common=R*(dr2+rho*rho*G2)+64*R**3*rho*rho*(rho*rho-1)
  gc=-self.H.T@(self.w*common)-R*np.einsum('pj,pbj,p->b',dr,self.dH,self.w*rho)
  A=self.w@(.5*dr2+.5*rho*rho*G2);B=self.w@e4;C=self.w@(.25*(rho*rho-1)**2)
  return np.r_[R*A-B/R+192*R**3*C,gc,ga.ravel()]

def embed(p,old,new):
 q=np.zeros(new.size);q[:6]=p[:6]
 for i,e in enumerate(old.amp_exps):
  if e in new.amp_exps:q[6+new.amp_exps.index(e)]=p[6+i]
 for i,e in enumerate(old.exps):q[new.offset+3*new.exps.index(e):new.offset+3*new.exps.index(e)+3]=p[old.offset+3*i:old.offset+3*i+3]
 return q

def checks():
 m=Model(40,16,32,1);p=np.zeros(m.size);A,B,C=m.components(p)
 er=max(abs(A/(16*np.pi**2)-1),abs(B/(16*np.pi**2)-1),abs(C))
 assert er<2e-7,(A,B,C,er)
 rng=np.random.default_rng(89);p[1:6]=[.03,.02,.01,0,0];p[6:]=.05*rng.normal(size=m.size-6)
 r,dr,n,G=m.fields(p)
 unit=float(np.max(abs(np.sum(n*n,axis=1)-1)));tangent=float(np.max(abs(np.sum(n[:,None,:]*G,axis=2))))
 assert unit<2e-14 and tangent<2e-13
 g=m.jac(p);gcs=m.jac_cs(p);analytic_cs_error=float(np.max(abs(g-gcs)));assert analytic_cs_error<1e-8;d=rng.normal(size=len(p));d/=np.linalg.norm(d)
 fd=(m.energy(p+1e-5*d)-m.energy(p-1e-5*d))/2e-5
 graderr=abs(fd-g@d)/max(1,abs(fd));assert graderr<1e-7
 # At frozen undeformed amplitude, E(R)=16pi^2(R+1/R).
 radius=[]
 for R in [.5,1,2]:
  z=np.zeros(m.size);z[0]=np.log(R);err=abs(m.energy(z)/(16*np.pi**2*(R+1/R))-1);assert err<2e-7;radius.append(err)
 res=dict(frozen_hopf_integral_relative_error=float(er),unit_vector_max_error=unit,tangent_derivative_max_error=tangent,analytic_directional_gradient_relative_error=float(graderr),analytic_vs_complex_step_max_error=analytic_cs_error,radius_controls=radius)
 (OUT/'implementation_checks.json').write_text(json.dumps(res,indent=2));print(json.dumps(res),flush=True)

def optimize(m,p,iterations,label):
 start=time.time();history=[];tick=time.time();it=0
 def cb(x):
  nonlocal tick,it
  it+=1
  if it%10==0:history.append(dict(iteration=it,energy=float(m.energy(x))))
  if time.time()-tick>15:
   print(json.dumps(dict(label=label,iteration=it,energy=float(m.energy(x)),seconds=time.time()-start)),flush=True);tick=time.time()
 bounds=[(np.log(.2),np.log(5))]+[(-2.,2.)]*(m.size-1)
 opt=minimize(m.energy,p,jac=m.jac,bounds=bounds,method='L-BFGS-B',callback=cb,options=dict(maxiter=iterations,maxcor=15,ftol=1e-13,gtol=1e-5,maxls=25))
 p=opt.x;rho,_,_,_=m.fields(p);com=m.components(p);R=np.exp(p[0]);g=m.jac(p)
 res=dict(label=label,order=m.order,parameters=p.tolist(),grid=m.grid,parameter_count=m.size,energy=float(m.energy(p)),R=float(R),components=dict(quadratic=float(R*com[0]),quartic=float(com[1]/R),potential=float(64*R**3*com[2])),min_rho=float(rho.min()),max_rho=float(rho.max()),max_parameter_gradient=float(abs(g).max()),success=bool(opt.success),message=str(opt.message),iterations=opt.nit,history=history,elapsed_seconds=time.time()-start,bound_active=bool(any(abs(x-lo)<1e-5 or abs(x-hi)<1e-5 for x,(lo,hi) in zip(p,bounds))))
 np.savez_compressed(OUT/(label+'_parameters.npz'),p=p)
 (OUT/(label+'.json')).write_text(json.dumps(res,indent=2));print(json.dumps({k:v for k,v in res.items() if k not in ['parameters','history','components']}),flush=True)
 return p,res

def main():
 checks()
 # Radial-amplitude/radius-only control, with orientation frozen.
 m1=Model(order=1);p=np.zeros(m1.size)
 def basefun(z):
  q=np.zeros(m1.size);q[:6]=z;return m1.energy(q)
 def basejac(z):
  q=np.zeros(m1.size);q[:6]=z;return m1.jac(q)[:6]
 ob=minimize(basefun,p[:6],jac=basejac,method='L-BFGS-B',bounds=[(np.log(.2),np.log(5))]+[(-2,2)]*5,options=dict(maxiter=150,gtol=1e-6,ftol=1e-14))
 p[:6]=ob.x;bg=m1.jac(p)
 base=dict(energy=float(m1.energy(p)),R=float(np.exp(p[0])),max_amplitude_radius_gradient=float(abs(bg[:6]).max()),max_omitted_orientation_gradient=float(abs(bg[6:]).max()),success=bool(ob.success),parameters=p.tolist())
 (OUT/'radial_control.json').write_text(json.dumps(base,indent=2));print(json.dumps(base),flush=True)
 p,r1=optimize(m1,p,160,'order1')
 m2=Model(order=2);p2=embed(p,m1,m2);p2,r2=optimize(m2,p2,220,'order2')
 # Extra smooth modes test whether the optimized family is stationary in a larger family.
 m3=Model(36,16,32,3);p3=embed(p2,m2,m3);g3=m3.jac(p3)
 extra=[6+3*i+j for i,e in enumerate(m3.exps) if sum(e)==3 for j in range(3)]
 refinements=[]
 for grid in [(28,12,24),(40,18,36),(56,24,48)]:
  mf=Model(*grid,order=2);refinements.append(dict(grid=list(grid),energy=float(mf.energy(p2)),max_gradient=float(abs(mf.jac(p2)).max())))
 result=dict(radial_control_energy=base['energy'],order1_energy=r1['energy'],order2_energy=r2['energy'],energy_reduction_from_radial=base['energy']-r2['energy'],larger_family_energy=float(m3.energy(p3)),max_new_cubic_mode_gradient=float(abs(g3[extra]).max()),quadrature_refinement=refinements,topology='Exactly inherited charge-one class for all finite parameters: smooth globally nonsingular Cayley rotations tend to identity at infinity; rho is strictly positive. This is an analytic homotopy statement, not an approximate lattice charge estimate.',limitation='Finite variational family. Parameter stationarity does not establish full PDE stationarity or dynamical stability. Quadrature and basis refinement are independent gates.')
 (OUT/'summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
if __name__=='__main__':main()

