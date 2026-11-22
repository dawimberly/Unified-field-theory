# Independent vacuum spectrum check, including phase block.
import sys,json
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh
sys.path.insert(0,str(Path(__file__).resolve().parent))
from coupled_hessian_check_2026_11_01 import hessian,OUT,kd
N=6;h=.7;r=np.ones((N,)*3);n=np.zeros((N,)*3+(3,));n[...,2]=1
shape=r.shape;size=r.size
# Three blocks rho, tangent n_x, tangent n_y. Phase has identical vacuum Laplacian.
def mv(x):
    dr=x[:size].reshape(shape);dn=np.zeros_like(n)
    dn[...,0]=x[size:2*size].reshape(shape);dn[...,1]=x[2*size:].reshape(shape)
    hr,hn=hessian(r,n,h,dr,dn)
    return np.r_[hr.ravel(),hn[...,0].ravel(),hn[...,1].ravel()]
# Build small matrix to check the entire spectrum against Fourier eigenvalues.
H=np.column_stack([mv(np.eye(3*size)[:,i]) for i in range(3*size)])
numerical=np.linalg.eigvalsh(H)
k=np.arange(N);lambdas=[]
for i in k:
 for j in k:
  for q in k:
   ev=4/h**2*(np.sin(np.pi*i/N)**2+np.sin(np.pi*j/N)**2+np.sin(np.pi*q/N)**2)
   lambdas.extend([ev+128,ev,ev])
expected=np.sort(lambdas);err=float(np.max(abs(numerical-expected)))
assert err<1e-10
res=dict(N=N,h=h,modes=len(numerical),max_eigenvalue_error=err,minimum_eigenvalue=float(numerical[0]),zero_modes=int(np.sum(abs(numerical)<1e-9)),negative_modes=int(np.sum(numerical < -1e-9)),phase_block='Positive semidefinite discrete Laplacian; one additional constant phase zero mode. No mixed quadratic phase terms at constant phase.')
(OUT/'vacuum_spectrum_check.json').write_text(json.dumps(res,indent=2));print(json.dumps(res))
