# FFff(B) restored vertex anchor

Follows FFFF_B_MIDPOINT_ANCHOR_AUDIT_2026-11-22.md. Write c for the overall arithmetic midpoint coordinate and F=F₁F₂. At first order in κ, using T for the spectator tensor including its exponential and vertex factors, the gradient part after the two integrations by parts with the original vertex-0 anchor is exactly

```
∫ du dv dw dc δ(c+(u−w)/2) [
  F ∂ui ∂vj T_ij − (1/2) ∂ci(F ∂vj T_ij)
].
```

This follows by applying ∂ui to the anchor delta and then moving ∂ci off the delta. The v derivative does not hit the anchor. Thus evaluate the first spectator derivative using x_k−x0 and the second via the total u derivative of those anchored spectator arguments. The leftover product term is −(1/2)(∂ci F) ∂vj T. At order κ only the scalar mass perturbation in F contributes to ∂ci F. With the paired endpoints already expressed relative to x0 and the flat Gaussian divided out, its κz multiplier is −∑_paired lengths s (x_p+x_q)_i. The E split still has no explicit dependence on c at this order.

I checked the restored operator by direct differentiation and polynomial Gaussian integration in d=2. For the two length tuples from the preceding checkpoint, the flat, mass-independent and mass-dependent expectations are respectively (0.07308345, −1.05620157, 2.15881834) and (0.47781325, −2.83722821, 2.99005423), agreeing with the original anchored covariance after E is removed within 1.4e−15. The anchor-only discrepancy is therefore explained at the kernel level. When translated back to the original anchoring, the midpoint method is an identity for the same analytically continued integration and does not add a separate pole. An independent near-four-dimensional Laurent extraction using only the differentiated polynomials remains open. The vertex-anchored −1/144 double residual and +L/48 simple-pole log are the two checks for subsequent bare graph and local counterterm audits.

Reproduce with the original FFFF_B_PROJECTIVE.md Python block saved as fffb_checkpoint.py in the same directory.

```python
import itertools
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from fffb_checkpoint import kernel

def midpoint(lengths, dim=2):
    # Normalized Gaussian expectation at fixed overall midpoint.
    a,b,c,d,e,f = lengths
    M=np.array([[1,0,0],[.5,.5,1],[0,1,0],[.5,.5,-1]])
    m=M[[0,1,2,3,0,2]]
    lower=np.linalg.cholesky(2*np.linalg.inv(m.T@np.diag(1/lengths)@m))
    nodes,weights=hermegauss(5)
    ids=np.array(list(itertools.product(range(5),repeat=3*dim)))
    wt=np.prod(weights[ids],axis=1)/(2*np.pi)**(3*dim/2)
    coords=np.einsum("ij,njk->nik",lower,nodes[ids].reshape(-1,3,dim))
    u,v,w=coords.transpose(1,0,2)
    g=lambda x,y: np.sum(x*y,axis=-1)
    cp=dim*(4-dim)/3
    qa=(1/a+1/e)/2
    qc=(1/c+1/f)/2
    qs=(1/b+1/d)/8
    # The spectator polynomial is degree at most four in translations.
    # Roots-of-unity coefficient extraction differentiates it without
    # finite-difference truncation. Differentiate the exponential analytically.
    phase=np.exp(2j*np.pi*np.arange(5)/5)
    def spectator(i,j,z1,z2):
        shift=np.eye(dim)[i]*z1+np.eye(dim)[j]*z2
        r=(u+v+shift)/2+w
        t=(u+v+shift)/2-w
        U=u+np.eye(dim)[i]*z1
        V=v+np.eye(dim)[j]*z2
        x0=(U-w)/2; x2=(-U-w)/2
        x1=(V+w)/2; x3=(-V+w)/2
        x1,x2,x3,x0=x1-x0,x2-x0,x3-x0,x0-x0
        H=lambda x,y,s: (g(x,x)+g(y,y)+g(x,y)+(dim-1)*s)/3
        T=(r[:,i]*t[:,j]+t[:,i]*r[:,j]-(i==j)*g(r,t))/(4*b*d)
        return np.array([T, -(4-dim)/2*sum(g(x,x) for x in (x0,x1,x2,x3))*T,
                         2*(b*H(x1,x2,b)+d*H(x0,x3,d))*T])
    def hs(r,w,s):
        return (3*g(w,w)/4+g(r,r)/4+dim*s)/3
    pairB=2*sum(s*hs(r,w,s) for r,s in ((u,a),(u,e),(v,c),(v,f)))
    x0=(u-w)/2
    pairB=pairB+2*((a+e)*(g(w,x0)+g(x0,x0))+(c+f)*(-g(w,x0)+g(x0,x0)))
    pair_derivative=-(a+e)*(-w-2*x0)-(c+f)*(w-2*x0)
    pairA=-cp*(a+e+c+f)
    vals=np.zeros(3)
    r=(u+v)/2+w; t=(u+v)/2-w
    q=(r/b+t/d)/4
    for i,j in itertools.product(range(dim),repeat=2):
        P=np.array([[spectator(i,j,p,k) for k in phase] for p in phase])
        coeff=lambda n,m: np.einsum("a,b,abkn->kn",phase**(-n),phase**(-m),P).real/25
        out=coeff(1,1)-q[:,i]*coeff(0,1)-q[:,j]*coeff(1,0)+(q[:,i]*q[:,j]-qs*(i==j))*coeff(0,0)
        vals += np.array([wt@out[0],wt@(out[1]+pairA*out[0]),wt@(out[2]+pairB*out[0]+pair_derivative[:,i]*(coeff(0,1)[0]-q[:,j]*coeff(0,0)[0]))])
    return vals/4

for lengths in (np.array([1.1,.8,1.3,.7,.5,1.9]),np.array([.03,.4,.9,.2,1.1,.7])):
    cp=2*2/3
    def direct(s):
        n,A,B=kernel(s[None,:],2,0)
        R=np.eye(4)[[0,1,1,0,0,1],1:]-np.eye(4)[[2,2,3,3,2,3],1:]
        U=np.prod(s)*np.linalg.det(R.T@np.diag(1/s)@R)
        return np.array([n[0],A[0]-cp*(s[0]+s[2])*n[0],B[0]])*U
    ds=np.mean([direct(lengths[list(p)]) for p in ((0,1,2,3,4,5),(4,1,2,3,0,5),(0,1,5,3,4,2),(4,1,5,3,0,2))],axis=0)
    result=midpoint(lengths)
    print(lengths,result,ds, np.max(abs(result-ds)))
    np.testing.assert_allclose(result,ds,atol=2e-12,rtol=2e-12)

```
