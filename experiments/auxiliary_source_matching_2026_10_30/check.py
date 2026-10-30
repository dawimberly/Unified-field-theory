"""Algebraic auxiliary-metric Hessian audit on a definite real branch.
Run: python experiments/auxiliary_source_matching_2026_10_30/check.py
Standard library only. This does not define a quantum integration cycle.
"""
from pathlib import Path
from decimal import Decimal,localcontext
from fractions import Fraction as F
import json,platform

def det(a):
    a=[[F(v) for v in row] for row in a];n=len(a);ans=F(1)
    for k in range(n):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:return F(0)
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];ans=-ans
        p=a[k][k];ans*=p
        for i in range(k+1,n):
            f=a[i][k]/p
            for j in range(k+1,n):a[i][j]-=f*a[k][j]
    return ans

def inv(a):
    n=len(a);b=[[F(v) for v in row]+[F(int(i==j)) for j in range(n)] for i,row in enumerate(a)]
    for k in range(n):
        pivot=next(i for i in range(k,n) if b[i][k]);b[k],b[pivot]=b[pivot],b[k]
        p=b[k][k];b[k]=[v/p for v in b[k]]
        for i in range(n):
            if i==k:continue
            f=b[i][k];b[i]=[a-f*c for a,c in zip(b[i],b[k])]
    return [row[n:] for row in b]

def dec(x):return Decimal(x.numerator)/Decimal(x.denominator)
def action_path(eps,h,V):
    g=[[F(int(i==j))+eps*h[i][j] for j in range(4)] for i in range(4)]
    # Fixed P=V I; the branch equation gives G=I in D=4.
    gi=inv(g);return dec(det(g)).sqrt()*(dec(F(V,2)*sum(gi[i][i] for i in range(4)))-Decimal(V))

def main():
    polarizations=[]
    polarizations.append(('trace',[[F(int(i==j)) for j in range(4)] for i in range(4)]))
    polarizations.append(('traceless_diagonal',[[F((1 if i==0 else -1 if i==1 else 0) if i==j else 0) for j in range(4)] for i in range(4)]))
    polarizations.append(('off_diagonal',[[F(int((i,j) in [(0,1),(1,0)])) for j in range(4)] for i in range(4)]))
    rows=[];determinants=[]
    with localcontext() as ctx:
        ctx.prec=70;eps=F(1,10**12)
        for V in [1,4,9]:
            for name,h in polarizations:
                trace=sum(h[i][i] for i in range(4));tr2=sum(h[i][j]*h[j][i] for i in range(4) for j in range(4))
                predicted=F(V,2)*(tr2-trace*trace/2)
                numeric=(action_path(eps,h,V)+action_path(-eps,h,V)-2*action_path(F(0),h,V))/dec(eps*eps)
                error=abs(numeric-dec(predicted));assert error<Decimal('1e-19')
                rows.append(dict(V=V,polarization=name,predicted_second_derivative=str(predicted),numeric_second_derivative=str(numeric),absolute_error=str(error)))
            # Orthonormal symmetric basis: 4 diagonals followed by 6 off-diagonals.
            v=[F(1)]*4+[F(0)]*6
            H=[[F(V,2)*(F(int(i==j))-v[i]*v[j]/2) for j in range(10)] for i in range(10)]
            J=[[-2*x for x in row] for row in H]
            hd=det(H);jd=det(J)
            assert hd==-F(V,2)**10 and jd==-F(V)**10
            determinants.append(dict(V=V,hessian_det=str(hd),constraint_jacobian_det=str(jd),hessian_traceless_eigenvalue=str(F(V,2)),hessian_trace_eigenvalue=str(-F(V,2)),stationary_phase_absolute_factor=str(F(2,V)**5)))
    result=dict(scope='D=4 algebraic S0=sqrt(det g)[tr(g^-1 P)/2-V], P=V I, definite G=I. No gauge, continuum, quantum vacuum or Lorentzian contour claim.',versions=dict(python=platform.python_version()),directional_hessians=rows,determinants=determinants,exact_general_results=dict(hessian='sqrt(det G)*V/(D-2)*[tr(G^-1 h G^-1 j)-tr(G^-1 h)tr(G^-1 j)/2]',metric_shift_jacobian='1 for the finite-dimensional triangular map (phi,h)->(phi,g=G(phi)+h)',raw_Gram_signature='P=B^T B is positive semidefinite for real positive-flavor gradients; full-rank P is positive definite'))
    (Path(__file__).resolve().parent/'results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(directional_checks=len(rows),max_absolute_error=max(float(r['absolute_error']) for r in rows),determinants=determinants),indent=2))

if __name__=='__main__':main()
