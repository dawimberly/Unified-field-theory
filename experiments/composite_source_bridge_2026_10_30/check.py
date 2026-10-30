"""Finite-dimensional exact Schur/source algebra, not a convergent path integral.
Run python experiments/composite_source_bridge_2026_10_30/check.py.
"""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util,json,platform
p=Path(__file__).resolve().parents[1]/'auxiliary_source_matching_2026_10_30/check.py'
s=importlib.util.spec_from_file_location('rational_helpers',p)
h=importlib.util.module_from_spec(s);s.loader.exec_module(h)
def eye(n):return [[Q(i==j) for j in range(n)] for i in range(n)]
def add(A,B):return [[a+b for a,b in zip(x,y)] for x,y in zip(A,B)]
def scale(A,c):return [[c*a for a in x] for x in A]
def mul(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def err(A,B):return max(abs(a-b) for x,y in zip(A,B) for a,b in zip(x,y))
def main():
 rows=[]
 for D in [3,4,5,6]:
  n=D*(D+1)//2;I=eye(n);v=[Q(i<D) for i in range(n)]
  P=[[Q(D-2,2)*I[i][j]-v[i]*v[j]/2 for j in range(n)] for i in range(n)]
  Pi=[[Q(2,D-2)*I[i][j]-v[i]*v[j]/(D-2) for j in range(n)] for i in range(n)]
  assert mul(P,Pi)==I
  # Symmetric positive B with diagonal and rank-one parts; it does not commute with P.
  u=[Q((i%3)-1,10) for i in range(n)]
  B=[[Q(10+i,10)*I[i][j]+u[i]*u[j] for j in range(n)] for i in range(n)]
  comm=err(mul(P,B),mul(B,P));assert comm>0
  for V in [11,17]:
   g=scale(P,Q(1,2*V));H=scale(Pi,Q(V,2))
   reduced=h.inv(add(h.inv(B),scale(P,Q(-1,2*V))))
   composite=scale(mul(mul(P,reduced),P),Q(1,V*V))
   auxiliary=h.inv(add(H,scale(B,Q(-1,4))))
   contact=scale(P,Q(2,V))
   # g+g B g+... has inverse g^-1-B. All identities retain noncommuting order.
   amplitude=h.inv(add(h.inv(g),scale(B,Q(-1))))
   e1=err(add(auxiliary,scale(composite,Q(-1))),contact)
   e2=err(auxiliary,scale(amplitude,Q(4)))
   expectedH=[[Q(V,D-2)*(I[i][j]-v[i]*v[j]/2) for j in range(n)] for i in range(n)]
   e3=err(H,expectedH)
   assert e1==e2==e3==0
   # Source quadratic forms for three rational nonconserved test vectors.
   forms=[]
   for shift in [0,1,2]:
    j=[[Q(((i+shift)%5)-2,7)] for i in range(n)]
    jt=list(map(list,zip(*j)))
    diff=mul(mul(jt,add(auxiliary,scale(composite,Q(-1)))),j)[0][0]
    exp=mul(mul(jt,contact),j)[0][0];assert diff==exp
    forms.append(str(diff))
   rows.append(dict(D=D,symmetric_components=n,V=V,noncommuting_test=str(comm),auxiliary_minus_composite_error=str(e1),auxiliary_vs_four_amplitude_error=str(e2),hessian_normalization_error=str(e3),source_contact_quadratic_forms=forms))
 result=dict(scope='Exact rational matrix Schur identities for formal Euclidean-sign quadratic algebra. No auxiliary integration cycle, positive measure, constrained CEV equivalence or source-gauge prescription is defined.',versions=dict(python=platform.python_version()),cases=rows,identities=dict(P='(D-2) I_sym/2 - vv^T/2',H='V P^-1/2',reduced_covariance='(B^-1-P/(2V))^-1',composite_covariance='P reduced_covariance P/V^2',auxiliary_covariance='(H-B/4)^-1',source_contact='auxiliary_covariance-composite_covariance=2P/V',chain_amplitude='(2V P^-1-B)^-1=auxiliary_covariance/4'))
 (Path(__file__).resolve().parent/'results.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(dict(exact_cases=len(rows),source_forms=3*len(rows),maximum_errors='0 (rational arithmetic)'),indent=2))
if __name__=='__main__':main()
