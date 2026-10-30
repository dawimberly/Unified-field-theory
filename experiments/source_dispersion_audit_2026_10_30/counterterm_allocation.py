"""Exact block-Gaussian illustration of source ambiguity at fixed metric Hessian.
Not a new physical action, convergent path integral, or actual PV matching model.
"""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,platform
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('bridge',ROOT/'experiments/composite_source_bridge_2026_10_30/check.py');b=importlib.util.module_from_spec(s);s.loader.exec_module(b)
add=b.add;mul=b.mul;scale=b.scale;inv=b.h.inv
I=b.eye(2);H=[[F(11,2),F(0)],[F(0),F(-11,2)]];B=[[F(2),F(1,5)],[F(1,5),F(3)]]
R=[[F(1,3),F(1,7)],[F(1,7),F(2,5)]];S=[[F(1,11),F(0)],[F(0),F(1,13)]]
def difference(A,B):return add(A,scale(B,F(-1)))
def covariance(HH,BB):
 # Direct full 4x4 block inverse; compare its hh block to a separate Schur formula.
 Binv=inv(BB)
 full=[Binv[i]+[-I[i][j]/2 for j in range(2)] for i in range(2)]
 full += [[-I[i][j]/2 for j in range(2)]+HH[i] for i in range(2)]
 cov=inv(full);CT=[row[:2] for row in cov[:2]];Dh=[row[2:] for row in cov[2:]]
 Hinv=inv(HH);CO=scale(mul(mul(Hinv,CT),Hinv),F(1,4))
 assert Dh==inv(add(HH,scale(BB,F(-1,4))))
 assert difference(Dh,CO)==Hinv
 return Dh,CO,Hinv

def main():
 rows=[]
 for q in map(F,[0,1,2,3]):
  delta=add(scale(R,q),scale(S,q*q))
  DhA,COA,contactA=covariance(add(H,delta),B)
  DhB,COB,contactB=covariance(H,add(B,scale(delta,F(-4))))
  assert DhA==DhB
  assert difference(COA,COB)==difference(contactB,contactA)
  if q:assert COA!=COB and contactA!=contactB
  else:assert COA==COB
  forms=[]
  for j in [[[F(1,3)],[F(-2,7)]],[[F(1)],[F(0)]],[[F(0)],[F(1)]]]:
   jt=list(map(list,zip(*j)))
   form=mul(mul(jt,difference(COA,COB)),j)[0][0]
   contact=mul(mul(jt,difference(contactB,contactA)),j)[0][0];assert form==contact
   forms.append(str(form))
  rows.append(dict(q=str(q),same_metric_covariance=True,composite_covariances_equal=COA==COB,metric_term_assignment_contact=[[str(v) for v in row] for row in contactA],bubble_contact_assignment_contact=[[str(v) for v in row] for row in contactB],source_forms_of_composite_difference=forms))
 result=dict(scope='Exact finite-dimensional formal Schur algebra only. Two allocations with the same H-B/4 need not define the same composite/source functional. No actual regulated composite model or integration cycle is selected.',python=platform.python_version(),cases=rows,maximum_algebraic_residual='0 (Fraction arithmetic)')
 (HERE/'counterterm_results.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(dict(block_inverse_cases=8,source_forms=12,maximum_residual='0',nonzero_q_distinguishes_composites=True)))
if __name__=='__main__':main()
