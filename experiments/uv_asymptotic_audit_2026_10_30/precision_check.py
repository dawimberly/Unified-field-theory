"""Expose cancellation in the compact large-q formula; verify the UV alternative.
Run after check.py: python experiments/uv_asymptotic_audit_2026_10_30/precision_check.py.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import localcontext
import importlib.util,json
s=importlib.util.spec_from_file_location('uv',Path(__file__).with_name('check.py'))
u=importlib.util.module_from_spec(s);s.loader.exec_module(u)
rows=[]
with localcontext() as ctx:
 ctx.prec=160
 for L in [3,5,10,30,100]:
  z,_=u.m.spec(L)
  for multiple in [100,10000,1000000,100000000]:
   R=multiple*z[-1];ref=u.reference(u.dec(R),L,'A')[0];series=u.approximation(u.dec(R),L,'A')[0]
   bound=u.remainder(L,'A',R)+u.coefficient_log_uncertainty(L,'A',R)
   assert abs(series-ref)<u.dec(bound)
   naive=float((u.m.base.base.A(int(R),L)/u.m.base.C).real)
   rows.append(dict(L=L,q=int(R),decimal_reference=float(ref),compact_longdouble=naive,compact_relative_error=float(abs(u.D(str(naive))-ref)/abs(ref)),uv_relative_error=str(abs((series-ref)/ref)),uv_total_error_bound=float(bound)))
Path(__file__).with_name('precision_results.json').write_text(json.dumps(dict(scope='Numerical conditioning audit only, not new model behavior.',cases=rows),indent=2))
print(json.dumps(dict(cases=len(rows),largest_compact_relative_error=max(x['compact_relative_error'] for x in rows),largest_uv_relative_error=max(float(x['uv_relative_error']) for x in rows)),indent=2))
# Audit only the old asymptotic-value rows affected by this new conditioning concern.
old=json.loads((u.ROOT/'experiments/inverse_scalar_dispersion_2026_10_30/results.json').read_text())
legacy=[]
with localcontext() as ctx:
 ctx.prec=160;pi=16*u.atan(u.D(1)/5)-4*u.atan(u.D(1)/239);C=1/(4*pi)**2
 for row in old['asymptotic']:
  correct=u.reference(u.D(str(row['q'])),row['L'],'A')[0]*C
  error=abs((u.D(str(row['A']))-correct)/correct)
  legacy.append(dict(L=row['L'],q=row['q'],archived_A=row['A'],high_precision_A=float(correct),relative_error=float(error)))
p=Path(__file__).with_name('precision_results.json');result=json.loads(p.read_text());result['legacy_asymptotic_value_audit']=legacy;p.write_text(json.dumps(result,indent=2))
print(json.dumps(dict(legacy_rows=len(legacy),maximum_legacy_relative_error=max(x['relative_error'] for x in legacy)),indent=2))
