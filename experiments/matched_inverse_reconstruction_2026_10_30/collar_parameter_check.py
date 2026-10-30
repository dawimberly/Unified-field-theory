"""Defining principal-log parameter quadrature at the independently certified pair.
Run after local_pole_certificate.py and check.py.
"""
from pathlib import Path
import importlib.util,json
import numpy as np
s=importlib.util.spec_from_file_location('reconstruction',Path(__file__).with_name('check.py'))
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
certificate=json.loads(Path(__file__).with_name('local_pole_certificate.json').read_text());q=complex(*certificate['center_float']);model=m.Model(30);rows=[]
for point in [q,q.conjugate()]:
 raw,error=m.th.principal_parameter(point,30)
 qq=m.CD(point);value=m.CD(raw)+model.d1*qq+model.d2*qq*qq
 assert abs(value)<1e-9
 rows.append(dict(q=[point.real,point.imag],parameter_kernel_value=[float(value.real),float(value.imag)],absolute_parameter_kernel_residual=float(abs(value)),quadrature_error_estimate=float(error),compact_kernel_residual=float(abs(model.A(qq))),independent_derivative_relative_error=float(abs(model.derivative_quad(qq)/model.derivative(qq)-1))))
Path(__file__).with_name('collar_parameter_results.json').write_text(json.dumps(dict(scope='Principal-log integrals on both conjugate first-sheet sides; no continuation through the cut.',cases=rows),indent=2));print(json.dumps(rows,indent=2))
