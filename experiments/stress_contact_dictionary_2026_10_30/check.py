"""Independent radial PV tadpoles and zero-momentum Wick stress contacts.
Run python experiments/stress_contact_dictionary_2026_10_30/check.py.
"""
from pathlib import Path
import importlib.util,json,platform
import numpy as np
import scipy
from scipy.integrate import quad
p=Path(__file__).resolve().parents[1]/'spin2_stieltjes_2026_10_30/check.py'
s=importlib.util.spec_from_file_location('spin_helpers',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def main():
 rows=[]
 for L in [3,5,10,30,100]:
  z,c,w0=m.spectrum(L);z=np.array(z,dtype=float);w0=float(w0)
  # Exact barycentric identity avoids subtracting large signed propagators.
  f=lambda u:-1/(w0*np.prod(u+z))
  I,e=quad(lambda u:u*u*f(u),0,np.inf,epsabs=1e-9,epsrel=1e-11,limit=300)
  J,ej=quad(lambda u:u**3*f(u)*np.sum(1/(u+z)),0,np.inf,epsabs=1e-9,epsrel=1e-11,limit=300)
  p_numeric=float(m.C)*I/4
  p_exact=float(-m.C*np.dot(c,m.spectrum(L)[0]**2*np.log(m.spectrum(L)[0]))/4)
  # In an orthonormal symmetric basis: vv^T trace form, identity Frobenius form.
  v=np.r_[np.ones(4),np.zeros(6)];P=np.eye(10)-np.outer(v,v)/2
  pi_numeric=2*float(m.C)*((np.outer(v,v)/24+np.eye(10)/12)*J-np.outer(v,v)*I/4)
  seagull=p_exact*P;cosmological=-p_exact*P/2;stationary=seagull+cosmological
  pi_exact=2*p_exact*P
  rows.append(dict(L=L,p_radial=p_numeric,p_log=p_exact,radial_relative_error=abs(p_numeric/p_exact-1),integration_by_parts_relative_error=abs(J/(3*I)-1),zero_momentum_stress_relative_error=float(np.max(abs(pi_numeric-pi_exact))/p_exact),scalar_tadpole_coefficient=-p_exact/2,cosmological_tadpole_coefficient=p_exact/2,seagull_P_coefficient=p_exact,cosmological_P_coefficient=-p_exact/2,stationary_P_coefficient=p_exact/2,connected_stress_P_coefficient=2*p_exact,zero_momentum_total_absolute=float(np.max(abs(stationary-pi_numeric/4))),radial_quad_error=e,radial_derivative_quad_error=ej))
 result=dict(scope='Free Euclidean minimally coupled PV scalar determinant, D=4, per physical species. Signed PV sums are not positive covariances. No constrained/composite measure is assumed.',versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),cases=rows,identities=dict(first_variation='S_h=-h:T/2',effective_hessian='Gamma_hh=<S_hh>-Pi/4+lambda volume_hh',p='-C sum c z^2 log(z)/4',stationarity='lambda=p',fixed_field_seagull='p P',cosmological='-p P/2',connected_stress_at_zero='Pi(0)=2p P',full_zero_momentum_kernel='p P/2-Pi(0)/4=0'))
 (Path(__file__).resolve().parent/'results.json').write_text(json.dumps(result,indent=2))
 assert max(r['radial_relative_error'] for r in rows)<1e-10
 assert max(r['integration_by_parts_relative_error'] for r in rows)<1e-10
 assert max(r['zero_momentum_stress_relative_error'] for r in rows)<1e-9
 print(json.dumps(dict(cases=len(rows),max_radial_relative=max(r['radial_relative_error'] for r in rows),max_stress_relative=max(r['zero_momentum_stress_relative_error'] for r in rows)),indent=2))
if __name__=='__main__':main()
