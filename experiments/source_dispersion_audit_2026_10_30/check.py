"""Source-projected dispersion and contact/subtraction audit, no new physical model.
Run python experiments/source_dispersion_audit_2026_10_30/check.py.
"""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util,json,platform
import numpy as np
import scipy
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('inverse',ROOT/'experiments/infinite_scalar_dispersion_2026_10_30/check.py');d=importlib.util.module_from_spec(s);s.loader.exec_module(d)
LD=d.LD;CD=d.CD;PI=d.PI;C=d.C
prior=json.loads((ROOT/'experiments/infinite_scalar_dispersion_2026_10_30/results.json').read_text())

def factor(q):return (q-2)**2/16
def W(response,q):return factor(q)/response.A(q)
def Wprime(response,q):return (q-2)/(8*response.A(q))-factor(q)*response.derivative(q)/response.A(q)**2
class Source:
 def __init__(self,response):self.response=response;self.model=response.model;self.L=response.L
 def rho(self,s):return (factor(-s)/self.response.lower(s)).imag/PI

def cert(L,matched):
 old=next(x for x in prior['exact_tail_certificates'] if x['L']==L and x['matched_to_L10']==matched)
 R=old['R'];K=Q(old['rho_K_upper'])*(1+Q(2,R))**2/16
 return dict(L=L,matched_to_L10=matched,R=R,rho_K_upper=str(K),rho_K_decimal=float(K),rho_decay_power=old['rho_decay_power']-2,parent_coefficient_certificate=old)

def ps_list(response):
 row=next(x for x in prior['poles'] if x['L']==response.L and x['matched_to_L10']==response.matched)
 return [dict(q=CD(complex(*p['q'])),residue=CD(complex(*p['residue']))*factor(CD(complex(*p['q'])))) for p in row['poles']]
def qbound(q):return Q(str(float(abs(q.real)+abs(q.imag))))
def cauchy_derivative(response):
 values=[];q0=CD(1);r=LD(1)/4
 for n in [32,64]:
  phase=np.exp(2j*PI*np.arange(n,dtype=LD)/n)
  vals=np.array([W(response,q0+r*x) for x in phase],dtype=CD)
  values.append(np.sum(vals/phase)/(n*r))
 exact=Wprime(response,q0)
 error=d.normerr(exact,values[-1]);assert error<1e-10
 return dict(L=response.L,matched_to_L10=response.matched,analytic=d.pair(exact),cauchy_32=d.pair(values[0]),cauchy_64=d.pair(values[1]),normalized_error=error)

def twice_subtracted(response,ps,b,q):
 q=CD(q);q0=CD(1);h=q-q0;source=Source(response);R=b['R'];K=Q(b['rho_K_upper'])
 # Factored second-subtraction kernel prevents loss from subtracting three UV terms.
 kernel=lambda ss:h*h/((ss+q)*(ss+q0)**2)
 cut,ce=d.cut_integral(source,R,kernel,f'{response.L} {response.matched} source 2sub {q}')
 tail,te=d.tail_integral(source,R,kernel,f'{response.L} {response.matched} source 2sub tail {q}')
 pole=sum(p['residue']*h*h/((q-p['q'])*(q0-p['q'])**2) for p in ps)
 delta=W(response,q)-W(response,q0)-h*Wprime(response,q0)
 finf=1/(C*response.model.constant);polynomial=finf*h*h/16
 bound=K*qbound(h)**2/((1-max(qbound(q),Q(1))/R)**3*R)
 assert abs(tail)<float(bound)
 reconstructed=polynomial+pole+cut+tail
 return dict(L=response.L,matched_to_L10=response.matched,q=d.pair(q),q0=d.pair(q0),direct_second_difference=d.pair(delta),known_quadratic_remainder=d.pair(polynomial),pole_second_difference=d.pair(pole),finite_cut=d.pair(cut),numerical_UV_tail=d.pair(tail),UV_tail_abs_bound=str(bound),UV_tail_abs_bound_decimal=float(bound),quadrature_error_estimate=ce+te,normalized_closure_error=d.normerr(delta,reconstructed),error_omitting_known_quadratic=d.normerr(delta,reconstructed-polynomial))

def matched_absolute(response,ps,b,q):
 q=CD(q);source=Source(response);R=b['R'];K=Q(b['rho_K_upper'])
 kernel=lambda ss:1/(ss+q)
 cut,ce=d.cut_integral(source,R,kernel,f'{response.L} source absolute {q}')
 tail,te=d.tail_integral(source,R,kernel,f'{response.L} source absolute tail {q}')
 pole=sum(p['residue']/(q-p['q']) for p in ps);contact=1/(16*response.model.d2)
 value=W(response,q);recon=contact+pole+cut+tail
 bound=K/(3*(1-qbound(q)/R)*R**3);assert abs(tail)<float(bound)
 return dict(L=response.L,q=d.pair(q),direct_W=d.pair(value),source_contact=float(contact),pole_sum=d.pair(pole),finite_cut=d.pair(cut),numerical_UV_tail=d.pair(tail),UV_tail_abs_bound=str(bound),UV_tail_abs_bound_decimal=float(bound),quadrature_error_estimate=ce+te,normalized_closure_error=d.normerr(value,recon),error_omitting_source_contact=d.normerr(value,recon-contact))

def matched_moments(response,ps,b):
 source=Source(response);R=b['R'];K=Q(b['rho_K_upper']);a1=response.model.d1;a2=response.model.d2;ainf=C*response.model.constant
 targets=[(-4/a2-a1/a2**2)/16,(4/a2+4*a1/a2**2+a1*a1/a2**3-ainf/a2**2)/16];rows=[]
 for n in [0,1]:
  kernel=lambda ss:ss**n
  cut,ce=d.cut_integral(source,R,kernel,f'{response.L} source moment {n}')
  tail,te=d.tail_integral(source,R,kernel,f'{response.L} source moment tail {n}')
  pole=sum(p['residue']*p['q']**n for p in ps);value=pole+(-1)**n*(cut+tail)
  scale=max(LD(1),abs(pole)+abs(cut)+abs(tail)+abs(targets[n]));bound=K/Q((2-n)*R**(2-n));assert abs(tail)<float(bound)
  rows.append(dict(L=response.L,n=n,pole_moment=d.pair(pole),finite_cut_moment=d.pair(cut),numerical_tail_moment=d.pair(tail),expected_coefficient=float(targets[n]),reconstructed=d.pair(value),absolute_error=float(abs(value-targets[n])),component_normalized_error=d.normerr(value,targets[n],scale),UV_tail_bound=str(bound),UV_tail_bound_decimal=float(bound),quadrature_error_estimate=ce+te))
 return rows

def main():
 certificates=[];derivatives=[];raw=[];matched=[];moments=[];source_poles=[]
 for L in [3,5,10,30,100]:
  model=d.m.Model(L)
  for ismatched in [False,True]:
   if ismatched and L==10:continue # Same prescription, not a second independent case.
   response=d.Response(model,ismatched);b=cert(L,ismatched);certificates.append(b);ps=ps_list(response)
   derivatives.append(cauchy_derivative(response))
   source_poles.append(dict(L=L,matched_to_L10=ismatched,source_residues=[dict(q=d.pair(p['q']),residue=d.pair(p['residue'])) for p in ps]))
   if b['rho_decay_power']==3:
    for q in [.5,2.,1+2j]:matched.append(matched_absolute(response,ps,b,q))
    moments+=matched_moments(response,ps,b)
   else:
    for q in [.5,2.,1+2j]:raw.append(twice_subtracted(response,ps,b,q))
 result=dict(scope='Fixed identical minimally coupled constituent source quadratic connected coefficient W=(q-2)^2/(16A). Not an S-matrix, arbitrary-source scalar kernel or constrained CEV correlator.',baseline='5a6a498fb88c02e837e755d58a5ca38ce625bc47',versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,longdouble_bits=np.finfo(LD).nmant),exact_source_tail_certificates=certificates,source_poles=source_poles,new_source_derivative_checks=derivatives,raw_reference_twice_subtracted=raw,nonreference_matched_absolute=matched,matched_source_moments=moments,quadrature_warnings=d.WARNINGS)
 (HERE/'results.json').write_text(json.dumps(result,indent=2))
 summary=dict(raw_checks=len(raw),matched_checks=len(matched),source_moment_checks=len(moments),derivative_checks=len(derivatives),max_raw_closure=max(x['normalized_closure_error'] for x in raw),max_matched_closure=max(x['normalized_closure_error'] for x in matched),max_source_moment_error=max(x['component_normalized_error'] for x in moments),warnings=len(d.WARNINGS))
 print(json.dumps(summary,indent=2))
 assert summary['max_raw_closure']<1e-8 and summary['max_matched_closure']<1e-8 and summary['max_source_moment_error']<1e-8
if __name__=='__main__':main()
