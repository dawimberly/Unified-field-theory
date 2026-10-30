"""Scalar divided-difference identity and formal low-energy matching audit.
Run: python experiments/scalar_identity_matching_2026_10_30/check.py
No physical counterterm choice or constrained-theory equivalence is made.
"""
from pathlib import Path
import importlib.util,json,math,platform
from fractions import Fraction
import numpy as np
import scipy
from scipy.special import roots_legendre
ROOT=Path(__file__).resolve().parents[2]
def load(relative,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/relative)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
scalar=load('experiments/inverse_scalar_dispersion_2026_10_30/check.py','scalar_reference')
spin=load('experiments/spin2_stieltjes_2026_10_30/check.py','spin_reference')
C=scalar.C

def scalar_simplex(q,L,n):
    z,_,w0=spin.spectrum(L);x,w=roots_legendre(n);x=(x+1)/2;w=w/2
    a,b,d=np.meshgrid(x,x,x,indexing='ij');wa,wb,wd=np.meshgrid(w,w,w,indexing='ij')
    Z=((1-a)*z[0]+a*(1-b)*z[1]+a*b*(1-d)*z[2]+a*b*d*z[3]).ravel()
    weight=(wa*wb*wd*a*a*b).ravel();total=0j
    for xx,ww in zip(x,w):
        t=xx*(1-xx);y=Z+t*q;delta=q*(t+.5)
        third=2/y+2*delta/y**2+2*delta**2/y**3-2/Z-q/(2*Z**2)
        total+=ww*np.sum(weight*third)
    return C/(4*w0)*total

def moment(n):return Fraction(math.factorial(n)**2,math.factorial(2*n+1))
def rational_coeff(n,channel):
    if channel=='A':return Fraction((-1)**(n+1),4)*(moment(n)/n+moment(n-1)/(n-1)+moment(n-2)/(4*(n-2)))
    return Fraction((-1)**(n-3),4*n*(n-1)*(n-2))*moment(n)
def coefficients(L,channel):
    z,c=scalar.spectrum(L);G=C*np.dot(c,z*np.log(z));S=np.dot(c,np.log(z))
    vals={1:-G/8 if channel=='A' else G/24,2:C*S/(16 if channel=='A' else 240)}
    for n in range(3,9):
        r=rational_coeff(n,channel)
        vals[n]=C*np.longdouble(r.numerator)/r.denominator*np.dot(c,z**(2-n))
    return vals

def main():
    checks=[];matching=[];series=[];reference=10
    for L in [3,10,100]:
        for q in [1+2j,-2+.3j]:
            nodes=[24,48,96] if L==100 else [24,48]
            exact=complex(scalar.A(q,L));vals=[complex(scalar_simplex(q,L,n)) for n in nodes]
            checks.append(dict(L=L,q=[q.real,q.imag],exact=[exact.real,exact.imag],nodes=nodes,relative_errors=[abs(v-exact)/abs(exact) for v in vals]))
    for L in [3,5,10,30,100]:
        for channel in ['A','K2']:
            co=coefficients(L,channel);ref=coefficients(reference,channel)
            kernel=scalar.A if channel=='A' else spin.kernel
            for q in [.05,.1,.2]:
                counter=(ref[1]-co[1])*q+(ref[2]-co[2])*q*q
                exact=float((kernel(q,L)+counter-kernel(q,reference)).real)
                predicted=sum((co[n]-ref[n])*q**n for n in range(3,9))
                matching.append(dict(L=L,channel=channel,q=q,reference_L=reference,formal_counterterm=float(counter),matched_difference=exact,series_difference=float(predicted),absolute_error=abs(exact-float(predicted)),difference_over_q3=exact/q**3,cubic_coefficient_difference=float(co[3]-ref[3])))
            series.append(dict(L=L,channel=channel,coefficients={str(n):float(v) for n,v in co.items()}))
    rational={ch:{str(n):str(rational_coeff(n,ch)) for n in range(3,9)} for ch in ['A','K2']}
    result=dict(scope='Algebraic scalar identity and bookkeeping comparison through q^2. No global scalar pole theorem or adopted bare action.',versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),scalar_simplex_checks=checks,matched_comparisons=matching,low_energy_series=series,rational_factors=rational)
    (Path(__file__).resolve().parent/'results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(max_scalar_simplex_relative=max(x['relative_errors'][-1] for x in checks),max_matching_absolute_error=max(x['absolute_error'] for x in matching),comparisons=len(matching),rational_factors=rational),indent=2))
    assert rational_coeff(3,'A')==Fraction(17,1120)
    assert rational_coeff(3,'K2')==Fraction(1,3360)
    assert max(x['relative_errors'][-1] for x in checks)<1e-5
    assert max(x['absolute_error'] for x in matching)<1e-10

if __name__=='__main__':main()
