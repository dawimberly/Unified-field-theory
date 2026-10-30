"""Independent local threshold/sheet verification and explicit kernel matching.
Run: python experiments/threshold_independent_2026_10_30/check.py
Uses a uniformizer and Decimal threshold calculation, not the provisional script.
"""
from pathlib import Path
from decimal import Decimal,localcontext
import importlib.util,json,platform
import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('base',ROOT/'experiments/inverse_scalar_dispersion_2026_10_30/check.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
C=base.C;PI=np.arccos(np.longdouble(-1));LD=np.longdouble;CD=np.clongdouble

def coeff(L):
    z,c=base.spectrum(L)
    return -C*np.dot(c,z*np.log(z))/8,C*np.dot(c,np.log(z))/16

def uniform(u,L,matched=False):
    u=CD(u);s=4/(1+u*u);q=-s;z,c=base.spectrum(L)
    logs=np.empty(4,dtype=CD)
    # Analytic near u=0, on both local sheets. Physical subthreshold u>0.
    logs[0]=-2+PI*u-2*u*np.arctan(u)
    rr=np.sqrt(q/(q+4*z[1:]));logs[1:]=np.log(z[1:])-2+2*np.arctanh(rr)/rr
    U=z+s/2;A=C/4*np.sum(c*(U*U*logs-U*z*np.log(z)))
    if matched:
        r=coeff(10);v=coeff(L);A+=(r[0]-v[0])*q+(r[1]-v[1])*q*q
    return A

def atan_decimal(x):
    term=x;total=x;n=1
    while True:
        term*=-x*x;new=total+term/Decimal(2*n+1)
        if new==total:return new
        total=new;n+=1

def threshold_decimal(L):
    z=[Decimal(1),L,2*L,5*L+1]
    w=[1/np_product([z[i]-z[j] for j in range(4) if j!=i]) for i in range(4)]
    c=[v/w[0] for v in w];total=Decimal(0)
    for m,ci in zip(z,c):
        if m==1:J=Decimal(-2)
        else:
            a=(m-1).sqrt();J=m.ln()-2+2*a*atan_decimal(1/a)
        U=m+2;total+=ci*(U*U*J-U*m*m.ln())
    # Positive C/4 is irrelevant to the threshold zero.
    return total

def np_product(xs):
    out=Decimal(1)
    for x in xs:out*=x
    return out

def critical():
    with localcontext() as ctx:
        ctx.prec=70;lo=Decimal(14);hi=Decimal(15)
        assert threshold_decimal(lo)<0<threshold_decimal(hi)
        for _ in range(190):
            mid=(lo+hi)/2
            if threshold_decimal(mid)>0:hi=mid
            else:lo=mid
        root=(lo+hi)/2
        return str(root),str(threshold_decimal(root))

def pole(L,matched=False):
    f=lambda u:float(uniform(u,L,matched).real)
    # A local root in the uniformizing variable. No global uniqueness assertion.
    u=brentq(f,-.7,.7,xtol=1e-14);s=4/(1+u*u);gap=4*u*u/(1+u*u)
    derivatives=[float(uniform(u+1j*h,L,matched).imag/h) for h in [1e-5,1e-6,1e-7]]
    aq=derivatives[-1]/(8*u/(1+u*u)**2)
    residue=(2+s)**2/(16*aq)
    k=9*C*PI/8;at=float(uniform(0,L,matched).real)
    return dict(L=L,matched_to_L10=matched,u=u,s=s,gap=gap,sheet='I' if u>0 else 'II',threshold_A=at,gap_ratio=gap/float((at/k)**2),source_residue=residue,residue_ratio=residue/float(4.5*np.sqrt(gap)/k),Au_relative_change=abs(derivatives[-1]/derivatives[-2]-1))

def principal_parameter(q,L):
    z,c=base.spectrum(L);total=0j;err=0.
    for m,ci in zip(z,c):
        points=[]
        if q.real<-4*m:
            beta=np.sqrt(1+4*float(m)/q.real);roots=[(1-beta)/2,(1+beta)/2]
            # Resolve the narrow complex-log layer, which generic adaptivity
            # can miss while returning an optimistic error estimate.
            width=abs(q.imag)*float(m)/(q.real*q.real*beta)
            points=sorted(set(p for root in roots for k in [-1000,-100,-10,-1,0,1,10,100,1000] if 0<(p:=root+k*width)<1))
        def f(x):return np.log(complex(m)+q*x*(1-x))
        re,e1=quad(lambda x:f(x).real,0,1,points=points,epsabs=1e-12,epsrel=1e-11,limit=400)
        im,e2=quad(lambda x:f(x).imag,0,1,points=points,epsabs=1e-12,epsrel=1e-11,limit=400)
        U=complex(m)-q/2;total+=complex(ci)*(U*U*complex(re,im)-U*float(m)*float(np.log(m)))
        err+=abs(complex(ci)*U*U)*(e1+e2)
    return complex(C/4)*total,float(C/4)*err

def main():
    lc_text,res=critical();lc=float(lc_text);near=[];branches=[];cuts=[];matched=[];flips=[]
    for offset in [-.03,-.003,-.0003,-.00003,.00003,.0003,.003,.03]:
        row=pole(lc*(1+offset));row['relative_offset']=offset;near.append(row)
    for L in [3.,lc]:
        for u in [.0001,.01,.1]:
            s=4/(1+u*u);first=complex(uniform(u,L));second=complex(uniform(-u,L))
            expected=first-float(C*PI/2)*(1+s/2)**2*u
            flips.append(dict(L=L,u=u,s=s,sheet_flip_error=abs(second-expected)))
        for q in [-3.9+.01j,-3.9-.01j,-4.1+.01j,-4.1-.01j]:
            u=np.sqrt(CD(-4)/q-1);direct=complex(base.A(q,L));other=complex(uniform(u,L))
            branches.append(dict(L=L,q=[q.real,q.imag],physical_overlap_error=abs(direct-other)))
        for s in [4.01,4.3,6.]:
            beta=np.sqrt(1-4/s);lower=complex(uniform(-1j*beta,L));upper=complex(uniform(1j*beta,L))
            expected=-float(C)*float(PI)/4*(1+s/2)**2*beta
            finite=[]
            for eta in [1e-4,1e-6]:
                q=-s-1j*eta;val,e=principal_parameter(q,L)
                u=np.sqrt(CD(-4)/q-1);exact=complex(uniform(u,L))
                finite.append(dict(eta=eta,principal_quad_error=abs(val-exact),quad_error_estimate=e,boundary_difference=abs(val-lower)))
            # Adjacent sheet in the upper-q lip is u=-i beta and equals lower physical lip.
            cuts.append(dict(L=L,s=s,lower_imag=lower.imag,expected_imag=expected,discontinuity_error=abs((lower-upper)-2j*expected),upper_lower_conjugation_error=abs(upper-lower.conjugate()),finite_eta_checks=finite))
    for L in [3,5,10,lc,30,100]:
        row=pole(L,matched=True);row['unmatched_threshold_A']=float(uniform(0,L).real);matched.append(row)
    ref=coeff(10);maps=[]
    for L in [3,5,10,lc,30,100]:
        v=coeff(L);d1=ref[0]-v[0];d2=ref[1]-v[1]
        maps.append(dict(L=L,reference_L=10,delta_A_linear=float(d1),delta_A_quadratic=float(d2),delta_K2_linear=float(-d1/3),delta_K2_quadratic=float(d2/15)))
    result=dict(scope='Local physical/adjacent-sheet threshold audit, free PV determinant only. Reference matching is explicit bookkeeping, not a physical UV completion.',versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),critical_L_decimal=lc_text,decimal_threshold_residual=res,provisional_L=14.557206455684199,critical_L_difference_from_provisional=lc-14.557206455684199,near_threshold=near,physical_overlap=branches,physical_cut=cuts,below_threshold_sheet_flips=flips,matched_kernel_map=maps,matched_comparison=matched)
    (Path(__file__).resolve().parent/'results.json').write_text(json.dumps(result,indent=2))
    assert abs(lc-14.557206455684199)<1e-10
    assert all(r['sheet']==('I' if r['relative_offset']<0 else 'II') for r in near)
    assert max(r['physical_overlap_error'] for r in branches)<1e-10
    assert max(r['discontinuity_error'] for r in cuts)<1e-10
    assert max(r['sheet_flip_error'] for r in flips)<1e-10
    assert max(x['principal_quad_error'] for r in cuts for x in r['finite_eta_checks'])<1e-8
    assert all(r['sheet']=='I' for r in matched)
    assert max(abs(r['gap_ratio']-1) for r in near if abs(r['relative_offset'])<1e-4)<1e-3
    print(json.dumps(dict(critical_L=lc_text,provisional_difference=lc-14.557206455684199,max_overlap_error=max(r['physical_overlap_error'] for r in branches),max_cut_error=max(r['discontinuity_error'] for r in cuts),max_quad_error=max(x['principal_quad_error'] for r in cuts for x in r['finite_eta_checks']),near_threshold=near,matched=matched),indent=2))

if __name__=='__main__':main()
