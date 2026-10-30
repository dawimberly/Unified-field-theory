"""Exact signs and local existence completing the matched spin-2 global count.
Run python experiments/spin2_global_count_2026_10_30/check.py.
Global exclusion uses the stated positive Stieltjes measure, not a numerical scan.
"""
from pathlib import Path
from fractions import Fraction as F
import importlib.util,json,platform
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent

def load(path,name):
 s=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
v=load('experiments/matched_inverse_reconstruction_2026_10_30/local_pole_certificate.py','intervals')
uv=load('experiments/uv_asymptotic_audit_2026_10_30/check.py','uv')
Z=v.Z

def coeff(L):
 z,c=v.m.spec(L);k1=k2=Z
 for mass,ci in zip(z,c):
  ln=v.log_positive(mass);k1=v.add(k1,v.scale(ln,ci*mass/24));k2=v.add(k2,v.scale(ln,ci/240))
 return k1,k2
REF=coeff(10)
def shifts(L):
 k1,k2=coeff(L);return v.sub(REF[0],k1),v.sub(REF[1],k2)
def ser(x):return dict(exact=[str(a) for a in x],decimal=[float(a) for a in x])

def point(L,center):
 """Enclose K2/C and its q derivative on upper half-plane first sheet."""
 q=v.complex_iv(*center);one=v.complex_iv(F(1));out=(Z,Z);der=(Z,Z);z,c=v.m.spec(L)
 for mass,ci in zip(z,c):
  beta=v.csqrt(v.ca(one,v.cscale(v.ci(q),4*mass)))
  lnratio=v.clog(v.cd(v.ca(beta,one),v.cs(beta,one)))
  I=v.ca(v.cm(beta,lnratio),(v.add(v.log_positive(mass),(F(-2),F(-2))),Z))
  Ip=v.cs(v.ci(q),v.cscale(v.cd(lnratio,v.cm(v.cm(q,q),beta)),2*mass))
  U=v.ca(v.complex_iv(mass),v.cscale(q,F(1,4)));lnz=(v.log_positive(mass),Z)
  subterm=v.cm(v.cs(v.complex_iv(mass*mass),v.cscale(q,mass/8)),lnz)
  out=v.ca(out,v.cscale(v.cs(v.cm(v.cm(U,U),I),subterm),ci/15))
  deriv=v.ca(v.ca(v.cscale(v.cm(U,I),F(1,2)),v.cm(v.cm(U,U),Ip)),v.cscale(lnz,mass/8))
  der=v.ca(der,v.cscale(deriv,ci/15))
 a,b=shifts(L)
 out=v.ca(out,v.ca(v.cm((a,Z),q),v.cm((b,Z),v.cm(q,q))))
 der=v.ca(der,v.ca((a,Z),v.cscale(v.cm((b,Z),q),2)))
 return out,der

def realpoint(L,q):
 z,c=v.m.spec(L);out=Z
 for mass,ci in zip(z,c):
  beta=v.sqrt_iv(v.add((F(1),F(1)),v.scale(v.reciprocal((q,q)),4*mass)))
  ratio=v.mul(v.add(beta,(F(1),F(1))),v.reciprocal(v.sub(beta,(F(1),F(1)))))
  I=v.add(v.add(v.log_positive(mass),(F(-2),F(-2))),v.mul(beta,v.log_iv(ratio)))
  U=mass+q/4
  term=v.sub(v.scale(I,U*U),v.scale(v.log_positive(mass),mass*mass-mass*q/8))
  out=v.add(out,v.scale(term,ci/15))
 a,b=shifts(L);return v.add(out,v.add(v.scale(a,q),v.scale(b,q*q)))

def complex_certificate(L,seed):
 center=tuple(F(str(x)) for x in seed);trace=[]
 for _ in range(4):
  value,der=point(L,center);a,b=[(x[0]+x[1])/2 for x in value];c,d=[(x[0]+x[1])/2 for x in der];den=c*c+d*d
  trace.append(dict(center=[str(x) for x in center],residual_L1_upper=str(v.abs_upper(value[0])+v.abs_upper(value[1]))))
  re=center[0]-(a*c+b*d)/den;im=center[1]-(b*c-a*d)/den
  center=tuple(F((x*v.GRID).__floor__(),v.GRID) for x in [re,im])
 r=F(1,10**9);value,der=point(L,center)
 Q=abs(center[0])+abs(center[1])+r;y=center[1]-r;assert y>0
 z,c=v.m.spec(L);M2=F(0)
 for mass,ci in zip(z,c):
  lower=mass*y/Q;upper=mass+Q/4
  logbound=max(v.abs_upper(v.log_positive(lower)),v.abs_upper(v.log_positive(upper)))+4
  U=mass+Q/4
  # (U^2 I)''=I/8+U I'+U^2 I'', |I'|<=1/y, |I''|<=1/y^2.
  M2+=abs(ci)*(logbound/8+U/y+U*U/(y*y))/15
 M2+=2*v.abs_upper(shifts(L)[1])
 error=v.abs_upper(value[0])+v.abs_upper(value[1])+M2*r*r/2
 linear=max(v.distance(der[0]),v.distance(der[1]))*r;assert error<linear
 return dict(L=L,center_exact=[str(x) for x in center],center_decimal=[float(x) for x in center],radius=str(r),scaled_residual=[ser(x) for x in value],scaled_derivative=[ser(x) for x in der],second_derivative_bound=str(M2),Rouche_error_upper=str(error),Rouche_linear_lower=str(linear),strict_margin=str(linear-error),strict_margin_decimal=float(linear-error),midpoint_Newton_locator_trace=trace,conclusion='Exactly one simple upper-half-plane zero in this disk; Schwarz reflection gives its conjugate. Certificate uses rational inequalities, not the midpoint iteration.')

def main():
 coefficients=[];real_bounds=[];real_roots=[];complex_roots=[];exterior=[];counts=[]
 prior=json.loads((ROOT/'experiments/matched_spin2_2026_10_30/results.json').read_text())
 for L in [3,5,10,30,100]:
  k1,k2=coeff(L);a,b=shifts(L);assert k1[0]>0 and k2[1]<0 and REF[0][0]>0 and REF[1][1]<0
  coefficients.append(dict(L=L,k1_over_C=ser(k1),k2_over_C=ser(k2),a_over_C=ser(a),b_over_C=ser(b)))
  if L<10:
   assert b[1]<0
   seed=next(x['q'][0] for x in prior['roots'] if x['L']==L);center=F(str(seed));r=F(1,10**8)
   lo=realpoint(L,center-r);hi=realpoint(L,center+r);assert lo[0]>0 and hi[1]<0
   real_roots.append(dict(L=L,bracket_exact=[str(center-r),str(center+r)],bracket_decimal=[float(center-r),float(center+r)],K_over_C_left=ser(lo),K_over_C_right=ser(hi),conclusion='Unique positive real root by strict monotonicity of K_matched/q, isolated by exact endpoint signs.'))
   number=1
  elif L>10:
   assert b[0]>0
   aa=v.neg(v.mul(b,k2));bb=v.sub(v.mul(b,k1),v.mul(a,k2));cc=v.mul(k1,REF[0])
   discriminant=v.sub(v.square(bb),v.scale(v.mul(aa,cc),4))
   assert aa[0]>0 and cc[0]>0 and discriminant[1]<0
   real_bounds.append(dict(L=L,quadratic=ser(aa),linear=ser(bb),constant=ser(cc),discriminant=ser(discriminant),conclusion='Cauchy-Schwarz lower numerator strictly positive on all real q; F(q)>0 for q>=0. Convexity and matched slope exclude -4<=q<=0.'))
   seed=next(x['q'] for x in prior['roots'] if x['L']==L);complex_roots.append(complex_certificate(L,seed));number=2
  else:number=0
  # Independent uniform UV dominance using prior exact asymptotic bounds.
  R=F(256*(5*L+1));ext=uv.exterior(L,'K2',R)
  if L==10:margin=ext['raw_no_zero_margin']
  else:margin=v.distance(b)-v.abs_upper(a)/R-(uv.absolute(ext['constant'])+ext['raw_difference_bound'])/R**2
  assert margin>0
  exterior.append(dict(L=L,R=str(R),positive_margin_summary=ser(v.iv(margin)),comparison='constant for unchanged L10; quadratic polynomial dominance for nonreference matched kernels',conclusion='K2_matched nonzero for all complex first-sheet |q|>=R.'))
  counts.append(dict(L=L,additional_first_sheet_zeros=number,massless_zero_multiplicity=1,cut_boundary_zeros=0,basis='Analytic Stieltjes exclusion plus exact real signs/local existence; no finite scan supplies global completeness.'))
 result=dict(scope='Five explicitly Lref10 matched auxiliary spin-2 kernels. Not a constrained-CEV equivalence, physical Lorentzian contour or UV-completion claim.',baseline='649a2d4b287407342ab5fb9e7fd2b4c8a4591d3c',python=platform.python_version(),interval_grid=str(v.GRID),coefficients=coefficients,real_axis_Cauchy_Schwarz_certificates=real_bounds,positive_real_root_brackets=real_roots,complex_local_certificates=complex_roots,uniform_exterior_certificates=exterior,global_counts=counts)
 (HERE/'results.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(dict(global_additional_counts=[x['additional_first_sheet_zeros'] for x in counts],exact_real_root_brackets=len(real_roots),exact_complex_disks=len(complex_roots),exact_real_exclusion_bounds=len(real_bounds),uniform_exterior_bounds=len(exterior),minimum_Rouche_margin=min(x['strict_margin_decimal'] for x in complex_roots)),indent=2))
if __name__=='__main__':main()
