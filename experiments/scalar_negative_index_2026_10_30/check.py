"""Applicability audit of finite-index generalized Nevanlinna methods.
Run: python experiments/scalar_negative_index_2026_10_30/check.py
The infinite-index conclusion uses the continuity/diagonal argument in REPORT.md;
finite matrices are numerical illustrations, not that proof.
"""
from pathlib import Path
from fractions import Fraction
import importlib.util,json,platform
import numpy as np
import scipy
root=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('inverse_scalar',root/'experiments/inverse_scalar_dispersion_2026_10_30/check.py')
scalar=importlib.util.module_from_spec(spec);spec.loader.exec_module(scalar)

def h(z,L):return np.asarray(1/scalar.A(-z,L),dtype=complex)

def exact_cut_sign(L):
    # At s=6L only physical and first PV species are open.
    c1=-Fraction(5*(2*L-1),4*L+1)
    positive_squared=Fraction((1+3*L)**4)*(1-Fraction(2,3*L))
    negative_squared=c1*c1*Fraction((4*L)**4,3)
    margin=negative_squared-positive_squared
    assert c1<0 and margin>0
    return dict(L=L,s=6*L,c1=str(c1),squared_margin=str(margin),rho_inverse=float(scalar.rho(6*L,L)))

def matrix_check(L,n,relative_eta):
    eta=L*relative_eta;x=np.linspace(5.8*L,6.2*L,n)
    z=x+1j*eta;f=h(z,L)
    K=(f[:,None]-f.conj()[None,:])/(z[:,None]-z.conj()[None,:])
    normalization=max(abs(h(6*L+1j*L*1e-8,L).imag),1e-20)
    scaled=eta*K/normalization
    error=float(np.max(abs(scaled-scaled.conj().T)))
    eig=np.linalg.eigvalsh((scaled+scaled.conj().T)/2)
    diagonal_limit=np.array([np.pi*scalar.rho(float(xx),L) for xx in x])/normalization
    off=scaled.copy();np.fill_diagonal(off,0)
    return dict(L=L,n=n,relative_eta=relative_eta,eta=eta,negative_eigenvalues=int(sum(eig<-1e-10)),min_eigenvalue=float(eig[0]),max_eigenvalue=float(eig[-1]),hermitian_error=error,max_scaled_off_diagonal=float(np.max(abs(off))),diagonal_limit_error=float(max(abs(np.diag(scaled).real-diagonal_limit))))

def main():
    signs=[];rows=[]
    for L in [3,5,10,30,100]:
        signs.append(exact_cut_sign(L))
        for n in [2,4,8,16,32]:
            for eta in [1e-2,1e-4,1e-6]:rows.append(matrix_check(L,n,eta))
    result=dict(scope='Fixed auxiliary scalar inverse, h(z)=1/A(-z). Finite-index theorem applicability audit, not a new physical ghost or complex-pole count.',versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),exact_sign_points=signs,gram_matrices=rows)
    assert all(r['rho_inverse']<0 for r in signs)
    assert all(r['negative_eigenvalues']==r['n'] for r in rows if r['relative_eta']==1e-6)
    assert max(r['hermitian_error'] for r in rows)<1e-12
    (Path(__file__).resolve().parent/'results.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(exact_sign_points=signs,matrix_cases=len(rows),resolved_counts=[dict(L=r['L'],n=r['n'],count=r['negative_eigenvalues']) for r in rows if r['relative_eta']==1e-6],max_hermitian_error=max(r['hermitian_error'] for r in rows)),indent=2))

if __name__=='__main__':main()
