"""Checks for the October 30 cross-reference. No physical unification claimed.

Run with Python + NumPy/SciPy. Existing research files are read only.
"""
from pathlib import Path
import json, math, sys, os
import numpy as np
from scipy.optimize import minimize_scalar

OUT = Path(__file__).resolve().parent
ROOT = Path(os.environ.get('PHYSICS_RESEARCH_ROOT', r'C:\workspace\physics'))
if not (ROOT/'scripts'/'hopf_stability_experiment.py').exists():
    ROOT = OUT/'reference_inputs'

def hopf(x):
    a,b,c,d = np.moveaxis(x, -1, 0)
    return np.stack([2*(a*c+b*d), 2*(b*c-a*d), a*a+b*b-c*c-d*d],axis=-1)

def dhopf(x):
    a,b,c,d=x
    return np.array([[2*c,2*d,2*a,2*b],[-2*d,2*c,2*b,-2*a],
                     [2*a,2*b,-2*c,-2*d]])

def main():
    rng=np.random.default_rng(20261030)
    geometric=[]
    for _ in range(32):
        x=rng.normal(size=4);x/=np.linalg.norm(x)
        # Independent tangent-frame construction by SVD.
        frame=np.linalg.svd(x[None,:],full_matrices=True)[2][1:].T
        n=hopf(x);D=dhopf(x)@frame
        eps=1e-6
        numeric=np.column_stack([(hopf((x+eps*e)/np.linalg.norm(x+eps*e))-
                                  hopf((x-eps*e)/np.linalg.norm(x-eps*e)))/(2*eps)
                                 for e in frame.T])
        F=np.array([[n@np.cross(D[:,i],D[:,j]) for j in range(3)] for i in range(3)])
        error=float(np.max(abs(D-numeric)))
        assert error<1e-8
        assert abs(np.sum(D*D)-8)<1e-12
        assert abs(np.sum(F*F)-32)<1e-11
        assert np.max(abs(np.linalg.svd(D,compute_uv=False)-[2,2,0]))<1e-12
        geometric.append(error)

    # Independent finite-element evaluation of the undilated exact continuum map.
    sys.path.insert(0,str(ROOT/'scripts'))
    import hopf_stability_experiment as fem
    mesh_controls=[]
    for level in range(3):
        saved=np.load(ROOT/'experiments'/'hopf_first_run'/f'level{level}_hopf.npz')
        x=saved['vertices'];t,edges,vol,basis=fem.mesh(x)
        E,_,E2,E4=fem.energy(hopf(x),t,vol,basis)
        mesh_controls.append(dict(level=level,vertices=len(x),energy=E,
                                  E2=E2,E4=E4,relative_to_continuum=E/(24*math.pi**2),
                                  chordal_volume=float(sum(vol))))

    # Frozen-amplitude energy matches the published FS_rho normalization.
    normalizations=[]
    for v,R,kappa in [(1,1,1),(2,.2,1),(.3,4,2),(1,2,1)]:
        rho=v*R/math.sqrt(kappa)
        own=8*math.pi**2*v*v*R+16*math.pi**2*kappa/R
        theorem=(v*v*R/2)*16*math.pi**2*(1+2/rho**2)
        assert abs(own-theorem)<1e-10
        normalizations.append(dict(v=v,R=R,kappa=kappa,theorem_parameter=rho,
                                    within_explicit_global_minimum_range=rho<=1,
                                    energy=own))

    # The theorem cannot simply be imported into the variable-amplitude theory.
    amplitudes=[]
    for lam in [1,4,8,16,64,256]:
        v=R=kappa=1.
        volume=2*math.pi**2*R**3
        def E(r):
            return volume*(4*r*r/R**2+lam/4*(r*r-v*v)**2+8*kappa/R**4)
        r=math.sqrt(max(0,v*v-8/(lam*R*R)))
        sol=minimize_scalar(E,bounds=(0,2),method='bounded',options={'xatol':1e-12})
        assert abs(sol.fun-E(r))<1e-8
        assert E(r)<E(v)
        amplitudes.append(dict(lam=lam,r_stationary=r,energy=E(r),
                               frozen_energy=E(v),numerical_minimizer=float(sol.x),
                               scope='Amplitude minimization with n fixed to the round-sphere Hopf map only.'))

    # Metric (Levi-Civita) Higgs-curvature model; v, MP are reference inputs.
    v=246.;MP=2.435e18;xi_all=(MP/v)**2
    higgs=[]
    for xi in [0,1,1e12,1e15,2.6e15,xi_all]:
        f=xi*v*v/MP**2
        K=1+6*xi*xi*v*v/MP**2
        normalization=1/math.sqrt(K)
        # Exact tree-level massive-fermion/vector linear coupling after rescaling
        # m_E(r) proportional to r/sqrt(F(r)), F=M0^2+xi*r^2.
        coupling=(1-f)*normalization
        higgs.append(dict(xi=xi,fraction_Einstein_coefficient=f,kinetic_K=K,
                          radial_normalization=normalization,
                          relative_tree_massive_matter_coupling=coupling))
    assert xi_all>1e31
    assert higgs[-1]['radial_normalization']<1e-15
    assert abs(higgs[-1]['relative_tree_massive_matter_coupling'])<1e-30

    # Independent finite differences in dimensionless r/v, with finite f controls.
    coupling_controls=[]
    ratio=v/MP
    for f in [0,.1,.5,.9,1]:
        xi=f/ratio**2
        def m(s):return s/math.sqrt(1-f+f*s*s)
        eps=1e-5
        derivative=(m(1+eps)-m(1-eps))/(2*eps)
        assert abs(derivative-(1-f))<1e-8
        coupling_controls.append(dict(f=f,numeric_derivative=derivative,analytic=1-f))

    # Independent honeycomb tight-binding expansion at an unstrained Dirac point.
    delta=np.array([[0.,-1.],[math.sqrt(3)/2,.5],[-math.sqrt(3)/2,.5]])
    Kpoint=np.array([4*math.pi/(3*math.sqrt(3)),0.])
    def band(k,t=1):return abs(t*np.exp(1j*(delta@k)).sum())
    graphene=[]
    assert band(Kpoint)<1e-14
    for t in [.5,1,2]:
        for angle in [0,.3,1.1,2.3]:
            direction=np.array([math.cos(angle),math.sin(angle)])
            eps=1e-6
            speed=(band(Kpoint+eps*direction,t)+band(Kpoint-eps*direction,t))/(2*eps)
            expected=1.5*t
            assert abs(speed-expected)<1e-8
            graphene.append(dict(t=t,angle=angle,speed=speed,expected=expected))

    result=dict(scope='Model identity checks, not experiments or a unified theory.',
        hopf_derivative_cases=len(geometric),max_derivative_error=max(geometric),
        hopf_invariants=dict(gradient_squared=8,F_ij_squared=32,
                            continuum_reference_energy=24*math.pi**2),
        independent_mesh_controls=mesh_controls,theorem_normalizations=normalizations,
        amplitude_controls=amplitudes,higgs_curvature_controls=higgs,
        higgs_exact_coupling_derivative_controls=coupling_controls,
        Higgs_mass_GeV=125,Higgs_Compton_range_m=.1973269804e-15/125,
        graphene_Dirac_controls=graphene,
        exclusions=['No transfer of compact fixed-amplitude stability to the full field model.',
                    'Higgs metric-coupling action already assumes spacetime curvature.',
                    'Graphene cone is not a 3+1-dimensional universal spin-2 mode.'])
    (OUT/'gravity_cross_reference_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
