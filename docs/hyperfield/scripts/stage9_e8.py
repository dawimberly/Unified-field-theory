"""Stage 9: does the E8-lattice torus R^8/Gamma_8 give the Standard Model on the brane?"""
import json, os, itertools, numpy as np
from math import pi, sqrt
CK = 'results/stage9.json'
out = json.load(open(CK)) if os.path.exists(CK) else {}
def save(): json.dump(out, open(CK, 'w'), indent=1, default=float)
HBARC = 1.973269804e-16; MBAR = 2.435e18

# E8 roots (even coordinate system: D8 roots + half-integer vectors with even # of minus signs)
I = np.eye(8)
R = [s1*I[i]+s2*I[j] for i in range(8) for j in range(i+1, 8) for s1 in (1, -1) for s2 in (1, -1)]
for k in range(256):
    s = np.array([1 if (k >> b) & 1 else -1 for b in range(8)])
    if np.sum(s < 0) % 2 == 0: R.append(s/2)
R = np.array(R); assert len(R) == 240

# ---- 1. KK reduction: isometries of flat torus; point group; root-state masses
def R_of(n, MD): return (MBAR**2/(MD*1e3)**(n+2))**(1/n)*HBARC
kk = {}
for MD in (6.0, 10.0):
    Rm = R_of(8, MD)                       # m;  equal volume (2 pi R)^8, E8 det 1, self-dual
    M1 = sqrt(2)*HBARC/Rm                  # GeV: m = 2 pi |k|/L, L = 2 pi R, |k|^2 = 2
    kk[str(MD)] = dict(R_m=Rm, M_root_KK_GeV=M1, lambda_m=Rm/sqrt(2))
out['kk_reduction'] = dict(
    isometry=('DERIVED: Killing vectors of a flat metric on R^8/Gamma are the 8 constant translations -> continuous isometry U(1)^8; '
              'KK reduction of pure gravity gives 8 abelian graviphotons g_{mu i} (+ 36 metric moduli). Rotations are not isometries; '
              'only the finite point group Aut(E8)=W(E8), order 696,729,600, survives as a DISCRETE global symmetry. No non-abelian gauge group.'),
    root_states_point_particle=('the 240 states with momentum k = root (|k|^2=2) are the first KK shell: massive spin-2 (graviton KK) states carrying U(1)^8 charges; '
                                'mass M = sqrt2/R at equal volume, never massless'),
    MD_rows=kk,
    string_enhancement=('closed string on T^8 with lattice Gamma_8 at the self-dual point (G_ij = C_ij/2, B_ij = antisymmetrised Cartan matrix, alpha\'=1; '
                        'Elitzur-Gross-Rabinovici-Seiberg 1987, from memory): winding+momentum states with p_R=0, p_L^2=2 are massless -> E8 gauge bosons. '
                        'Near the point, along an overall scale t: M_root ~ (sqrt2/sqrt(alpha\')) |t - 1/t| (circle-analog, ASSUMED form). '
                        'Requires R ~ sqrt(alpha\') ~ 1/M_string, i.e. R^-1 ~ M_s ~ TeV-Planck, NOT the ADD radius R^-1 ~ GeV.'),
    heterotic=('heterotic E8xE8: the E8xE8 roots live in the 16-dimensional INTERNAL left-moving lattice Gamma_8+Gamma_8, massless at any radius; '
               'that lattice is not a geometric spacetime torus. Dan\'s geometric E8 torus is a different object.'))
save()

# ---- 2. branchings, verified on the E8 root system
def subsystem_centralizer(sub):
    orth = [r for r in R if np.all(np.abs(sub @ r) < 1e-9)]
    return np.array(orth)
def find_chain(k):
    """simple roots of A_k chain: consecutive inner product -1, others 0"""
    chain = [R[0]]
    while len(chain) < k:
        for r in R:
            if abs(r @ chain[-1] + 1) < 1e-9 and all(abs(r @ c) < 1e-9 for c in chain[:-1]):
                chain.append(r); break
    return np.array(chain)
def roots_of_span(simple):
    P = simple.T @ np.linalg.pinv(simple.T)          # projector onto span
    return np.array([r for r in R if np.allclose(P @ r, r)])
br = {}
for name, k, cent_name, expect in [('SU(3)', 2, 'E6', 72), ('SU(5)', 4, 'SU(5)', 20), ('SU(4)', 3, 'SO(10)', 40)]:
    S = find_chain(k); sub_roots = roots_of_span(S); C = subsystem_centralizer(sub_roots)
    rest = [r for r in R if not any(np.allclose(r, x) for x in sub_roots) and not any(np.allclose(r, x) for x in C)]
    # project rest onto the subalgebra plane -> its weights; count multiplicities
    Pm = S.T @ np.linalg.pinv(S.T)
    w = np.round(np.array([Pm @ r for r in rest]), 6)
    u, c = np.unique(w, axis=0, return_counts=True)
    br[f'E8 > {cent_name} x {name}'] = dict(sub_roots=len(sub_roots), centralizer_roots=len(C), expected_centralizer_roots=expect,
                                           remaining_roots=len(rest), distinct_subalgebra_weights=len(u), multiplicities=sorted(set(c.tolist())))
br['dimension_checks'] = {
    'E8 > E6 x SU(3)': '248 = (78,1)+(1,8)+(27,3)+(27b,3b) = 78+8+81+81 = %d' % (78+8+81+81),
    'E8 > SU(5) x SU(5)': '248 = (24,1)+(1,24)+(10,5)+(10b,5b)+(5,10b)+(5b,10) = %d' % (24+24+50*4),
    'E8 > SO(10) x SU(4)': '248 = (45,1)+(1,15)+(10,6)+(16,4)+(16b,4b) = %d' % (45+15+60+64+64),
    'E8 > SO(16)': '248 = 120 + 128 = %d' % (120+128),
    'E6 > SO(10) x U(1)': '27 = 16 + 10 + 1 = %d' % (16+10+1),
    'SO(10) 16 > SM': 'Q(3,2)+u^c(3b,1)+d^c(3b,1)+L(1,2)+e^c(1,1)+nu^c(1,1) = 6+3+3+2+1+1 = %d (one generation + nu_R)' % (6+3+3+2+1+1)}
out['branchings'] = br
out['chirality'] = dict(
    flat_torus='DERIVED (index theorem): Dirac index on flat T^n vanishes -> zero net chiral zero modes; KK reduction on R^8/Gamma_8 gives a NON-CHIRAL 4D spectrum (Witten 1983 no-go for smooth/flat). Lattice shape (E8 vs Z8) is irrelevant.',
    fixes='chirality needs orbifold fixed points, singularities, fluxes or a domain wall (Stage 5 kink, Q=1)',
    heterotic_Z3='heterotic E8xE8 on T^6/Z3, standard embedding: E6xSU(3)xE8, net 36 generations of 27 = 9 untwisted [(27,3): 3 x 3] + 27 twisted (one per fixed point) (Dixon-Harvey-Vafa-Witten 1985; from memory)',
    calabi_yau='standard embedding on CY3: N_gen = |chi|/2; three generations need |chi| = 6 (e.g. Tian-Yau/Z3) (from memory)',
    three_flag=('COINCIDENCE/FEATURE: the "3" of E8 > E6 x SU(3) is the SU(3) holonomy of a 6-dim compact space, not the generation number; '
                'it multiplies untwisted 27s (giving 9 on T^6/Z3), and N_gen is set by topology (chi or fixed points). Not a prediction of 3.'),
    lisi=('Distler & Garibaldi, Commun. Math. Phys. 298, 419 (2010), arXiv:0905.2658 (verified): no subgroup SL(2,C).G of any real or complex form of E8 '
          'gives a chiral gauge theory; 3 generations need >= 180 states in the fermionic (-1) eigenspace, real E8 forms provide only 112 or 128.'),
    dims_note='heterotic uses 6 compact dims (10D string); Dan\'s 8 compact dims (12D total) has no known consistent string/M-theory with an E8-torus; E8 here is a lattice, not a gauge group.')
save()

# ---- 3. couplings and KK spectrum
cp = []
for MD in (6.0, 10.0):
    V8M8 = (2*pi)**8*(MBAR/(MD*1e3))**2       # V8 M_D^8 = (2pi)^8 Mbar^2/M_D^2 (PDG convention Mbar^2 = M_D^{n+2} R^n)
    for gD2 in (1.0, 4*pi):
        g4sq = gD2/V8M8
        cp.append(dict(M_D_TeV=MD, gD2_hat=gD2, V8_MD8=V8M8, g4_squared=g4sq, alpha4=g4sq/(4*pi), alpha_em_over_alpha4=(1/137.036)/(g4sq/(4*pi)),
                       gD2_hat_needed_for_alpha_em=(4*pi/137.036)*V8M8))
out['couplings'] = dict(rows=cp, formula='1/g4^2 = V8/g_D^2, g_D^2 = gD2_hat / M_D^8 (dimensionless bulk coupling ASSUMED O(1)-4pi)',
    verdict='bulk gauge fields give alpha_4 ~ 1e-37: excluded by ~35 orders unless g_D is non-perturbatively huge (gD2_hat ~ 1e33: tuning). '
            'Bulk SM gauge fields would also have KK copies at ~2-4 GeV, excluded (UED-type limits >~ TeV, from memory). => SM gauge fields must be brane-confined; then the bulk lattice does not set the gauge group.')
# theta series
def zn_shells(n, Nmax):
    one = np.zeros(Nmax+1); k = np.arange(0, int(np.sqrt(Nmax))+1); one[k**2] += 1; one[k[1:]**2] += 1
    c = one.copy()
    for _ in range(n-1): c = np.convolve(c, one)[:Nmax+1]
    return c
def s3(k): return sum(d**3 for d in range(1, k+1) if k % d == 0)
e8 = [(2*k, 240*s3(k)) for k in range(1, 7)]
z8c = zn_shells(8, 13); z8 = [(k, int(z8c[k])) for k in range(1, 13)]
# verify E8 shells by brute counting of lattice vectors with norm 2,4
def count_e8(norm):
    cnt = 0
    for x in itertools.product(range(-2, 3), repeat=8):
        v = np.array(x); 
        if v.sum() % 2 == 0 and v @ v == norm: cnt += 1
    for x in itertools.product((-1.5, -0.5, 0.5, 1.5), repeat=8):
        v = np.array(x)
        if (v.sum()) % 2 == 0 and abs(v @ v - norm) < 1e-9: cnt += 1
    return cnt
ver = {2: count_e8(2), 4: count_e8(4)}
out['kk_spectrum'] = dict(E8_shells_norm_count=e8, E8_bruteforce_check=ver, Z8_shells_norm_count=z8,
    E8_masses_over_M1=[sqrt(n/2) for n, c in e8], Z8_masses_over_M1z=[sqrt(n) for n, c in z8],
    ratio_2160_240=2160/240, ratio_note='E8 shell degeneracies 240 sigma3(k): 240, 2160, 6720, 17520 (ratio 9 = sigma3(2)); Z8: 16, 112, 448, 1136',
    M1_E8_GeV={k: v['M_root_KK_GeV'] for k, v in kk.items()},
    weyl_shape_effect_1TeV={k: (v['M_root_KK_GeV']/1000)**2 for k, v in kk.items()},
    observability=('each KK graviton couples as 1/Mbar (2.4e18 GeV): individual shells unresolvable; summed monojet/diphoton rates are shape-independent '
                   'up to (M1/E)^2 ~ 4e-6 at 1 TeV (Weyl law). Gravity deviation at r~lambda=1e-16 m unmeasurable (Stage 7). n=8 collider bound UNVERIFIED (CMS covers n<=7).'))
save()
print(json.dumps(out['branchings'], indent=1)); print(out['couplings']['rows']); print(ver, e8[:3], z8[:4]); print(kk)
