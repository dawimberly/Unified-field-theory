"""Signed Weyl traces, multiplicities, couplings and Martin-formula checks
for the Yukawa-containing three-loop vacuum classes QFf, CcFf, FFff(A), FFff(B).
Python 3 with numpy and sympy."""
import itertools, math
import numpy as np
import sympy as sp

rng = np.random.default_rng(7)

# ---------------------------------------------------------------- Euclidean sigma matrices
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]], complex)
s3 = np.array([[1, 0], [0, -1]], complex); one = np.eye(2, dtype=complex)
sig  = [-1j*s1, -1j*s2, -1j*s3, one]      # sigma_mu
sigb = [ 1j*s1,  1j*s2,  1j*s3, one]      # sigmabar_mu
eps = 1j*s2                               # eps_{ab} = [[0,1],[-1,0]]
sl  = lambda k: sum(k[m]*sig[m]  for m in range(4))
slb = lambda k: sum(k[m]*sigb[m] for m in range(4))
for m in range(4):
    assert np.allclose(sigb[m], eps @ sig[m].T @ np.linalg.inv(eps))        # sigmabar = eps sigma^T eps^-1
    for n in range(4):
        assert np.isclose(np.trace(sig[m] @ sigb[n]), 2*(m == n))            # tr sigma sigmabar = 2 delta
c = 1j                                    # kinetic propagator <psi_a psibar_adot>(k) = c (sigma.k)_{a adot}/(k^2+m^2)
S = lambda k: c*sl(k)

def perm_sign(seq):
    inv = sum(1 for i in range(len(seq)) for j in range(i+1, len(seq)) if seq[i] > seq[j])
    return -1 if inv % 2 else 1

def amp(pi, ks):
    """Full Wick contraction of F^n f^n slot operators: F = psi_a eps_ab psi_b, f = psibar eps psibar.
    pi[i] = index of the psibar slot contracted with psi slot i; ks[i] = momentum of that line (F->f)."""
    n = len(pi)
    order = []
    for i in range(n):
        order += [i, n + pi[i]]                       # reorder to (psi_i psibar_pi(i)) pairs
    sgn = perm_sign(order)
    A = "abcdefgh"[:n]; B = "ijklmnop"[:n]
    terms = [A[2*v:2*v+2] for v in range(n//2)] + [B[2*v:2*v+2] for v in range(n//2)]
    ops = [eps]*len(terms)
    for i in range(n):
        terms.append(A[i] + B[pi[i]]); ops.append(S(ks[i]))
    return sgn*np.einsum(",".join(terms) + "->", *ops)

dot = lambda a, b: float(np.dot(a, b))

# ---------------------------------------------------------------- 1. two-vertex trace (Ff)
k, q = rng.normal(size=4), rng.normal(size=4)
a_id, a_sw = amp((0, 1), (k, q)), amp((1, 0), (q, k))     # swapping slots swaps which momentum sits on which slot
r2 = a_id / (-2*dot(k, q))
assert np.isclose(a_id, a_sw) and np.isclose(abs(r2), 1)
print("two-vertex: amp/(-2 k.q) =", np.round(r2, 12), "  (both slot pairings equal)")
kappa = r2.real                                           # overall sign of one F f pair; enters squared at 4 vertices
N2 = lambda k, q: kappa*(-2*dot(k, q))                    # signed single-loop numerator

# ---------------------------------------------------------------- 2. four-vertex FFff: all 24 fermion matchings
trchain = lambda a, b, cc, d: np.trace(sl(a) @ slb(b) @ sl(cc) @ slb(d))
def levi(a, b, cc, d): return np.linalg.det(np.array([a, b, cc, d]))
ratioA, ratioB = set(), set()
for pi in itertools.permutations(range(4)):
    ks = [rng.normal(size=4) for _ in range(4)]
    val = amp(pi, ks)
    f_of = [p // 2 for p in pi]
    if f_of[0] == f_of[1] and f_of[2] == f_of[3]:                        # two 2-cycles
        pa, pb = (pi[0] % 2, pi[1] % 2), (pi[2] % 2, pi[3] % 2)
        ref = amp(pa, (ks[0], ks[1])) * amp(pb, (ks[2], ks[3]))
        ratioA.add(round((val/ref).real, 9))
    else:                                                                 # one 4-cycle F0 f0 F1 f1
        line = {(i//2, pi[i]//2): ks[i] for i in range(4)}
        ref = c**4 * trchain(line[0, 0], line[1, 0], line[1, 1], line[0, 1])
        ratioB.add(round((val/ref).real, 9))
        # trace = 2(dd-dd+dd) + epsilon term
        a_, b_, c_, d_ = line[0, 0], line[1, 0], line[1, 1], line[0, 1]
        sym = 2*(dot(a_, b_)*dot(c_, d_) - dot(a_, c_)*dot(b_, d_) + dot(a_, d_)*dot(b_, c_))
        assert np.isclose(abs(trchain(a_, b_, c_, d_) - sym), 2*abs(levi(a_, b_, c_, d_)))
print("class A (2x2-cycle): amp / [product of two-vertex amps] =", ratioA)
print("class B (4-cycle)  : amp / [c^4 tr(s sb s sb)]            =", ratioB)
assert ratioA == {1.0} and len(ratioB) == 1
rB = ratioB.pop()
# four-cycle with c^4=1:  amp = rB * tr[...]; momentum conservation sets two opposite lines equal => eps term vanishes
u_, v_, w_ = (rng.normal(size=4) for _ in range(3))
for chord_parallel in ("k00,k11", "k01,k10"):
    if chord_parallel == "k00,k11": line = dict(k00=u_, k10=v_, k11=w_, k01=v_)   # kb = kd
    else:                           line = dict(k00=u_, k10=v_, k11=u_, k01=w_)   # k00 = k11
    assert abs(levi(line["k00"], line["k10"], line["k11"], line["k01"])) < 1e-12
    t = trchain(line["k00"], line["k10"], line["k11"], line["k01"])
    assert abs(t.imag) < 1e-12
print("momentum conservation: chord pair", chord_parallel, "-> eps term = 0; no gamma5/eps ambiguity at flat level")
a_, v2 = u_, v_
t = trchain(u_, v_, w_, v_).real
assert np.isclose(t, 2*(2*dot(u_, v_)*dot(v_, w_) - dot(v_, v_)*dot(u_, w_)))
print("class B numerator: tr = 2[2(ka.v)(v.kc) - v^2 ka.kc]  ->  N_B = %+d * that" % (kappa**2*rB))
sB = kappa**2*rB                                          # overall sign of N_B relative to +tr
# ---------------------------------------------------------------- 3. multiplicities from slot counting
def count(vertices):
    slots = {"Q": "xxXX", "C": "xxY", "c": "XXy", "F": "xff", "f": "XFF"}
    fam = []
    for sp_ in "xyf":
        src = [i for i, v in enumerate(vertices) for _ in range(slots[v].count(sp_))]
        snk = [i for i, v in enumerate(vertices) for _ in range(slots[v].count(sp_.upper()))]
        fam.append([list(zip(src, p)) for p in itertools.permutations(snk)] if len(src) == len(snk) else [])
    n, total = len(vertices), 0
    for e in itertools.product(*fam):
        edges = sum(e, [])
        par = list(range(n))
        def find(a):
            while par[a] != a: a = par[a]
            return a
        for u, v in edges: par[find(u)] = find(v)
        if len({find(i) for i in range(n)}) == 1:
            total += 1
    return total
mult = {v: count(tuple(v)) for v in ("QFf", "CcFf", "FFff")}
print("connected slot matchings:", mult)
assert mult == {"QFf": 8, "CcFf": 8, "FFff": 40}
nA = sum(1 for pi in itertools.permutations(range(4)) if pi[0]//2 == pi[1]//2 and pi[2]//2 == pi[3]//2)   # fermion matchings
nB = 24 - nA
assert (nA*1, nB*2) == (8, 32)                            # class A: 1 connected scalar choice; class B: both scalar choices

# ---------------------------------------------------------------- 4. component couplings (real-scalar basis, Martin's Y^{jIJ}, lambda^{jkl})
R1, I1, R2, I2 = sp.symbols("R1 I1 R2 I2", real=True)
ym, mm = 1.3*np.exp(0.7j), 0.9*np.exp(0.4j)
X = (R1 + sp.I*I1)/sp.sqrt(2); Yf = (R2 + sp.I*I2)/sp.sqrt(2)
WX, WY = mm*Yf + ym*X**2/2, mm*X
V = sp.expand(WX*sp.conjugate(WX) + WY*sp.conjugate(WY))
fl = [R1, I1, R2, I2]
d3 = lambda a, b, c_: complex(sp.diff(V, a, b, c_).subs({f: 0 for f in fl}))
d4 = lambda a, b, c_, d: complex(sp.diff(V, a, b, c_, d).subs({f: 0 for f in fl}))
lam3 = {(a, b, c_): d3(fl[a], fl[b], fl[c_]) for a in range(4) for b in range(4) for c_ in range(4)}
lam4 = lambda a, b, c_, d: d4(fl[a], fl[b], fl[c_], fl[d])
a_j = np.array([ym/np.sqrt(2), 1j*ym/np.sqrt(2), 0, 0])                 # Y^{j X X}
Mjk = np.outer(a_j, a_j.conj())                                          # Y^{jXX} Ybar_{kXX}
y4, z = abs(ym)**4, abs(mm)**2
SigB = sum(abs(a_j[j])**2*abs(a_j[k])**2 for j in range(4) for k in range(4))
SigA, SigAcc = np.sum(abs(Mjk)**2), np.sum(Mjk*Mjk)
SigQ = sum(lam4(j, k, l, l)*Mjk[j, k] for j in range(4) for k in range(4) for l in range(4))
SigC = sum(lam3[j, l, m]*lam3[k, l, m]*Mjk[j, k] for j in range(4) for k in range(4) for l in range(4) for m in range(4))
cp = {"B (1/2)S": SigB/2, "A (1/8)(S+cc)": (SigA + SigAcc.conjugate())/8, "QFf (1/4)S": SigQ/4, "CcFf (1/4)S": SigC/4}
print({k_: np.round(v_/y4, 10) for k_, v_ in cp.items()}, " [units |y|^4; CcFf also x z]")
assert np.isclose(cp["B (1/2)S"], y4/2) and np.isclose(cp["A (1/8)(S+cc)"], y4/8)
assert np.isclose(cp["QFf (1/4)S"], y4/2) and np.isclose(cp["CcFf (1/4)S"], z*y4/2) and abs(SigAcc) < 1e-12
# helicity-violating (mass-insertion) contractions vanish: M_XX = 0
Mm = np.array([[0, mm], [mm, 0]]); Yt = np.zeros((4, 2, 2), complex); Yt[:, 0, 0] = a_j
flip2 = sum(Yt[j, I, J]*Yt[j, Ip, Jp]*np.conj(Mm[I, Ip])*np.conj(Mm[J, Jp])
            for j in range(4) for I in range(2) for J in range(2) for Ip in range(2) for Jp in range(2))
assert abs(flip2) < 1e-12
print("final prefactors (Martin coefficient x |y|^4 units) match; F-bar F-bar contraction =", abs(flip2))

# ---------------------------------------------------------------- 5. generic-mass reduction vs Martin 1709.02397
x, w, u, zz, y, v, P = sp.symbols("x w u z y v P")
A = {m_: sp.Symbol("A_" + str(m_)) for m_ in (x, w, u, zz, y, v)}
I3 = {(m_, bub): sp.Symbol("I3_%s_%d" % (m_, bub)) for m_ in (x, w) for bub in (1, 2)}
Gx, Gw, E = sp.symbols("Gx Gw E")                       # G(x;(u,z),(y,v)), G(w;...), basketball E
A2 = (A[w] - A[x])/(x - w); I4 = {b: (I3[w, b] - I3[x, b])/(x - w) for b in (1, 2)}; K = (Gw - Gx)/(x - w)
bub0 = {1: A[u]*A[zz], 2: A[y]*A[v]}
def rule(n, i, j):                         # int G_x G_w p^{2n} b1^i b2^j, doubled lines x,w carry the same momentum
    if (i, j) == (0, 0): return [A2, A[x] - w*A2, -x*A[x] - w*A[w] + x*w*A2][n]
    if (i, j) in ((1, 0), (0, 1)):
        b = 1 if i else 2
        return [I4[b], I3[x, b] - w*I4[b], bub0[b] - x*I3[x, b] - w*I3[w, b] + x*w*I4[b]][n]
    return [K, Gx - w*K, E - x*Gx - w*Gw + x*w*K][n]
b1, b2 = sp.symbols("b1 b2")
def integrate(expr):
    poly = sp.Poly(sp.expand(expr), P, b1, b2)
    return sum(co*rule(*mon) for mon, co in poly.terms())
def Bferm(m1, m2, b): return (A[m1] + A[m2]) - (m1 + m2 + P)*b            # fermion bubble numerator k^2+q^2-p^2
cancel = lambda e_: sp.cancel(sp.together(e_))

def martin(A2_, I41_, I42_, K_):           # Martin's (3.27), (3.17), (3.24), (3.19) with the doubled-propagator functions as arguments
    J   = A[u]*((w - v - y)*I42_ - I3[x, 2] + A2_*(A[v] + A[y]))
    Kc  = (v + y - (w + x)/2)*K_ + (Gw + Gx)/2 - (A[v] + A[y])*I41_
    KA  = (((w + x)*(u + v + y + zz) - 2*(u + zz)*(v + y) - w**2 - x**2)*K_/2 - E
           + (w + x - u - v - y - zz)*(Gw + Gx)/2 + (u - w + zz)*(A[v] + A[y])*I41_
           + (v - w + y)*(A[u] + A[zz])*I42_ + (A[v] + A[y])*I3[x, 1] + (A[u] + A[zz])*I3[x, 2]
           - A2_*(A[u] + A[zz])*(A[v] + A[y]))
    KB  = ((x**2 + w**2 + 2*u*v - 2*u*y - 2*v*zz + 2*y*zz + u*w + v*w + u*x + v*x - w*y - x*y - w*zz - x*zz)*K_/4
           + (y + zz - u - v - w - x)*(Gw + Gx)/4 + (u + w - zz)*(A[v] - A[y])*I41_/2
           + (v + w - y)*(A[u] - A[zz])*I42_/2 + E/2 + (A[y] - A[v])*I3[x, 1]/2 + (A[zz] - A[u])*I3[x, 2]/2
           + A2_*(A[zz] - A[u])*(A[y] - A[v])/2)
    return dict(J=J, Kc=Kc, KA=KA, KB=KB)
M = martin(A2, I4[1], I4[2], K)
f1 = (A[u] - A[zz]) - (P + zz - u)*b1                       # 2p^2 x coefficient of v^mu in  int ka^mu G_u G_z  (fermion u, scalar z)
f2 = (A[v] - A[y]) - (P + y - v)*b2                         # same for bubble 2 (fermion v, scalar y)
mine_generic = dict(
    J  =  A[u]*integrate(Bferm(y, v, b2)),                  # QFf
    Kc = -integrate(b1*Bferm(y, v, b2)),                    # CcFf
    KA = -integrate(Bferm(u, zz, b1)*Bferm(y, v, b2)),      # FFff class A
    KB =  integrate(sp.Rational(1, 2)*f1*f2))               # FFff class B
for name in ("J", "Kc", "KA", "KB"):
    assert cancel(mine_generic[name] - M[name]) == 0, name
print("generic masses: our reductions == Martin (3.27) J_SSFFS, (3.17) K_SSSSFF, (3.24) K_SSFFFF, (3.19) K_FFFSSF (exact)")

# ---------------------------------------------------------------- 6. equal-mass closed forms
Az, T2, I21, I3z, Ez, Gz, Kz, zs = sp.symbols("Az T2 I21 I3z Ez Gz Kz zs")
Me = martin(T2, I21, I21, Kz)
def equal_mass(e_):
    e_ = e_.subs({Gx: Gz, Gw: Gz, E: Ez}, simultaneous=True)
    e_ = e_.subs({A[m_]: Az for m_ in A}, simultaneous=True).subs({I3[k_]: I3z for k_ in I3}, simultaneous=True)
    return sp.expand(e_.subs({x: zs, w: zs, u: zs, zz: zs, y: zs, v: zs}, simultaneous=True))
ours = {   # Gamma/|y|^4 from the Wick prefactors: +8/16 (QFf), -8 z/16 (CcFf), -8/64 (A), -32/64*(-1/2) (B)
  "QFf":  sp.Rational(1, 2)*Az*(2*Az*T2 - I3z - zs*I21),
  "CcFf": sp.Rational(1, 2)*zs*(zs*Kz + Gz - 2*Az*I21),
  "A":   -sp.Rational(1, 8)*(Ez + 2*zs*Gz + zs**2*Kz + 4*Az**2*T2 - 4*Az*I3z - 4*zs*Az*I21),
  "B":    sp.Rational(1, 4)*(Ez - 2*zs*Gz + zs**2*Kz)}
martin_eq = {"QFf": sp.Rational(1, 2)*equal_mass(Me["J"]), "CcFf": sp.Rational(1, 2)*zs*equal_mass(Me["Kc"])*(-1)*(-1),
             "A": sp.Rational(1, 8)*equal_mass(Me["KA"]), "B": sp.Rational(1, 2)*equal_mass(Me["KB"])}
for name in ours:
    assert sp.expand(martin_eq[name] - ours[name]) == 0, name
print("equal masses: (Martin coefficient) x (Martin function) == our Gamma/|y|^4 for QFf, CcFf, A, B")
print("ALL CHECKS PASSED  [N_B sign s_B = %+d]" % sB)
