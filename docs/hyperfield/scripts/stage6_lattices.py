"""Stage 6A: Dan's 'diagonal slide' = glue the half body-diagonal h=(1/2,...,1/2) onto Z^d / D_d, d=2..13."""
import json, os, numpy as np
from fractions import Fraction
from math import gamma, pi, factorial
CK = 'results/stage6.json'
out = json.load(open(CK)) if os.path.exists(CK) else {}
def save(): json.dump(out, open(CK, 'w'), indent=1, default=float)

def hnf_basis(gens2):
    """integer row-HNF of generator matrix (entries = 2*coords); returns basis (in 2*coords)"""
    A = [list(map(int, r)) for r in gens2]
    m, n = len(A), len(A[0]); row = 0
    for col in range(n):
        piv = [i for i in range(row, m) if A[i][col] != 0]
        if not piv: continue
        while True:
            piv = [i for i in range(row, m) if A[i][col] != 0]
            p = min(piv, key=lambda i: abs(A[i][col]))
            A[row], A[p] = A[p], A[row]
            done = True
            for i in range(row+1, m):
                if A[i][col]:
                    q = A[i][col] // A[row][col]
                    A[i] = [a - q*b for a, b in zip(A[i], A[row])]
                    if A[i][col]: done = False
            if done: break
        row += 1
    return np.array(A[:row], dtype=float)

def lll(B, delta=0.99):
    B = B.copy(); n = len(B)
    def gs(B):
        Bs = np.zeros_like(B); mu = np.zeros((n, n))
        for i in range(n):
            Bs[i] = B[i]
            for j in range(i):
                mu[i, j] = B[i] @ Bs[j] / (Bs[j] @ Bs[j]); Bs[i] -= mu[i, j]*Bs[j]
        return Bs, mu
    Bs, mu = gs(B); k = 1
    while k < n:
        for j in range(k-1, -1, -1):
            q = round(mu[k, j])
            if q: B[k] -= q*B[j]; Bs, mu = gs(B)
        if Bs[k] @ Bs[k] >= (delta - mu[k, k-1]**2)*(Bs[k-1] @ Bs[k-1]): k += 1
        else: B[[k, k-1]] = B[[k-1, k]]; Bs, mu = gs(B); k = max(k-1, 1)
    return B

def enum(B, C):
    G = B @ B.T; n = len(G); Q = np.zeros((n, n))
    for i in range(n):
        Q[i, i] = G[i, i] - sum(Q[k, k]*Q[k, i]**2 for k in range(i))
        for j in range(i+1, n):
            Q[i, j] = (G[i, j] - sum(Q[k, k]*Q[k, i]*Q[k, j] for k in range(i)))/Q[i, i]
    res = []; x = np.zeros(n, dtype=int)
    def rec(i, rem):
        if i < 0:
            if np.any(x): res.append(x.copy())
            return
        c = -sum(Q[i, j]*x[j] for j in range(i+1, n)); w = np.sqrt(max(rem, 0)/Q[i, i])
        for xi in range(int(np.ceil(c-w-1e-9)), int(np.floor(c+w+1e-9))+1):
            x[i] = xi; t = Q[i, i]*(xi-c)**2
            if t <= rem+1e-9: rec(i-1, rem-t)
        x[i] = 0
    rec(n-1, C)
    return np.array(res) @ B if res else np.zeros((0, n))

def minvecs(B):
    m = min(np.sum(B**2, 1))
    V = enum(B, m*1.0000001)
    nr = np.sum(V**2, 1); mn = nr.min()
    return mn, V[np.abs(nr-mn) < 1e-7]

def aut_order_bruteforce(B, cap=10**6):
    """count isometries: images of an LLL basis among lattice vectors of equal norm preserving Gram"""
    G = B @ B.T; n = len(B)
    cands = {}
    for i in range(n):
        V = enum(B, G[i, i]*1.0000001); cands[i] = V[np.abs(np.sum(V**2, 1)-G[i, i]) < 1e-7]
    cnt = 0; img = []
    def rec(i):
        nonlocal cnt
        if cnt > cap: return
        if i == n: cnt += 1; return
        for v in cands[i]:
            if all(abs(v @ img[j] - G[i, j]) < 1e-7 for j in range(i)):
                img.append(v); rec(i+1); img.pop()
    rec(0); return cnt

BEST_LAT_KISS = {2: 6, 3: 12, 4: 24, 5: 40, 6: 72, 7: 126, 8: 240, 9: 272, 10: 336, 11: 438, 12: 756, 13: 918}
BEST_LAT_DELTA = {2: 0.28868, 3: 0.17678, 4: 0.125, 5: 0.08839, 6: 0.07217, 7: 0.0625, 8: 0.0625,
                  9: 0.04419, 10: 0.03608, 11: 0.03208, 12: 0.03704, 13: 0.03208}
BEST_ANY_KISS_LB = {2: 6, 3: 12, 4: 24, 5: 40, 6: 72, 7: 126, 8: 240, 9: 306, 10: 510, 11: 604, 12: 841, 13: 1154}

def aut_formula(name, d):
    if name == 'Z^d': return 2**d*factorial(d)
    if name == 'D_d': return 1152 if d == 4 else (8 if d == 2 else 2**d*factorial(d))
    if name == 'D_d^+':
        return {4: 384, 8: 696729600}.get(d, 2**(d-1)*factorial(d) if d > 2 else 4)
    if name == 'BCC Z^d u (Z^d+h)':
        return {4: 1152, 2: 8, 3: 48}.get(d, 2**d*factorial(d))
def analyse(name, gens, d, h_in_base_norm=None):
    gens2 = np.round(2*np.array(gens)).astype(int)
    B = hnf_basis(gens2)/2
    if len(B) != d: return dict(lattice=False)
    B = lll(B)
    G = B @ B.T; det = float(np.linalg.det(G))
    integral = bool(np.allclose(G, np.round(G), atol=1e-9))
    even = integral and bool(np.allclose(np.diag(G) % 2, 0, atol=1e-9))
    G2 = 2*G; integral_sqrt2 = bool(np.allclose(G2, np.round(G2), atol=1e-9))
    even_sqrt2 = integral_sqrt2 and bool(np.allclose(np.diag(G2) % 2, 0, atol=1e-9))
    mn, MV = minvecs(B)
    delta = (np.sqrt(mn)/2)**d/np.sqrt(det)
    Vd = pi**(d/2)/gamma(d/2+1)
    r = dict(lattice=True, det_gram=round(det, 8), integral=integral, even=even, unimodular=integral and abs(det-1) < 1e-8,
             even_after_sqrt2_scaling=even_sqrt2, det_after_sqrt2=round(det*2**d, 6),
             min_norm=round(float(mn), 8), kissing=int(len(MV)), center_density=float(delta), packing_density=float(delta*Vd),
             delta_over_best_lattice=float(delta/BEST_LAT_DELTA[d]), kiss_over_best_lattice=len(MV)/BEST_LAT_KISS[d],
             aut_order=aut_formula(name, d), aut_source='standard (Conway-Sloane ch.4/8); brute-force checked for d<=5')
    if d <= 5:
        r['aut_bruteforce'] = aut_order_bruteforce(B)
    return r, MV

if 'partA' not in out:
    out['partA'] = {}
for d in range(2, 14):
    key = str(d)
    if key in out['partA']: continue
    h = np.full(d, 0.5); I = np.eye(d)
    Dg = [I[i]-I[i+1] for i in range(d-1)] + [I[d-2]+I[d-1]]
    res = dict(h_norm=d/4)
    for name, gens, closure in [('Z^d', list(I), True),
                                ('D_d', Dg, True),
                                ('D_d^+', Dg + [h], d % 2 == 0),        # 2h=(1..1) in D_d iff d even
                                ('BCC Z^d u (Z^d+h)', list(I) + [h], True)]:
        if not closure:
            res[name] = dict(lattice=False, reason='2h=(1,...,1) has odd coordinate sum, not in D_d: D_d u (D_d+h) not closed'); continue
        a, MV = analyse(name, gens, d)
        res[name] = a
    # slide-rule criteria
    res['criterion_h_on_first_shell_of_Zd'] = abs(d/4 - 1) < 1e-12
    res['criterion_h_on_first_shell_of_Dd'] = abs(d/4 - 2) < 1e-12
    res['criterion_even_unimodular'] = bool(res['D_d^+'].get('even') and res['D_d^+'].get('unimodular'))
    out['partA'][key] = res; save(); print(d, {k: (v.get('kissing'), v.get('det_gram'), v.get('aut_order'), v.get('aut_bruteforce')) for k, v in res.items() if isinstance(v, dict)}, flush=True)

# selection summary + the d=9 / d=13 statement
sel = dict(h_norm_equals_Zd_min=[int(k) for k, v in out['partA'].items() if v['criterion_h_on_first_shell_of_Zd']],
           h_norm_equals_Dd_min=[int(k) for k, v in out['partA'].items() if v['criterion_h_on_first_shell_of_Dd']],
           even_unimodular=[int(k) for k, v in out['partA'].items() if v['criterion_even_unimodular']],
           best_lattice_kissing_reached=[(int(k), n) for k, v in out['partA'].items() for n in ['Z^d', 'D_d', 'D_d^+', 'BCC Z^d u (Z^d+h)']
                                         if v[n].get('lattice') and abs(v[n]['kiss_over_best_lattice']-1) < 1e-12])
out['partA_selection'] = sel; save(); print(sel)
