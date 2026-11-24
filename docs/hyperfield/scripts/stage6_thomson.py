"""Stage 6B(2): does energy minimization of N points on S^{d-1} ('churn crystallizing') produce D4 / E8 / Lambda9 kissing configurations?"""
import json, os, sys, time, numpy as np
from scipy.optimize import minimize
CK = 'results/stage6_thomson_' + (sys.argv[1] if len(sys.argv) > 1 else 'all') + '.json'
out = json.load(open(CK)) if os.path.exists(CK) else {}
def save(): json.dump(out, open(CK, 'w'), indent=1, default=float)
exec(open('scripts/stage6_lattices.py').read().split("BEST_LAT_KISS")[0].split("CK = ")[0])  # imports only
src = open('scripts/stage6_lattices.py').read()
ns = {}; exec(src.split("BEST_LAT_KISS")[0].replace("out = json.load(open(CK)) if os.path.exists(CK) else {}", "out={}"), ns)
def ref_config(name):
    if name == 'D4':
        I = np.eye(4); V = [s1*I[i]+s2*I[j] for i in range(4) for j in range(i+1, 4) for s1 in (1, -1) for s2 in (1, -1)]
        return np.array(V)/np.sqrt(2)
    if name == 'E8':
        I = np.eye(8); V = [s1*I[i]+s2*I[j] for i in range(8) for j in range(i+1, 8) for s1 in (1, -1) for s2 in (1, -1)]
        for k in range(256):
            s = np.array([1 if (k >> b) & 1 else -1 for b in range(8)])
            if np.sum(s < 0) % 2 == 0: V.append(s/2)
        return np.array(V)/np.sqrt(2)
    if name == 'L9':
        E8 = np.array([[2,0,0,0,0,0,0,0],[-1,1,0,0,0,0,0,0],[0,-1,1,0,0,0,0,0],[0,0,-1,1,0,0,0,0],
                       [0,0,0,-1,1,0,0,0],[0,0,0,0,-1,1,0,0],[0,0,0,0,0,-1,1,0],[.5]*8])
        B = np.zeros((9, 9)); B[:8, :8] = E8; B[8, 0] = 1; B[8, 8] = 1
        V = ns['enum'](B, 2.0000001); V = V[np.abs(np.sum(V**2, 1)-2) < 1e-7]
        return V/np.sqrt(2)
POTS = {'riesz1': lambda t: (t**-0.5, -0.5*t**-1.5),
        'harmonic': None,   # s = d-2 set per d
        'gauss2': lambda t: (np.exp(-2*t), -2*np.exp(-2*t)),
        'gauss6': lambda t: (np.exp(-6*t), -6*np.exp(-6*t)),
        'log': lambda t: (-0.5*np.log(t), -0.5/t)}
def pot(name, d):
    if name == 'harmonic':
        s = d-2; return lambda t: (t**(-s/2), -s/2*t**(-s/2-1))
    return POTS[name]
def energy(X, f):
    G = X @ X.T; t = np.maximum(2-2*G, 1e-12); iu = np.triu_indices(len(X), 1)
    return np.sum(f(t[iu])[0])
def solve(N, d, f, seed):
    rng = np.random.default_rng(seed)
    def fg(z):
        Y = z.reshape(N, d); nr = np.linalg.norm(Y, axis=1, keepdims=True); X = Y/nr
        G = X @ X.T; t = np.maximum(2-2*G, 1e-12); np.fill_diagonal(t, 1.0)
        v, dv = f(t); np.fill_diagonal(v, 0); np.fill_diagonal(dv, 0)
        E = v.sum()/2
        gX = -2*(dv @ X)                       # dE/dX (t=2-2x.y)
        gY = (gX - np.sum(gX*X, 1, keepdims=True)*X)/nr
        return E, gY.ravel()
    r = minimize(fg, rng.normal(size=N*d), jac=True, method='L-BFGS-B', options=dict(maxiter=20000, gtol=1e-11, ftol=1e-15))
    Y = r.x.reshape(N, d); return Y/np.linalg.norm(Y, axis=1, keepdims=True), r.fun
def ip_signature(X):
    G = np.round(X @ X.T, 4); u, c = np.unique(G, return_counts=True)
    return u, c
cases = [('D4', 24, 4, 12), ('E8', 240, 8, 4), ('L9', 272, 9, 4)]
which = sys.argv[1:] or [c[0] for c in cases]
for name, N, d, nstart in cases:
    if name not in which: continue
    R = ref_config(name); assert len(R) == N, (name, len(R))
    for pn in ['riesz1', 'harmonic', 'gauss2', 'gauss6', 'log']:
        key = f'{name}_{pn}'
        if key in out: continue
        f = pot(pn, d); Eref = energy(R, f); best = None; t0 = time.time(); Es = []
        for s in range(nstart):
            X, E = solve(N, d, f, s); Es.append(float(E))
            if best is None or E < best[1]: best = (X, E)
        X, E = best
        u, c = ip_signature(X); ur, cr = ip_signature(R)
        same = (len(u) == len(ur)) and np.allclose(u, ur, atol=2e-3) and np.array_equal(c, cr)
        out[key] = dict(N=N, d=d, E_ref=float(Eref), E_min_found=float(E), rel_diff=float((E-Eref)/abs(Eref)),
                        all_start_energies=Es, n_starts_hitting_ref=int(sum(abs(e-Eref) < 1e-8*abs(Eref) for e in Es)),
                        matches_ref_innerproducts=bool(same), max_ip_found=float(np.max(X@X.T-2*np.eye(N))),
                        max_ip_ref=float(np.max(R@R.T-2*np.eye(N))), n_distinct_ip_found=int(len(u)), seconds=time.time()-t0)
        save(); print(key, out[key], flush=True)
