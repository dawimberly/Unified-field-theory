"""Is the reference (24-cell / E8 / Lambda9-kissing) configuration a local minimum? perturb by eps and re-minimize."""
import json, sys, numpy as np
sys.argv = ['x', 'none']
exec(open('scripts/stage6_thomson.py').read().split("cases = [")[0])
from scipy.optimize import minimize
which = sys.argv_cases if hasattr(sys, 'argv_cases') else None
import os
CK2 = 'results/stage6_thomson_stability_' + os.environ.get('CASE', 'D4') + '.json'
res = json.load(open(CK2)) if os.path.exists(CK2) else {}
name = os.environ.get('CASE', 'D4'); R = ref_config(name); N, d = R.shape
for pn in ['riesz1', 'harmonic', 'gauss2', 'gauss6', 'log']:
    if pn in res: continue
    f = pot(pn, d); Eref = energy(R, f); rng = np.random.default_rng(7)
    def fg(z):
        Y = z.reshape(N, d); nr = np.linalg.norm(Y, axis=1, keepdims=True); X = Y/nr
        G = X @ X.T; t = np.maximum(2-2*G, 1e-12); np.fill_diagonal(t, 1.0)
        v, dv = f(t); np.fill_diagonal(v, 0); np.fill_diagonal(dv, 0)
        gX = -2*(dv @ X); return v.sum()/2, ((gX - np.sum(gX*X, 1, keepdims=True)*X)/nr).ravel()
    Es = []
    for eps in (1e-3, 1e-2, 5e-2):
        z = (R + eps*rng.normal(size=R.shape)).ravel()
        r = minimize(fg, z, jac=True, method='L-BFGS-B', options=dict(maxiter=20000, gtol=1e-12, ftol=1e-16))
        Es.append(float((r.fun-Eref)/abs(Eref)))
    res[pn] = dict(E_ref=float(Eref), rel_dE_after_perturb=Es, local_min=bool(all(abs(e) < 1e-9 for e in Es[:2])))
    json.dump(res, open(CK2, 'w'), indent=1); print(name, pn, res[pn], flush=True)
