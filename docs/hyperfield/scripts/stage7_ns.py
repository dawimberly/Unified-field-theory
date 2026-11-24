"""Stage 7.4: neutron stars with a first-order (lock-freezing / deconfinement-analog) transition: twin stars?
Hadronic EOS = SLy piecewise polytrope (Read et al. 2009, parameters from memory, validated by R1.4/Mmax output);
quark phase = constant-speed-of-sound (CSS) with jump Delta eps at p_t (parameters scanned: NOT predicted)."""
import json, os, numpy as np
from scipy.integrate import solve_ivp
CK = 'results/stage7.json'
out = json.load(open(CK))
def save(): json.dump(out, open(CK, 'w'), indent=1, default=float)
G, c, Msun = 6.67430e-8, 2.99792458e10, 1.98847e33
geo = G/c**4*1e10          # erg/cm^3 -> 1/km^2
Mkm = G*Msun/c**2/1e5      # 1 Msun in km
G0, K0 = 1.35692, 3.59389e13
r1, r2, p1 = 10**14.7, 10**15.0, 10**34.384
G1, G2, G3 = 3.005, 2.988, 2.851
K1 = p1/r1**G1; K2 = p1/r1**G2; K3 = K2*r2**(G2-G3)
r0 = (K0/K1)**(1/(G1-G0))
pieces = [(0, r0, K0, G0), (r0, r1, K1, G1), (r1, r2, K2, G2), (r2, 1e17, K3, G3)]
a = [0.0]
for i in range(1, 4):
    rb = pieces[i][0]; Kp, Gp = pieces[i-1][2], pieces[i-1][3]; Ki, Gi = pieces[i][2], pieces[i][3]
    a.append(a[-1] + Kp*rb**(Gp-1)/((Gp-1)*c**2) - Ki*rb**(Gi-1)/((Gi-1)*c**2))
def had(rho):
    for i, (lo, hi, K, Gm) in enumerate(pieces):
        if rho < hi or i == 3:
            p = K*rho**Gm; e = rho*c**2*(1+a[i]) + K*rho**Gm/(Gm-1); return p, e
rho_grid = np.logspace(3, 16.2, 6000)
P = np.array([had(x)[0] for x in rho_grid])*geo; E = np.array([had(x)[1] for x in rho_grid])*geo; RHO = rho_grid*c**2*geo
lP = np.log(P)
def eos_factory(pt=None, de_frac=0.0, cs2=1.0):
    if pt is not None:
        et = np.exp(np.interp(np.log(pt), lP, np.log(E)))
    def eps(p):
        if pt is not None and p > pt: return et*(1+de_frac) + (p-pt)/cs2
        return np.exp(np.interp(np.log(p), lP, np.log(E)))
    def cs2f(p):
        if pt is not None and p > pt: return cs2
        h = 1e-4; return (2*h*p)/(eps(p*(1+h))-eps(p*(1-h)))
    def rhob(p):
        return np.exp(np.interp(np.log(p), lP, np.log(RHO))) if (pt is None or p <= pt) else np.nan
    return eps, cs2f, rhob, (et if pt is not None else None)
def tov(pc, eps, cs2f, rhob, pt=None, et=None, de_frac=0.0):
    psurf = P[0]*10
    def rhs(r, s):
        m, p, y, mb = s
        if p <= psurf: return [0, 0, 0, 0]
        e = eps(p); f = 1-2*m/r
        dm = 4*np.pi*r**2*e
        dp = -(e+p)*(m+4*np.pi*r**3*p)/(r**2*f)
        F = (1-4*np.pi*r**2*(e-p))/f
        Q = 4*np.pi*(5*e+9*p+(e+p)/cs2f(p))/f - 6/(r**2*f) - 4*(m+4*np.pi*r**3*p)**2/(r**4*f**2)
        dy = (-y**2 - y*F - r**2*Q)/r
        rb = rhob(p); dmb = 4*np.pi*r**2*(rb if np.isfinite(rb) else 0)/np.sqrt(f)
        return [dm, dp, dy, dmb]
    r0_ = 1e-6; e0 = eps(pc)
    s0 = [4/3*np.pi*r0_**3*e0, pc, 2.0, 0.0]
    ev_s = lambda r, s: s[1]-psurf; ev_s.terminal = True; ev_s.direction = -1
    if pt is not None and pc > pt:
        ev_t = lambda r, s: s[1]-pt; ev_t.terminal = True; ev_t.direction = -1
        sol = solve_ivp(rhs, (r0_, 50), s0, events=ev_t, rtol=1e-8, atol=1e-12, max_step=0.05)
        r_ = sol.t[-1]; m, p, y, mb = sol.y[:, -1]
        e_in = et*(1+de_frac); e_out = et
        y = y - 4*np.pi*r_**3*(e_in-e_out)/(m+4*np.pi*r_**3*p)
        sol2 = solve_ivp(rhs, (r_, 50), [m, pt*(1-1e-9), y, mb], events=ev_s, rtol=1e-8, atol=1e-12, max_step=0.05)
        R = sol2.t[-1]; M, _, yR, MB = sol2.y[:, -1]
    else:
        sol = solve_ivp(rhs, (r0_, 50), s0, events=ev_s, rtol=1e-8, atol=1e-12, max_step=0.05)
        R = sol.t[-1]; M, _, yR, MB = sol.y[:, -1]
    C = M/R; yv = yR
    k2 = (8*C**5/5)*(1-2*C)**2*(2+2*C*(yv-1)-yv)/(2*C*(6-3*yv+3*C*(5*yv-8)) + 4*C**3*(13-11*yv+C*(3*yv-2)+2*C**2*(1+yv)) + 3*(1-2*C)**2*(2-yv+2*C*(yv-1))*np.log(1-2*C))
    Lam = 2*k2/(3*C**5)
    return M/Mkm, R, Lam, MB/Mkm
def sequence(pt=None, de=0.0, cs2=1.0, n=70):
    eps, cs2f, rhob, et = eos_factory(pt, de, cs2)
    pcs = np.logspace(np.log10(P[np.searchsorted(RHO, 2.0e14*c**2*geo)]), np.log10(P[-1]*0.9 if pt is None else max(P[-1]*0.9, pt*80)), n)
    rows = []
    for pc in pcs:
        try: rows.append((pc,)+tov(pc, eps, cs2f, rhob, pt, et, de))
        except Exception as ex: pass
    return np.array(rows)
def analyse(seq, tol=5e-3):
    pc, M, R, L, MB = seq.T
    # extrema with hysteresis tolerance: alternate running max / running min
    ext = [0]; direction = 1; best = 0
    for i in range(1, len(M)):
        if direction == 1:
            if M[i] > M[best]: best = i
            elif M[i] < M[best] - tol: ext.append(best); direction = -1; best = i
        else:
            if M[i] < M[best]: best = i
            elif M[i] > M[best] + tol: ext.append(best); direction = 1; best = i
    ext.append(best if direction == 1 else len(M)-1)
    br = [(ext[k], ext[k+1]) for k in range(0, len(ext)-1, 2)]
    br = [(a, b) for a, b in br if M[b] - M[a] > tol]
    Mmax = float(M.max())
    res = dict(M_max=Mmax, n_stable_branches=len(br), twin=len(br) >= 2, branch_Mmax=[float(M[b]) for a, b in br])
    for a, b in br:
        idx = np.arange(a, b+1)
        if M[idx].min() <= 1.4 <= M[idx].max():
            o = np.argsort(M[idx]); bb = idx[o]
            res.setdefault('R14_km', []).append(float(np.interp(1.4, M[bb], R[bb])))
            res.setdefault('Lam14', []).append(float(np.interp(1.4, M[bb], L[bb])))
    return res
if 'ns' not in out: out['ns'] = {}
ns = out['ns']
if 'SLy' not in ns:
    s = sequence(None); A = analyse(s)
    pc, M, R, L, MB = s.T
    o = np.argsort(M[:np.argmax(M)+1])
    A['binding_fraction_1.4'] = float((np.interp(1.4, M[o], MB[o])-1.4)/1.4)
    A['validation'] = 'literature SLy: Mmax~2.05, R1.4~11.7 km, Lambda1.4~300 (from memory)'
    ns['SLy'] = A; ns['SLy_curve'] = s[:, 1:4].tolist(); save(); print('SLy', A, flush=True)
rho_sat = 2.7e14
for nfac in [1.5, 2.0, 3.0, 4.0]:
    pt = had(nfac*rho_sat)[0]*geo
    for cs2 in [1/3, 1.0]:
        for de in [0.2, 0.4, 0.6, 0.8, 1.0]:
            key = f'n{nfac}_cs{cs2:.2f}_de{de}'
            if key in ns: continue
            s = sequence(pt, de, cs2); A = analyse(s)
            A.update(n_t_over_nsat=nfac, p_t_MeVfm3=pt/geo/1.602176634e33, cs2=cs2, de_over_et=de,
                     seidov_unstable=bool(de > 0.5+1.5*pt/np.exp(np.interp(np.log(pt), lP, np.log(E)))))
            r14 = A.get('R14_km', []); l14 = A.get('Lam14', [])
            A['passes'] = bool(A['M_max'] >= 2.01 and any(11.0 <= x <= 13.5 for x in r14) and any(x <= 800 for x in l14))
            ns[key] = A; ns[key+'_curve'] = s[:, 1:4].tolist(); save(); print(key, {k: v for k, v in A.items()}, flush=True)
