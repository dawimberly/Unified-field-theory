"""Stage 10: detectability of the gravitational-wave background from the first-order 'lock' transition."""
import json, os, numpy as np
from math import pi, log, sqrt
CK = 'results/stage10.json'
out = json.load(open(CK)) if os.path.exists(CK) else {}
def save(): json.dump(out, open(CK, 'w'), indent=1, default=float)
H100 = 3.2407793e-18          # s^-1  (100 km/s/Mpc)
YR = 3.15576e7
MBAR = 2.435e18

# ---------------- templates: Caprini et al. 2016 (arXiv:1512.06239, eqs 7-20, verified)
def kappa_v(a):  return a/(0.73+0.083*np.sqrt(a)+a)
def spectra(f, T, alpha, bH, vw=0.95, gs=106.75, eps_turb=0.05, kappa_phi=0.0, suppress=True):
    kv = kappa_v(alpha); K = alpha/(1+alpha); gfac = (100/gs)**(1/3); g6 = (gs/100)**(1/6); T2 = T/100.0
    fsw = 1.9e-5/vw*bH*T2*g6
    Ssw = (f/fsw)**3*(7/(4+3*(f/fsw)**2))**3.5
    sw = 2.65e-6/bH*(kv*K)**2*gfac*vw*Ssw
    if suppress:   # finite sound-wave lifetime (Caprini et al. 2020 / Guo et al. 2020)
        HR = (8*pi)**(1/3)*vw/bH; Uf = np.sqrt(0.75*kv*K); tau = HR/Uf
        sw = sw*(1-1/np.sqrt(1+2*tau))
    hs = 1.65e-5*T2*g6
    fturb = 2.7e-5/vw*bH*T2*g6
    Sturb = (f/fturb)**3/((1+f/fturb)**(11/3)*(1+8*pi*f/hs))
    turb = 3.35e-4/bH*(eps_turb*kv*K)**1.5*gfac*vw*Sturb
    env = 0*f
    if kappa_phi > 0:
        fenv = 1.65e-5*(0.62/(1.8-0.1*vw+vw**2))*bH*T2*g6
        Senv = 3.8*(f/fenv)**2.8/(1+2.8*(f/fenv)**3.8)
        env = 1.67e-5/bH**2*(kappa_phi*K)**2*gfac*(0.11*vw**3/(0.42+vw**2))*Senv
    return sw, turb, env, fsw

# ---------------- detectors (sky-averaged strain sensitivity S_n) -> h^2 Omega_n = 2 pi^2 f^3 S_n / (3 H100^2)
def lisa_Sn(f, galactic=True):    # Robson, Cornish & Liu 2019 (arXiv:1803.01944), from memory
    L = 2.5e9; fs = 19.09e-3
    Poms = (1.5e-11)**2*(1+(2e-3/f)**4)
    Pacc = (3e-15)**2*(1+(0.4e-3/f)**2)*(1+(f/8e-3)**4)
    Sn = 10/(3*L**2)*(Poms+2*(1+np.cos(f/fs)**2)*Pacc/(2*pi*f)**4)*(1+0.6*(f/fs)**2)
    if galactic:
        A, al, be, ka, ga, fk = 9e-45, 0.138, -221, 521, 1680, 0.00113
        Sn = Sn + A*f**(-7/3)*np.exp(-f**al+be*f*np.sin(ka*f))*(1+np.tanh(ga*(fk-f)))
    return Sn
def bbo_Sn(f):     # Yagi & Seto 2011 (from memory)
    return 2.00e-49*f**2 + 4.58e-49 + 1.26e-51*f**-4
def decigo_Sn(f):  # Yagi & Seto 2011 (from memory)
    fp = 7.36
    return 7.05e-48*(1+(f/fp)**2) + 4.8e-51*f**-4/(1+(f/fp)**2) + 5.33e-52*f**-4
def ground_Sn(f, floor, fl, fh):   # crude analytic ET/CE: seismic wall ~f^-8 below fl, flat floor, shot noise above fh (ASSUMED, factor ~3)
    return floor**2*((fl/f)**8 + 1 + (f/fh)**2)
DET = {'LISA': (lambda f: lisa_Sn(f), (1e-5, 1.0), 3.0),
       'LISA (no galactic foreground)': (lambda f: lisa_Sn(f, False), (1e-5, 1.0), 3.0),
       'DECIGO': (decigo_Sn, (1e-3, 100), 4.0), 'BBO': (bbo_Sn, (1e-3, 100), 4.0),
       'Einstein Telescope (crude)': (lambda f: ground_Sn(f, 3e-25, 8, 300), (1, 5000), 4.0),
       'Cosmic Explorer (crude)': (lambda f: ground_Sn(f, 2e-25, 6, 500), (1, 5000), 4.0)}
def omega_n(Sn, f): return 2*pi**2*f**3*Sn(f)/(3*H100**2)
def snr(f, h2om, Sn, Tyr):
    on = omega_n(Sn, f); return np.sqrt(Tyr*YR*np.trapezoid((h2om/on)**2, f))
FG = np.logspace(-6, 4, 6000)
def total(f, *a, **k):
    sw, tu, en, fsw = spectra(f, *a, **k); return sw+tu+en, fsw
def all_snr(T, alpha, bH, **k):
    res = {}
    for name, (Sn, (lo, hi), Ty) in DET.items():
        f = FG[(FG >= lo) & (FG <= hi)]; h2, _ = total(f, T, alpha, bH, **k); res[name] = float(snr(f, h2, Sn, Ty))
    return res
# LVK O3: Omega_GW <= 5.8e-9 flat, 20-76.6 Hz (95%), H0=67.9 -> h^2 Omega <= 2.67e-9 (verified PRD 104 022004)
LVK_h2 = 5.8e-9*0.679**2

# ---------------- 1. lock-model derivation (thin wall, phi^6 of Stages 5/6)
def lock_params(lam, a1, v_over_Tc, Tc=6000.0, gs=106.75):
    """V = lam v^4 [(p^2-1)^2 (p^2+a)/4], a(T) = a1*(1 - T/Tc) (ASSUMED linear), v T-independent (ASSUMED).
    sigma(a=0) = v^3 sqrt(lam)/(4 sqrt2), dV = lam v^4 a/4 (DERIVED). S3 = 16 pi sigma^3/(3 dV^2)."""
    v = v_over_Tc*Tc
    sig = v**3*sqrt(lam)/(4*sqrt(2))
    Sc = 146 - 2*log(gs/100) - 4*log(Tc/100)          # S3/T at nucleation (radiation era, beta/H~100; standard estimate)
    # S3/T = 16 pi sig^3 / (3 T (lam v^4 a1 d/4)^2), T = Tc(1-d) -> solve for d
    def S(d): return 16*pi*sig**3/(3*Tc*(1-d)*(lam*v**4*a1*d/4)**2)
    d_spin = 0.5/a1                                       # symmetric-phase spinodal a = 1/2 (Stage 5 hysteresis window, DERIVED)
    if d_spin >= 1:   # a(T=0)=a1 < 1/2: symmetric phase stays metastable down to T=0 -> no thermal completion
        return dict(lam=lam, a1=a1, v_over_Tc=v_over_Tc, delta=None, Tn=None, S3_over_T=Sc, beta_over_H=None, alpha=None,
                    regime='no thermal completion (a1<1/2: spinodal never reached; vacuum-dominated, templates invalid)', a_n=None)
    lo, hi = 1e-8, d_spin*0.999
    if S(hi) > Sc:     # barrier never low enough before spinodal: transition at the spinodal (thin wall invalid)
        d = hi; regime = 'spinodal (strong supercooling; thin-wall invalid)'
    else:
        for _ in range(200):
            m = sqrt(lo*hi)
            if S(m) > Sc: lo = m
            else: hi = m
        d = sqrt(lo*hi); regime = 'thin-wall' if a1*d < 0.15 else 'thin-wall marginal'
    Tn = Tc*(1-d)
    bH = Sc*(2*(1-d)/d - 1)
    dtheta = lam*v**4*a1/4*(d + (1-d)/4)
    alpha = dtheta/(pi**2*gs*Tn**4/30)
    return dict(lam=lam, a1=a1, v_over_Tc=v_over_Tc, delta=d, Tn=Tn, S3_over_T=Sc, beta_over_H=bH, alpha=alpha, regime=regime, a_n=a1*d)
if 'lock_scan' not in out:
    rows = []
    for lam in (0.1, 0.3, 1.0):
        for a1 in (0.3, 1.0, 3.0):
            for vt in (0.5, 1.0, 2.0):
                p = lock_params(lam, a1, vt); p['snr'] = all_snr(6000.0, p['alpha'], p['beta_over_H']) if p['alpha'] else None; rows.append(p)
    out['lock_scan'] = rows
    out['lock_notes'] = dict(
        derived=['sigma = v^3 sqrt(lam)/(4 sqrt2)', 'Delta V = lam v^4 a/4', 'thin-wall S3', 'beta/H = (S3/T)[2(1-delta)/delta - 1]',
                 'max supercooling a_n <= 1/2 (spinodal of the symmetric phase; Stage-5 hysteresis window -1<a<1/2)'],
        assumed=['a(T) = a1 (1 - T/Tc) linear', 'v independent of T', 'lam, a1, v/Tc in natural O(1) box (0.1-1, 0.3-3, 0.5-2)', 'Tc = 6 TeV (Stage 5)',
                 'S3/T_n ~ 146 - 4 ln(T/100 GeV) nucleation criterion', 'v_w = 0.95, g* = 106.75, eps_turb = 0.05, no runaway (kappa_phi = 0)'],
        note_9_8='The Stage-6 9/8 ratio belongs to the cubic-quartic Landau form; for the phi^6 lock the equivalent statement is the window 0 < a_n < 1/2 used here.')
    save()

# ---------------- 2/3. benchmark spectra, T scan, map
if 'benchmarks' not in out:
    B = {}
    for name, (T, a, b) in {'stage5_illustrative (alpha=0.1, beta/H=100, T=6 TeV)': (6000, 0.1, 100),
                            'strong (alpha=0.5, beta/H=10, T=6 TeV)': (6000, 0.5, 10),
                            'lock natural median': (6000, None, None)}.items():
        if a is None:
            V=[r for r in out['lock_scan'] if r['alpha']]; p = sorted(V, key=lambda r: r['alpha'])[len(V)//2]; a, b = p['alpha'], p['beta_over_H']
        f = FG; h2, fsw = total(f, T, a, b)
        B[name] = dict(T=T, alpha=a, beta_over_H=b, f_sw_Hz=fsw, peak_h2Omega=float(h2.max()), f_peak_Hz=float(f[np.argmax(h2)]),
                       snr=all_snr(T, a, b), h2Omega_25Hz=float(np.interp(25, f, h2)), LVK_O3_h2_limit=LVK_h2)
    out['benchmarks'] = B; save()
if 'T_scan' not in out:
    Ts = np.logspace(2, 5, 25); res = []
    for T in Ts:
        for a, b in [(0.1, 100), (0.5, 10), (0.01, 1000)]:
            h2, fsw = total(FG, T, a, b)
            res.append(dict(T=float(T), alpha=a, beta_over_H=b, f_sw=fsw, snr=all_snr(T, a, b), LVK_excluded=bool(np.interp(25, FG, h2) > LVK_h2)))
    out['T_scan'] = res; save()
if 'map' not in out:
    al = np.logspace(-3, 1, 41); bh = np.logspace(0, 4, 41); M = {k: np.zeros((41, 41)) for k in ('LISA', 'DECIGO', 'BBO')}
    for i, a in enumerate(al):
        for j, b in enumerate(bh):
            s = all_snr(6000.0, a, b)
            for k in M: M[k][i, j] = s[k]
    out['map'] = dict(alpha=al.tolist(), beta_over_H=bh.tolist(), **{k: v.tolist() for k, v in M.items()}); save()
# minimal alpha for LISA SNR>10 at given beta/H (T=6 TeV)
m = out['map']; al = np.array(m['alpha']); bh = np.array(m['beta_over_H']); L = np.array(m['LISA'])
thr = {}
for j, b in enumerate(bh):
    idx = np.where(L[:, j] > 10)[0]; thr[f'{b:.3g}'] = float(al[idx[0]]) if len(idx) else None
out['lisa_threshold_alpha_vs_betaH'] = thr; save()
print(json.dumps(out['benchmarks'], indent=1))
for r in out['lock_scan']: print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k != 'snr'}, r['snr'] and {k: round(v, 3) for k, v in r['snr'].items()})
for n,(Sn,(lo,hi),T) in DET.items():
    f=FG[(FG>=lo)&(FG<=hi)]; on=omega_n(Sn,f); print(n,'min h2Omega_n %.2e at %.3g Hz'%(on.min(), f[on.argmin()]))
print({k:round(v['snr']['LISA'],2) for k,v in out['benchmarks'].items()})
for r in out['T_scan'][::6]: print(r['T'], r['alpha'], r['beta_over_H'], '%.3g'%r['f_sw'], {k: '%.2g'%v for k,v in r['snr'].items()}, r['LVK_excluded'])
print({k: v for k, v in list(thr.items())[::4]})
