import json, numpy as np, importlib.util
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import minimize
OUT = '/workspace/uft-stage33/docs/hyperfield/'
R = json.load(open(OUT+'results/stage51.json'))
d = np.genfromtxt(OUT+'data/stage31/grundmann2020_Fig3_data.txt', delimiter=',', skip_header=2)
ca = np.unique(d[:, 2]); cb = np.unique(d[:, 3]); H = np.zeros((len(ca), len(cb)))
for r in d: H[np.searchsorted(ca, r[2]), np.searchsorted(cb, r[3])] = r[4]
fw = H[:, cb > 0].sum(1); bw = H[:, cb < 0].sum(1)
def model(q, x):
    N, s, V, k, x0 = q; return np.clip(N*(1+s*x)*(1+V*np.cos(k*(x-x0))), 1e-9, None)
def fit(y, W=0.45):
    m = np.abs(ca) <= W+1e-9; x = ca[m]; yy = y[m]
    nll = lambda q: np.sum(model(q, x)-yy*np.log(model(q, x)))
    return min((minimize(nll, [yy.mean(), 0, 0.5, k0, 0], method='Nelder-Mead', options=dict(maxiter=20000)) for k0 in [10.5, 11.4]), key=lambda r: r.fun).x
fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
a = ax[0]
xx = np.linspace(-0.45, 0.45, 200)
for y, lab, col, mir in [(fw, 'forward half, cos(beta) > 0', 'tab:blue', False), (bw, 'backward half, cos(beta) < 0, mirrored in cos(alpha)', 'tab:red', True)]:
    yy = y[::-1] if mir else y; q = fit(yy); n = yy.sum()
    a.errorbar(ca, yy/n, np.sqrt(yy)/n, fmt='o', ms=3, color=col, label=lab)
    a.plot(xx, model(q, xx)/n, '-', color=col, label='fit: V = %.3f' % q[2])
a.axvspan(-0.45, 0.45, color='grey', alpha=0.08); a.set_xlabel('cos(alpha) (fast electron vs molecular axis)'); a.set_ylabel('fraction of half per bin')
a.set_title('Grundmann 2020 public histogram: central fringe by half'); a.legend(fontsize=7, loc='upper center'); a.grid(alpha=0.3)
b = ax[1]
P = R['pairs']; c = [p['cos_beta'] for p in P]
b.errorbar(np.array(c)-0.01, [p['Vf'] for p in P], [p['sVf'] for p in P], fmt='o', color='tab:blue', label='V at +|cos beta|')
b.errorbar(np.array(c)+0.01, [p['Vb'] for p in P], [p['sVb'] for p in P], fmt='s', color='tab:red', label='V at -|cos beta|')
b.axhline(R['kirpitch']['Vf'], color='tab:blue', ls=':', label='Kirpitch forward 0.579'); b.axhline(R['kirpitch']['Vb'], color='tab:red', ls=':', label='Kirpitch backward 0.508')
b.set_xlabel('|cos(beta)|'); b.set_ylabel('central-fringe visibility V (1 sigma bootstrap)'); b.set_ylim(0.3, 0.75); b.legend(fontsize=7, ncol=2); b.grid(alpha=0.3)
b.set_title('pairwise V(+cos beta) vs V(-cos beta): dV = %.3f +- %.3f' % (R['pairs_mean']['dV'], R['pairs_mean']['sdV']))
e = ax[2]
labs = ['fit W=0.40', 'fit W=0.45', 'fit W=0.50', 'raw max/min', 'pairwise', 'Kirpitch\n(quoted)']
vals = [R['W0.40']['dV'], R['W0.45']['dV'], R['W0.50']['dV'], R['raw']['dV'], R['pairs_mean']['dV'], 0.071]
errs = [R['W0.40']['sdV'], R['W0.45']['sdV'], R['W0.50']['sdV'], np.hypot(R['raw']['sF'], R['raw']['sB']), R['pairs_mean']['sdV'], 0.071/4.2]
e.errorbar(range(6), vals, errs, fmt='o', color='k', capsize=4); e.axhline(0, color='tab:green', lw=2, label='every physical model (rotation symmetry about k): dV = 0')
e.set_xticks(range(6)); e.set_xticklabels(labs, fontsize=8); e.set_ylabel('V_forward - V_backward'); e.set_ylim(-0.03, 0.14); e.grid(alpha=0.3)
e.text(0.05, 0.125, 'event-count asymmetry of the halves: %.1f%% (%.0f sigma)\nmirror-shape chi2 = %.0f / %d  -> instrument asymmetry' % (100*R['counts']['asym'], R['counts']['asym_sigma'], R['mirror_shape_chi2']['chi2'], R['mirror_shape_chi2']['dof']), fontsize=8)
e.legend(fontsize=7, loc='lower right'); e.set_title('forward-backward visibility difference')
fig.suptitle('Stage 51 - Kirpitch forward/backward visibility difference REPRODUCES (dV = %.3f +- %.3f, %.1f sigma) but symmetry says it must be instrumental, not physics' % (R['W0.45']['dV'], R['W0.45']['sdV'], R['W0.45']['z']), fontsize=11)
fig.tight_layout(); fig.savefig(OUT+'figures/stage51_fb_visibility.png', dpi=110)
