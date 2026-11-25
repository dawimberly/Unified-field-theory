"""Stage 51: forward/backward fringe visibility in the public Grundmann 2020 H2 data (Zenodo Fig3_data.txt).
Check of Kirpitch 2026 (Zenodo preprint): "the visibility differs between forward and backward light directions (0.579 against 0.508, 4.2 standard deviations)".
Halves: cos(beta) > 0 vs < 0 (beta = photon direction vs molecular axis; the axis sign is a labelling convention: first-detected -> second-detected proton, CITED).
Visibility: fit of the central (zeroth-order) fringe, I(x) = N (1 + s x) (1 + V cos(k (x - x0))), |cos(alpha)| <= W, Poisson likelihood.
Errors: parametric Poisson bootstrap of the histogram (equivalent to an event bootstrap for binned data).
"""
import json, csv, numpy as np
from scipy.optimize import minimize
OUT = '/workspace/uft-stage33/docs/hyperfield/'
d = np.genfromtxt(OUT+'data/stage31/grundmann2020_Fig3_data.txt', delimiter=',', skip_header=2)
ca = np.unique(d[:, 2]); cb = np.unique(d[:, 3]); H = np.zeros((len(ca), len(cb)))
for r in d: H[np.searchsorted(ca, r[2]), np.searchsorted(cb, r[3])] = r[4]
rng = np.random.default_rng(51)
def model(q, x):
    N, s, V, k, x0 = q
    return np.clip(N*(1+s*x)*(1+V*np.cos(k*(x-x0))), 1e-9, None)
def fitV(y, W=0.45):
    m = np.abs(ca) <= W + 1e-9; x = ca[m]; yy = y[m]
    q0 = [yy.mean(), 0.0, 0.5, 2*np.pi/0.6, 0.0]
    nll = lambda q: np.sum(model(q, x) - yy*np.log(model(q, x)))
    best = None
    for k0 in [2*np.pi/0.55, 2*np.pi/0.6, 2*np.pi/0.65]:
        r = minimize(nll, q0[:3]+[k0, 0.0], method='Nelder-Mead', options=dict(maxiter=20000, xatol=1e-7, fatol=1e-9))
        if best is None or r.fun < best.fun: best = r
    return best.x
def boot(y, W=0.45, n=400):
    q = fitV(y, W); m = np.abs(ca) <= W + 1e-9
    lam = y.copy()
    vs = [fitV(rng.poisson(lam).astype(float), W)[2] for _ in range(n)]
    return q, float(np.std(vs))
R = {}
fw = H[:, cb > 0].sum(1); bw = H[:, cb < 0].sum(1); al = H.sum(1)
R['counts'] = dict(forward=float(fw.sum()), backward=float(bw.sum()), cosb0=float(H[:, cb == 0].sum()))
R['counts']['asym'] = (R['counts']['forward']-R['counts']['backward'])/(R['counts']['forward']+R['counts']['backward'])
R['counts']['asym_sigma'] = R['counts']['asym']/np.sqrt(1/(R['counts']['forward']+R['counts']['backward']))
for W in [0.40, 0.45, 0.50]:
    res = {}
    for name, y in [('all', al), ('forward', fw), ('backward', bw)]:
        q, s = boot(y, W); res[name] = dict(V=float(q[2]), sV=s, x0=float(q[4]), period=float(2*np.pi/q[3]))
    dV = res['forward']['V']-res['backward']['V']; sdV = np.hypot(res['forward']['sV'], res['backward']['sV'])
    res['dV'] = dV; res['sdV'] = sdV; res['z'] = dV/sdV
    R['W%.2f' % W] = res
# raw max/min estimator (no fit), central max vs mean of the two adjacent minima
def rawV(y):
    c = np.argmin(np.abs(ca)); i0 = c-4+np.argmax(y[c-4:c+5])
    lo = y[(ca > -0.45) & (ca < -0.15)].min(); hi = y[(ca > 0.15) & (ca < 0.45)].min(); mn = 0.5*(lo+hi); mx = y[i0]
    return (mx-mn)/(mx+mn)
rb = [(rawV(rng.poisson(fw).astype(float)), rawV(rng.poisson(bw).astype(float)), rawV(rng.poisson(al).astype(float))) for _ in range(2000)]
rb = np.array(rb)
R['raw'] = dict(forward=rawV(fw), backward=rawV(bw), all=rawV(al), sF=float(rb[:, 0].std()), sB=float(rb[:, 1].std()), sA=float(rb[:, 2].std()))
R['raw']['dV'] = R['raw']['forward']-R['raw']['backward']; R['raw']['z'] = R['raw']['dV']/np.hypot(R['raw']['sF'], R['raw']['sB'])
# per |cos beta| pairs: V(+cb) - V(-cb), controls for the cos(beta) mix (fringe shift ~ cos beta smears halves differently if their |cos beta| mix differs)
pairs = []
for j in range(len(cb)):
    if cb[j] <= 0: continue
    jm = np.argmin(np.abs(cb + cb[j]))
    qf, sf = boot(H[:, j], 0.45, 200); qb, sb = boot(H[:, jm], 0.45, 200)
    pairs.append(dict(cos_beta=float(cb[j]), Vf=float(qf[2]), sVf=sf, Vb=float(qb[2]), sVb=sb, nf=float(H[:, j].sum()), nb=float(H[:, jm].sum())))
dv = np.array([p['Vf']-p['Vb'] for p in pairs]); sv = np.array([np.hypot(p['sVf'], p['sVb']) for p in pairs])
w = 1/sv**2; R['pairs'] = pairs
R['pairs_mean'] = dict(dV=float(np.sum(w*dv)/np.sum(w)), sdV=float(1/np.sqrt(np.sum(w))), chi2=float(np.sum(w*(dv-np.sum(w*dv)/np.sum(w))**2)), dof=len(dv)-1)
R['pairs_mean']['z'] = R['pairs_mean']['dV']/R['pairs_mean']['sdV']
# mirror (relabelling) symmetry test: forward half vs backward half mirrored in cos(alpha), shape only
fn = fw/fw.sum(); bn = bw[::-1]/bw.sum()
var = fw/fw.sum()**2 + bw[::-1]/bw.sum()**2
R['mirror_shape_chi2'] = dict(chi2=float(np.sum((fn-bn)**2/var)), dof=len(ca)-1)
R['kirpitch'] = dict(Vf=0.579, Vb=0.508, z=4.2, V_all=0.543, V_all_err=0.013)
R['symmetry_argument'] = ('The axis sign is the proton detection order (CITED). A rotation by pi about the photon direction k (an exact symmetry of an isotropic gas in a '
    'circularly polarized beam) maps every cos(beta)>0 configuration onto a cos(beta)<0 configuration with cos(alpha) -> -cos(alpha); only the spectrometer is not invariant. '
    'So any physical mechanism - standard nondipole physics, a do-si-do d1 term, or the preprint\'s coherence law - predicts V(+cos beta) = V(-cos beta). A difference measures instrument acceptance/resolution/background.')
json.dump(R, open(OUT+'results/stage51.json', 'w'), indent=1, default=float)
print(json.dumps({k: v for k, v in R.items() if k not in ('pairs',)}, indent=1, default=float))
for p in pairs: print(p)
S = [dict(key='grundmann_data', ref='Grundmann et al. 2020, Zenodo 10.5281/zenodo.3899920 Fig3_data.txt (cos(alpha) x cos(beta) histogram, 45 x 15 bins), copy in data/stage31/',
          url='https://zenodo.org/api/records/3899920/files/Fig3_data.txt/content', location='file columns "Bin#x, Bin#y, cos(alpha), cos(beta), count"', quote='"Variables = Bin#x, Bin#y, cos(alpha), cos(beta), count"', values='{}'),
     dict(key='grundmann_axis', ref='Grundmann et al., Science 370, 339 (2020), arXiv:2010.08298 (re-fetched)', url='https://arxiv.org/pdf/2010.08298', location='Materials and Methods',
          quote='"To cover the full range of -1 < cos(beta) < 1, we label the two detected protons by their sequence of impact on the detector. Thus, the molecule direction ... points from the proton detected first to the one detected second." ; "800 eV left-handed circularly polarized photons at right angle" ; "randomly oriented with respect to the light propagation"', values='{}'),
     dict(key='grundmann_nondipole', ref='same', url='https://arxiv.org/pdf/2010.08298', location='main text (red line) and Fig. S3 caption',
          quote='"the atomic nondipole effect tilts the electron angular distributions from each center slightly in forward direction with respect to the photon propagation" ; "The differences are most likely due to the integration over cos(beta) necessary for the display of the experimental results."', values='{}'),
     dict(key='kirpitch2026', ref='Kirpitch, Zenodo 10.5281/zenodo.23372925 (preprint v2), abstract as fetched in Stage 50 (PDF download blocked 403 this session)', url='https://zenodo.org/records/23370670', location='abstract',
          quote='"The fast-electron fringe visibility is V = 0.543 +- 0.013" ; "the visibility differs between forward and backward light directions (0.579 against 0.508, 4.2 standard deviations) while the product of visibility and delay is the same in both halves"', values='{"Vf":0.579,"Vb":0.508}')]
with open(OUT+'results/stage51_sources.csv', 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(S[0].keys())); wr.writeheader(); [wr.writerow(s) for s in S]
