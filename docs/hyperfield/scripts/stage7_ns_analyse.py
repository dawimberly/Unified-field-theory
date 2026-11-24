"""Re-analyse stored M-R-Lambda sequences with the Bardeen-Thorne-Meltzer stability rule, make figure."""
import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
out = json.load(open('results/stage7.json')); ns = out['ns']
def btm(curve, tol=3e-3):
    M, R, L = np.array(curve).T
    # extrema with tolerance
    ext = []; direction = 1; best = 0
    for i in range(1, len(M)):
        if direction == 1:
            if M[i] > M[best]: best = i
            elif M[i] < M[best]-tol: ext.append(('max', best)); direction = -1; best = i
        else:
            if M[i] < M[best]: best = i
            elif M[i] > M[best]+tol: ext.append(('min', best)); direction = 1; best = i
    unstable = 0; seg_start = 0; branches = []
    for kind, i in ext:
        if unstable == 0: branches.append((seg_start, i))
        j0, j1 = max(i-2, 0), min(i+2, len(R)-1)
        Rdec = R[j1] < R[j0]
        if kind == 'max': unstable += 1 if Rdec else -1
        else: unstable += -1 if Rdec else 1
        unstable = max(unstable, 0); seg_start = i
    if unstable == 0 and (not ext or ext[-1][0] == 'min'): branches.append((seg_start, len(M)-1))
    branches = [(a, b) for a, b in branches if M[b]-M[a] > tol]
    res = dict(M_max_stable=float(max(M[b] for a, b in branches)), n_stable_branches=len(branches), twin=len(branches) >= 2,
               branch_ranges=[(float(M[a]), float(M[b])) for a, b in branches])
    for a, b in branches:
        idx = np.arange(a, b+1)
        if M[idx].min() <= 1.4 <= M[idx].max():
            o = np.argsort(M[idx]); bb = idx[o]
            res.setdefault('R14_km', []).append(float(np.interp(1.4, M[bb], R[bb])))
            res.setdefault('Lam14', []).append(float(np.interp(1.4, M[bb], L[bb])))
    # twin at observable masses: second branch overlapping first in mass range above 1.0 Msun
    if len(branches) >= 2:
        (a1, b1), (a2, b2) = branches[0], branches[1]
        res['twin_mass_window'] = (float(max(M[a2], M[a1])), float(min(M[b1], M[b2])))
    return res
summ = {}
for k, v in list(ns.items()):
    if k.endswith('_curve'): continue
    A = btm(ns[k+'_curve']); A.update({x: v[x] for x in ('n_t_over_nsat', 'p_t_MeVfm3', 'cs2', 'de_over_et', 'seidov_unstable', 'binding_fraction_1.4') if x in v})
    r14 = A.get('R14_km', []); l14 = A.get('Lam14', [])
    A['passes_Mmax2.01_R11-13.5_Lam800'] = bool(A['M_max_stable'] >= 2.01 and any(11.0 <= x <= 13.5 for x in r14) and any(x <= 800 for x in l14))
    summ[k] = A
out['ns_summary'] = summ
out['ns_meta'] = dict(constraints='M_max>=2.01 (PSR J0348+0432; J0740+6620 2.08+-0.07), R1.4 in 11-13.5 km (NICER, broad), Lambda1.4<=800 (GW170817 low-spin 90%; refined 190+390-120) -- from memory, UNVERIFIED digits',
    hysteresis_note=('Stage-6 9/8 ratio is in the Landau control parameter r; mapping r <-> baryon chemical potential needs an unknown slope and offset: '
                     'NOT predictive for NS. A cold NS sits on the Maxwell construction; hysteresis would only appear as metastable (supercooled) hadronic cores before nucleation.'))
json.dump(out, open('results/stage7.json', 'w'), indent=1, default=float)
fig, ax = plt.subplots(figsize=(8, 6))
M, R, L = np.array(ns['SLy_curve']).T; ax.plot(R, M, 'k', lw=2, label='SLy (hadronic)')
for k in ['n2.0_cs1.00_de0.6', 'n2.0_cs1.00_de1.0', 'n3.0_cs1.00_de0.6', 'n1.5_cs1.00_de0.8', 'n4.0_cs1.00_de0.4']:
    if k+'_curve' in ns:
        M, R, L = np.array(ns[k+'_curve']).T; ax.plot(R, M, label=k+(' TWIN' if summ[k]['twin'] else ''))
ax.axhspan(2.01, 2.15, color='gray', alpha=.2); ax.axvspan(11, 13.5, color='b', alpha=.05)
ax.set_xlim(8, 16); ax.set_ylim(0.5, 3.2); ax.set_xlabel('R [km]'); ax.set_ylabel('M [Msun]'); ax.legend(fontsize=7)
ax.set_title('Stage 7: first-order transition (CSS) on SLy: twin stars need chosen (p_t, Delta eps)')
fig.tight_layout(); fig.savefig('figures/stage7_ns_MR.png', dpi=110)
for k, A in summ.items(): print(k, A.get('M_max_stable'), A['n_stable_branches'], A['twin'], A.get('twin_mass_window'), A.get('R14_km'), A.get('Lam14'), A['passes_Mmax2.01_R11-13.5_Lam800'])
