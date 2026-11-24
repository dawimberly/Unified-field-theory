import json, glob, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
A = json.load(open('results/stage6.json'))['partA']
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
ds = sorted(map(int, A))
for n, mk in [('Z^d', 'o'), ('D_d', 's'), ('D_d^+', '^'), ('BCC Z^d u (Z^d+h)', 'D')]:
    x = [d for d in ds if A[str(d)][n].get('lattice')]
    ax[0].plot(x, [A[str(d)][n]['delta_over_best_lattice'] for d in x], mk+'-', label=n)
    ax[1].plot(x, [A[str(d)][n]['kiss_over_best_lattice'] for d in x], mk+'-', label=n)
for a, t in zip(ax, ['center density / best known lattice', 'kissing / best known lattice']):
    a.axvline(4, ls=':', c='gray'); a.axvline(8, ls=':', c='gray'); a.axvline(9, ls='--', c='r'); a.set_xlabel('d'); a.set_title(t); a.set_yscale('log')
ax[0].legend(fontsize=8); fig.suptitle("Stage 6A: half-diagonal 'slide' lattices; only d=4 (BCC->D4) and d=8 (D8+->E8) hit the best lattice; d=9 (red) nothing")
fig.tight_layout(); fig.savefig('figures/stage6_slide_lattices.png', dpi=110)
T = {}
for f in glob.glob('results/stage6_thomson_*.json'):
    if 'stability' in f: continue
    T.update(json.load(open(f)))
if T:
    fig, ax = plt.subplots(figsize=(9, 4.5)); labs = sorted(T)
    ax.bar(range(len(labs)), [max(abs(T[k]['rel_diff']), 1e-16) for k in labs], color=['g' if T[k]['rel_diff'] < 1e-9 else ('r' if T[k]['rel_diff'] < 0 else 'orange') for k in labs])
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, rotation=60, fontsize=7); ax.set_yscale('log')
    ax.set_ylabel('abs(E_best_random - E_lattice)/E_lattice (red: random beats lattice)'); ax.set_title('Thomson-type minimization: does random start reach the lattice config? (green=yes)')
    fig.tight_layout(); fig.savefig('figures/stage6_thomson.png', dpi=110)
