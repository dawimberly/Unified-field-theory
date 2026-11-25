import json, os, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
H = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
d = json.load(open(os.path.join(H, 'results', 'stage53.json')))
fig, ax = plt.subplots(figsize=(8, 5))
for e, c, lab in ((d['envelope'], 'C0', 'LPF Δg L4, 17 segments: tightest 95% band bound (valid)'),
                  (d['indicative_above_cut'], 'C7', 'above common −3 dB cut (indicative only)')):
    f = [x for r in e for x in (r['f_lo'], r['f_hi'])]; l = [x for r in e for x in (r['lam_band_95'],) * 2]
    ax.plot(f, l, c=c, lw=2, label=lab)
fl = [x for r in d['envelope'] for x in (r['f_lo'], r['f_hi'])]; ll = [x for r in d['envelope'] for x in (r['lam_line_95'],) * 2]
ax.plot(fl, ll, c='C1', lw=1, ls='--', label='narrow line (single 0.5 mHz bin), 95%')
ax.hlines(2.96e-8, 0.7e-3, 20e-3, colors='k', label='Helou 2017, 0.7–20 mHz (CITED)')
ax.axvspan(0.19, 9.3, color='C3', alpha=0.08, label='empty window 0.19–9.3 Hz (Stage 36)')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(5e-4, 10); ax.set_xlabel('frequency (Hz)')
ax.set_ylabel('λ (s⁻¹), CSL-equivalent, r_C = 100 nm'); ax.set_title('Stage 53: LISA Pathfinder Δg above 30 mHz → λ(f)')
ax.legend(fontsize=7, loc='upper left'); ax.grid(alpha=0.3, which='both')
fig.tight_layout(); os.makedirs(os.path.join(H, 'figures'), exist_ok=True)
fig.savefig(os.path.join(H, 'figures', 'stage53_lpf_lambda_vs_f.png'), dpi=130)
