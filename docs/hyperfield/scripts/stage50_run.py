"""Stage 50: angle-dependent release time d1 (per unit cos beta) - search for independent data and combination.
Parametrization (Stage 31/32): birth delay t(beta) = R cos(beta)/c + nondipole(authors' model) + d1 cos(beta) + d2 P2(cos beta).
Result: no independent dataset could be converted; a sensitivity floor for the best candidate (He non-dipole RABBITT) is DERIVED.
"""
import json, csv, numpy as np
OUT = '/workspace/uft-stage33/docs/hyperfield/'
R = {}
# --- Grundmann 2020 (Stage 31/32), residual over authors' red line (travel + nondipole) ---
R['grundmann'] = dict(d1=235.0, sig=128.0, d1_95=485.0, src='Stage 31 (CITED/DIGITIZED Grundmann 2020 Fig. 3C); Stage 32 red line 271.7 zs')
# --- Kirpitch 2026 re-analysis of the SAME Grundmann data (not independent) ---
slope = 391.0; s_slope = float(np.sqrt(39**2 + 52**2 + 50**2))   # CITED components (stat, method, zero point), added in quadrature (DERIVED)
red = 271.7                                                       # authors' travel + nondipole model, Stage 32 (DIGITIZED)
R['kirpitch_same_data'] = dict(slope=slope, sig=s_slope, d1=slope-red, d1_sig=s_slope, independent=False,
    note='same events as Grundmann 2020; unrefereed preprint; shown as a cross-check only, never combined')
# --- Liang 2024 He non-dipole RABBITT: sensitivity floor (DERIVED, optimistic) ---
# a d1 cos(beta_k) term gives a forward/backward delay difference 2 d1 cos(beta_k); detection pair theta = 70/110 deg (CITED),
# ASSUMED mirror images about the plane normal to k -> |cos beta_k| = cos 70 deg.
cb = float(np.cos(np.radians(70.0))); sig_fb_as = 10.0   # CITED +-10 as, used as an optimistic forward/backward-delay precision (ASSUMED)
sig_d1 = sig_fb_as*1e3/(2*cb)                            # zs
R['liang2024_floor'] = dict(cos_beta=cb, sig_FB_as=sig_fb_as, sig_d1_zs=sig_d1, d1_95_floor_zs=1.96*sig_d1, vs_grundmann=1.96*sig_d1/485.0,
    note='full text and data (ETH Research Collection) not reachable (Nature paywall; ETH 429/500); the paper reports a dipole-quadrupole transition delay, which cannot be mapped to d1 without its formalism. Best case only: sigma(d1) >= %.0f zs.' % sig_d1)
# --- combination (only independent sets): Grundmann + Liang floor (central 0 ASSUMED) ---
w = np.array([1/128.0**2, 1/sig_d1**2]); m = np.array([235.0, 0.0])
comb = float(np.sum(w*m)/np.sum(w)); sc = float(1/np.sqrt(np.sum(w)))
R['combined_best_case'] = dict(d1=comb, sig=sc, d1_95=abs(comb)+1.96*sc, weight_liang=float(w[1]/w.sum()),
    note='even granting the He floor with central value 0, the combination is unchanged: |d1| < 485 zs stands')
R['dropped'] = [
  'Rezvan et al. 2024 J. Phys. B 57 145101 (N2 1s, 880 eV): forward/backward MFPAD intensity asymmetry, no timing; large statistical uncertainties at beta = 45 deg',
  'CO 1s at 905 eV, PRL 129, 253201 (2022): intensity asymmetry, no timing',
  'Heuser 2016 He RABBITT (Stage 38): hemispheres along k not separated, so d1 not accessible',
  'Molecular-frame RABBITT (e.g. Nat. Commun. 2021): angle to polarization / molecular axis, not to photon direction; dipole symmetric along k',
  'Chiral forward/backward delays (Nature 2025): dominated by PECD (tens of as); standard-theory term far larger than d1']
json.dump(R, open(OUT+'results/stage50.json', 'w'), indent=1, default=float)
print(json.dumps(R, indent=1, default=float))
S = [dict(key='stage31_grundmann', ref='This repo, Stage 31/32 (Grundmann et al., Science 370, 339 (2020))', url='docs/hyperfield/HYPERFIELD_RESULTS.md', location='Stage 31 sec. 31.2, Stage 32',
          quote='"|d1| < 485 zs per unit cos beta" ; "residual +235 +- 128 zs (1.8 sigma)" ; "red (nondipole + travel) 271.7"', values='{"d1":235,"sig":128}'),
     dict(key='kirpitch2026', ref='Kirpitch, R., Zenodo 10.5281/zenodo.23372925 (preprint v2, 14 Oct 2026)', url='https://zenodo.org/records/23370670', location='abstract',
          quote='"We re-analyse the published molecular-frame data of the zeptosecond birth-delay experiment ... [Grundmann et al., Science 370, 339 (2020)]. Four estimators ... give a delay slope of 391 +- 39 (stat) +- 52 (method) +- 50 (zero point) zs, (1.59 +- 0.33) R/c, linear in cos beta"',
          values='{"slope":391,"stat":39,"method":52,"zero":50}'),
     dict(key='liang2024', ref='Liang, Han, Liao, ... Woerner, Nature Photonics 18, 311 (2024), doi 10.1038/s41566-023-01349-z (preview page)', url='https://www.nature.com/articles/s41566-023-01349-z',
          location='abstract; Extended Data Fig. 2 caption; Data availability',
          quote='"we measure a time delay of 15 +- 10 as between the electric-dipole and electric-quadrupole transitions by resolving the asymmetry of the photoelectron forward-backward yields" ; "Measured photoelectron spectra at theta = 70 deg (blue solid line) and theta = 110 deg (purple line)" ; "The data ... are available via the ETH Zurich Research Collection at https://doi.org/10.3929/ethz-b-000638773"',
          values='{"sig_as":10,"theta":[70,110]}'),
     dict(key='rezvan2024', ref='Rezvan et al., J. Phys. B 57, 145101 (2024)', url='https://doi.org/10.3204/pubdb-2024-05394', location='abstract; results text (search snippet)',
          quote='"Including nondipole (retardation) effects in the calculations results in a simultaneous increase and decrease of the scattering peaks towards the nitrogen atoms pointing in the forward and backward directions" ; "owing to the large statistical uncertainties, the distribution in figure 2(d), measured for beta = 45 deg, cannot be considered as an unambiguous proof of the effect"', values='{}'),
     dict(key='mst2024', ref='arXiv:2402.02435v3 (multiple-scattering theory of core photoemission delays, N2)', url='https://arxiv.org/pdf/2402.02435v3', location='abstract; refs [45],[46]',
          quote='"the incoherent sums of the photoemission time delays show a slight variation of a few hundred zeptoseconds" ; "when considering the contributions from the gerade and ungerade delocalized states in incoherent sums, the two-center interference terms cancel each other"', values='{}')]
with open(OUT+'results/stage50_sources.csv', 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(S[0].keys())); wr.writeheader(); [wr.writerow(s) for s in S]
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
a = ax[0]
rows = [("Grundmann 2020 H2\n(Stage 31, authors' model subtracted)", 235, 128, 'tab:blue', 'o', False),
        ('Kirpitch 2026 re-analysis\n(SAME data - not combined)', R['kirpitch_same_data']['d1'], s_slope, 'grey', 'o', True),
        ('combined, independent only\n(best case incl. He floor)', comb, sc, 'k', 'D', False)]
for i, (lab, v, s, col, mk, hollow) in enumerate(rows):
    a.errorbar(v, i, xerr=1.96*s, fmt=mk, color=col, capsize=4, mfc='white' if hollow else col)
a.annotate('', xy=(1500, 3), xytext=(-1500, 3), arrowprops=dict(arrowstyle='<->', color='tab:red'))
a.plot(0, 3, 's', color='tab:red', mfc='white')
a.set_yticks(range(4)); a.set_yticklabels([r[0] for r in rows] + ['Liang 2024 He non-dipole RABBITT\n(sensitivity floor, centre 0 ASSUMED)'], fontsize=8)
a.axvline(0, color='k', lw=0.8); a.axvline(247, color='tab:green', ls='--', lw=1); a.text(255, 2.6, 'R/c = 247 zs', color='tab:green', fontsize=8)
a.set_ylim(-0.5, 3.6); a.set_xlim(-1500, 1500); a.set_xlabel('d1 (zs per unit cos beta), 95% bars'); a.invert_yaxis(); a.grid(alpha=0.3)
a.text(-1450, 3.35, 'He floor 95%%: +-%.0f zs (off scale)' % (1.96*sig_d1), fontsize=7, color='tab:red')
a.set_title('angle-dependent release term d1 beyond travel + nondipole')
b = ax[1]
labs = ['Grundmann H2\n(Stage 31)', 'He non-dipole\nRABBITT floor', 'combined']
vals = [485, 1.96*sig_d1, abs(comb)+1.96*sc]
b.bar(range(3), vals, color=['tab:blue', 'tab:red', 'k']); b.set_yscale('log'); b.set_xticks(range(3)); b.set_xticklabels(labs, fontsize=8)
for i, v in enumerate(vals): b.text(i, v*1.15, '%.0f zs' % v, ha='center', fontsize=8)
b.axhline(247, color='tab:green', ls='--', lw=1, label='light travel R/c (247 zs)'); b.legend(fontsize=8)
b.set_ylabel('|d1| 95% upper limit (zs)'); b.set_ylim(100, 1e5); b.set_title('no independent dataset tightens d1')
fig.suptitle('Stage 50 - angle-dependent release time: |d1| < 485 zs unchanged (Grundmann only); best independent candidate (He RABBITT) >= %.0f zs' % (1.96*sig_d1), fontsize=11)
fig.tight_layout(); fig.savefig(OUT+'figures/stage50_d1_search.png', dpi=110)
