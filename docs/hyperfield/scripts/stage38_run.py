#!/usr/bin/env python3
"""Stage 38 - other do-si-do elements with public data not used in Stages 31-34 (at most two).
Picked: element 1 (angle), the even d2*P2(cos beta) release-delay term, in a second, independent system:
angle-resolved RABBITT on the helium atom (Heuser et al. 2015/2016), vs Stage 31's H2 molecule (Grundmann 2020).
Others tried and dropped (reasons in results/stage38.json). Every input fetched this session (results/stage38_sources.csv).
"""
import json, os, csv
import numpy as np
from scipy.integrate import quad
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
S = {}
def src(key, ref, url, loc, quote, **vals): S[key] = dict(ref=ref, url=url, location=loc, quote=quote, values=vals)
src('heuser15', 'Heuser, Jimenez Galan, Cirelli, ..., Martin & Keller, "Time delay anisotropy in photoelectron emission from the isotropic ground state of helium", PRA 94, 063409 (2016), arXiv:1503.08966',
    'https://arxiv.org/pdf/1503.08966', 'Results (Fig. 2-3 text)',
    '"Panel (c) comprises electrons emitted within a cone of emission between 60 and 65" ; "a delay at 15 is understood as the delay between electrons emitted at angles between 0 and 30 (reference)" ; '
    '"shown with error bars in Fig. 3 for four consecutive sidebands, SB 18-SB 24. For all sidebands, the measurements deviate significantly from zero for angles larger than 50" ; '
    '"The results of the ab initio calculations are in quantitative agreement with the measurement" ; "for the considered angular range, the discrepancy between experimental data and the theory curves is not larger than 10-15 as"',
    ref_cone=(0.0, 30.0), last_bin=(60.0, 65.0), max_disc_as=15.0)
src('ivanov16', 'Ivanov, Dahlstrom, Lindroth & Kheifets, "On the angular dependence of the photoemission time delay in helium", arXiv:1605.04539', 'https://arxiv.org/pdf/1605.04539', 'abstract; Sec. 3',
    '"In the case of helium with only one 1s -> Ep photoemission channel, the Wigner time delay tau_W does not depend on the photoelectron detection angle" ; "In the experimentally accessible angular range of 0 to 65, where the RABBITT signal is sufficiently strong" ; '
    '"the experimental data (Heuser et al 2015) is insufficient to extract these parameters"')
src('stage31', 'Stage 31 (this repo): Grundmann 2020 H2 reanalysis', 'docs/hyperfield/HYPERFIELD_RESULTS.md', 'Stage 31, element 1',
    '"For release timing, the delay is t = R cos(beta)/c + d1 cos(beta) + d2 P2(cos(beta))" ; "The cos^2 (P2) term: |d2| < 102 zs"', d2_H2_zs=102.0)
h = S['heuser15']['values']
P2 = lambda c: 0.5 * (3 * c * c - 1)
w = lambda th: np.cos(th) ** 2 * np.sin(th)          # ASSUMED: dipole (1s -> p) angular distribution x solid angle
def avgP2(a, b):
    a, b = np.radians(a), np.radians(b)
    return quad(lambda t: P2(np.cos(t)) * w(t), a, b)[0] / quad(w, a, b)[0]
p_ref, p_last = avgP2(*h['ref_cone']), avgP2(*h['last_bin'])
lever = p_last - p_ref
d2_max_as = h['max_disc_as'] / abs(lever)             # ASSUMED: theory exact; the quoted max discrepancy is an envelope on any extra term
res = dict(P2_ref=p_ref, P2_60_65=p_last, lever=lever, d2_max_as=d2_max_as, d2_max_zs=d2_max_as * 1e3, ratio_to_H2=d2_max_as * 1e3 / S['stage31']['values']['d2_H2_zs'])
dropped = [
    dict(item='Element 4 spread / element 6 release time: non-exponential nuclear decay (Norman et al. PRL 60, 2246 (1988), 56Mn over 45 half-lives)', why='only the qualitative statement "measured up to 45 half lives ... none were found" could be fetched (Aston arXiv:1204.5953); no numerical limit on the non-exponential amplitude -> dropped.'),
    dict(item='Element 1 d1 (odd, cos beta) term in helium', why='the fetched text does not say whether the two emission hemispheres were separated; only the even P2 term is bounded.'),
    dict(item='Elements 2 (gate threshold) and 7 (hysteresis)', why='no new public dataset beyond Stage 31 (single-photon microwave counting; FIRAS mu) was found this stage.'),
]
with open(os.path.join(ROOT, 'results', 'stage38_sources.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh); wr.writerow(['key', 'ref', 'url', 'location', 'quote', 'values'])
    for kk, v in S.items(): wr.writerow([kk, v['ref'], v['url'], v['location'], v['quote'], json.dumps(v['values'])])
json.dump(dict(stage=38, results=res, dropped=dropped, sources=S), open(os.path.join(ROOT, 'results', 'stage38.json'), 'w'), indent=1, default=float)
print(json.dumps(res, indent=1))
