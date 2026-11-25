#!/usr/bin/env python3
"""Stage 53 - LISA Pathfinder Delta-g above 30 mHz -> frequency-resolved CSL-equivalent lambda(f) per nucleon, r_C = 100 nm.

DATA (fetched this session, anonymous HTTPS, no login needed):
  LPF Legacy Archive (ESAC) metadata API
    https://lpf.esac.esa.int/lpfsa-sl/metadata-action?RESOURCE_CLASS=ANALYSIS_OBJECT&SELECTED_FIELDS=...&RETURN_TYPE=JSON&QUERY=ANALYSIS_OBJECT.ORIGIN='D'
  -> 337 Delta-g analysis objects (levels L0..L4). Level 4 = "Delta_g_noise" (L3 cleaned of glitches, noise-only runs), 17 segments.
  download: https://lpf.esac.esa.int/lpfsa-sl/data-action?ANALYSIS_OBJECT.ANALYSIS_OBJECT_OID=<oid>  (LTPDA XML, fs = 10 Hz)
  Raw XML (~2 GB) is NOT committed; oid + sha256 per file are in results/stage53.json.

METHOD (Stage 33, unchanged): white CSL force on each cube, thick-slab limit of Eq. (13) of Carlesso-Ferialdi-Bassi 2018:
  S_F = 8 pi hbar^2 lambda rho^2 A r_C^2 / m0^2 (two-sided, per TM); one-sided Delta-g PSD = 2 x 2(TMs) x S_F / m^2.
  All measured Delta-g power is attributed to the effect (conservative). Check: 5.2 fm -> 2.986e-8 vs Helou 2.96e-8.
FREQUENCY VALIDITY (ASSUMED rule, fixed before looking at bounds): above 30 mHz Delta-g is readout-dominated, i.e.
  displacement-equivalent noise ASD/(2 pi f)^2 ~ flat. A band is kept while that quantity stays >= 1/sqrt(2) (-3 dB) of its
  median over 0.1-0.3 Hz; the first band below that marks the on-board low-pass (not modelled), and it and all higher bands are dropped.
"""
import re, sys, os, glob, json, csv, hashlib, numpy as np
from scipy.signal import welch
from scipy.stats import chi2

RAW = sys.argv[1] if len(sys.argv) > 1 else '/workspace/s53/raw'
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
pi = np.pi
HBAR = 1.054571817e-34      # CITED NIST CODATA (exact), fetched this session
M0 = 1.66053906892e-27      # CITED NIST CODATA atomic mass constant; m0 = 1 u ASSUMED (as Stage 33)
RC = 1e-7                   # ASSUMED convention r_C = 100 nm
L_TM, M_TM = 0.046, 1.928   # CITED Carlesso-Ferialdi-Bassi 2018 Sec. IV "length 4.6 cm ... mass 1.928 kg"
HELOU_LAM, HELOU_SA = 2.96e-8, 5.2e-15   # CITED Helou et al. arXiv:1606.03637 abstract
READOUT = 34.8e-15          # CITED Armano 2018 PRL 120 061101: S_x^1/2 = (34.8 +- 0.3) fm/rtHz, effective Delta-g ASD S_x^1/2 (2 pi f)^2
FLOOR18 = 1.74e-15          # CITED Armano 2018: (1.74 +- 0.05) fm s^-2/rtHz above 2 mHz
OIDS = {}                   # filled from metadata list

rho, A = M_TM / L_TM**3, L_TM**2
S_slab_per_lam = 8 * pi * HBAR**2 * rho**2 * A * RC**2 / M0**2

def lam_from_psd(S1):
    return S1 * M_TM**2 / (2 * 2 * S_slab_per_lam)

def load_ao(path):
    t = open(path).read()
    body = t[t.find('</historyRoot>'):]
    name = re.search(r'<name shape[^>]*>(.*?)</name>', body).group(1)
    ts = re.search(r'<tsdata.*?<realData shape="1x(\d+)"[^>]*>(.*?)</realData>', body, re.S)
    chunks = re.findall(r'<doubleVector type="double">(.*?)</doubleVector>', ts.group(2), re.S)   # 50000-sample chunks
    y = np.array(' '.join(chunks).split(), dtype=float)
    assert len(y) == int(ts.group(1))
    fs = float(re.findall(r'<fs shape="1x1" type="double">([^<]+)</fs>', body)[-1])
    t0 = float(re.findall(r'utc_epoch_milli="([^"]+)"', body)[0]) / 1000
    return dict(name=name, y=y, fs=fs, t0=t0)

def sha(path):
    hh = hashlib.sha256()
    with open(path, 'rb') as f:
        for ch in iter(lambda: f.read(1 << 22), b''): hh.update(ch)
    return hh.hexdigest()

def psd(y, fs, Tseg):
    n = int(Tseg * fs)
    f, P = welch(y - y.mean(), fs=fs, window='blackmanharris', nperseg=n, noverlap=n // 2, detrend='linear', scaling='density')
    return f, P, int((len(y) - n) // (n // 2) + 1)

EDGES = 10 ** np.arange(np.log10(0.02), np.log10(5.0) + 1e-9, 0.05)   # 1/20-decade bands (ASSUMED binning)

def bands(f, P, K):
    rows = []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        sel = (f >= lo) & (f < hi); nb = int(sel.sum())
        if nb < 2: continue
        Pm = P[sel].mean(); fc = np.sqrt(lo * hi)
        dof = 2 * K * nb / 1.9          # ~2K dof per bin, adjacent BH bins partly correlated (ASSUMED /1.9)
        Pup = Pm * dof / chi2.ppf(0.05, dof)
        Pline = P[sel].max() * (2 * K) / chi2.ppf(0.05, 2 * K)
        xeq = np.sqrt(Pm) / (2 * pi * fc)**2
        rows.append(dict(f_lo=lo, f_hi=hi, f_c=fc, nbins=nb, asd_fm=np.sqrt(Pm) * 1e15, xeq_fm=xeq * 1e15,
                         readout_model_fm=np.hypot(READOUT * (2 * pi * fc)**2, FLOOR18) * 1e15,
                         lam_band=lam_from_psd(Pm), lam_band_95=lam_from_psd(Pup), lam_line_95=lam_from_psd(Pline)))
    ref = np.median([r['xeq_fm'] for r in rows if 0.1 <= r['f_c'] <= 0.3])
    ok = True
    for r in rows:
        if r['f_lo'] >= 0.0299 and r['xeq_fm'] < ref / np.sqrt(2): ok = False
        r['valid'] = ok and r['f_lo'] >= 0.0299
        r['valid_or_lowf'] = ok
    return rows, ref

if __name__ == '__main__':
    calib = lam_from_psd(HELOU_SA**2)
    print(f'formula check: lambda(5.2 fm) = {calib:.4e} vs Helou {HELOU_LAM:.3e} ratio {calib/HELOU_LAM:.4f}')
    catf = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'stage53', 'lpfsa_deltag_catalogue.json')
    meta = {r['ANALYSIS_OBJECT__NAME']: r for r in json.load(open(catf))['data']}
    res = dict(calib=dict(lambda_5p2fm=calib, helou=HELOU_LAM, ratio=calib / HELOU_LAM), segments={})
    allrows = {}
    for p in sorted(glob.glob(os.path.join(RAW, 'Delta-g-x-L4-*.xml'))):
        key = os.path.basename(p)[:-4]
        try:
            ao = load_ao(p)
        except Exception as e:
            print('skip', key, e); continue
        y = ao['y'][np.isfinite(ao['y'])]; fs = ao['fs']; dur = len(y) / fs
        if dur < 4 * 2000: print('skip short', key, dur); continue
        fH, PH, KH = psd(y, fs, 2000.0)
        fL, PL, KL = psd(y, fs, min(20000.0, dur / 4))
        low = (fL >= 2e-3) & (fL <= 20e-3)
        asd_low = np.sqrt(PL[low].mean())
        rows, ref = bands(fH, PH, KH)
        allrows[key] = rows
        vmax = max([r['f_hi'] for r in rows if r['valid']], default=None)
        m = meta.get(key, {})
        res['segments'][key] = dict(oid=m.get('ANALYSIS_OBJECT__ANALYSIS_OBJECT_OID'), main_doy=m.get('ANALYSIS_OBJECT__MAIN_DOY'),
                                    ao_name=ao['name'], n=int(len(y)), fs=fs, duration_d=dur / 86400, t0_utc_epoch_s=ao['t0'],
                                    sha256=sha(p), asd_2_20mHz_fm=asd_low * 1e15, lam_2_20mHz=float(lam_from_psd(asd_low**2)),
                                    xeq_ref_0p1_0p3Hz_fm=ref, f_valid_max=vmax, K_2000s=KH,
                                    ratio_xeq_to_34p8=ref / (READOUT * 1e15))
        print(f"{key:18s} {ao['name']} {dur/86400:6.2f} d  ASD(2-20mHz) {asd_low*1e15:6.3f} fm  xeq(0.1-0.3Hz) {ref:6.1f} fm  valid<= {vmax}")
        if key == 'Delta-g-x-L4-66':
            with open(os.path.join(OUT, 'stage53_lpf_asd_L4-66.csv'), 'w') as fh:
                fh.write('f_Hz,asd_m_s2_rtHz  # Welch, 2000 s Blackman-Harris segments, 50% overlap\n')
                [fh.write(f'{a:.6e},{np.sqrt(b):.6e}\n') for a, b in zip(fH[1:], PH[1:]) if a >= 0.01]
            res['L4_66_peaks'] = []
            from scipy.signal import medfilt
            r_ = PH / medfilt(PH, 101); s_ = (fH > 0.03) & (fH < 2.5)
            for i in np.argsort(r_[s_])[::-1][:8]: res['L4_66_peaks'].append([float(fH[s_][i]), float(r_[s_][i])])
    # envelope: per band, the tightest valid 95% bound over segments
    os.makedirs(OUT, exist_ok=True)
    keys = list(allrows); nb = len(allrows[keys[0]])
    # common cutoff (ASSUMED rule): the on-board low-pass is the same instrument for every segment, so the robust
    # upper limit is the LOWEST per-segment -3 dB point; bands above it are kept only as 'indicative'
    FCUT = min(res['segments'][k]['f_valid_max'] for k in keys)
    res['f_cut_common'] = FCUT
    env, indicative = [], []
    for j in range(nb):
        r0 = allrows[keys[0]][j]
        tgt = env if r0['f_hi'] <= FCUT * (1 + 1e-9) else indicative
        cand = [(allrows[k][j]['lam_band_95'], k) for k in keys if allrows[k][j]['valid']]
        candl = [(allrows[k][j]['lam_line_95'], k) for k in keys if allrows[k][j]['valid']]
        if not cand: continue
        b, kb = min(cand); bl, kl = min(candl)
        tgt.append(dict(f_lo=r0['f_lo'], f_hi=r0['f_hi'], f_c=r0['f_c'], lam_band_95=b, seg_band=kb, lam_line_95=bl, seg_line=kl,
                        n_valid_segments=len(cand), lam_band_95_L4_66=allrows.get('Delta-g-x-L4-66', [{}]*nb)[j].get('lam_band_95')))
    with open(os.path.join(OUT, 'stage53_lpf_lambda_vs_f.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(env[0])); w.writeheader()
        [w.writerow({k: (f'{v:.4e}' if isinstance(v, float) else v) for k, v in r.items()}) for r in env]
    with open(os.path.join(OUT, 'stage53_lpf_bands_per_segment.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['segment'] + list(allrows[keys[0]][0]))
        for k in keys:
            for r in allrows[k]: w.writerow([k] + [f'{v:.4e}' if isinstance(v, float) else v for v in r.values()])
    def at(f0): return min(env, key=lambda r: abs(np.log(r['f_c'] / f0)))
    res['envelope'] = env; res['indicative_above_cut'] = indicative
    with open(os.path.join(OUT, 'stage53_lpf_lambda_indicative_above_cut.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(indicative[0])); w.writeheader()
        [w.writerow({k: (f'{v:.4e}' if isinstance(v, float) else v) for k, v in r.items()}) for r in indicative]
    res['best'] = min(env, key=lambda r: r['lam_band_95'])
    res['at'] = {str(f0): at(f0) for f0 in (0.03, 0.1, 0.19, 0.3, 0.5, 1.0, 1.5, 2.0)}
    win = [r for r in env if r['f_hi'] > 0.19]
    res['window_0p19_9p3'] = dict(f_max=max(r['f_hi'] for r in env), n_bands=len(win),
                                  lam_min=min(r['lam_band_95'] for r in win) if win else None,
                                  lam_max=max(r['lam_band_95'] for r in win) if win else None)
    json.dump(res, open(os.path.join(OUT, 'stage53.json'), 'w'), indent=1, default=float)
    print('\nENVELOPE  f_c   lam_band95 (seg)   lam_line95   n_seg')
    print('common cut', FCUT)
    for r in env + indicative: print(('  ' if r in env else '~ ') + f"{r['f_c']:.4f} {r['lam_band_95']:.3e} {r['seg_band'][10:]:6s} {r['lam_line_95']:.3e} {r['n_valid_segments']}")
    print('best', res['best']['f_c'], res['best']['lam_band_95']); print('window', res['window_0p19_9p3'])
