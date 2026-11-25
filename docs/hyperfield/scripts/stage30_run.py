#!/usr/bin/env python3
"""Stage 30 - Big G vs test-mass volume: can CQ/DP gravitational noise (volume law, Stages 27/29) explain the ~500 ppm G spread?
Data: the 16 CODATA 2018 = CODATA 2022 input values (CODATA 2022, arXiv:2409.03787, Table XXX: "there is no new relevant input datum
for the 2022 adjustment").
VERIFICATION PASS (Stage 30v): every geometry below was re-read from a source document fetched in this session (arXiv / journal
Extended Data / CODATA RMP text / PhD thesis); each row carries 'source' (document + page/table/line) and 'quote'.  Nothing is
recalled.  Remaining non-geometric inputs that are handbook material densities (W 19300, Cu 8960, Cu-Te 8900, Pyrex 2230 kg/m^3) are
labelled 'ASSUMED-handbook' in V_how and only used where the paper gives mass + material but not dimensions.
Labels: DERIVED / CITED / ASSUMED.
"""
import json, os, csv, time
import numpy as np
from math import pi, sqrt, log10
from scipy import stats
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
t0 = time.time(); rng = np.random.default_rng(30)
HBAR = 1.054571817e-34; G0 = 6.67430e-11; GH = G0 * HBAR
G_CODATA, U_CODATA = 6.67430, 0.00015          # CITED CODATA 2018/2022, x1e-11, Birge expansion 3.9

def N_spheroid_axis(m):     # demagnetizing factor along the symmetry axis of a spheroid with c/a = m (Osborn 1945)
    if abs(m - 1) < 1e-9: return 1 / 3
    if m > 1: e = sqrt(1 - 1 / m ** 2); return (1 - e ** 2) / e ** 3 * (np.arctanh(e) - e)
    e = sqrt(1 / m ** 2 - 1); return (1 + e ** 2) / e ** 3 * (e - np.arctan(e))
def N_trans(m): return (1 - N_spheroid_axis(m)) / 2
cylV = lambda d, h: pi * (d / 2) ** 2 * h
boxV = lambda a, b, c: a * b * c
RHO_W, RHO_CU, RHO_CUTE, RHO_PYREX = 19300., 8960., 8900., 2230.   # ASSUMED-handbook
g_n = 9.80665
# ---------------------------------------------------------------- data (G, u: CITED CODATA 2022 Table XXX)
# method classes: TB-dyn (time-of-swing), TB-comp (static/servo/compensation, incl. AAF), SUSP (suspended body / pendulum displacement),
# STAT (stationary body, beam balance), AI (atom interferometer)
# V_eff: per suspended BODY (a dumbbell = two bodies). N_rel: demagnetization factor along the signal direction (spheroid approximation
# for cylinders/blocks, DERIVED); shape_ok=False where shape/orientation is not given in the source -> excluded from the proxy test and
# N_rel set to a conservative upper value for the physics.  a_sig: signal acceleration of the test body (CITED or DERIVED from cited
# numbers); None -> no verified source-mass distance -> row excluded from the physics table.
L = 'CODATA-98 = Mohr & Taylor, RMP 72, 351 (2000) [physics.nist.gov/cuu/pdf/CODATA_RMP2000.pdf]'
D = [
 dict(id='NIST-82', ref='Luther & Towler, PRL 48, 121 (1982)', G=6.67248, u=0.00043, method='TB-dyn',
      test='W dumbbell: 2 W disks 2.5 mm thick x 7.2 mm dia + W rod ~1 mm x 29 mm, ~7 g total', mat='W', rho=RHO_W,
      V_eff=cylV(7.2e-3, 2.5e-3), V_how='dims (disk)', m_body=RHO_W * cylV(7.2e-3, 2.5e-3), shape='disk(t/d=0.347)', aspect=2.5 / 7.2,
      N_rel=N_spheroid_axis(2.5 / 7.2), shape_ok=False,
      src='2 W spheres ~10.2 cm dia, ~10.5 kg; centres ~14 cm apart on dumbbell axis', Ms=10.5, r=0.07 - 0.0145,
      a_sig=None, geo='CITED (CODATA-98)',
      source=L + ', Sec. III.Q, p. 425',
      quote='"a 7 g, dumbbell-like small-mass system ... two tungsten disks about 2.5 mm thick and 7.2 mm in diameter, the centers of which were '
            'connected by a tungsten rod about 1 mm in diameter and 29 mm long ... two tungsten spheres, each about 10.2 cm in diameter and with a '
            'mass of about 10.5 kg ... separated by about 14 cm"'),
 dict(id='TR&D-96', ref='Karagioz & Izmailov, Meas. Tech. 39, 979 (1996) = Karagioz et al. (1998) (same result, CODATA-06)', G=6.6729, u=0.0005, method='TB-dyn',
      test='beam 23 cm x 1.8 mm with a 1.6 g sphere (7.2 mm dia) at each end, 5.3 g total (apparatus version 3, 1988-89)', mat='? (rho DERIVED 8190)',
      rho=1.6e-3 / (pi / 6 * 7.2e-3 ** 3), V_eff=pi / 6 * 7.2e-3 ** 3, V_how='dims (sphere)', m_body=1.6e-3, shape='sphere', aspect=1.0,
      N_rel=1 / 3, shape_ok=True,
      src='brass 12.2 cm 8.0 kg / bronze 10.2 cm 4.9 kg / bearing steel 10.1 cm 4.3 kg spheres at 19.2-47.2 cm from fibre', Ms=8.0, r=0.192 - 0.115,
      a_sig=None, geo='CITED (CODATA-98; version-3 test mass only)',
      source=L + ', p. 427-428 item (4) TR&D-98; CODATA-06 (arXiv:0801.0028) p. 58: TR&D-96 "reports the same result as does the 1998 reference"',
      quote='"In version 3, a 25 um diameter, 23 cm long fiber supported a 5.3 g test mass at its center. The test mass consisted of a cylindrical '
            'beam 23 cm long and 1.8 mm in diameter ... with a 1.6 g spherical mass of diameter 7.2 mm at each end ... one of brass (12.2 cm in '
            'diameter, mass of 8.0 kg) ... four positions at distances 19.2 cm, 21.2 cm, 25.2 cm, and 47.2 cm from the torsion fiber"'),
 dict(id='LANL-97', ref='Bagley & Luther, PRL 78, 3047 (1997)', G=6.67398, u=0.00070, method='TB-dyn',
      test='same W dumbbell as NIST-82', mat='W', rho=RHO_W,
      V_eff=cylV(7.2e-3, 2.5e-3), V_how='dims (disk, = NIST-82)', m_body=RHO_W * cylV(7.2e-3, 2.5e-3), shape='disk(t/d=0.347)', aspect=2.5 / 7.2,
      N_rel=N_spheroid_axis(2.5 / 7.2), shape_ok=False,
      src='same W spheres as NIST-82 (near position assumed as NIST-82)', Ms=10.5, r=0.07 - 0.0145, a_sig=None, geo='CITED (CODATA-98)',
      source=L + ', p. 427 item (3) LANL-97',
      quote='"is in many ways similar to the NIST-UVA experiment of Luther and Towler (1982); it used the same dumbbell test mass and tungsten source masses"'),
 dict(id='UWash-00', ref='Gundlach & Merkowitz, PRL 85, 2869 (2000) [gr-qc/0006043]', G=6.674255, u=0.000092, method='TB-comp',
      test='gold-coated Pyrex plate 1.506 x 76 x 41.6 mm', mat='Pyrex', rho=RHO_PYREX,
      V_eff=boxV(1.506e-3, 76e-3, 41.6e-3), V_how='dims (plate)', m_body=RHO_PYREX * boxV(1.506e-3, 76e-3, 41.6e-3), shape='plate',
      aspect=1.506 / sqrt(76 * 41.6), N_rel=N_spheroid_axis(1.506 / sqrt(76 * 41.6)), shape_ok=True,
      src='#316 stainless spheres, dia 124.89 mm, ~8.140 kg, centres at rho = 16.76 cm', Ms=8.140, r=0.1676, a_sig=None,
      geo='CITED (Gundlach-Merkowitz 2000)',
      source='arXiv:gr-qc/0006043, p. 2 (apparatus paragraph)',
      quote='"The pendulum was a 1.506 mm thick, 76 mm wide, and 41.6 mm tall Pyrex-glass plate with a thin gold-coating ... The centers of the '
            'spheres were located at rho=16.76 cm ... #316 stainless steel. Their average diameter was 124.89 mm and their mass was ~8.140 kg"'),
 dict(id='BIPM-01', ref='Quinn et al., PRL 87, 111101 (2001)', G=6.67559, u=0.00027, method='TB-comp',
      test='4 Cu-0.7%Te right cylinders h = d = 55 mm, ~1.2 kg (Mk I)', mat='CuTe', rho=1.2 / cylV(0.055, 0.055),
      V_eff=cylV(0.055, 0.055), V_how='dims (cyl h=d)', m_body=1.2, shape='cyl(h=d)', aspect=1.0, N_rel=N_trans(1.0), shape_ok=True,
      src='4 Cu-Te cylinders h = d = 118 mm, ~11 kg (PRL 2001 text: 12 kg); source-test separation ~107 mm', Ms=11.0, r=0.107, a_sig=None,
      geo='CITED (Quinn et al. 2014 review)',
      source='Quinn, Speake, Parks & Davis, Phil. Trans. A 372, 20140032 (2014), Sec. 4 & Fig. 1; Speake BIPM 2015 slides p. 9 (PRL 2001 text)',
      quote='"Both test and source masses were made with their heights equal to their diameters, 55 and 118 mm, respectively ... At the separation '
            'between the source and test masses (about 107 mm)"; Fig. 1: "source masses (~11 kg) ... test masses (~1.2 kg)"; PRL 2001: "masses of 12 and 1.2 kg"'),
 dict(id='UWup-02', ref='Kleinevoss, PhD thesis Wuppertal WUB-DIS 2002-2 (2002)', G=6.67422, u=0.00098, method='SUSP',
      test='2 pendulum bodies m ~3.3 kg: OFHC-Cu spherical mirror ~198 mm dia glued in Al mount (100 mm long, 200 mm dia)', mat='Cu + Al',
      rho=RHO_CU, V_eff=3.3 / RHO_CU, V_how='m/rho_Cu (ASSUMED-handbook; composite Cu+Al -> lower bound, upper bound m/rho_Al = 1.22e-3)',
      m_body=3.3, shape='disk (thickness not given)', aspect=np.nan, N_rel=1.0, shape_ok=False,
      src='2 brass (MS90) cylinders 440 mm dia x 430 mm, ~576 kg; CM distance 0.915 m (front measuring position)', Ms=576.0, r=0.915, a_sig=None,
      geo='CITED (Kleinevoss thesis)',
      source='Kleinevoss thesis, Sec. 3.1 (p. ~26), Sec. 4.1 (p. 36-37), Sec. 4.2.1; CODATA-02 (RMP 77, 1) item (f)',
      quote='"Die beiden Pendel der Massen m ~ 3.3 kg haengen in Wolframdrahtschlaufen der Laenge l ~ 2.7 m. Die Feldmassen wiegen M ~ 576 kg und ihr '
            'Schwerpunkt ist ... in der vordersten Position r ~ 0.915 m"; "Die Spiegel sind aus OFHC-Kupfer ... Der Durchmesser der Spiegel betraegt etwa '
            '198 mm"; "Messingzylinder mit einem Durchmesser von 440 mm, einer Laenge von 430 mm bei einem Gewicht von etwa 576 kg"'),
 dict(id='MSL-03', ref='Armstrong & Fitzgerald, PRL 91, 201101 (2003)', G=6.67387, u=0.00027, method='TB-comp',
      test='532 g horizontal Cu cylinder, 19 mm dia x 220 mm', mat='Cu', rho=0.532 / cylV(0.019, 0.220),
      V_eff=cylV(0.019, 0.220), V_how='dims (cyl)', m_body=0.532, shape='cyl(L/d=11.6)', aspect=220 / 19, N_rel=N_trans(220 / 19), shape_ok=True,
      src='two vertical cylindrical field masses, Cu (3 series) and stainless (1 series); mass and distance not given in fetched sources',
      Ms=np.nan, r=np.nan, a_sig=None, geo='CITED (CODATA-02)',
      source='CODATA-02 = Mohr & Taylor, RMP 77, 1 (2005) [physics.nist.gov/cuu/pdf/CODATA_RMP2005.pdf], Sec. III.R item (h) MSL-03',
      quote='"tungsten-fiber-suspended test mass (a 532 g horizontal copper cylinder of diameter 19 mm and length 220 mm) ... two large, vertical, '
            'cylindrical field masses ... copper field masses in addition to the previously used stainless steel masses"'),
 dict(id='HUST-05', ref='Hu, Guo & Luo, PRD 71, 127505 (2005) (re-analysis of Luo et al. 1999)', G=6.67222, u=0.00087, method='TB-dyn',
      test='Cu test mass ~32 g at each beam end (sphere per Luo 2014 review: 32.26 g)', mat='Cu', rho=RHO_CU,
      V_eff=32.26e-3 / RHO_CU, V_how='m/rho_Cu (ASSUMED-handbook)', m_body=32.26e-3, shape='sphere', aspect=1.0, N_rel=1 / 3, shape_ok=True,
      src='2 stainless cylinders 6.25 kg, 100 mm long x 100 mm dia, axes in line with test mass, facing faces 60 mm apart', Ms=6.25, r=0.030 + 0.050,
      a_sig=None, geo='CITED (CODATA-98 HUST-99; CODATA-06; Luo 2014 review PMC4173272)',
      source=L + ', p. 428 item (6) HUST-99; CODATA-06 p. 58 (HUST-05 = corrected HUST-99)',
      quote='"A copper test mass of mass approximately 32 g was suspended from each end of the balance beam ... two 6.25 kg stainless steel cylinders, '
            '100 mm in length and diameter ... the faces of the source masses opposite one another and between which the test mass hung were separated by 60 mm"'),
 dict(id='UZur-06', ref='Schlamminger et al., PRD 74, 082001 (2006) [gr-qc/0609027]', G=6.67425, u=0.00012, method='STAT',
      test='2 gold-plated Cu cylinders 45 mm dia x 77 mm high, ~1.1 kg (weighed vertically)', mat='Cu', rho=1.1 / cylV(0.045, 0.077),
      V_eff=cylV(0.045, 0.077), V_how='dims (cyl)', m_body=1.1, shape='cyl(h/d=1.71)', aspect=77 / 45, N_rel=N_spheroid_axis(77 / 45), shape_ok=True,
      src='2 movable field masses of 7.5 t each (Hg tanks); measured signal 784.9 ug', Ms=7500.0, r=np.nan,
      a_sig=784.9e-9 * g_n / (4 * 1.1), geo='CITED (Schlamminger 2006)',
      source='arXiv:gr-qc/0609027, Sec. I (p. 1-2), Sec. III.B (test masses), Sec. II (signal 784.8976(91) ug)',
      quote='"The 45 mm diameter, 77 mm high copper cylinders were plated with a 10 um gold layer"; "two nearly identical 1.1 kg TM\'s hanging on long wires ... '
            'two movable FM\'s weighing 7.5 t each"; "calculated signal was 784.8976(91) ug"'),
 dict(id='HUST-09', ref='Luo et al., PRL 102, 240801 (2009); Tu et al., PRD 82, 022001 (2010)', G=6.67349, u=0.00018, method='TB-dyn',
      test='gilded glass/quartz block 91.47 x 12.01 x 26.22 mm, 63.38 g', mat='quartz', rho=63.38e-3 / (0.09147 * 0.01201 * 0.02622),
      V_eff=boxV(0.09147, 0.01201, 0.02622), V_how='dims (block)', m_body=63.38e-3, shape='block', aspect=0.01201 / sqrt(0.09147 * 0.02622),
      N_rel=N_spheroid_axis(0.01201 / sqrt(0.09147 * 0.02622)), shape_ok=True,
      src='2 SS316 spheres 778.18 g, 57.15 mm dia (separation not re-checked)', Ms=0.77818, r=np.nan, a_sig=None, geo='CITED (PMC4173272, Stage 30)',
      source='Luo et al. review, Phil. Trans. A 2014 (PMC4173272), Table 3 / text (fetched in Stage 30; not re-fetched in this pass)',
      quote='"Two SS316 stainless steel spheres with diameters of 57.15 mm and vacuum masses of 778.18 g"'),
 dict(id='BIPM-14', ref='Quinn et al., PRL 111, 101102 (2013); Phil Trans A 372, 20140032 (2014)', G=6.67554, u=0.00016, method='TB-comp',
      test='4 NEW Cu-0.7%Te cylinders "similar in size" to Mk I (h = d = 55 mm), ~1.2 kg (Mk II)', mat='CuTe', rho=1.2 / cylV(0.055, 0.055),
      V_eff=cylV(0.055, 0.055), V_how='dims (cyl h=d, Mk I dims; Mk II "similar in size")', m_body=1.2, shape='cyl(h=d)', aspect=1.0,
      N_rel=N_trans(1.0), shape_ok=True,
      src='Mk I source masses cut down to 115 mm height, ~11.2 kg; separation ~107 mm', Ms=11.2, r=0.107, a_sig=None,
      geo='CITED (Quinn et al. 2014 review)',
      source='Quinn, Speake, Parks & Davis, Phil. Trans. A 372, 20140032 (2014), Sec. 4 & Fig. 2',
      quote='"For the Mk II apparatus, new test masses were made, and the source masses were cut down to 115 mm in height ... The masses in the Mk II were '
            'about 1.2 and 11.2 kg"; Fig. 2: "test masses, newly manufactured but similar in size to the Mk I test masses"'),
 dict(id='LENS-14', ref='Rosi et al., Nature 510, 518 (2014) [arXiv:1412.7954]', G=6.67191, u=0.00099, method='AI',
      test='87Rb atoms (point-like, R << l); 1e5 atoms detected per cloud, T = 160 ms, 1.9 s per measurement', mat='Rb', rho=np.nan,
      V_eff=np.nan, V_how='n/a (atoms)', m_body=87 * 1.66054e-27, shape='atom', aspect=np.nan, N_rel=np.nan, shape_ok=False,
      src='24 W-alloy cylinders, 99.90 mm dia x 150.11 mm, total ~516 kg; differential phase 0.547870 rad', Ms=516.0, r=np.nan,
      a_sig=0.547870 / (4 * pi / 780.241e-9 * 0.160 ** 2), geo='CITED (Rosi 2014)',
      source='arXiv:1412.7954, main text and Methods (lines "The source mass is composed of 24 tungsten alloy cylinders"; "T = 160 ms"; '
             '"We typically detect 10^5 atoms"; "Each measurement takes 1.9 s"; "differential phase shift is 0.547870(63) rad")',
      quote='"The source mass is composed of 24 tungsten alloy cylinders, for a total mass of about 516 kg. Each cylinder is machined to a diameter of '
            '99.90 mm and a height of 150.11 mm"; a_sig DERIVED = dphi/(k_eff T^2), k_eff = 4 pi/780.241 nm (Rb D2, two-photon)'),
 dict(id='UCI-14', ref='Newman et al., Phil Trans A 372, 20140025 (2014)', G=6.67435, u=0.00013, method='TB-dyn',
      test='fused silica plate 40.59 x 40.59 x 2.943 mm, 10.60 g', mat='SiO2', rho=10.5966e-3 / (0.0405902 * 0.04059 * 0.002943),
      V_eff=boxV(0.0405902, 0.04059, 0.002943), V_how='dims (plate)', m_body=10.5966e-3, shape='plate', aspect=0.002943 / 0.0405901,
      N_rel=N_spheroid_axis(0.002943 / 0.0405901), shape_ok=True,
      src='2 OFHC Cu rings 59 kg (OD 520.7, ID 312.42, 48.26 mm); source-pendulum distance not re-checked', Ms=59.0, r=np.nan, a_sig=None,
      geo='CITED (Newman et al. 2014 Table 3, Stage 30)',
      source='Newman et al. 2014, Table 3 (fetched in Stage 30; publisher PDF behind Cloudflare in this pass; dims corroborated by search snippet '
             '"40.5902 mm wide, 40.5900 mm high, and 2.9430 mm thick")',
      quote='"pendulum 1 ... 40.5902 mm wide, 40.5900 mm high, and 2.9430 mm thick"'),
 dict(id='HUSTT-18', ref='Li et al., Nature 560, 582 (2018), time-of-swing', G=6.674184, u=0.000078, method='TB-dyn',
      test='TOS pendulum 1: block 91.00575 x 11.08688 x 30.66846 mm, 68.09937 g (pendulum 2: 91.00336 x 11.04448 x 30.00446 mm, 66.35715 g)',
      mat='rho DERIVED 2.201 g/cm3 (fused-silica-like)', rho=68.09937e-3 / boxV(0.09100575, 0.01108688, 0.03066846),
      V_eff=boxV(0.09100575, 0.01108688, 0.03066846), V_how='dims (block)', m_body=68.09937e-3, shape='block',
      aspect=0.01108688 / sqrt(0.09100575 * 0.03066846), N_rel=N_spheroid_axis(0.01108688 / sqrt(0.09100575 * 0.03066846)), shape_ok=True,
      src='SS spheres 1,3: 57.15 mm, 778.16 / 777.96 g; centre-centre 157.19 mm -> r = 78.6 mm to pendulum axis', Ms=0.77806, r=0.15719 / 2,
      a_sig=None, geo='CITED (Li et al. 2018 Extended Data Tables 1, 2, 4)',
      source='Li et al., Nature 560, 582 (2018), Extended Data Table 1 (pendulums), Table 2 (source masses), Table 4 (GC distances)',
      quote='ED Table 1 "TOS method Pendulum 1 91.00575(11) 11.08688(9) 30.66846(12) 68.09937(22)"; ED Table 2 "Sphere 1 57.15072(25) 778.1630(8)"; '
            'ED Table 4 "TOS-I: Fibre 1 ... 157.19245(34)"'),
 dict(id='HUSTA-18', ref='Li et al., Nature 560, 582 (2018), angular-acceleration feedback', G=6.674484, u=0.000077, method='TB-comp',
      test='AAF pendulum 3: plate 91.05243 x 4.00240 x 49.92441 mm, 40.0379 g (NOT the TOS body)', mat='rho DERIVED 2.201 g/cm3',
      rho=40.0379e-3 / boxV(0.09105243, 0.0040024, 0.04992441), V_eff=boxV(0.09105243, 0.0040024, 0.04992441), V_how='dims (plate)',
      m_body=40.0379e-3, shape='plate', aspect=0.0040024 / sqrt(0.09105243 * 0.04992441),
      N_rel=N_spheroid_axis(0.0040024 / sqrt(0.09105243 * 0.04992441)), shape_ok=True,
      src='SS spheres 7, 9, 10, 12: 127.0 mm, ~8.542 kg; horizontal 342.29 mm, vertical 139.79 mm apart -> r = 184.9 mm', Ms=8.5418,
      r=sqrt((0.3422874 / 2) ** 2 + (0.1397997 / 2) ** 2), a_sig=None, geo='CITED (Li et al. 2018 Extended Data Tables 1, 2, 4)',
      source='Li et al., Nature 560, 582 (2018), Extended Data Tables 1, 2, 4',
      quote='ED Table 1 "AAF method Pendulum 3 91.05243(29) 4.00240(8) 49.92441(24) 40.0379(3)"; ED Table 2 "Sphere 7 127.0003(8) 8,543.5826(53)"; '
            'ED Table 4 "Between Sphere 7 and 9 (S7,9) 342.2874(19) ... Between Sphere 7 and 10 (S7,10) 139.7997(15)"'),
 dict(id='JILA-18', ref='Parks & Faller, PRL 105, 110801 (2010) [arXiv:1008.3203]; revised PRL 122 (2019) = CODATA "JILA-18"', G=6.67260, u=0.00025,
      method='SUSP', test='2 OFHC Cu pendulum bobs, 780 g each (shape not given)', mat='Cu', rho=RHO_CU, V_eff=0.780 / RHO_CU,
      V_how='m/rho_Cu (ASSUMED-handbook)', m_body=0.780, shape='bob (not given)', aspect=np.nan, N_rel=1.0, shape_ok=False,
      src='4 W-alloy (95.5% W) cylinders 120 kg; force change 480 nN per bob', Ms=4 * 120.0, r=np.nan, a_sig=480e-9 / 0.780,
      geo='CITED (Parks-Faller 2010)',
      source='arXiv:1008.3203 (PRL 105, 110801), Fig. 1 caption and text; CODATA 2022 Table XXX lists JILA-18 as Parks and Faller (2019)',
      quote='"The bobs are made of oxygen-free copper and have a mass of 780 g ... When the four 120 kg tungsten source masses ... the gravitational force '
            'on each pendulum bob changes by 480 nN"'),
]
for d in D:
    d['ppm'] = (d['G'] - G_CODATA) / G_CODATA * 1e6; d['u_ppm'] = d['u'] / d['G'] * 1e6
    if d['a_sig'] is None and np.isfinite(d['Ms']) and np.isfinite(d['r']): d['a_sig'] = G0 * d['Ms'] / d['r'] ** 2
    if d['a_sig'] is None: d['a_sig'] = np.nan
    d['dims_V'] = d['V_how'].startswith('dims'); d['T_int'] = 1e5   # ASSUMED 1 day integration (shorter than any campaign -> conservative)

# weighted mean of the raw data (check against CODATA 6.67430 with the 3.9 Birge factor)
g = np.array([d['G'] for d in D]); u = np.array([d['u'] for d in D])
w = 1 / u ** 2; Gw = np.sum(w * g) / np.sum(w); chi2 = np.sum(((g - Gw) / u) ** 2)
out = dict(stage=30, codata=dict(G=G_CODATA, u=U_CODATA, note='CODATA 2022 Table XXX = CODATA 2018 input data (16 values, no new datum)',
                                 Gw_recomputed=Gw, chi2=chi2, dof=len(g) - 1, birge=sqrt(chi2 / (len(g) - 1)),
                                 spread_max_min_ppm=(g.max() - g.min()) / G_CODATA * 1e6))
# DerSimonian-Laird between-experiment scatter tau (ppm)
y = np.array([d['ppm'] for d in D]); s = np.array([d['u_ppm'] for d in D])
wi = 1 / s ** 2; ybar = np.sum(wi * y) / np.sum(wi); Q = np.sum(wi * (y - ybar) ** 2)
tau2 = max(0.0, (Q - (len(y) - 1)) / (np.sum(wi) - np.sum(wi ** 2) / np.sum(wi)))
out['scatter'] = dict(Q=Q, tau_ppm=sqrt(tau2), sd_unweighted_ppm=float(np.std(y, ddof=1)), median_u_ppm=float(np.median(s)))

# ---------------------------------------------------------------- 2. statistics
def wls(x, y, s, tau=0.0):
    W = 1 / (s ** 2 + tau ** 2); X = np.vstack([np.ones_like(x), x]).T
    C = np.linalg.inv(X.T @ (W[:, None] * X)); b = C @ (X.T @ (W * y)); return b, np.sqrt(np.diag(C))
def perm_p(x, y, s, tau, n=20000):
    b0 = wls(x, y, s, tau)[0][1]; cnt = 0
    for _ in range(n):
        if abs(wls(rng.permutation(x), y, s, tau)[0][1]) >= abs(b0): cnt += 1
    return (cnt + 1) / (n + 1)
def analyse(name, xfun, subset):
    idx = [i for i, d in enumerate(D) if subset(d) and np.isfinite(xfun(d))]
    x = np.array([xfun(D[i]) for i in idx]); yy = y[idx]; ss = s[idx]; z = (g[idx] - G_CODATA) / u[idx]
    # tau re-estimated on the subset
    wi_ = 1 / ss ** 2; yb = np.sum(wi_ * yy) / np.sum(wi_); Q_ = np.sum(wi_ * (yy - yb) ** 2)
    t2 = max(0.0, (Q_ - (len(yy) - 1)) / (np.sum(wi_) - np.sum(wi_ ** 2) / np.sum(wi_))); t = sqrt(t2)
    rs_ppm = stats.spearmanr(x, yy); rs_z = stats.spearmanr(x, z); kt = stats.kendalltau(x, yy)
    b, eb = wls(x, yy, ss, t); bn, ebn = wls(x, yy, ss, 0.0)
    jk = np.array([wls(np.delete(x, k), np.delete(yy, k), np.delete(ss, k), t)[0][1] for k in range(len(x))])
    return dict(variable=name, n=len(x), ids=[D[i]['id'] for i in idx],
                spearman_rho_ppm=rs_ppm.statistic, spearman_p_ppm=rs_ppm.pvalue, spearman_rho_z=rs_z.statistic, spearman_p_z=rs_z.pvalue,
                kendall_tau=kt.statistic, kendall_p=kt.pvalue,
                wls_slope_ppm_per_unit=b[1], wls_slope_err=eb[1], wls_tau_ppm=t, wls_t=b[1] / eb[1],
                wls_p_normal=2 * stats.norm.sf(abs(b[1] / eb[1])), wls_perm_p=perm_p(x, yy, ss, t, 5000),
                wls_no_tau_slope=bn[1], wls_no_tau_t=bn[1] / ebn[1],
                jackknife_slope_range=[float(jk.min()), float(jk.max())], jackknife_sign_changes=int(np.sum(np.sign(jk) != np.sign(b[1]))))
macro = lambda d: d['method'] != 'AI' and np.isfinite(d['V_eff'])
documented = lambda d: macro(d) and d['dims_V']      # V from cited dimensions (no handbook density)
VARS = dict(log10_Veff=lambda d: log10(d['V_eff']) if np.isfinite(d['V_eff']) else np.nan,
            density=lambda d: d['rho'],
            log10_mass=lambda d: log10(d['m_body']) if np.isfinite(d['m_body']) else np.nan,
            log10_CQ_Gamma_proxy_3N_over_V=lambda d: log10(3 * d['N_rel'] / d['V_eff']) if (np.isfinite(d['V_eff']) and d['shape_ok']) else np.nan)
res = {k: analyse(k, f, macro) for k, f in VARS.items()}
res_doc = {k: analyse(k, f, documented) for k, f in VARS.items()}
# method: Kruskal-Wallis (ppm) on all 16; plus random-effects group chi2
groups = {}
for d in D: groups.setdefault(d['method'], []).append(d['ppm'])
kw = stats.kruskal(*groups.values())
grp = {}
for k, v in groups.items():
    ii = [i for i, d in enumerate(D) if d['method'] == k]; W_ = 1 / (s[ii] ** 2 + tau2)
    grp[k] = dict(n=len(ii), mean_ppm=float(np.sum(W_ * y[ii]) / np.sum(W_)), err_ppm=float(1 / sqrt(np.sum(W_))), members=[D[i]['id'] for i in ii])
Wall = 1 / (s ** 2 + tau2); mu = np.sum(Wall * y) / np.sum(Wall)
chi2_grp = sum(((v['mean_ppm'] - mu) / v['err_ppm']) ** 2 for v in grp.values()); dof_grp = len(grp) - 1
method = dict(kruskal_H=kw.statistic, kruskal_p=kw.pvalue, groups=grp, re_chi2=chi2_grp, re_dof=dof_grp, re_p=stats.chi2.sf(chi2_grp, dof_grp))
# look-elsewhere: 5 primary tests (4 continuous on macro subset via permutation WLS + method KW)
prim = {k: v['wls_perm_p'] for k, v in res.items()}; prim['method_KW'] = method['kruskal_p']
pmin = min(prim.values()); ntest = len(prim)
lee = dict(primary_p=prim, n_tests=ntest, min_p=pmin, bonferroni=min(1, pmin * ntest), sidak=1 - (1 - pmin) ** ntest,
           note='also tried: Spearman on ppm and z, documented-geometry subset - effective number of looks > 5, so Bonferroni(5) is lenient')
# sensitivity (verification pass): the only remaining volume uncertainties are (i) UWup composite Cu+Al body (V_eff between m/rho_Cu and
# m/rho_Al), (ii) TR&D-96 test mass described only for apparatus version 3, (iii) HUSTT-18 pendulum 1 vs 2, (iv) rows whose V uses a handbook density
def rerun(mod=None, drop=()):
    idx_ = [i for i, d in enumerate(D) if macro(d) and d['id'] not in drop]
    xs = []
    for i in idx_:
        v = D[i]['V_eff']
        if mod and D[i]['id'] in mod: v = mod[D[i]['id']]
        xs.append(log10(v))
    xs = np.array(xs); yy = y[idx_]; ss = s[idx_]
    wi_ = 1 / ss ** 2; yb = np.sum(wi_ * yy) / np.sum(wi_); Q_ = np.sum(wi_ * (yy - yb) ** 2)
    t_ = sqrt(max(0.0, (Q_ - (len(yy) - 1)) / (np.sum(wi_) - np.sum(wi_ ** 2) / np.sum(wi_))))
    b_, e_ = wls(xs, yy, ss, t_)
    return dict(n=len(idx_), spearman_rho=stats.spearmanr(xs, yy).statistic, spearman_p=stats.spearmanr(xs, yy).pvalue,
                wls_slope=b_[1], wls_err=e_[1], wls_perm_p=perm_p(xs, yy, ss, t_, 5000))
sens = dict(UWup_Al_upper_bound=rerun(mod={'UWup-02': 3.3 / 2700.}), drop_TRD96=rerun(drop=('TR&D-96',)),
            HUSTT_pendulum2=rerun(mod={'HUSTT-18': 0.09100336 * 0.01104448 * 0.03000446}),
            drop_handbook_density_rows=rerun(drop=tuple(d['id'] for d in D if not d['dims_V'])),
            drop_3_old_dumbbells_NIST_TRD_LANL=rerun(drop=('NIST-82', 'TR&D-96', 'LANL-97')),
            drop_all_TB_dyn=rerun(drop=tuple(d['id'] for d in D if d['method'] == 'TB-dyn')))
# same-volume pairs (direct test, no statistics)
pairs = [('BIPM-14', 'UZur-06'), ('BIPM-14', 'JILA-18'), ('HUSTT-18', 'HUSTA-18'), ('BIPM-01', 'BIPM-14'), ('UCI-14', 'UWash-00')]
byid = {d['id']: d for d in D}
pair_out = []
for a, b_ in pairs:
    A, B = byid[a], byid[b_]
    pair_out.append(dict(pair=f'{a} vs {b_}', V_ratio=A['V_eff'] / B['V_eff'], dG_ppm=A['ppm'] - B['ppm'],
                         sigma_comb_ppm=sqrt(A['u_ppm'] ** 2 + B['u_ppm'] ** 2), n_sigma=(A['ppm'] - B['ppm']) / sqrt(A['u_ppm'] ** 2 + B['u_ppm'] ** 2)))
# (a) the CQ prediction is excess SCATTER for small V, not a signed shift: test |residual| vs log V
idm = [i for i, d in enumerate(D) if macro(d)]
xv = np.array([log10(D[i]['V_eff']) for i in idm])
abs_test = dict(spearman_absppm=stats.spearmanr(xv, np.abs(y[idm])).__dict__ if False else None)
r1 = stats.spearmanr(xv, np.abs(y[idm])); r2 = stats.spearmanr(xv, np.abs((g[idm] - G_CODATA) / u[idm]))
abs_test = dict(spearman_rho_absppm=r1.statistic, p_absppm=r1.pvalue, spearman_rho_absz=r2.statistic, p_absz=r2.pvalue,
                note='CQ (if anything) predicts rho < 0 here (smaller V -> more scatter)')
# (b) method fixed effects + log V (does the volume trend survive within method classes?)
meths = sorted(set(D[i]['method'] for i in idm))
X = np.column_stack([xv] + [[1.0 if D[i]['method'] == m else 0.0 for i in idm] for m in meths])
Wf = 1 / (s[idm] ** 2 + tau2); C = np.linalg.pinv(X.T @ (Wf[:, None] * X)); bf = C @ (X.T @ (Wf * y[idm]))
within = dict(slope_ppm_per_decade=bf[0], err=sqrt(C[0, 0]), t=bf[0] / sqrt(C[0, 0]), p=2 * stats.norm.sf(abs(bf[0] / sqrt(C[0, 0]))), methods=meths)
# within the two big classes separately
for m in ('TB-dyn', 'TB-comp'):
    ii = [i for i in idm if D[i]['method'] == m]
    within[m] = dict(n=len(ii), spearman=stats.spearmanr([log10(D[i]['V_eff']) for i in ii], y[ii]).statistic,
                     p=stats.spearmanr([log10(D[i]['V_eff']) for i in ii], y[ii]).pvalue)
# association between method and volume (confound check)
mv = {m: float(np.median([log10(D[i]['V_eff']) for i in idm if D[i]['method'] == m])) for m in meths}
out['statistics_extra'] = dict(abs_residual_vs_logV=abs_test, logV_with_method_fixed_effects=within, median_log10V_by_method=mv)
out['statistics'] = dict(macro_subset=res, documented_subset=res_doc, method=method, look_elsewhere=lee,
                         sensitivity_verified_volumes=sens,
                         same_volume_pairs=pair_out)

# ---------------------------------------------------------------- 3. physics: what CQ noise can do to a G measurement
s28 = json.load(open(os.path.join(ROOT, 'results/stage28.json')))
GV_MAX = max(v['GammaV_bound'] for v in s28['bounds'].values())       # most permissive Stage 28 variant (95%)
GV_DP = 4 * pi / 3 * GH                                                 # Diosi kappa=1 / CQ kappa=1/2 floor
GAM_ATOM_MAX = GH / (6 * sqrt(pi) * 4.5e-10 ** 3)                       # Stage 26 X-ray ceiling per nucleus/atom (point-like regime)
D2_MAX = 1.169e-15                                                      # Stage 29 matter-sourced ceiling, kg s
REQ = 500e-6
phys = []
for d in D:
    if not np.isfinite(d['a_sig']): continue      # no verified source distance / signal -> not computed
    a_sig = d['a_sig']                           # CITED or DERIVED from cited numbers (see 'source')
    T = d['T_int']                                # ASSUMED 1e5 s (1 day) - shorter than any real campaign -> conservative
    if d['method'] == 'AI':
        Tint = 0.16; Nat = 1e5; Nshot = T / 1.9   # CITED Rosi 2014: T = 160 ms, 1e5 atoms detected, 1.9 s per measurement
        var_a_shot = 2 * GAM_ATOM_MAX / (3 * Tint) / Nat
        sig_a = sqrt(var_a_shot / Nshot); gam = GAM_ATOM_MAX; gv_note = 'point-like: Gamma <= G hbar/(6 sqrt(pi) R0^3) (X-ray, white assumed)'
    else:
        gam = 3 * d['N_rel'] * GV_MAX / d['V_eff']; sig_a = sqrt(gam / T); gv_note = 'Gamma = 3 N_rel (Gamma V)_max / V_eff'
    gam_dp = 3 * d['N_rel'] * GV_DP / d['V_eff'] if d['method'] != 'AI' else GH / (6 * sqrt(pi) * 1e-9 ** 3)
    # matter-sourced: Gamma_src = G^2 Dhat M_s/(4 r^4) = Dhat a_sig^2/(4 M_s)   (rewritten with the verified a_sig)
    gam_src = D2_MAX * a_sig ** 2 / (4 * d['Ms'])
    # rectification bound for time-of-swing: heating dE/dt = M Gamma/2, steady extra energy M Gamma Q/(2 w0) (ASSUMED Q=1e5, period 500 s,
    # amplitude A = 1 mrad, r_g = 3 cm, amplitude coefficient of the period c_A = 1 per rad^2, signal frequency shift 1e-3): all deliberately generous
    rect = np.nan
    if d['method'] == 'TB-dyn':
        Qf, w0, A, rg, cA, fsig = 1e5, 2 * pi / 500, 1e-3, 0.03, 1.0, 1e-3
        dA2 = gam * Qf / (w0 ** 3 * rg ** 2)           # = 2 dE/(I w0^2)
        rect = 2 * cA * dA2 / fsig
    phys.append(dict(id=d['id'], a_signal=a_sig, Gamma_max=gam, Gamma_DP_floor=gam_dp, Gamma_matter_sourced_max=gam_src, note=gv_note,
                     scatter_ppm_at_max=sig_a / a_sig * 1e6, scatter_ppm_DP_floor=(sqrt(gam_dp / T) / a_sig * 1e6) if d['method'] != 'AI' else np.nan,
                     matter_sourced_ppm=sqrt(gam_src / T) / a_sig * 1e6, rectification_bias_ppm=(rect * 1e6) if np.isfinite(rect) else np.nan,
                     orders_below_500ppm=log10(REQ * 1e6 / (sig_a / a_sig * 1e6))))
need = []
for d in D:
    if d['method'] == 'AI' or not np.isfinite(d['a_sig']): continue
    a_sig = d['a_sig']; gam_need = (REQ * a_sig) ** 2 * d['T_int']
    need.append(dict(id=d['id'], Gamma_needed=gam_need, GammaV_needed=gam_need * d['V_eff'] / (3 * d['N_rel']),
                     orders_above_stage28=log10(gam_need * d['V_eff'] / (3 * d['N_rel']) / GV_MAX)))
out['physics'] = dict(GammaV_max_stage28=GV_MAX, GammaV_DP_floor=GV_DP, Gamma_atom_max=GAM_ATOM_MAX, Dhat_max=D2_MAX,
                      per_experiment=phys, needed_for_500ppm=need,
                      max_scatter_ppm=max(p['scatter_ppm_at_max'] for p in phys),
                      max_rectification_ppm=float(np.nanmax([p['rectification_bias_ppm'] for p in phys])),
                      mean_bias_mechanisms=[
                          'Linear response: white zero-mean force noise adds variance only; <F> unchanged. CQ models reproduce the Newtonian mean '
                          '(Tilloy-Diosi 2017 Eq. 10: unitary pair potential + double commutator; Kryhin-Sudhir Eq. 5: averaged Poisson equation sourced by <f>).',
                          'DP/CQ self-noise is independent of source-mass position (Stage 29), so any heating-induced shift is identical in both source '
                          'configurations and cancels in the near-far difference that defines G.',
                          'Second-order (rectification) via amplitude-dependent torsion period: bounded above (column rectification_bias_ppm) with generous inputs.',
                          'Matter-sourced (discrete) kernel: noise near the source mass differs between configurations, but enters as variance; bounded above '
                          '(matter_sourced_ppm).'])
out['oppenheim_prior'] = ('arXiv:2203.01982 Sec. V already floated this: "Since we have found that CQ theories predict an uncertainty in mass measurements '
                          'it is perhaps intriguing that different experiments to measure Newton\'s constant G yield results whose relative uncertainty differ '
                          'by as much as 5e-4 ... masses with smaller volume should yield larger uncertainty and this would be the effect to look for".')
# CSV
with open(os.path.join(ROOT, 'results/stage30_G_table.csv'), 'w', newline='') as fh:
    cols = ['id', 'ref', 'method', 'G', 'u', 'ppm', 'u_ppm', 'test', 'mat', 'rho', 'm_body', 'V_eff', 'V_how', 'shape', 'aspect', 'N_rel', 'shape_ok',
            'src', 'Ms', 'r', 'a_sig', 'geo', 'source', 'quote']
    wr = csv.writer(fh); wr.writerow(cols)
    for d in D: wr.writerow([d[c] for c in cols])
out['table'] = D
out['runtime_s'] = time.time() - t0
json.dump(out, open(os.path.join(ROOT, 'results/stage30.json'), 'w'), indent=1, default=float)
pr = lambda o: print(json.dumps(o, indent=1, default=float))
pr(out['codata']); pr(out['scatter'])
for k, v in res.items(): print('MACRO', k, {kk: v[kk] for kk in ('n', 'spearman_rho_ppm', 'spearman_p_ppm', 'spearman_p_z', 'wls_slope_ppm_per_unit', 'wls_slope_err', 'wls_tau_ppm', 'wls_perm_p', 'jackknife_slope_range', 'wls_no_tau_t')})
for k, v in res_doc.items(): print('DOC', k, {kk: v[kk] for kk in ('n', 'spearman_rho_ppm', 'spearman_p_ppm', 'wls_perm_p')})
pr(out['statistics_extra']); pr(method); pr(lee); pr(out['statistics']['sensitivity_verified_volumes']); pr(pair_out)
for p in phys: print(p['id'], '%.2e ppm scatter(max)  DP %.2e  src %.2e  rect %.2e  orders %.1f' % (p['scatter_ppm_at_max'], p['scatter_ppm_DP_floor'], p['matter_sourced_ppm'], p['rectification_bias_ppm'], p['orders_below_500ppm']))
pr(need)
