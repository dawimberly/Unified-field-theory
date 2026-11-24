"""Stage 3 SCOPING ONLY (nothing is run): compute cost of adding more fields/particles at D=13.
All numbers are order-of-magnitude estimates; assumptions labeled."""
import json, math
out = {}
FLOPS_EXA = 1e18          # ~1 exaflop/s machine (Frontier-class), ASSUMED sustained
BYTES_PER = 16            # complex double

# (1) Helium, brute-force grid in 2*D configuration dims
he = []
for D in [3, 13]:
    for npts in [32, 100]:
        pts = float(npts) ** (2 * D)
        he.append(dict(D=D, pts_per_axis=npts, grid_points=pts, memory_bytes=pts * BYTES_PER))
out['helium_bruteforce_grid'] = he
# (2) Helium, symmetry-reduced: for L=0 states the internal coords are (r1, r2, r12) for ANY D>=2
#     (nucleus + 2 electrons span a plane). Hylleraas basis of size K -> dense eig K^3.
out['helium_hylleraas'] = [dict(D=D, internal_coords=3, basis_K=K, flops=K**3, seconds_on_laptop=K**3/1e10)
                           for D in [3, 13] for K in [500, 2000]]
# (3) N electrons: internal coordinates N(N+1)/2 when D >= N+1 (else DN - D(D+1)/2 ... general count)
def internal(N, D):
    # positions of N electrons relative to a fixed nucleus in R^D, modulo SO(D)
    k = min(N, D)
    return N * D - (D * (D - 1) // 2 - (D - k) * (D - k - 1) // 2)
out['N_electron_internal_coords'] = [dict(N=N, D3=internal(N, 3), D13=internal(N, 13)) for N in [2, 3, 6, 10, 26]]
# (4) Lattice field theory with full SM field content in 13 space + 1 time = 14 spacetime dims
#     ASSUMED: gauge group dim 12, 3 generations x 16 Weyl fermions (+nu_R) = 48 Weyl; Higgs 4 real.
d = 14
spinor_complex = 2 ** (d // 2 - 1)       # Weyl spinor complex components in even d=14: 2^(d/2-1)=64
gauge_link_reals = 12 * d                # one Lie-algebra element per link direction (lower bound)
fermion_reals = 48 * spinor_complex * 2
higgs_reals = 4
per_site = gauge_link_reals + fermion_reals + higgs_reals
lat = []
for L in [8, 16, 32]:
    sites = float(L) ** d
    lat.append(dict(L=L, sites=sites, reals_per_site=per_site, memory_bytes=sites * per_site * 8,
                    one_sweep_flops=sites * per_site * 1e3,     # ASSUMED ~1e3 flop per dof per update
                    seconds_per_sweep_exa=sites * per_site * 1e3 / FLOPS_EXA))
out['lattice_SM_d14'] = dict(per_site_reals=per_site, spinor_complex_components=spinor_complex, cases=lat)
# reference: state-of-art 4D lattice QCD ~ 144^3 x 288 sites
out['reference_4d_lattice_sites'] = 144**3 * 288
# power counting: gauge coupling dimension in d spacetime dims = (4-d)/2 -> non-renormalizable for d>4
out['gauge_coupling_mass_dim_d14'] = (4 - d) / 2
json.dump(out, open('results/stage3.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
