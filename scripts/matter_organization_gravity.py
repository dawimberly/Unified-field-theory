"""Matched matter arrangements and their Newtonian pull.

Reproducible calculation for docs/MATTER_ORGANIZATION_GRAVITY_RESULTS.md.
The Candidate-0 action in this repository does not define an organization-dependent
correction on top of stress-energy. This script computes the ordinary baseline and
one explicitly labeled toy coupling. It does not use an Earth or Moon interior model.

Run from anywhere:
    python scripts/matter_organization_gravity.py
"""

from __future__ import annotations

import json
import platform
from pathlib import Path

import numpy as np
from scipy import constants
import scipy

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "experiments" / "matter_organization_gravity_2026_11_19"
C_LIGHT = float(constants.c)
G_NEWTON = float(constants.G)
# Repo benchmark (dimensionless). Not a fit to a particle mass.
V_VEV = 1.0
LAMBDA_POT = 64.0
KAPPA_SKYRME = 1.0


def wrap_delta(delta: np.ndarray | float) -> np.ndarray | float:
    """Principal value in (-pi, pi]."""
    return (np.asarray(delta, dtype=float) + np.pi) % (2.0 * np.pi) - np.pi


def potential(rho: float, lam: float = LAMBDA_POT, v: float = V_VEV) -> float:
    return 0.25 * lam * (rho * rho - v * v) ** 2


def potential_second_derivative(rho: float, lam: float = LAMBDA_POT, v: float = V_VEV) -> float:
    # d²/dρ² [λ (ρ²-v²)² / 4] = λ (3ρ² - v²).
    return lam * (3.0 * rho * rho - v * v)


def numerical_second_derivative(func, x: float, step: float = 1e-5) -> float:
    return (func(x + step) - 2.0 * func(x) + func(x - step)) / step**2


def cube_vertices() -> np.ndarray:
    """Eight vertices of the unit cube, then centered. Shape (8, 3)."""
    pts = np.array([[x, y, z] for x in (0.0, 1.0) for y in (0.0, 1.0) for z in (0.0, 1.0)])
    return pts - 0.5


def cube_edges(points: np.ndarray) -> list[tuple[int, int]]:
    """Geometric edges of the embedded cube. Same graph as the periodic 2×2×2 lattice."""
    edges = []
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            if np.count_nonzero(np.abs(points[i] - points[j]) > 0.5) == 1 and np.isclose(
                np.linalg.norm(points[i] - points[j]), 1.0
            ):
                edges.append((i, j))
    if len(edges) != 12:
        raise RuntimeError(f"expected 12 cube edges, found {len(edges)}")
    return edges


def cycle_laplacian(n: int) -> np.ndarray:
    lap = np.zeros((n, n))
    for i in range(n):
        lap[i, i] = 2.0
        lap[i, (i + 1) % n] = -1.0
        lap[i, (i - 1) % n] = -1.0
    return lap


def hypercube_laplacian_eigenvalues(dim: int = 3) -> np.ndarray:
    """Laplacian of the dim-dimensional hypercube: eigenvalue 2k with multiplicity C(dim, k)."""
    from math import comb

    vals = []
    for k in range(dim + 1):
        vals.extend([2.0 * k] * comb(dim, k))
    return np.array(vals, dtype=float)


def phase_energy(theta: np.ndarray, edges: list[tuple[int, int]]) -> float:
    acc = 0.0
    for i, j in edges:
        d = float(wrap_delta(theta[i] - theta[j]))
        acc += 0.5 * d * d
    return acc


def orientation_energy(spins: np.ndarray, edges: list[tuple[int, int]]) -> float:
    acc = 0.0
    for i, j in edges:
        acc += 1.0 - float(np.dot(spins[i], spins[j]))
    return acc


def pair_cosine(theta: np.ndarray, pairs: list[tuple[int, int]]) -> float:
    if not pairs:
        return float("nan")
    acc = 0.0
    for i, j in pairs:
        acc += float(np.cos(wrap_delta(theta[i] - theta[j])))
    return acc / len(pairs)


def ring_second_neighbor_pairs(n: int) -> list[tuple[int, int]]:
    return [(i, (i + 2) % n) for i in range(n)]


def ring_edges(n: int) -> list[tuple[int, int]]:
    return [(i, (i + 1) % n) for i in range(n)]


def winding(theta: np.ndarray) -> float:
    edges = ring_edges(len(theta))
    total = 0.0
    for i, j in edges:
        total += float(wrap_delta(theta[j] - theta[i]))
    return total / (2.0 * np.pi)


def plaquette_skyrme_density(spins: np.ndarray, points: np.ndarray) -> float:
    """Sum of [n · (Δ_a n × Δ_b n)]² on the six cube faces.

    Any discretization built from n · (Δn × Δn) vanishes when every spin is
    collinear, so this is a structural zero rather than a tuned cancellation.
    """
    # Faces: constant coordinate x, y, or z.
    acc = 0.0
    coords = np.round(points + 0.5).astype(int)  # back to {0,1}^3
    index = {tuple(coords[i]): i for i in range(len(coords))}
    for axis in range(3):
        for level in (0, 1):
            axes = [a for a in range(3) if a != axis]
            corners = []
            for u in (0, 1):
                for v in (0, 1):
                    key = [0, 0, 0]
                    key[axis] = level
                    key[axes[0]] = u
                    key[axes[1]] = v
                    corners.append(index[tuple(key)])
            # corners ordered (u,v) = (0,0), (0,1), (1,0), (1,1) if u is inner.
            # Rebuild explicitly.
            c00 = index[tuple(0 if a != axis else level for a in range(3))]
            # The comprehension above is wrong for off-axis bits. Use direct keys.
            def at(u: int, v: int) -> int:
                key = [0, 0, 0]
                key[axis] = level
                key[axes[0]] = u
                key[axes[1]] = v
                return index[tuple(key)]

            i00, i10, i01 = at(0, 0), at(1, 0), at(0, 1)
            n0 = spins[i00]
            d1 = spins[i10] - n0
            d2 = spins[i01] - n0
            f = float(np.dot(n0, np.cross(d1, d2)))
            acc += f * f
    return acc


def rotate_all(spins: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    c, s = np.cos(angle), np.sin(angle)
    out = []
    for n in spins:
        out.append(n * c + np.cross(axis, n) * s + axis * np.dot(axis, n) * (1.0 - c))
    return np.array(out)


def multipoles(masses: np.ndarray, positions: np.ndarray) -> dict:
    masses = np.asarray(masses, dtype=float)
    positions = np.asarray(positions, dtype=float)
    total = float(masses.sum())
    if total == 0.0:
        raise RuntimeError("massless arrangement")
    dipole_origin = masses @ positions
    center = dipole_origin / total
    rel = positions - center
    rr = np.einsum("a,ai,aj->ij", masses, rel, rel)
    r2 = float(np.einsum("a,ai,ai->", masses, rel, rel))
    quad = 3.0 * rr - r2 * np.eye(3)
    return {
        "mass": total,
        "center_of_mass": center,
        "dipole_about_origin": dipole_origin,
        "quadrupole_about_cm": quad,
        "quadrupole_frobenius": float(np.linalg.norm(quad)),
        "quadrupole_trace": float(np.trace(quad)),
    }


def potential_multipole(masses, positions, point, g_newton: float = 1.0) -> float:
    point = np.asarray(point, dtype=float)
    acc = 0.0
    for m, r in zip(masses, positions):
        acc += -g_newton * float(m) / np.linalg.norm(point - r)
    return acc


def potential_from_multipoles(multi: dict, point, g_newton: float = 1.0) -> float:
    point = np.asarray(point, dtype=float)
    # Expansion about the center of mass.
    rel = point - multi["center_of_mass"]
    radius = np.linalg.norm(rel)
    n = rel / radius
    quad = multi["quadrupole_about_cm"]
    return -g_newton * multi["mass"] / radius - g_newton * float(n @ quad @ n) / (2.0 * radius**3)


def mass_distribution(site_masses, site_positions, bond_masses, bond_positions):
    masses = np.concatenate([np.asarray(site_masses, dtype=float), np.asarray(bond_masses, dtype=float)])
    positions = np.concatenate(
        [np.asarray(site_positions, dtype=float), np.asarray(bond_positions, dtype=float)], axis=0
    )
    return masses, positions


def bond_midpoints(points: np.ndarray, edges: list[tuple[int, int]]) -> np.ndarray:
    return np.array([(points[i] + points[j]) / 2.0 for i, j in edges])


def two_mass_multipole_check() -> dict:
    """Analytic control: two equal masses on the z axis. Quadrupole formula through O(1/R³)."""
    mass = 1.7
    separation = 1.0
    positions = np.array([[0.0, 0.0, -separation], [0.0, 0.0, separation]])
    masses = np.array([mass, mass])
    multi = multipoles(masses, positions)
    point = np.array([0.0, 0.0, 100.0])
    direct = potential_multipole(masses, positions, point)
    approx = potential_from_multipoles(multi, point)
    q_zz_expected = 4.0 * mass * separation**2
    return {
        "quadrupole_zz": float(multi["quadrupole_about_cm"][2, 2]),
        "quadrupole_zz_expected": q_zz_expected,
        "trace": multi["quadrupole_trace"],
        "relative_error_at_R100": abs(direct - approx) / abs(direct),
        "direct": direct,
        "approx": approx,
    }


def ring_arrangements() -> dict:
    n = 4
    edges = ring_edges(n)
    seconds = ring_second_neighbor_pairs(n)
    locked = np.zeros(n)
    circulating = np.array([0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi])
    anticirculating = -circulating
    # Same nearest-neighbor |Δθ| as the circulation, but successive steps reverse.
    folded = np.array([0.0, 0.5 * np.pi, 0.0, -0.5 * np.pi])
    configs = {
        "locked": locked,
        "circulating": circulating,
        "anticirculating": anticirculating,
        "folded": folded,
    }
    lap = cycle_laplacian(n)
    eigs = np.linalg.eigvalsh(lap)
    # Hessian of h = (1/2) θᵀ L θ is L. Zeros are global phase rotations.
    out = {"n_sites": n, "laplacian_eigenvalues": eigs.tolist(), "configurations": {}}
    for name, theta in configs.items():
        nn = []
        for i, j in edges:
            nn.append(abs(float(wrap_delta(theta[i] - theta[j]))))
        grad = lap @ theta if name != "folded" else None
        # Folded state: derivative of (1/2) Σ wrap(Δ)². At |Δ|=π/2 < π the wrap is inactive.
        forces = []
        for i in range(n):
            d_prev = float(wrap_delta(theta[i] - theta[(i - 1) % n]))
            d_next = float(wrap_delta(theta[i] - theta[(i + 1) % n]))
            forces.append(d_prev + d_next)
        out["configurations"][name] = {
            "theta": theta.tolist(),
            "nearest_neighbor_abs_delta": nn,
            "nearest_neighbor_energy": phase_energy(theta, edges),
            "C2": pair_cosine(theta, seconds),
            "winding": winding(theta),
            "max_abs_force": float(np.max(np.abs(forces))),
            "energy_under_global_phase_shift": phase_energy(theta + 1.234, edges),
            "C2_under_global_phase_shift": pair_cosine(theta + 1.234, seconds),
        }
    # Long-wave relation C2 ≈ cos(2Δ) and 1 - 4*(bond energy) for a common small step.
    small = 0.01
    slow = np.array([0.0, small, 2 * small, 3 * small])
    # On a 4-ring a uniform small step is not closed. Use an open 3-bond comparison instead.
    delta = 0.01
    bond_energy = 0.5 * delta * delta
    c2_slow = float(np.cos(2.0 * delta))
    out["long_wave_check"] = {
        "delta": delta,
        "bond_energy": bond_energy,
        "C2": c2_slow,
        "C2_plus_4_bond_energy": c2_slow + 4.0 * bond_energy,
        "deviation_from_1": c2_slow + 4.0 * bond_energy - 1.0,
    }
    out["stable_critical_windings"] = {
        "locked_max_force": out["configurations"]["locked"]["max_abs_force"],
        "circulating_max_force": out["configurations"]["circulating"]["max_abs_force"],
        "anticirculating_max_force": out["configurations"]["anticirculating"]["max_abs_force"],
        "folded_max_force": out["configurations"]["folded"]["max_abs_force"],
        "hessian_eigenvalues": eigs.tolist(),
        "n_negative_hessian_eigenvalues": int(np.sum(eigs < -1e-10)),
        "n_zero_hessian_eigenvalues": int(np.sum(np.abs(eigs) <= 1e-10)),
    }
    return out


def orientation_sector() -> dict:
    points = cube_vertices()
    edges = cube_edges(points)
    aligned = np.tile(np.array([0.0, 0.0, 1.0]), (8, 1))
    parity = np.array([1 if int(round(x + y + z)) % 2 == 0 else -1 for x, y, z in np.round(points + 0.5)])
    neel = parity[:, None] * aligned
    flipped = aligned.copy()
    flipped[-1] *= -1.0

    def pack(spins: np.ndarray) -> dict:
        bond_e = [1.0 - float(np.dot(spins[i], spins[j])) for i, j in edges]
        return {
            "energy": float(sum(bond_e)),
            "bond_energies": bond_e,
            "skyrme_density": plaquette_skyrme_density(spins, points),
            "mean_neighbor_dot": float(np.mean([np.dot(spins[i], spins[j]) for i, j in edges])),
        }

    # Global rotation must leave both stationary collinear energies invariant.
    rotated_aligned = rotate_all(aligned, np.array([1.0, 1.0, 0.3]), 0.4)
    rotated_neel = rotate_all(neel, np.array([0.2, -0.5, 0.7]), -0.7)

    # Unstable direction of the Néel point: identical lab-frame transverse tilt, then renormalize.
    def tilted_neel(eps: float) -> np.ndarray:
        raw = np.stack([np.full(8, eps), np.zeros(8), parity.astype(float)], axis=1)
        return raw / np.linalg.norm(raw, axis=1, keepdims=True)

    eps = 1e-4
    energy_neel = orientation_energy(neel, edges)
    energy_tilt = orientation_energy(tilted_neel(eps), edges)
    # h(ε) = 24 / (1+ε²) for this normalization, so [h(ε)-h(0)]/ε² → -24.
    curvature_estimate = (energy_tilt - energy_neel) / eps**2

    # Vacuum transverse Hessian: each of 2 lab directions has the hypercube Laplacian.
    # h ≈ Σ_bonds (1/2) |π_i - π_j|² = (1/2) Σ_a π^a · L · π^a, Hessian per component = L.
    vac_eigs = hypercube_laplacian_eigenvalues(3)
    # Numerical check of one Laplacian entry via finite differences on a single transverse coordinate.
    def energy_one_transverse(pi: np.ndarray) -> float:
        spins = aligned.copy()
        # pi is the x-component; z = sqrt(1-pi²), y=0. Valid for small pi.
        spins[:, 0] = pi
        spins[:, 2] = np.sqrt(np.maximum(0.0, 1.0 - pi * pi))
        return orientation_energy(spins, edges)

    basis = np.eye(8)
    hess = np.zeros((8, 8))
    # Step is large enough that the O(step²) energy is well above roundoff, and small
    # enough that the quartic remainder in the unit-vector constraint is negligible.
    step = 1e-4
    for a in range(8):
        for b in range(a, 8):
            epp = energy_one_transverse(step * basis[a] + step * basis[b])
            epm = energy_one_transverse(step * basis[a] - step * basis[b])
            emp = energy_one_transverse(-step * basis[a] + step * basis[b])
            emm = energy_one_transverse(-step * basis[a] - step * basis[b])
            val = (epp - epm - emp + emm) / (4.0 * step * step)
            hess[a, b] = hess[b, a] = val
    num_eigs = np.linalg.eigvalsh(hess)
    return {
        "n_sites": 8,
        "n_edges": len(edges),
        "aligned": pack(aligned),
        "neel": pack(neel),
        "single_flip": pack(flipped),
        "rotation_invariance_aligned_abs_diff": abs(pack(rotated_aligned)["energy"] - pack(aligned)["energy"]),
        "rotation_invariance_neel_abs_diff": abs(pack(rotated_neel)["energy"] - pack(neel)["energy"]),
        "neel_tilt_curvature_over_eps2": curvature_estimate,
        "neel_tilt_curvature_expected": -24.0,
        "vacuum_orientation_eigenvalues_analytic_per_transverse_component": vac_eigs.tolist(),
        "vacuum_orientation_eigenvalues_numerical_x": num_eigs.tolist(),
        "vacuum_orientation_eigenvalue_abs_error": float(np.max(np.abs(np.sort(num_eigs) - np.sort(vac_eigs)))),
        "points": points,
        "edges": edges,
        "parity": parity,
        "spins": {"aligned": aligned, "neel": neel, "single_flip": flipped},
    }


def newtonian_baseline(orient: dict, j_over_rest: float) -> dict:
    """Rest mass on sites, gradient energy on bond midpoints. Model units c = 1, ε_0 = 1, a = 1.

    j_over_rest is J/ε_0 with J = v² a. It is an illustrative ratio, not a measured binding fraction.
    """
    points = orient["points"]
    edges = orient["edges"]
    mids = bond_midpoints(points, edges)
    rest = np.ones(8)
    j = j_over_rest  # because ε_0 = 1

    def assemble(bond_dimensionless: list[float]) -> dict:
        bond_mass = j * np.asarray(bond_dimensionless, dtype=float)  # / c² with c=1
        masses, pos = mass_distribution(rest, points, bond_mass, mids)
        multi = multipoles(masses, pos)
        far = np.array([0.0, 0.0, 80.0])
        direct = potential_multipole(masses, pos, far)
        approx = potential_from_multipoles(multi, far)
        return {
            "mass": multi["mass"],
            "center_of_mass": multi["center_of_mass"].tolist(),
            "dipole_about_origin": multi["dipole_about_origin"].tolist(),
            "quadrupole_frobenius": multi["quadrupole_frobenius"],
            "quadrupole_trace": multi["quadrupole_trace"],
            "quadrupole_eigenvalues": np.linalg.eigvalsh(multi["quadrupole_about_cm"]).tolist(),
            "potential_direct_R80": direct,
            "potential_multipole_R80": approx,
            "multipole_relative_error": abs(direct - approx) / abs(direct),
            "masses": masses,
            "positions": pos,
            "multi": multi,
        }

    aligned_bonds = orient["aligned"]["bond_energies"]
    neel_bonds = orient["neel"]["bond_energies"]
    flip_bonds = orient["single_flip"]["bond_energies"]
    aligned = assemble(aligned_bonds)
    neel = assemble(neel_bonds)
    flip = assemble(flip_bonds)
    # Same total binding as the flip, spread uniformly over the eight sites.
    extra = j * float(np.sum(flip_bonds))
    smeared_site = rest + extra / 8.0
    smeared_masses, smeared_pos = mass_distribution(smeared_site, points, np.zeros(12), mids)
    smeared_multi = multipoles(smeared_masses, smeared_pos)
    # The flip quadrupole is traceless with Q_xx=Q_yy=Q_zz=0 in the cube axes, so a
    # coordinate-axis ray does not see it. The body diagonal is the nonzero eigendirection.
    diagonal = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    radii = (10.0, 40.0, 80.0)
    # Common geometric origin: a center-of-mass shift is a dipole. After each
    # distribution is translated so its own center of mass sits at the origin,
    # that dipole is matched and the residual is the quadrupole.
    flip_centered = flip["positions"] - flip["multi"]["center_of_mass"]
    smear_centered = smeared_pos - smeared_multi["center_of_mass"]
    falloff = []
    for radius in radii:
        far = diagonal * radius
        d_flip = potential_multipole(flip["masses"], flip["positions"], far)
        d_smear = potential_multipole(smeared_masses, smeared_pos, far)
        d_flip_cm = potential_multipole(flip["masses"], flip_centered, far)
        d_smear_cm = potential_multipole(smeared_masses, smear_centered, far)
        falloff.append(
            {
                "R": radius,
                "direction": "body diagonal",
                "fractional_difference_common_origin": (d_flip - d_smear) / d_smear,
                "fractional_difference_cm_matched": (d_flip_cm - d_smear_cm) / d_smear_cm,
            }
        )
    # Drop non-serializable objects for the returned public part.
    def pub(block: dict) -> dict:
        return {k: v for k, v in block.items() if k not in ("masses", "positions", "multi")}

    delta_e = (neel["mass"] - aligned["mass"])  # c=1, this is ΔE
    return {
        "units": "c=1, a=1, rest energy per site ε_0=1, J/ε_0 = j_over_rest. Illustrative, not a material prediction.",
        "j_over_rest": j_over_rest,
        "delta_E_neel_minus_aligned": delta_e,
        "delta_M_from_delta_E_over_c2": delta_e,  # c=1
        "fractional_monopole_neel_vs_aligned": (neel["mass"] - aligned["mass"]) / aligned["mass"],
        "aligned": pub(aligned),
        "neel": pub(neel),
        "single_flip": pub(flip),
        "flip_smeared_same_total_mass": {
            "mass": smeared_multi["mass"],
            "quadrupole_frobenius": smeared_multi["quadrupole_frobenius"],
            "mass_minus_flip": smeared_multi["mass"] - flip["mass"],
        },
        "flip_versus_smeared_fractional_potential": falloff,
        "_arrays": {"aligned": aligned, "neel": neel, "flip": flip},
    }


def phase_order_statistics(n_samples: int = 20000, n_sites: int = 8, seed: int = 0) -> dict:
    """Exact and Monte-Carlo moments. Mean shifts are not assigned a 1/√N factor."""
    rng = np.random.default_rng(seed)
    # Sample-spanning order on N independent uniform phases.
    # E[Φ] = 1/N, Var(Φ) = (N-1)/N³, Φ = |N⁻¹ Σ exp(iθ)|².
    exact = {}
    mc = {}
    for n in (4, 8, 64, 512):
        exact[str(n)] = {
            "mean_Phi": 1.0 / n,
            "variance_Phi": (n - 1) / n**3,
            "rms_Phi": float(np.sqrt((n - 1) / n**3)),
        }
    phases = rng.uniform(0.0, 2.0 * np.pi, size=(n_samples, n_sites))
    amps = np.exp(1j * phases)
    phi = np.abs(amps.mean(axis=1)) ** 2
    # Nearest-neighbor mean on the cube graph.
    points = cube_vertices()
    edges = cube_edges(points)
    cosines = []
    for i, j in edges:
        cosines.append(np.cos(phases[:, i] - phases[:, j]))
    cbar = np.mean(np.stack(cosines, axis=1), axis=1)
    mc = {
        "n_sites": n_sites,
        "n_samples": n_samples,
        "seed": seed,
        "Phi_mean": float(phi.mean()),
        "Phi_mean_exact": 1.0 / n_sites,
        "Phi_rms": float(phi.std(ddof=0)),
        "Phi_rms_exact": float(np.sqrt((n_sites - 1) / n_sites**3)),
        "short_range_mean": float(cbar.mean()),
        "short_range_mean_exact": 0.0,
        "short_range_rms": float(cbar.std(ddof=0)),
        "iid_bond_rms_estimate": float(np.sqrt(0.5 / len(edges))),
        "n_bonds": len(edges),
    }
    return {"sample_spanning_exact": exact, "monte_carlo": mc}


def scaling_table(beta: float = 1.0) -> list[dict]:
    """Fractional source-swap pull at fixed density and fixed lattice spacing.

    β here is the dimensionless ratio α / (ρ c²), held fixed while N grows.
    The mean shift is β times the mean correlator. The 1/√N column is the rms of
    a sample-spanning fluctuation and is not added to the mean.
    """
    rows = []
    for n_side in (2, 4, 8, 16):
        n = n_side**3
        rows.append(
            {
                "N": n,
                "short_range_mean_shift_if_DeltaC_is_1": beta * 1.0,
                "sample_spanning_locked_mean_shift": beta * 1.0,
                "sample_spanning_random_mean_shift": beta / n,
                "sample_spanning_random_rms": beta * float(np.sqrt((n - 1) / n**3)),
            }
        )
    return rows


def wafer_grid(n_lat: int, n_thick: int, z0: float, thickness: float, lx: float, ly: float) -> np.ndarray:
    xs = (np.arange(n_lat) + 0.5) / n_lat - 0.5
    ys = (np.arange(n_lat) + 0.5) / n_lat - 0.5
    zs = (np.arange(n_thick) + 0.5) / n_thick
    x, y, z = np.meshgrid(xs * lx, ys * ly, z0 + zs * thickness, indexing="ij")
    return np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)


def pairwise_fz(points_a: np.ndarray, mass_a: float, points_b: np.ndarray, mass_b: float) -> float:
    """z-component of the Newtonian force on body b due to body a. Negative if b is pulled toward smaller z."""
    if mass_a == 0.0 or mass_b == 0.0 or len(points_a) == 0 or len(points_b) == 0:
        return 0.0
    dm_a = mass_a / len(points_a)
    dm_b = mass_b / len(points_b)
    fz = 0.0
    chunk = 300
    for i0 in range(0, len(points_b), chunk):
        sl = points_b[i0 : i0 + chunk]
        d = sl[:, None, :] - points_a[None, :, :]
        r2 = np.einsum("ijk,ijk->ij", d, d)
        r3 = r2 * np.sqrt(r2)
        fz += float(np.sum(-G_NEWTON * dm_b * dm_a * d[:, :, 2] / r3))
    return fz


def qin_geometry() -> dict:
    """Wafer geometry from Qin et al., Nat. Commun. (2026), 10.1038/s41467-026-78187-4.

    The published 471 aN figure is the 16ω_d torque divided by the 23.694 mm arm, not the
    static attraction of one face-on pair. Geometry integrals below are the static pair.
    """
    lx, ly, lz = 2.598e-3, 2.597e-3, 0.487e-3
    mass = 7.6e-6
    volume = lx * ly * lz
    density = mass / volume
    gaps = (0.654e-3, 1.554e-3, 50.0e-3)
    skin = 50.0e-6
    resolutions = ((8, 4), (12, 6), (16, 8))
    finest = resolutions[-1]
    rows = []
    for gap in gaps:
        center = gap + lz
        point = G_NEWTON * mass * mass / center**2
        forces = []
        skin_forces = []
        for n_lat, n_thick in resolutions:
            test = wafer_grid(n_lat, n_thick, 0.0, lz, lx, ly)
            source = wafer_grid(n_lat, n_thick, lz + gap, lz, lx, ly)
            fz = pairwise_fz(test, mass, source, mass)
            attractive = -fz
            skin_pts = wafer_grid(n_lat, max(3, n_thick // 2), lz + gap, skin, lx, ly)
            skin_mass = mass * (skin / lz)
            fz_skin = pairwise_fz(test, mass, skin_pts, skin_mass)
            forces.append(attractive)
            skin_forces.append(-fz_skin)
        # Last interval as a plain convergence estimate. Not a Richardson proof.
        rel = abs(forces[-1] - forces[-2]) / abs(forces[-1])
        rel_skin = abs(skin_forces[-1] - skin_forces[-2]) / abs(skin_forces[-1])
        rows.append(
            {
                "gap_m": gap,
                "center_separation_m": center,
                "point_force_N": point,
                "extended_force_N": forces[-1],
                "extended_over_point": forces[-1] / point,
                "resolution_forces_N": forces,
                "relative_change_last_two_grids": rel,
                "skin_thickness_m": skin,
                "skin_mass_fraction": skin / lz,
                "skin_force_N": skin_forces[-1],
                "skin_force_over_extended": skin_forces[-1] / forces[-1],
                "skin_resolution_forces_N": skin_forces,
                "skin_relative_change_last_two_grids": rel_skin,
            }
        )
    # Linearity: a uniform density rescaling of the source must rescale the force by the same factor.
    gap = 0.654e-3
    n_lat, n_thick = finest
    test = wafer_grid(n_lat, n_thick, 0.0, lz, lx, ly)
    source = wafer_grid(n_lat, n_thick, lz + gap, lz, lx, ly)
    f1 = -pairwise_fz(test, mass, source, mass)
    f2 = -pairwise_fz(test, mass, source, mass * 1.02)
    # Distance ratio of the extended force versus a skin-only extra piece.
    f_near = rows[0]["extended_force_N"]
    f_far = rows[1]["extended_force_N"]
    s_near = rows[0]["skin_force_N"]
    s_far = rows[1]["skin_force_N"]
    torque_total_aNm = 11.17
    spoke_aNm = 0.53
    lever_m = 23.694e-3
    equivalent_force_N = torque_total_aNm * 1e-18 / lever_m
    return {
        "source": "Qin et al., Nature Communications (2026), DOI 10.1038/s41467-026-78187-4",
        "wafer_m": [lx, ly, lz],
        "mass_kg": mass,
        "volume_m3": volume,
        "density_kg_m3": density,
        "published_equivalent_force_aN": 471.0,
        "published_force_uncertainty_aN_1sigma": 23.0,
        "published_G_fractional_uncertainty": 0.053,
        "torque_cross_check_aNm": torque_total_aNm,
        "lever_arm_m": lever_m,
        "torque_over_lever_aN": equivalent_force_N * 1e18,
        "spoke_torque_fraction": spoke_aNm / torque_total_aNm,
        "gaps": rows,
        "uniform_source_rescaling": {
            "factor": 1.02,
            "force_ratio": f2 / f1,
        },
        "distance_ratios_0.654_to_1.554": {
            "extended": f_far / f_near,
            "point": rows[1]["point_force_N"] / rows[0]["point_force_N"],
            "skin_over_its_near_value": s_far / s_near,
            "skin_fraction_of_extended_near": s_near / f_near,
            "skin_fraction_of_extended_far": s_far / f_far,
        },
        "grids": resolutions,
    }


def qin_bound(density: float, delta_c: float = 1.0) -> dict:
    force = 471.0e-18
    sigma = 23.0e-18
    fractional = sigma / force
    # Protocol: swap the source arrangement, leave the test wafer unchanged.
    # Equilibrium reading of an ordinary mass shift would be matched away.
    # The number below bounds the assumed active-mass coefficient that is NOT included in the matched inertial mass.
    alpha_over_rho_c2 = fractional / abs(delta_c)
    alpha_j_per_m3 = alpha_over_rho_c2 * density * C_LIGHT**2
    both_wafers = alpha_over_rho_c2 / 2.0
    # External chemical scale, for comparison only. Not an output of the lattice model.
    cohesive_eV = 4.63
    silicon_u = 28.085
    joule_per_eV = float(constants.e)
    atom_mass = silicon_u * float(constants.atomic_mass)
    fractional_chemical = (cohesive_eV * joule_per_eV) / (atom_mass * C_LIGHT**2)
    return {
        "protocol": "source wafer swapped, test wafer unchanged, geometric kernel identical to the measured channel",
        "delta_C2": delta_c,
        "fractional_force_1sigma": fractional,
        "alpha_over_rho_c2_1sigma": alpha_over_rho_c2,
        "alpha_J_per_m3_1sigma": alpha_j_per_m3,
        "alpha_over_rho_c2_if_both_wafers_change": both_wafers,
        "G_inference_fractional_uncertainty": 0.053,
        "chemical_binding_fraction_Si_external_reference": fractional_chemical,
        "chemical_reference": "4.63 eV/atom cohesive energy over 28.085 u c²; comparison scale only, not a model input",
        "force_resolution_aN_for_chemical_fraction": fractional_chemical * 471.0,
        "assumption": (
            "Bound applies only if |ΔC2| is maintained at matched ordinary inertial mass and the toy "
            "coefficient contributes to active gravitational mass. That split is not a consequence of "
            "the Candidate-0 action and does not embed in one conserved minimally coupled stress tensor."
        ),
    }


def time_average_note(variance: float, tau_s: float, duration_s: float) -> dict:
    """What averaging does IF the correlator is a stationary exponential process.

    The Hamiltonian calculation is static. This factor is not applied to the mean shift.
    For an exponential autocorrelation exp(-|t|/τ), the variance of the average over T≫τ
    is approximately (2τ/T) times the instantaneous variance.
    """
    factor = 2.0 * tau_s / duration_s
    return {
        "assumed_correlation_time_s": tau_s,
        "assumed_duration_s": duration_s,
        "variance_reduction_factor_2tau_over_T": factor,
        "rms_reduction_factor": float(np.sqrt(factor)),
        "applied_to_mean_shift": False,
        "justification": "Not implied by the static Hamiltonian. Recorded only so the fluctuation formula is explicit.",
    }


def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items() if not str(k).startswith("_")}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, float):
        if not np.isfinite(value):
            raise RuntimeError("non-finite result")
        return value
    return value


def assert_close(name: str, got: float, expected: float, tol: float) -> None:
    if abs(got - expected) > tol:
        raise AssertionError(f"{name}: got {got}, expected {expected}, tol {tol}")


def main() -> None:
    amp_second = potential_second_derivative(V_VEV)
    amp_second_num = numerical_second_derivative(lambda r: potential(r), V_VEV)
    amp_at_zero = potential_second_derivative(0.0)
    assert_close("U''(v)", amp_second, 2.0 * LAMBDA_POT * V_VEV**2, 1e-12)
    assert_close("U''(v) numerical", amp_second_num, amp_second, 1e-6)
    if amp_at_zero >= 0.0 or amp_second <= 0.0:
        raise AssertionError("amplitude curvature at the vacuum or at zero has the wrong sign")

    ring = ring_arrangements()
    circ = ring["configurations"]["circulating"]
    anti = ring["configurations"]["anticirculating"]
    folded = ring["configurations"]["folded"]
    locked = ring["configurations"]["locked"]
    assert_close("circulating force", circ["max_abs_force"], 0.0, 1e-12)
    assert_close("anticirculating force", anti["max_abs_force"], 0.0, 1e-12)
    assert_close("locked force", locked["max_abs_force"], 0.0, 1e-12)
    if folded["max_abs_force"] < 1.0:
        raise AssertionError("folded arrangement should not be stationary")
    assert_close(
        "matched NN energy",
        circ["nearest_neighbor_energy"],
        folded["nearest_neighbor_energy"],
        1e-12,
    )
    if abs(circ["C2"] - folded["C2"]) < 0.5:
        raise AssertionError("C2 should distinguish the matched-energy pair")
    assert_close("opposite windings C2", circ["C2"], anti["C2"], 1e-12)
    assert_close("opposite windings energy", circ["nearest_neighbor_energy"], anti["nearest_neighbor_energy"], 1e-12)
    assert_close("winding +1", circ["winding"], 1.0, 1e-12)
    assert_close("winding -1", anti["winding"], -1.0, 1e-12)
    assert_close("phase-shift energy", circ["energy_under_global_phase_shift"], circ["nearest_neighbor_energy"], 1e-12)
    if ring["stable_critical_windings"]["n_negative_hessian_eigenvalues"] != 0:
        raise AssertionError("ring Hessian should be positive semidefinite")
    if ring["stable_critical_windings"]["n_zero_hessian_eigenvalues"] != 1:
        raise AssertionError("ring Hessian should have one global-phase zero")
    if abs(ring["long_wave_check"]["deviation_from_1"]) > 1e-6:
        raise AssertionError("long-wave C2 identity failed")

    orient = orientation_sector()
    assert_close("aligned energy", orient["aligned"]["energy"], 0.0, 1e-12)
    assert_close("neel energy", orient["neel"]["energy"], 24.0, 1e-12)
    assert_close("flip energy", orient["single_flip"]["energy"], 6.0, 1e-12)
    assert_close("aligned skyrme", orient["aligned"]["skyrme_density"], 0.0, 1e-12)
    assert_close("neel skyrme", orient["neel"]["skyrme_density"], 0.0, 1e-12)
    assert_close("rotation aligned", orient["rotation_invariance_aligned_abs_diff"], 0.0, 1e-12)
    assert_close("rotation neel", orient["rotation_invariance_neel_abs_diff"], 0.0, 1e-12)
    assert_close("neel curvature", orient["neel_tilt_curvature_over_eps2"], -24.0, 1e-4)
    if orient["vacuum_orientation_eigenvalue_abs_error"] > 1e-6:
        raise AssertionError("vacuum orientation Hessian does not match the hypercube Laplacian")
    vac_eigs = np.array(orient["vacuum_orientation_eigenvalues_analytic_per_transverse_component"])
    if np.sum(np.abs(vac_eigs) <= 1e-12) != 1 or np.any(vac_eigs < -1e-12):
        raise AssertionError("vacuum orientation spectrum is not a single zero plus positive modes")

    control = two_mass_multipole_check()
    assert_close("two-mass Qzz", control["quadrupole_zz"], control["quadrupole_zz_expected"], 1e-10)
    if abs(control["trace"]) > 1e-10 or control["relative_error_at_R100"] > 1e-7:
        raise AssertionError("multipole control failed")

    # Illustrative J comparable to rest energy, so the bookkeeping is visible. Not a silicon estimate.
    baseline = newtonian_baseline(orient, j_over_rest=1.0)
    if baseline["aligned"]["quadrupole_frobenius"] > 1e-10:
        raise AssertionError("aligned quadrupole should vanish")
    if baseline["neel"]["quadrupole_frobenius"] > 1e-10:
        raise AssertionError("Néel excess is uniform on the cube and should have no quadrupole")
    if baseline["single_flip"]["quadrupole_frobenius"] <= 1e-6:
        raise AssertionError("localized flip should carry a quadrupole about the center of mass")
    if abs(baseline["flip_smeared_same_total_mass"]["mass_minus_flip"]) > 1e-10:
        raise AssertionError("smeared comparison is not mass-matched")
    fr = baseline["flip_versus_smeared_fractional_potential"]
    # Common-origin residual is a dipole: fractional potential ~ 1/R, so R=40 over R=10 is ~1/4.
    dipole_ratio = abs(fr[1]["fractional_difference_common_origin"] / fr[0]["fractional_difference_common_origin"])
    # After the centers of mass coincide, the residual is a quadrupole: fractional ~ 1/R², ratio ~1/16.
    quad_ratio = abs(fr[1]["fractional_difference_cm_matched"] / fr[0]["fractional_difference_cm_matched"])
    if not 0.20 < dipole_ratio < 0.30:
        raise AssertionError(f"dipole falloff ratio {dipole_ratio} is not near 1/4")
    if not 0.04 < quad_ratio < 0.09:
        raise AssertionError(f"quadrupole falloff ratio {quad_ratio} is not near 1/16")

    stats = phase_order_statistics()
    mc = stats["monte_carlo"]
    if abs(mc["Phi_mean"] - mc["Phi_mean_exact"]) > 0.01:
        raise AssertionError("sample-spanning Monte Carlo mean disagrees with 1/N")
    if abs(mc["short_range_mean"]) > 0.02:
        raise AssertionError("short-range Monte Carlo mean should vanish")

    geometry = qin_geometry()
    if abs(geometry["uniform_source_rescaling"]["force_ratio"] - 1.02) > 1e-9:
        raise AssertionError("uniform density rescaling did not pass through the force linearly")
    near = geometry["gaps"][0]
    far = geometry["gaps"][1]
    wide = geometry["gaps"][2]
    if wide["extended_over_point"] < 0.98 or wide["extended_over_point"] > 1.02:
        raise AssertionError("extended wafers should approach the point-mass force at 50 mm")
    if abs(near["skin_force_over_extended"] - near["skin_mass_fraction"]) < 1e-3:
        raise AssertionError("skin support should not reproduce the bulk mass fraction at 0.654 mm")
    ratios = geometry["distance_ratios_0.654_to_1.554"]
    if abs(ratios["skin_fraction_of_extended_near"] - ratios["skin_fraction_of_extended_far"]) < 1e-3:
        raise AssertionError("skin fraction of the force should change with gap")

    bound = qin_bound(geometry["density_kg_m3"])
    # Fluctuation formula only. τ and T are assumptions, not outputs of the Hamiltonian.
    fluctuation = time_average_note(mc["short_range_rms"] ** 2, tau_s=1.0, duration_s=100.0 * 3600.0)
    scale = scaling_table(beta=1.0)

    flags = {
        "action_defines_extra_organization_coupling": False,
        "toy_coefficient_is_assumed": True,
        "stable_matched_pair_split_by_toy_C2": False,
        "unstable_folded_pair_split_by_toy_C2": True,
        "C2_double_counts_nearest_neighbor_energy_on_matched_pair": False,
        "C2_repeats_gradient_energy_in_long_wave_limit": True,
        "nonlocal_C2_embeds_in_one_conserved_stress_tensor": False,
        "signal_survives_matching_of_full_active_mass": False,
        "one_over_sqrt_N_used_as_mean_shift": False,
        "earth_or_moon_interior_used": False,
        "rho_plus_3p_used_as_equilibrium_wafer_prediction": False,
    }

    conclusion = {
        "stable_arrangements_differ_only_by_ordinary_mass_energy": True,
        "additional_signal_among_stable_states_of_this_hamiltonian": False,
        "additional_signal_if_the_assumed_coupling_acts_on_a_held_C2_difference": True,
        "plain_statement": (
            "On the lattices computed here, arrangements that are stable under the repository "
            "Hamiltonian differ gravitationally only through ordinary stress-energy. Opposite "
            "windings match in energy, multipoles, and the even pair correlator. The Néel "
            "arrangement is stationary but unstable, and its far field differs from the vacuum "
            "by ΔM = ΔE/c² with vanishing extra multipoles. The one toy coupling can split a "
            "same-stress pair, but that pair is not stable, and the coupling is not part of the "
            "action. Matched full active mass leaves no residual."
        ),
    }

    ring_public = jsonable(ring)
    orient_public = {
        "n_sites": orient["n_sites"],
        "n_edges": orient["n_edges"],
        "aligned_energy": orient["aligned"]["energy"],
        "neel_energy": orient["neel"]["energy"],
        "single_flip_energy": orient["single_flip"]["energy"],
        "aligned_skyrme": orient["aligned"]["skyrme_density"],
        "neel_skyrme": orient["neel"]["skyrme_density"],
        "mean_neighbor_dot_aligned": orient["aligned"]["mean_neighbor_dot"],
        "mean_neighbor_dot_neel": orient["neel"]["mean_neighbor_dot"],
        "rotation_invariance_aligned_abs_diff": orient["rotation_invariance_aligned_abs_diff"],
        "rotation_invariance_neel_abs_diff": orient["rotation_invariance_neel_abs_diff"],
        "neel_tilt_curvature_over_eps2": orient["neel_tilt_curvature_over_eps2"],
        "vacuum_orientation_eigenvalues": orient["vacuum_orientation_eigenvalues_analytic_per_transverse_component"],
        "vacuum_orientation_eigenvalue_abs_error": orient["vacuum_orientation_eigenvalue_abs_error"],
        "amplitude_U_second_at_v": amp_second,
        "amplitude_U_second_at_0": amp_at_zero,
        "amplitude_gap_omega2_minus_k2": amp_second,
        "kappa_skyrme_recorded": KAPPA_SKYRME,
        "note": (
            "Skyrme density built from n·(Δn×Δn) is identically zero on both collinear "
            "arrangements, so it does not change ΔE. Vacuum stability of the two-derivative "
            "term is the hypercube Laplacian. The Néel curvature -24 is the two-derivative term."
        ),
    }

    payload = {
        "model": {
            "name": "Candidate 0 static energy, lattice transcription",
            "potential": "λ(ρ²-v²)²/4",
            "v": V_VEV,
            "lambda": LAMBDA_POT,
            "kappa": KAPPA_SKYRME,
            "gradient_normalization": (
                "E = J h with J = v² a, h_θ = Σ_bonds (1/2) wrap(Δθ)², "
                "h_n = Σ_bonds (1 - n_i·n_j). This matches ∫(1/2)ρ²(|∇θ|²+|∇n|²) "
                "on a cubic lattice of spacing a at ρ=v."
            ),
            "metric": "imposed flat background; no dynamical metric in the matter action",
            "observable_phase": "cos(θ_i-θ_j), not the global phase",
            "observable_orientation": "n_i·n_j, not a global rotation",
        },
        "ring": ring_public,
        "orientation": orient_public,
        "newtonian_baseline": jsonable(baseline),
        "multipole_control": jsonable(control),
        "statistics": jsonable(stats),
        "scaling_at_fixed_density_beta_1": scale,
        "fluctuation_time_average_assumption": fluctuation,
        "qin_geometry": jsonable(geometry),
        "qin_bound": jsonable(bound),
        "flags": flags,
        "conclusion": conclusion,
        "constants": {"c_m_per_s": C_LIGHT, "G_SI": G_NEWTON},
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    out_path = OUT / "results.json"
    out_path.write_text(json.dumps(payload, indent=2) + "\n")

    print("Matter organization and Newtonian pull")
    print(f"python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}")
    print(f"amplitude U''(v) = {amp_second:.12g}; U''(0) = {amp_at_zero:.12g}")
    print(
        "ring: locked E,C2 = "
        f"{locked['nearest_neighbor_energy']:.12g}, {locked['C2']:.12g}; "
        f"circulating E,C2,w = {circ['nearest_neighbor_energy']:.12g}, {circ['C2']:.12g}, {circ['winding']:.12g}; "
        f"anticirculating E,C2,w = {anti['nearest_neighbor_energy']:.12g}, {anti['C2']:.12g}, {anti['winding']:.12g}; "
        f"folded E,C2,max|F| = {folded['nearest_neighbor_energy']:.12g}, {folded['C2']:.12g}, {folded['max_abs_force']:.12g}"
    )
    print(
        "ring Hessian eigenvalues "
        f"{ring['stable_critical_windings']['hessian_eigenvalues']}; "
        f"negative {ring['stable_critical_windings']['n_negative_hessian_eigenvalues']}, "
        f"zero {ring['stable_critical_windings']['n_zero_hessian_eigenvalues']}"
    )
    print(
        "cube orientation energies aligned, Neel, flip = "
        f"{orient['aligned']['energy']:.12g}, {orient['neel']['energy']:.12g}, {orient['single_flip']['energy']:.12g}"
    )
    print(
        "Neel tilt curvature / eps^2 = "
        f"{orient['neel_tilt_curvature_over_eps2']:.12g} (expected -24); "
        f"vacuum Hessian max abs error {orient['vacuum_orientation_eigenvalue_abs_error']:.3e}"
    )
    print(
        "baseline masses aligned, Neel, flip = "
        f"{baseline['aligned']['mass']:.12g}, {baseline['neel']['mass']:.12g}, {baseline['single_flip']['mass']:.12g}"
    )
    print(
        "baseline quadrupoles aligned, Neel, flip = "
        f"{baseline['aligned']['quadrupole_frobenius']:.3e}, "
        f"{baseline['neel']['quadrupole_frobenius']:.3e}, "
        f"{baseline['single_flip']['quadrupole_frobenius']:.12g}"
    )
    print(
        "ΔM(Neel-aligned) = "
        f"{baseline['delta_M_from_delta_E_over_c2']:.12g}; "
        f"fractional monopole {baseline['fractional_monopole_neel_vs_aligned']:.12g}"
    )
    for row in fr:
        print(
            f"flip vs smeared at R={row['R']}: common-origin {row['fractional_difference_common_origin']:.12g}, "
            f"CM-matched {row['fractional_difference_cm_matched']:.12g}"
        )
    print(
        "random-phase Monte Carlo: "
        f"Phi mean {mc['Phi_mean']:.6g} (exact {mc['Phi_mean_exact']:.6g}), "
        f"rms {mc['Phi_rms']:.6g} (exact {mc['Phi_rms_exact']:.6g}); "
        f"short-range mean {mc['short_range_mean']:.6g}, rms {mc['short_range_rms']:.6g}, "
        f"iid-bond rms estimate {mc['iid_bond_rms_estimate']:.6g}"
    )
    print(
        f"Qin density {geometry['density_kg_m3']:.6g} kg/m^3; "
        f"torque/lever {geometry['torque_over_lever_aN']:.6g} aN; "
        f"spoke fraction {geometry['spoke_torque_fraction']:.6g}"
    )
    for row in geometry["gaps"]:
        print(
            f"gap {row['gap_m']*1e3:.6g} mm: point {row['point_force_N']*1e18:.6g} aN, "
            f"extended {row['extended_force_N']*1e18:.6g} aN, "
            f"extended/point {row['extended_over_point']:.6g}, "
            f"grid change {row['relative_change_last_two_grids']:.3e}, "
            f"skin/extended {row['skin_force_over_extended']:.6g} "
            f"(mass fraction {row['skin_mass_fraction']:.6g})"
        )
    print(
        "distance ratios 1.554/0.654 mm: "
        f"extended {ratios['extended']:.6g}, point {ratios['point']:.6g}, "
        f"skin fraction near {ratios['skin_fraction_of_extended_near']:.6g}, "
        f"far {ratios['skin_fraction_of_extended_far']:.6g}"
    )
    print(
        "uniform source rescaling 1.02 -> force ratio "
        f"{geometry['uniform_source_rescaling']['force_ratio']:.12g}"
    )
    print(
        "toy bound |alpha|/(rho c^2) at |ΔC2|=1, source swap, 1σ: "
        f"{bound['alpha_over_rho_c2_1sigma']:.6g}; "
        f"|alpha| {bound['alpha_J_per_m3_1sigma']:.6g} J/m^3; "
        f"both-wafers reading {bound['alpha_over_rho_c2_if_both_wafers_change']:.6g}"
    )
    print(
        "external Si cohesive fraction "
        f"{bound['chemical_binding_fraction_Si_external_reference']:.6g}; "
        f"force resolution for that fraction {bound['force_resolution_aN_for_chemical_fraction']:.6g} aN"
    )
    for row in scale:
        print(
            f"N={row['N']}: short-range mean {row['short_range_mean_shift_if_DeltaC_is_1']:.6g}, "
            f"locked spanning {row['sample_spanning_locked_mean_shift']:.6g}, "
            f"random spanning mean {row['sample_spanning_random_mean_shift']:.6g}, "
            f"random spanning rms {row['sample_spanning_random_rms']:.6g}"
        )
    print(
        "assumed fluctuation average: variance factor "
        f"{fluctuation['variance_reduction_factor_2tau_over_T']:.6g} "
        "(not applied to the mean)"
    )
    print(conclusion["plain_statement"])
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
