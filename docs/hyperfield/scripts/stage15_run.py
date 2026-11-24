#!/usr/bin/env python3
"""
Stage 15 — Oppenheim classical-gravity decoherence–diffusion bound tracker.

Reads (or, if missing, writes) results/stage15_datapoints.csv, recomputes the
three-model-class tradeoff of Oppenheim, Sparaciari, Soda, Weller-Davies,
Nat. Commun. 14, 7910 (2023), with the Janse et al. Phys. Rev. Research 6,
033076 (2024) force-noise compilation and later coherence data, and writes
results/stage15.json.

Usage:
    python scripts/stage15_run.py              # read CSV if present, else seed it
    python scripts/stage15_run.py --write-csv  # refresh CSV from stage15_data.py
    python scripts/stage15_fig.py              # after a run, replot

To add a new experiment later: append a row to results/stage15_datapoints.csv
(see the header and HYPERFIELD_RESULTS.md §15.6) and re-run both scripts.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from stage15_data import DATAPOINTS  # noqa: E402

CSV_PATH = os.path.join(ROOT, "results/stage15_datapoints.csv")
JSON_PATH = os.path.join(ROOT, "results/stage15.json")

# ---- constants (SI) ----
G = 6.67430e-11          # CODATA
U = 1.66053906660e-27    # atomic mass unit
M_N = 1.67262192369e-27  # proton mass; Oppenheim uses ~1.67e-27
R_N = 1.0e-15            # Oppenheim 2023 "very conservative" nuclear radius, m
R_E = 6.371e6            # Earth radius, m
H_ATM = 1.0e4            # ASSUMED "atmospheric height" ~ tropopause, m
V_B_ATM = math.pi * R_E**2 * H_ATM
# V_b that reproduces the published ultra-local-continuous upper bound
# D2 <= 1e-41 kg^2 s m^{-3} from FOM=1e14 and r_N=1e-15 (see writeup).
V_B_TABLE1 = (1.0e14 * R_N**4) / (1.0e-41 * G**2)

# Conversion FOM -> D2 (Oppenheim Eqs. 44, 46, 47; Janse Eqs. 1–3)
# FOM = sigma_a^2 * N * Delta T = S_a * N, units m^2 s^{-3}
K_DISC = R_N**4 / (M_N * G**2)          # (l_P^3/m_P) D2  [kg]  = K_DISC * FOM
K_NL = R_N**3 / G**2                    # l_P^2 D2 [kg s m^{-1}] = K_NL * FOM
K_CONT_ATM = R_N**4 / (V_B_ATM * G**2)  # D2 [kg^2 s m^{-3}]
K_CONT_TAB = R_N**4 / (V_B_TABLE1 * G**2)

CSV_FIELDS = [
    "id", "role", "year", "short_name", "citation", "experiment_type",
    "element", "mass_kg", "N", "n_constituents", "mass_amu",
    "f0_Hz", "sqrt_SF", "sqrt_Sa", "FOM", "tau_s", "lambda_Hz",
    "delocalization_m", "applicability", "caveat", "source",
]


def _fmt(x):
    if x is None or x == "":
        return ""
    if isinstance(x, float):
        return f"{x:.6e}"
    return str(x)


def write_csv(path=CSV_PATH, rows=DATAPOINTS):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            out = {k: "" for k in CSV_FIELDS}
            out.update({k: _fmt(r.get(k, "")) for k in CSV_FIELDS})
            # lambda from tau if missing
            if not out["lambda_Hz"] and r.get("tau_s"):
                out["lambda_Hz"] = _fmt(1.0 / float(r["tau_s"]))
            if not out["N"] and r.get("n_constituents"):
                out["N"] = _fmt(r["n_constituents"])
            if not out["n_constituents"] and r.get("N"):
                out["n_constituents"] = _fmt(r["N"])
            w.writerow(out)
    return path


def _num(x):
    if x is None or x == "":
        return None
    if isinstance(x, (int, float)):
        return float(x)
    return float(x)


def load_csv(path=CSV_PATH):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            for k in ("year", "mass_kg", "N", "n_constituents", "mass_amu",
                      "f0_Hz", "sqrt_SF", "sqrt_Sa", "FOM", "tau_s",
                      "lambda_Hz", "delocalization_m"):
                r[k] = _num(r.get(k, ""))
            r["year"] = int(r["year"]) if r["year"] is not None else None
            if r.get("lambda_Hz") is None and r.get("tau_s"):
                r["lambda_Hz"] = 1.0 / r["tau_s"]
            if r.get("tau_s") is None and r.get("lambda_Hz"):
                r["tau_s"] = 1.0 / r["lambda_Hz"]
            if r.get("n_constituents") is None and r.get("N") is not None:
                r["n_constituents"] = r["N"]
            rows.append(r)
    return rows


def volume_and_R(row):
    """Oppenheim convention: V = r_N^2 * Delta x, R = V^{1/3}."""
    dx = row.get("delocalization_m") or 1.0e-7
    V = R_N**2 * dx
    return V, V ** (1.0 / 3.0), dx


def lower_bounds(row):
    """D2 lower bounds from a coherence experiment (Oppenheim 44, 46, 47)."""
    M = row["mass_kg"]
    lam = row["lambda_Hz"]
    N = row.get("n_constituents") or row.get("N") or 1.0
    V, R, dx = volume_and_R(row)
    disc = M / lam                              # (l_P^3/m_P) D2, kg
    nl = N * M**2 / (R * lam)                   # l_P^2 D2, kg s / m
    cont = N * M**2 / (V * lam)                 # D2, kg^2 s / m^3
    return {
        "ultralocal_discrete_lP3_over_mP_D2_kg": disc,
        "nonlocal_lP2_D2_kg_s_per_m": nl,
        "ultralocal_continuous_D2_kg2_s_per_m3": cont,
        "V_m3": V, "R_m": R, "delocalization_m": dx,
        "M_kg": M, "lambda_Hz": lam, "N_constituents": N,
    }


def upper_bounds(row):
    """D2 upper bounds from a force-noise experiment (Janse 1–3 / Opp. 44–47)."""
    FOM = row["FOM"]
    return {
        "FOM_m2_per_s3": FOM,
        "ultralocal_discrete_lP3_over_mP_D2_kg": K_DISC * FOM,
        "nonlocal_lP2_D2_kg_s_per_m": K_NL * FOM,
        "ultralocal_continuous_D2_kg2_s_per_m3_atm": K_CONT_ATM * FOM,
        "ultralocal_continuous_D2_kg2_s_per_m3_table1_Vb": K_CONT_TAB * FOM,
    }


def running_best(points, year_key, value_key, want="min"):
    """Step series: at each year, the best value so far (min upper / max lower)."""
    ordered = sorted(points, key=lambda p: (p[year_key], p[value_key]))
    best = None
    series = []
    for p in ordered:
        v = p[value_key]
        if best is None:
            best = v
        elif want == "min":
            best = min(best, v)
        else:
            best = max(best, v)
        series.append({
            "year": p[year_key],
            "id": p["id"],
            "short_name": p.get("short_name", p["id"]),
            "value": v,
            "best_so_far": best,
            "applicability": p.get("applicability"),
        })
    return series


def gap_orders(lower, upper):
    if lower is None or upper is None or lower <= 0 or upper <= 0:
        return None
    if upper < lower:
        return {"orders": 0.0, "ratio_upper_over_lower": upper / lower, "status": "ruled_out"}
    return {"orders": math.log10(upper / lower), "ratio_upper_over_lower": upper / lower, "status": "open"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-csv", action="store_true",
                    help="refresh results/stage15_datapoints.csv from stage15_data.py")
    args = ap.parse_args()

    if args.write_csv or not os.path.exists(CSV_PATH):
        write_csv()
        print("wrote", CSV_PATH)
    rows = load_csv()

    computed = []
    for r in rows:
        item = {k: r[k] for k in r}
        if r["role"] in ("lower_coherence", "projection_lower"):
            item["bounds"] = lower_bounds(r)
            item["side"] = "lower"
        elif r["role"] in ("upper_diffusion", "projection_upper"):
            item["bounds"] = upper_bounds(r)
            item["side"] = "upper"
        else:
            continue
        computed.append(item)

    def sel(role=None, appl=None, year_max=None, ids=None):
        out = []
        for p in computed:
            if role and p["role"] != role:
                continue
            if appl:
                if isinstance(appl, (list, tuple, set)):
                    if p.get("applicability") not in appl:
                        continue
                elif p.get("applicability") != appl:
                    continue
            if year_max is not None and p["year"] > year_max:
                continue
            if ids is not None and p["id"] not in ids:
                continue
            out.append(p)
        return out

    # Oppenheim 2023 Table 1 reproduction from their Cavendish + Gerlich numbers
    cav = next(p for p in computed if p["id"] == "cavendish1798")
    ger = next(p for p in computed if p["id"] == "gerlich2011")
    table1_check = {
        "published": {
            "ultralocal_continuous_D2": {"upper": 1e-41, "lower": 1e-9, "units": "kg^2 s m^{-3}", "status": "ruled_out"},
            "ultralocal_discrete_lP3_over_mP_D2": {"upper": 1e-1, "lower": 1e-25, "units": "kg", "status": "open"},
            "nonlocal_lP2_D2": {"upper": 1e-9, "lower": 1e-35, "units": "kg s m^{-1}", "status": "open"},
        },
        "recomputed_from_stated_Cavendish_and_Gerlich": {
            "FOM_Cavendish": cav["bounds"]["FOM_m2_per_s3"],
            "ultralocal_continuous_D2_atm_Vb": {
                "upper": cav["bounds"]["ultralocal_continuous_D2_kg2_s_per_m3_atm"],
                "lower": ger["bounds"]["ultralocal_continuous_D2_kg2_s_per_m3"],
                "note": "V_b = pi R_E^2 * 10 km (ASSUMED tropopause). Upper is a few orders below the published 1e-41; either way the window is empty.",
            },
            "ultralocal_continuous_D2_table1_Vb": {
                "upper": cav["bounds"]["ultralocal_continuous_D2_kg2_s_per_m3_table1_Vb"],
                "lower": ger["bounds"]["ultralocal_continuous_D2_kg2_s_per_m3"],
                "V_b_m3": V_B_TABLE1,
                "implied_h_m": V_B_TABLE1 / (math.pi * R_E**2),
                "note": "V_b chosen to reproduce the published 1e-41 upper bound; the implied height is ~18 m, not an atmospheric scale. We keep the published Table 1 number as CITED.",
            },
            "ultralocal_discrete": {
                "upper": cav["bounds"]["ultralocal_discrete_lP3_over_mP_D2_kg"],
                "lower": ger["bounds"]["ultralocal_discrete_lP3_over_mP_D2_kg"],
                "note": "Formula+Cavendish numbers give upper ~10^{+1} kg; Table 1 quotes 10^{-1} (about 2 orders more conservative).",
            },
            "nonlocal": {
                "upper": cav["bounds"]["nonlocal_lP2_D2_kg_s_per_m"],
                "lower": ger["bounds"]["nonlocal_lP2_D2_kg_s_per_m"],
                "note": "Formula+Cavendish numbers give upper ~10^{-11}; Table 1 quotes 10^{-9} (about 2 orders more conservative).",
            },
        },
        "janse_published_2024_conservative_rounded": {
            "ultralocal_discrete": {"upper": 1e-16, "lower": 1e-25, "units": "kg", "note": "Janse Eqs. (4)–(5); 15-order FOM improvement applied to Table 1, excluding open circles."},
            "nonlocal": {"upper": 1e-24, "lower": 1e-35, "units": "kg s m^{-1}"},
        },
    }

    # Current bests
    lowers = sel(role="lower_coherence")
    uppers_cons = sel(role="upper_diffusion", appl="absolute")
    uppers_aggr = sel(role="upper_diffusion",
                      appl={"absolute", "differential", "differential_L1"})

    def best_lower(model_key):
        best = None
        for p in lowers:
            v = p["bounds"][model_key]
            if best is None or v > best["value"]:
                best = {"id": p["id"], "short_name": p["short_name"], "year": p["year"],
                        "value": v, "citation": p["citation"]}
        return best

    def best_upper(pool, model_key):
        best = None
        for p in pool:
            v = p["bounds"][model_key]
            if best is None or v < best["value"]:
                best = {"id": p["id"], "short_name": p["short_name"], "year": p["year"],
                        "value": v, "FOM": p["bounds"]["FOM_m2_per_s3"],
                        "applicability": p["applicability"], "citation": p["citation"]}
        return best

    models = {
        "ultralocal_continuous": {
            "parameter": "D2",
            "units": "kg^2 s m^{-3}",
            "lower_key": "ultralocal_continuous_D2_kg2_s_per_m3",
            "upper_key": "ultralocal_continuous_D2_kg2_s_per_m3_table1_Vb",
            "table1_upper": 1e-41, "table1_lower": 1e-9,
        },
        "ultralocal_discrete": {
            "parameter": "(l_P^3/m_P) D2",
            "units": "kg",
            "lower_key": "ultralocal_discrete_lP3_over_mP_D2_kg",
            "upper_key": "ultralocal_discrete_lP3_over_mP_D2_kg",
            "table1_upper": 1e-1, "table1_lower": 1e-25,
        },
        "nonlocal_DP_like": {
            "parameter": "l_P^2 D2",
            "units": "kg s m^{-1}",
            "lower_key": "nonlocal_lP2_D2_kg_s_per_m",
            "upper_key": "nonlocal_lP2_D2_kg_s_per_m",
            "table1_upper": 1e-9, "table1_lower": 1e-35,
        },
    }

    current = {}
    for name, m in models.items():
        L = best_lower(m["lower_key"])
        U_cons = best_upper(uppers_cons, m["upper_key"])
        U_aggr = best_upper(uppers_aggr, m["upper_key"])
        current[name] = {
            "parameter": m["parameter"],
            "units": m["units"],
            "table1_2023": {"upper": m["table1_upper"], "lower": m["table1_lower"],
                            "gap": gap_orders(m["table1_lower"], m["table1_upper"])},
            "lower_now": L,
            "upper_conservative_absolute": U_cons,
            "upper_aggressive_including_differential": U_aggr,
            "gap_conservative": gap_orders(L["value"], U_cons["value"]),
            "gap_aggressive": gap_orders(L["value"], U_aggr["value"]),
        }

    # Historical series for the figure (conservative = absolute uppers only)
    history = {}
    for name, m in models.items():
        Lpts = [{"id": p["id"], "year": p["year"], "short_name": p["short_name"],
                 "value": p["bounds"][m["lower_key"]], "applicability": p["applicability"]}
                for p in lowers]
        Upts = [{"id": p["id"], "year": p["year"], "short_name": p["short_name"],
                 "value": p["bounds"][m["upper_key"]], "applicability": p["applicability"]}
                for p in uppers_cons]
        Upts_all = [{"id": p["id"], "year": p["year"], "short_name": p["short_name"],
                     "value": p["bounds"][m["upper_key"]], "applicability": p["applicability"]}
                    for p in uppers_aggr]
        history[name] = {
            "lower_running_max": running_best(Lpts, "year", "value", want="max"),
            "upper_running_min_absolute": running_best(Upts, "year", "value", want="min"),
            "upper_running_min_including_differential": running_best(Upts_all, "year", "value", want="min"),
        }

    # Projections vs current conservative window
    projections = []
    for p in computed:
        if not p["role"].startswith("projection"):
            continue
        entry = {"id": p["id"], "year": p["year"], "short_name": p["short_name"],
                 "role": p["role"], "citation": p["citation"], "caveat": p.get("caveat", ""),
                 "bounds": p["bounds"]}
        for name, m in models.items():
            if p["role"] == "projection_lower":
                v = p["bounds"][m["lower_key"]]
                U = current[name]["upper_conservative_absolute"]["value"]
                entry[name] = {"projected_lower": v, "vs_current_conservative_upper": gap_orders(v, U)}
            else:
                v = p["bounds"][m["upper_key"]]
                L = current[name]["lower_now"]["value"]
                entry[name] = {"projected_upper": v, "vs_current_lower": gap_orders(L, v)}
        projections.append(entry)

    # What would close each window
    def close_from_below(name):
        U = current[name]["upper_conservative_absolute"]["value"]
        L = current[name]["lower_now"]["value"]
        return {
            "need_lower_at_least": U,
            "current_lower": L,
            "factor_needed": (U / L) if L else None,
            "orders_needed": math.log10(U / L) if L and U > L else 0.0,
        }

    def close_from_above(name):
        L = current[name]["lower_now"]["value"]
        U = current[name]["upper_conservative_absolute"]["value"]
        FOM_now = current[name]["upper_conservative_absolute"]["FOM"]
        # FOM_close = L / K  so that K*FOM = L
        if name == "ultralocal_discrete":
            FOM_close = L / K_DISC
        elif name == "nonlocal_DP_like":
            FOM_close = L / K_NL
        else:
            FOM_close = L / K_CONT_TAB
        return {
            "need_upper_at_most": L,
            "current_upper": U,
            "current_FOM": FOM_now,
            "FOM_needed": FOM_close,
            "FOM_improvement_orders": math.log10(FOM_now / FOM_close) if FOM_close and FOM_now > FOM_close else 0.0,
        }

    closing = {name: {"from_below_coherence": close_from_below(name),
                      "from_above_force_noise": close_from_above(name)}
               for name in models}

    # Follow-up literature (not re-parameterized into the 2023 D2 units unless stated)
    followups = {
        "janse2024": {
            "citation": "Janse, Uitenbroek, van Everdingen, Plugge, Hensen, Oosterkamp, Phys. Rev. Research 6, 033076 (2024), arXiv:2403.08912",
            "role": "Tightens the upper side of the 2023 Newtonian D2 window by compiling low force-noise experiments. Conservative (absolute, on Earth) improvement of at least 15 orders of FOM relative to Cavendish. Aggressive inclusion of Asenbaum 2017 would drop below the discrete lower bound and leave ~1 order (on the 2011 floor) for the non-local model, but those points are differential.",
        },
        "layton2023_2024": {
            "citation": "Layton, Oppenheim, Russo, Weller-Davies, JHEP 08, 163 (2023); Layton, Oppenheim, Weller-Davies, Quantum 8, 1565 (2024)",
            "role": "Theory of the weak-field / 'healthier' semi-classical CQ dynamics. No new experimental D2 number.",
        },
        "grudka2024": {
            "citation": "Grudka, Morris, Oppenheim, Russo, Sajjad, arXiv:2402.17844 (2024)",
            "role": "Formal renormalisability of pure postquantum-classical gravity via a quadratic-gravity path integral. Relativistic corrections alter the gravitational two-point function and give milder stochastic acceleration fluctuations than a naive Newtonian extrapolation. They argue current upper and lower bounds are already tight enough that tabletop tests of the quantum vs classical nature of spacetime may be possible in the near term — a caveat on the Newtonian D2 mapping used in this tracker.",
        },
        "angeli2025": {
            "citation": "Angeli, Donadi, Di Bartolomeo, Gaona-Reyes, Vinante, Bassi, arXiv:2501.13030 (2025)",
            "role": "If gravity is an LOCC channel it must inject a minimum diffusion into the motion of quantum systems. They propose a millikelvin torsion-pendulum protocol to look for that diffusion without preparing a macroscopic superposition. Related in spirit to the upper side of the tradeoff, but the bound is setup-dependent and is not the same parameter as Oppenheim's D2.",
        },
        "oppenheim_sajjad2026": {
            "citation": "Oppenheim and Sajjad, arXiv:2605.05375 (2026)",
            "role": "Linearised relativistic theory: dimensionless diffusion constants for scalar and tensor channels. Matter-wave floor still quoted as D ≳ 10^{-71} from the 2011 molecular-superposition analysis; LVK tensor ceiling D ≲ 3.3e-67; window 10^{-71} ≲ D_{2,tensor} ≲ 10^{-67} (~4.5 orders). Full-positivity scalar ceiling D ≲ 1.3e-66 (~5.1 orders). LISA Pathfinder constrains the scalar channel under transverse positivity only. These are dimensionless couplings of the relativistic theory, not the dimensionful Newtonian D2 of Table 1; we report them separately and do not mix units.",
        },
        "carney2026": {
            "citation": "Fabiano, Fujita, Matsumura, Carney, arXiv:2603.26075 (2026)",
            "role": "Classification of non-quantized, Galilei-invariant Newtonian models. Any non-entangling model must inject a quantifiable minimal noise; measuring below that threshold would be equivalent to showing that Newtonian gravity is entangling. Does not directly replace the 2023 D2 windows.",
        },
        "bild2023_not_used": {
            "citation": "Bild et al., Science 380, 274 (2023)",
            "role": "16 µg mechanical-oscillator cat, decoherence ~1 MHz (Janse). A phonon cat is not a centre-of-mass mass-density superposition of the kind entering Oppenheim's M/lambda bound, so we do not promote it to a D2 floor.",
        },
        "wang2021_not_used": {
            "citation": "Wang et al., Nat. Commun. 12, 233 (2021)",
            "role": "Hour-scale ion-qubit coherence is a spin superposition, not a mass-density superposition (Janse correctly set it aside).",
        },
    }

    constants = {
        "G": G, "u": U, "m_N_kg": M_N, "r_N_m": R_N,
        "R_E_m": R_E, "h_atm_m_ASSUMED": H_ATM,
        "V_b_atmosphere_m3": V_B_ATM,
        "V_b_reproduces_table1_1e-41_m3": V_B_TABLE1,
        "K_disc_kg_s3_per_m2": K_DISC,
        "K_nl_kg_s4_per_m3": K_NL,
        "K_cont_table1_kg2_s4_per_m5": K_CONT_TAB,
        "notes": {
            "r_N": "CITED Oppenheim 2023, r_N ~ 10^{-15} m (very conservative nuclear radius).",
            "V": "DERIVED Oppenheim convention V = r_N^2 * Delta x (nuclear size in y,z; delocalization in x). The smearing-length vs wavepacket-volume choice can move the continuous and non-local lower bounds by many orders; the discrete lower bound M/lambda does not depend on V.",
            "FOM": "CITED Janse: FOM_D2 = sigma_a^2 N Delta T = S_a N, units m^2 s^{-3}.",
            "Table1_vs_formula": "Plugging the stated Cavendish FOM into Eqs. (46)–(47) gives uppers ~2 orders weaker (larger) than the rounded Table 1 values. We use the formulas as DERIVED and quote Table 1 / Janse rounded windows as CITED.",
        },
    }

    verdict = (
        "The ultra-local continuous model is ruled out (Oppenheim 2023 Table 1: "
        "upper 10^{-41} < lower 10^{-9} kg^2 s m^{-3}) and stays ruled out. "
        "The ultra-local discrete window, using formula bounds, Pedalino 2025 for the floor "
        "and Gisler 2022 for a conservative absolute ceiling, is still open by about 10 orders "
        "of magnitude; Janse's rounded conservative window was 9 orders on the 2011 floor. "
        "The non-local / Diosi-Penrose-like window is much tighter: about 3–4 orders remain "
        "conservatively (Pedalino floor vs Gisler ceiling). Including Asenbaum 2017 as if it "
        "were an absolute acceleration measurement would already close both remaining windows, "
        "but that point is differential and its applicability is questioned (Janse 2024). "
        "A 10^6 amu, 0.1 s interferometer would essentially close the non-local conservative "
        "window from below; closing the discrete conservative window from below wants "
        "MAQRO-scale mass and time (several 10^{10} amu, ~100 s), or a ~10-order improvement "
        "in an absolute FOM. The 2026 Oppenheim–Sajjad dimensionless relativistic window is "
        "a separate accounting (~4.5 orders on the tensor channel) and still uses the 2011 "
        "matter-wave floor."
    )

    # JSON-friendly computed rows (drop huge caveats already stored)
    points_out = []
    for p in computed:
        points_out.append({
            "id": p["id"], "role": p["role"], "year": p["year"],
            "short_name": p["short_name"], "citation": p["citation"],
            "experiment_type": p["experiment_type"], "element": p.get("element"),
            "mass_kg": p.get("mass_kg"), "applicability": p.get("applicability"),
            "FOM": (p.get("bounds") or {}).get("FOM_m2_per_s3"),
            "bounds": p.get("bounds"),
            "caveat": p.get("caveat"),
            "source": p.get("source"),
        })

    out = {
        "constants": constants,
        "table1_reproduction": table1_check,
        "current_windows": current,
        "history": history,
        "projections": projections,
        "closing": closing,
        "followups_2024_2026": followups,
        "datapoints": points_out,
        "verdict_one_liner": verdict,
    }

    def default(o):
        if isinstance(o, float) and (math.isnan(o) or math.isinf(o)):
            return None
        raise TypeError(type(o))

    os.makedirs(os.path.dirname(JSON_PATH), exist_ok=True)
    with open(JSON_PATH, "w") as f:
        json.dump(out, f, indent=2, default=default)
    print("wrote", JSON_PATH)

    # stdout summary
    print("\n--- Table 1 reproduction ---")
    rec = table1_check["recomputed_from_stated_Cavendish_and_Gerlich"]
    print("Cavendish FOM", rec["FOM_Cavendish"])
    print("discrete upper formula", rec["ultralocal_discrete"]["upper"], "lower", rec["ultralocal_discrete"]["lower"])
    print("non-local upper formula", rec["nonlocal"]["upper"], "lower", rec["nonlocal"]["lower"])
    print("continuous upper atm", rec["ultralocal_continuous_D2_atm_Vb"]["upper"],
          "table1-Vb", rec["ultralocal_continuous_D2_table1_Vb"]["upper"],
          "lower", rec["ultralocal_continuous_D2_atm_Vb"]["lower"])
    print("\n--- Current windows ---")
    for name, w in current.items():
        print(name, "units", w["units"])
        print("  lower", w["lower_now"])
        print("  upper conservative", w["upper_conservative_absolute"])
        print("  gap conservative", w["gap_conservative"])
        print("  upper aggressive", w["upper_aggressive_including_differential"])
        print("  gap aggressive", w["gap_aggressive"])
    print("\n--- Closing (conservative) ---")
    for name, c in closing.items():
        print(name, c)
    print("\nVERDICT:", verdict)


if __name__ == "__main__":
    main()
