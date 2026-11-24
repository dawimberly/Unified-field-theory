#!/usr/bin/env python3
"""
Stage 17 — Survey of existing and planned matter-wave / levitated / force-noise
experiments against the Oppenheim (2023) classical-gravity windows.

Reuses the Stage-15 formulae (scripts/stage15_run.py: lower_bounds, upper_bounds,
load_csv, gap_orders) and the shared tracker CSV results/stage15_datapoints.csv.
Rows whose caveat starts with "Stage 17." were added in this stage.

Selection rules (no tuning, fixed before looking at the numbers):
  conservative floor : role == lower_coherence (centre-of-mass superposition
                       larger than the object, interference/coherence verified)
  conservative ceiling: role == upper_diffusion and applicability == absolute
  flagged            : role == lower_coherence_flagged (Delta x below object size
                       or quantum/classical not discriminated) -> reported, NOT used
Writes results/stage17.json.
"""
from __future__ import annotations
import json, math, os, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from stage15_run import load_csv, lower_bounds, upper_bounds, gap_orders, R_N  # noqa: E402

JSON_PATH = os.path.join(ROOT, "results/stage17.json")
U = 1.66053906660e-27

MODELS = {
    "ultralocal_discrete": {"key_lo": "ultralocal_discrete_lP3_over_mP_D2_kg",
                            "key_up": "ultralocal_discrete_lP3_over_mP_D2_kg", "units": "kg"},
    "nonlocal_DP_like": {"key_lo": "nonlocal_lP2_D2_kg_s_per_m",
                         "key_up": "nonlocal_lP2_D2_kg_s_per_m", "units": "kg s m^-1"},
    "ultralocal_continuous": {"key_lo": "ultralocal_continuous_D2_kg2_s_per_m3",
                              "key_up": "ultralocal_continuous_D2_kg2_s_per_m3_table1_Vb",
                              "units": "kg^2 s m^-3"},
}


def is_s17(r):
    return (r.get("caveat") or "").startswith("Stage 17.")


def compute(rows):
    out = []
    for r in rows:
        p = dict(r)
        if r["role"] in ("lower_coherence", "lower_coherence_flagged",
                         "projection_lower", "projection_lower_flagged"):
            p["bounds"] = lower_bounds(r); p["side"] = "lower"
        elif r["role"] in ("upper_diffusion", "projection_upper"):
            p["bounds"] = upper_bounds(r); p["side"] = "upper"
        else:
            continue
        p["stage17"] = is_s17(r)
        out.append(p)
    return out


def best(pool, key, want):
    b = None
    for p in pool:
        v = p["bounds"][key]
        if b is None or (v > b["value"] if want == "max" else v < b["value"]):
            b = {"id": p["id"], "short_name": p["short_name"], "year": p["year"], "value": v}
    return b


def windows(pts, include_s17=True, extra_lower=()):
    lo = [p for p in pts if p["role"] == "lower_coherence" and (include_s17 or not p["stage17"])]
    lo += list(extra_lower)
    up = [p for p in pts if p["role"] == "upper_diffusion" and p["applicability"] == "absolute"
          and (include_s17 or not p["stage17"])]
    res = {}
    for name, m in MODELS.items():
        L = best(lo, m["key_lo"], "max"); Uu = best(up, m["key_up"], "min")
        res[name] = {"units": m["units"], "lower": L, "upper": Uu, "gap": gap_orders(L["value"], Uu["value"])}
    return res


def main():
    rows = load_csv()
    pts = compute(rows)
    byid = {p["id"]: p for p in pts}

    base = windows(pts, include_s17=False)          # Stage-15 baseline, same code
    now = windows(pts, include_s17=True)            # Stage 17 conservative
    flagged = [p for p in pts if p["role"] == "lower_coherence_flagged"]
    with_flag = windows(pts, include_s17=True, extra_lower=flagged)

    change = {}
    for n in MODELS:
        g0 = base[n]["gap"]; g1 = now[n]["gap"]
        change[n] = {"stage15_gap_orders": g0["orders"], "stage17_gap_orders": g1["orders"],
                     "closed_by_stage17_orders": g0["orders"] - g1["orders"],
                     "table1_2023_gap_orders": {"ultralocal_discrete": 24.0, "nonlocal_DP_like": 26.0,
                                                "ultralocal_continuous": 0.0}[n]}

    # Each Stage-17 row, its value and where it sits relative to the current window
    survey = []
    for p in pts:
        if not p["stage17"]:
            continue
        e = {"id": p["id"], "role": p["role"], "year": p["year"], "short_name": p["short_name"],
             "citation": p["citation"], "source": p["source"], "applicability": p["applicability"],
             "inputs": {k: p.get(k) for k in ("mass_kg", "mass_amu", "N", "tau_s", "delocalization_m", "FOM", "sqrt_SF", "sqrt_Sa")},
             "caveat": p["caveat"]}
        for n, m in MODELS.items():
            if n == "ultralocal_continuous":
                continue
            if p["side"] == "lower":
                v = p["bounds"][m["key_lo"]]; Uc = now[n]["upper"]["value"]; Lc = now[n]["lower"]["value"]
                e[n] = {"lower_bound": v, "orders_below_current_ceiling": math.log10(Uc / v),
                        "orders_vs_current_floor": math.log10(v / Lc)}
            else:
                v = p["bounds"][m["key_up"]]; Lc = now[n]["lower"]["value"]; Uc = now[n]["upper"]["value"]
                e[n] = {"upper_bound": v, "orders_above_current_floor": math.log10(v / Lc),
                        "orders_vs_current_ceiling": math.log10(v / Uc)}
        survey.append(e)

    # Near-term (<= 2030) projections: what each leaves open against today's ceiling
    near = []
    for p in pts:
        if p["role"] not in ("projection_lower", "projection_lower_flagged") or p["year"] > 2030:
            continue
        e = {"id": p["id"], "year": p["year"], "short_name": p["short_name"], "role": p["role"]}
        for n in ("ultralocal_discrete", "nonlocal_DP_like"):
            v = max(p["bounds"][MODELS[n]["key_lo"]], now[n]["lower"]["value"])
            e[n + "_gap_left_orders"] = math.log10(now[n]["upper"]["value"] / v)
        near.append(e)
    near.sort(key=lambda e: e["nonlocal_DP_like_gap_left_orders"])

    # What closes the non-local window from below with a MUSCLE-type beam (DERIVED, no fit)
    ped = byid["pedalino2025"]["bounds"]; Unl = now["nonlocal_DP_like"]["upper"]["value"]
    need_factor = Unl / ped["nonlocal_lP2_D2_kg_s_per_m"]
    tau_ped = byid["pedalino2025"]["tau_s"]; tau25 = 0.983 / 25.0
    # nonlocal lower ~ N M^2 tau / R ; N proportional to M for Na clusters -> M^3 scaling
    M_close_25 = 172e3 * (need_factor / (tau25 / tau_ped)) ** (1 / 3)
    M_close_160 = 172e3 * need_factor ** (1 / 3)
    # Discrete from below needs M*tau >= ceiling
    Ud = now["ultralocal_discrete"]["upper"]["value"]
    closing = {
        "nonlocal_from_below": {"factor_in_N_M2_tau_over_R_needed": need_factor,
                                "Na_cluster_mass_Da_needed_at_25_m_per_s": M_close_25,
                                "Na_cluster_mass_Da_needed_at_160_m_per_s": M_close_160,
                                "note": "same L=0.983 m, d=133 nm (ASSUMED); N = M/23 u"},
        "nonlocal_from_above": {"FOM_needed_m2_s3": now["nonlocal_DP_like"]["lower"]["value"] / (Unl / byid["gisler2022"]["FOM"]),
                                "best_absolute_FOM_now": byid["gisler2022"]["FOM"]},
        "discrete_from_below": {"M_tau_needed_kg_s": Ud, "e.g._mass_amu_at_100_s": Ud / 100 / U,
                                "atoms_hold_time_needed_Cs_s": Ud / (132.905 * U)},
    }

    # Trend extrapolation (ASSUMED: log-linear progress continues). Two-point rates.
    nl_lo_rate = math.log10(byid["pedalino2025"]["bounds"]["nonlocal_lP2_D2_kg_s_per_m"] /
                            byid["gerlich2011"]["bounds"]["nonlocal_lP2_D2_kg_s_per_m"]) / (2025 - 2011)
    fom_rate = math.log10(byid["teufel2009"]["FOM"] / byid["gisler2022"]["FOM"]) / (2022 - 2009)
    gap_nl = now["nonlocal_DP_like"]["gap"]["orders"]
    gap_d = now["ultralocal_discrete"]["gap"]["orders"]
    d_lo_rate = math.log10(byid["panda2024"]["bounds"]["ultralocal_discrete_lP3_over_mP_D2_kg"] /
                           byid["gerlich2011"]["bounds"]["ultralocal_discrete_lP3_over_mP_D2_kg"]) / (2024 - 2011)
    trend = {
        "assumption": "ASSUMED: past log-linear rates continue; not a schedule.",
        "nonlocal_floor_rate_orders_per_yr_2011_2025": nl_lo_rate,
        "absolute_FOM_rate_orders_per_yr_2009_2022": fom_rate,
        "discrete_floor_rate_orders_per_yr_2011_2024": d_lo_rate,
        "nonlocal_close_year_floor_only": 2025 + gap_nl / nl_lo_rate,
        "nonlocal_close_year_ceiling_only": 2022 + gap_nl / fom_rate,
        "nonlocal_close_year_both": 2025 + gap_nl / (nl_lo_rate + fom_rate),
        "discrete_close_year_both": 2025 + gap_d / (d_lo_rate + fom_rate),
        "note": "No absolute FOM record since Gisler 2022 (Leiden 2025 is 0.88, worse), so the ceiling-side rate has stalled for 3+ years.",
    }

    flagged_out = {n: {"lower_if_flagged_used": with_flag[n]["lower"], "gap": with_flag[n]["gap"]}
                   for n in ("ultralocal_discrete", "nonlocal_DP_like")}

    out = {"stage": 17, "baseline_stage15": base, "current_stage17": now, "change": change,
           "survey_rows": survey, "near_term_projections_by_2030": near, "closing": closing,
           "trend_extrapolation": trend, "flagged_temptation": flagged_out,
           "r_N_m": R_N}
    with open(JSON_PATH, "w") as f:
        json.dump(out, f, indent=2, default=lambda o: None)
    print("wrote", JSON_PATH)
    for n in MODELS:
        c = change[n]; w = now[n]
        print(f"{n}: stage15 {c['stage15_gap_orders']:.2f} -> stage17 {c['stage17_gap_orders']:.2f} orders "
              f"(floor {w['lower']['id']} {w['lower']['value']:.3e}, ceiling {w['upper']['id']} {w['upper']['value']:.3e})")
    print("\nSurvey rows:")
    for e in survey:
        d = e.get("ultralocal_discrete", {}); nl = e.get("nonlocal_DP_like", {})
        k = "lower_bound" if "lower_bound" in d else "upper_bound"
        print(f"  {e['id']:20s} {e['role']:24s} disc {k} {d.get(k):.3e}  nl {nl.get(k):.3e}")
    print("\nNear-term projections:")
    for e in near:
        print("  ", e)
    print("\nClosing:", json.dumps(closing, indent=1))
    print("Trend:", json.dumps(trend, indent=1))
    print("Flagged:", json.dumps(flagged_out, indent=1, default=str))


if __name__ == "__main__":
    main()
