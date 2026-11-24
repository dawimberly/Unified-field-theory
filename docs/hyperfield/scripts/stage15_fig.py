#!/usr/bin/env python3
"""Stage 15 figure: Oppenheim classical-gravity window vs year, FOM, projections."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
D = json.load(open(os.path.join(ROOT, "results/stage15.json")))
PTS = D["datapoints"]

LOWER_KEY = {
    "ultralocal_discrete": "ultralocal_discrete_lP3_over_mP_D2_kg",
    "nonlocal_DP_like": "nonlocal_lP2_D2_kg_s_per_m",
    "ultralocal_continuous": "ultralocal_continuous_D2_kg2_s_per_m3",
}
UPPER_KEY = {
    "ultralocal_discrete": "ultralocal_discrete_lP3_over_mP_D2_kg",
    "nonlocal_DP_like": "nonlocal_lP2_D2_kg_s_per_m",
    "ultralocal_continuous": "ultralocal_continuous_D2_kg2_s_per_m3_table1_Vb",
}


def ev_step(X, Y, xq):
    out = np.empty_like(xq)
    j, last = 0, Y[0]
    for i, x in enumerate(xq):
        while j < len(X) and X[j] <= x:
            last = Y[j]
            j += 1
        out[i] = last
    return out


def series_xy(running, lo=1998, hi=2036):
    xs, ys, prev = [], [], None
    for e in running:
        y, v = e["year"], e["best_so_far"]
        if y < lo:
            prev = v
            continue
        if not xs:
            xs.append(lo)
            ys.append(prev if prev is not None else v)
        xs.append(y)
        ys.append(v)
    if xs:
        xs.append(hi)
        ys.append(ys[-1])
    return np.array(xs, float), np.array(ys, float)


def step_xy(xs, ys):
    X, Y = [], []
    for i, (x, y) in enumerate(zip(xs, ys)):
        if i:
            X.append(x)
            Y.append(Y[-1])
        X.append(x)
        Y.append(y)
    return np.array(X), np.array(Y)


fig, axes = plt.subplots(2, 2, figsize=(12.4, 9.1))

# -------------------- (a, b) windows vs year --------------------
panel = [
    (axes[0, 0], "ultralocal_discrete",
     r"(a)  Ultra-local discrete   $\ell_P^3/m_P\,D_2$  [kg]",
     1e-26, 3e1),
    (axes[0, 1], "nonlocal_DP_like",
     r"(b)  Non-local / Diósi–Penrose-like   $\ell_P^2 D_2$  [kg s m$^{-1}$]",
     1e-37, 3e-9),
]
labels_drawn = False
for ax, key, title, ymin, ymax in panel:
    hist = D["history"][key]
    xL, yL = series_xy(hist["lower_running_max"])
    xU, yU = series_xy(hist["upper_running_min_absolute"])
    xA, yA = series_xy(hist["upper_running_min_including_differential"])
    XL, YL = step_xy(xL, yL)
    XU, YU = step_xy(xU, yU)
    XA, YA = step_xy(xA, yA)
    ax.semilogy(XL, YL, color="C0", lw=2.2, label="lower bound (coherence, running best)")
    ax.semilogy(XU, YU, color="C3", lw=2.2, label="upper bound (absolute force noise)")
    ax.semilogy(XA, YA, color="C3", lw=1.25, ls="--",
                label="upper if differential points are allowed")
    years = np.linspace(2011, 2026.0, 300)
    lo = ev_step(XL, YL, years)
    hi = ev_step(XU, YU, years)
    ax.fill_between(years, lo, hi, where=(hi > lo), color="C2", alpha=0.22, lw=0,
                    label="conservative allowed window")
    lk, uk = LOWER_KEY[key], UPPER_KEY[key]
    first_abs = first_diff = first_low = True
    for p in PTS:
        if p["role"] == "lower_coherence":
            ax.scatter(p["year"], p["bounds"][lk], color="C0", s=42, zorder=5,
                       edgecolor="k", lw=0.35,
                       label="coherence experiment" if first_low else None)
            first_low = False
        elif p["role"] == "upper_diffusion" and p["year"] >= 1998:
            v = p["bounds"][uk]
            if p["applicability"] == "absolute":
                ax.scatter(p["year"], v, color="C3", s=16, zorder=4, alpha=0.4,
                           label="absolute force-noise point" if first_abs else None)
                first_abs = False
            else:
                ax.scatter(p["year"], v, facecolors="none", edgecolors="C3",
                           s=32, zorder=5, lw=0.9,
                           label="differential (applicability questioned)" if first_diff else None)
                first_diff = False
    first_pl = first_pu = True
    for p in D["projections"]:
        if p["role"] == "projection_lower":
            ax.scatter(p["year"], p[key]["projected_lower"], marker="*",
                       color="C0", s=100, zorder=6, edgecolor="k", lw=0.3,
                       label="projection (coherence)" if first_pl else None)
            first_pl = False
        else:
            ax.scatter(p["year"], p[key]["projected_upper"], marker="*",
                       color="C3", s=100, zorder=6, edgecolor="k", lw=0.3,
                       label="projection (force noise)" if first_pu else None)
            first_pu = False
    cav = next(x for x in PTS if x["id"] == "cavendish1798")
    ax.axhline(cav["bounds"][uk], color="0.4", ls=":", lw=0.8)
    ax.text(1999.1, cav["bounds"][uk] * (1.6 if key.endswith("discrete") else 1.6),
            "Cavendish 1798", fontsize=7, color="0.3")
    ax.axvline(2026, color="0.55", ls=":", lw=0.7)
    ax.set_xlim(1998, 2036)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("year")
    ax.set_title(title, loc="left", fontsize=10)
    ax.legend(fontsize=5.7, loc="upper right", framealpha=0.9, handlelength=1.6)

# annotate a few named points on (a)
ax = axes[0, 0]
ann = {"gerlich2011": (8, 8), "fein2019": (8, -14), "pedalino2025": (6, 10),
       "gisler2022": (-28, 10), "proj_1e6amu": (4, 8), "proj_maqro": (-18, 8)}
for p in PTS:
    if p["id"] in ("gerlich2011", "fein2019", "pedalino2025"):
        ax.annotate(p["short_name"].split("  ")[0],
                    (p["year"], p["bounds"][LOWER_KEY["ultralocal_discrete"]]),
                    textcoords="offset points", xytext=ann[p["id"]],
                    fontsize=6.5, color="C0")
for p in PTS:
    if p["id"] == "gisler2022":
        ax.annotate("Gisler 2022",
                    (p["year"], p["bounds"][UPPER_KEY["ultralocal_discrete"]]),
                    textcoords="offset points", xytext=ann[p["id"]],
                    fontsize=6.5, color="C3")
for p in D["projections"]:
    if p["id"] in ("proj_1e6amu", "proj_maqro"):
        ax.annotate(p["short_name"].split("  ")[0] if p["id"] != "proj_maqro" else "MAQRO",
                    (p["year"], p["ultralocal_discrete"]["projected_lower"]),
                    textcoords="offset points", xytext=ann[p["id"]],
                    fontsize=6.5, color="C0")

# -------------------- (c) FOM vs mass (Janse-style) --------------------
ax = axes[1, 0]
# lower-bound equivalent FOM: FOM = lower / K
K_DISC = D["constants"]["K_disc_kg_s3_per_m2"]
K_NL = D["constants"]["K_nl_kg_s4_per_m3"]
Ldisc = D["current_windows"]["ultralocal_discrete"]["lower_now"]["value"]
Lnl = D["current_windows"]["nonlocal_DP_like"]["lower_now"]["value"]
ax.axhline(Ldisc / K_DISC, color="C0", ls="--", lw=1.2,
           label=r"Pedalino floor, discrete  (FOM needed to close)")
ax.axhline(Lnl / K_NL, color="C2", ls="--", lw=1.2,
           label=r"Pedalino floor, non-local  (FOM needed to close)")
style = {
    "absolute": dict(c="C3", marker="o", s=22, alpha=0.85),
    "differential": dict(facecolors="none", edgecolors="C1", marker="o", s=40, lw=1.0),
    "differential_L1": dict(facecolors="none", edgecolors="C4", marker="s", s=40, lw=1.0),
}
seen = set()
for p in PTS:
    if p["role"] != "upper_diffusion" or p.get("FOM") is None or not p.get("mass_kg"):
        continue
    appl = p["applicability"]
    kw = dict(style[appl])
    lab = None
    if appl not in seen:
        lab = {"absolute": "absolute (Janse filled)",
               "differential": "differential (atom interferometer)",
               "differential_L1": "LISA Pathfinder (differential, L1)"}[appl]
        seen.add(appl)
        kw["label"] = lab
    ax.scatter(p["mass_kg"], p["FOM"], **kw)
    ax.set_xscale("log"); ax.set_yscale("log")
    if p["id"] in ("cavendish1798", "gisler2022", "asenbaum2017",
                   "armano2018", "martynov2016", "westphal2021", "fuchs2024"):
        ax.annotate(p["short_name"].replace(" 2018", "").replace(" 2022", "")
                    .replace(" 2017", "").replace(" 2016", "").replace(" 2021", "")
                    .replace(" 2024", "").replace(" 1798", ""),
                    (p["mass_kg"], p["FOM"]), textcoords="offset points",
                    xytext=(4, 3), fontsize=6)
# projection FOM
for p in D["projections"]:
    if p["role"] == "projection_upper" and p.get("bounds", {}).get("FOM_m2_per_s3"):
        ax.scatter(4e-7, p["bounds"]["FOM_m2_per_s3"], marker="*", c="C3",
                   s=90, edgecolor="k", lw=0.3, zorder=6, label="projected thermal-limited FOM")
ax.set_xlim(3e-27, 3e2)
ax.set_ylim(1e-16, 1e16)
ax.set_xlabel("test mass  [kg]")
ax.set_ylabel(r"FOM$_{D_2}=\sigma_a^2 N\Delta T$   [m$^2$ s$^{-3}$]")
ax.set_title(r"(c)  Force-noise figure of merit vs mass  (Janse 2024 compilation)",
             loc="left", fontsize=10)
ax.legend(fontsize=6.3, loc="lower left", framealpha=0.92)

# -------------------- (d) gaps and what closes them --------------------
ax = axes[1, 1]
ax.set_axis_off()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_title(r"(d)  Current conservative gaps, and what would close them",
             loc="left", fontsize=10)

# a compact text panel — clearer than a cramped bar chart
lines = [
    "Conservative window  =  Pedalino 2025 floor  vs  Gisler 2022 ceiling.",
    "Differential points (Asenbaum, LISA Pathfinder) are excluded here.",
    "",
    "Ultra-local continuous (Newtonian, D2 \u03b4(x-y)):  RULED OUT (2023 Table 1).",
    "    10\u207b\u2074\u00b9 \u2265 D2 \u2265 10\u207b\u2079 kg\u00b2 s m\u207b\u00b3, window empty.",
    "",
    "Ultra-local discrete   (\u2113_P\u00b3/m_P) D2:   OPEN by  10.1  orders.",
    "    lower  3.51\u00d710\u207b\u00b2\u2074 kg   (Pedalino, M/\u03bb).",
    "    upper  4.00\u00d710\u207b\u00b9\u2074 kg   (Gisler, formula).",
    "    To close from below: need M\u03c4 \u2273 4\u00d710\u207b\u00b9\u2074 kg s",
    "        \u2014  10\u2076 amu at 0.1 s leaves 8.4 orders;",
    "        \u2014  MAQRO 3\u00d710\u00b9\u2070 amu, 100 s leaves 0.9 orders;",
    "        \u2014  a slightly heavier/longer MAQRO-class run would finish it.",
    "    To close from above: need FOM \u2272 2.6\u00d710\u207b\u00b9\u00b9  (~10 more orders",
    "        than Gisler). That is Asenbaum-level, but as an absolute measurement.",
    "",
    "Non-local / DP-like   \u2113_P\u00b2 D2:   OPEN by  3.7  orders.",
    "    lower 1.47\u00d710\u207b\u00b2\u2079,  upper 6.69\u00d710\u207b\u00b2\u2076  (same pair).",
    "    10\u2076 amu at 0.1 s leaves only 0.15 orders  (essentially closes).",
    "    A ~4-order FOM improvement (thermal-limited mag-lev / SQUID) also closes it.",
    "",
    "If Asenbaum 2017 were applicable it would already close both remaining windows.",
    "Janse 2024 flag it as a differential / phase-shift measurement.  Do not use it",
    "as a claim that classical gravity is already ruled out.",
    "",
    "Separate 2026 accounting (Oppenheim\u2013Sajjad, dimensionless D2):",
    "    10\u207b\u2077\u00b9 \u2272 D2,tensor \u2272 10\u207b\u2076\u2077  (~4.5 orders), still with the 2011 floor.",
]
text = "\n".join(lines)
ax.text(0.02, 0.98, text, transform=ax.transAxes, va="top", ha="left",
        fontsize=7.15, family="DejaVu Sans", linespacing=1.28,
        bbox=dict(boxstyle="round,pad=0.45", facecolor="#f7f7f4", edgecolor="0.75", lw=0.6))

fig.suptitle("Stage 15: Oppenheim classical-gravity decoherence–diffusion window",
             fontsize=12, y=1.01)
fig.tight_layout(rect=(0,0,1,0.98))
fig.subplots_adjust(hspace=0.32, wspace=0.28)
out = os.path.join(ROOT, "figures/stage15_oppenheim.png")
fig.savefig(out, dpi=140, bbox_inches="tight")
print("wrote", out)
