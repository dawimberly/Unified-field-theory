#!/usr/bin/env python3
"""Stage 17 figure: every tracker datapoint vs year in the two open Oppenheim windows."""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from stage17_run import load_csv, compute, MODELS  # noqa: E402

D = json.load(open(os.path.join(ROOT, "results/stage17.json")))
pts = compute(load_csv())

fig, axes = plt.subplots(1, 2, figsize=(15, 6.8))
for ax, name, title in ((axes[0], "ultralocal_discrete", r"Ultra-local discrete: $(l_P^3/m_P)D_2$ [kg]"),
                        (axes[1], "nonlocal_DP_like", r"Non-local / DP-like: $l_P^2 D_2$ [kg s m$^{-1}$]")):
    m = MODELS[name]; cur = D["current_stage17"][name]; base = D["baseline_stage15"][name]
    L, Uu = cur["lower"]["value"], cur["upper"]["value"]
    ax.axhspan(L, Uu, color="tab:green", alpha=0.12, zorder=0)
    ax.axhline(L, color="tab:blue", lw=1.2, ls="-"); ax.axhline(Uu, color="tab:red", lw=1.2, ls="-")
    if base["lower"]["value"] != L:
        ax.axhline(base["lower"]["value"], color="tab:blue", lw=1, ls=":")
        ax.text(1999, base["lower"]["value"], " Stage-15 floor", color="tab:blue", fontsize=7, va="top")
    for p in pts:
        yr = p["year"]
        if p["side"] == "lower":
            v = p["bounds"][m["key_lo"]]
            if p["role"] == "lower_coherence":
                st = dict(marker="o", color="tab:blue", mfc="tab:blue")
            elif p["role"] == "lower_coherence_flagged":
                st = dict(marker="o", color="gray", mfc="none")
            elif p["role"] == "projection_lower":
                st = dict(marker="*", color="tab:purple", mfc="tab:purple", ms=11)
            else:
                st = dict(marker="*", color="gray", mfc="none", ms=11)
        else:
            if p["applicability"] not in ("absolute", "projection"):
                continue
            v = p["bounds"][m["key_up"]]
            st = dict(marker="v", color="tab:red", mfc="tab:red" if p["applicability"] == "absolute" else "none")
        if yr < 1999:
            continue
        ax.semilogy([yr], [v], ls="none", ms=st.pop("ms", 7), **st)
        if p["stage17"] or p["id"] in ("pedalino2025", "gisler2022", "gerlich2011", "proj_1e6amu", "proj_maqro"):
            ax.annotate(p["short_name"], (yr, v), fontsize=6.5, xytext=(4, 2), textcoords="offset points",
                        color="dimgray" if "flagged" in p["role"] else "black")
    ax.set_xlim(1998, 2038)
    lo = min(L, 1e-40) if name == "nonlocal_DP_like" else 1e-26
    ax.set_ylim(1e-41 if name == "nonlocal_DP_like" else 1e-26, 1e5 if name == "nonlocal_DP_like" else 1e-10)
    g = cur["gap"]["orders"]; g0 = base["gap"]["orders"]
    ax.set_title(f"{title}\nconservative window: {g0:.1f} → {g:.2f} orders", fontsize=10)
    ax.set_xlabel("year"); ax.grid(alpha=0.3, which="both", lw=0.4)
axes[0].set_ylabel("bound value")
from matplotlib.lines import Line2D
h = [Line2D([], [], marker="o", ls="none", color="tab:blue", label="coherence floor (used)"),
     Line2D([], [], marker="o", ls="none", color="gray", mfc="none", label="flagged: Δx < object size / not discriminating (NOT used)"),
     Line2D([], [], marker="*", ls="none", color="tab:purple", ms=10, label="projection (year ASSUMED unless stated)"),
     Line2D([], [], marker="*", ls="none", color="gray", mfc="none", ms=10, label="flagged projection"),
     Line2D([], [], marker="v", ls="none", color="tab:red", label="absolute force-noise ceiling"),
     Line2D([], [], color="tab:green", lw=8, alpha=0.3, label="open window (Stage 17, conservative)")]
fig.legend(handles=h, loc="lower center", ncol=3, fontsize=8, frameon=False)
fig.suptitle("Stage 17 — Oppenheim classical-gravity windows after the 2026 experiment survey", fontsize=12)
fig.tight_layout(rect=(0, 0.08, 1, 0.95))
out = os.path.join(ROOT, "figures/stage17_window.png")
fig.savefig(out, dpi=150); print("wrote", out)
