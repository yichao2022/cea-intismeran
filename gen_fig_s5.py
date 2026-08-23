"""Figure S5: External validation - pembro predictions vs trial landmarks.
Uses only genuinely published landmark values — no extrapolated estimates."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import ModelParams
from rerun_primary import CEAModelV2

plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "figure.dpi": 300, "savefig.dpi": 300})

p = ModelParams(constraint_general_pop=True)
m = CEAModelV2(p)
t = m.t
sur = m._survival()

# Published landmarks (only confirmed values, no extrapolations)
# KEYNOTE-054: 7y RFS 50% (Eggermont NEJM 2025/ESMO 2024), 7y DMFS 54% (ESMO 2024)
# CheckMate 238: 5y RFS 51%, 7y RFS 46%, 9y RFS 44% (Ascierto NEJM 2025)
#               5y DMFS 59%, 7y DMFS 55% (Ascierto NEJM 2025)
#               OS 5y 76% (Ascierto JCO 2023), 7y OS 71% (Ascierto NEJM 2025)
landmarks = {
    "RFS": {
        "KEYNOTE-054": [(7, 50.0)],
        "CheckMate 238": [(5, 51.0), (7, 46.0)],
    },
    "DMFS": {
        "KEYNOTE-054": [(7, 54.0)],
        "CheckMate 238": [(5, 59.0), (7, 55.0)],
    },
    "OS": {
        "KEYNOTE-054": [],
        "CheckMate 238": [(5, 76.0), (7, 71.0)],
    },
}
colors = {"KEYNOTE-942": "black", "KEYNOTE-054": "#2166ac", "CheckMate 238": "#b2182b"}
markers = {"KEYNOTE-942": "o", "KEYNOTE-054": "s", "CheckMate 238": "^"}
titles = {"RFS": "Recurrence-Free Survival", "DMFS": "Distant Metastasis-Free Survival", "OS": "Overall Survival"}

# KEYNOTE-942 observed data (from Table S2)
obs = {
    "RFS": [(5, 49.1)],
    "DMFS": [],
    "OS": [(5, 71.3)],
}

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
for idx, curve in enumerate(["RFS", "DMFS", "OS"]):
    ax = axes[idx]
    for arm, label, ls in [("combo", "Intismeran + Pembro", "-"), ("pembro", "Pembrolizumab", "--")]:
        s = sur[f"rfs_{arm}"] if curve == "RFS" else sur[f"dmfs_{arm}"] if curve == "DMFS" else sur[f"os_{arm}"]
        ax.plot(t, s * 100, color="gray" if arm == "pembro" else "#4393c3", ls=ls, lw=1.5, label=label, alpha=0.8)
    # Observed
    for yr, val in obs.get(curve, []):
        ax.scatter([yr], [val], color="black", marker="o", s=50, zorder=10, label="KEYNOTE-942 (obs)" if idx == 0 else "")
    # External landmarks
    for trial, pts in landmarks[curve].items():
        if pts:
            yrs, vals = zip(*pts)
            ax.scatter(yrs, vals, color=colors[trial], marker=markers[trial], s=50, zorder=10, label=trial)
    ax.set_title(titles[curve], fontsize=11)
    ax.set_xlabel("Years")
    ax.set_ylabel("Survival (%)")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=6.5, loc="lower left")

fig.text(0.5, 0.01, 
    "Model pembrolizumab arm predictions (dashed gray) vs. external adjuvant trial landmarks.\n"
    "KEYNOTE-054 (blue squares): pembrolizumab, stage III, Eggermont et al. (NEJM 2025).\n"
    "CheckMate 238 (red triangles): nivolumab, stage IIIB-IV, Ascierto et al. (NEJM 2025).\n"
    "Landmarks are published values only; no extrapolated estimates are used.",
    ha="center", fontsize=7.5, fontstyle="italic")
plt.tight_layout(rect=[0, 0.06, 1, 1])
plt.savefig("/tmp/figures/fig_s5_external_validation.png", dpi=300)
plt.savefig("/tmp/figures/fig_s5_external_validation.pdf")
print("Saved Figure S5 (clean landmarks)")