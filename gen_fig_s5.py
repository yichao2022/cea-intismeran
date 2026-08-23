"""Generate Figure S5: External validation - pembro predictions vs trial landmarks."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import ModelParams, CEAModel

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10,
    "axes.titlesize": 12, "axes.labelsize": 11,
    "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.1,
})

p = ModelParams(constraint_general_pop=True)
model = CEAModel(p)
sur = model._survival()
t = model.t

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4), sharey=False)

# External landmarks
# KEYNOTE-054: Eggermont et al NEJM 2025, 7y update
# CheckMate 238: Ascierto et al JCO 2024, 7y update
landmarks = {
    "RFS": {
        "KEYNOTE-942": [(5, 49.1)],
        "KEYNOTE-054": [(5, 49.1), (7, 50.0)],
        "CheckMate 238": [(5, 50.0), (7, 44.0)],
    },
    "DMFS": {
        "KEYNOTE-942": [(5, 65.4)],
        "KEYNOTE-054": [(5, 58.0), (7, 54.0)],
        "CheckMate 238": [(5, 58.0), (7, 50.0)],
    },
    "OS": {
        "KEYNOTE-942": [(5, 71.3)],
        "KEYNOTE-054": [],
        "CheckMate 238": [(5, 76.0), (7, 70.0)],
    },
}

colors = {"KEYNOTE-942": "black", "KEYNOTE-054": "#2166ac", "CheckMate 238": "#b2182b"}
markers = {"KEYNOTE-942": "o", "KEYNOTE-054": "s", "CheckMate 238": "^"}

for ax, curve, title in [
    (ax1, "RFS", "Recurrence-Free Survival"),
    (ax2, "DMFS", "Distant Metastasis-Free Survival"),
    (ax3, "OS", "Overall Survival"),
]:
    # Model prediction
    ax.plot(t, sur[f"{curve.lower()}_pembro"] * 100, "-", color="gray", lw=2, label="Model (pembro)")

    # External landmarks
    for trial, pts in landmarks[curve].items():
        if pts:
            yrs = [p[0] for p in pts]
            vals = [p[1] for p in pts]
            ax.scatter(yrs, vals, marker=markers[trial], s=60,
                       color=colors[trial], zorder=5, label=trial)

    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("Years")
    ax.set_ylabel("Survival Probability (%)")
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 100)
    ax.legend(fontsize=7, loc="lower left")
    ax.grid(True, alpha=0.3)

fig.text(0.5, 0.01,
         "Model pembrolizumab arm predictions (solid gray line) vs. external adjuvant trial landmarks.\n"
         "KEYNOTE-942 (black circles): 5-year data from Khattak et al. (2026).\n"
         "KEYNOTE-054 (blue squares): pembrolizumab adjuvant, stage III, Eggermont et al. (NEJM 2025).\n"
         "CheckMate 238 (red triangles): nivolumab adjuvant, stage IIIB-IV, Ascierto et al. (JCO 2024).\n"
         "Differences reflect trial population heterogeneity; see text for discussion.",
         ha="center", fontsize=7, fontstyle="italic")

plt.tight_layout(rect=[0, 0.06, 1, 1])
plt.savefig("/tmp/figures/fig_s5_external_validation.png", dpi=300)
plt.savefig("/tmp/figures/fig_s5_external_validation.pdf")
print("Saved Figure S5")