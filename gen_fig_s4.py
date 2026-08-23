"""Figure S4: OS vs general-population survival — hazard floor (CEAModelV2)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import ModelParams
from rerun_primary import CEAModelV2, LT_surv_monthly

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10,
    "axes.titlesize": 12, "axes.labelsize": 11,
    "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.1,
})

p_gp = ModelParams(constraint_general_pop=True)
model_gp = CEAModelV2(p_gp)
sur_gp = model_gp._survival()
t_yrs = model_gp.t

p_unc = ModelParams(constraint_general_pop=False)
model_unc = CEAModelV2(p_unc)
sur_unc = model_unc._survival()

# GP survival from life table
gp_surv = np.interp(t_yrs * 12, np.arange(480), LT_surv_monthly)

fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(t_yrs, sur_gp["os_combo"] * 100, "-", color="#2166ac", lw=2, label="Intismeran + Pembro (GP-constrained)")
ax.plot(t_yrs, sur_gp["os_pembro"] * 100, "-", color="#b2182b", lw=2, label="Pembrolizumab (GP-constrained)")
ax.plot(t_yrs, gp_surv * 100, "--", color="gray", lw=1.5, label="US general population (age 61, 64% male)")

ax.plot(t_yrs, sur_unc["os_combo"] * 100, ":", color="#2166ac", lw=1.2, alpha=0.5, label="Intismeran + Pembro (unconstrained)")
ax.plot(t_yrs, sur_unc["os_pembro"] * 100, ":", color="#b2182b", lw=1.2, alpha=0.5, label="Pembrolizumab (unconstrained)")

ax.axvline(5, color="black", ls=":", alpha=0.4, label="GP hazard floor onset (60 mo)")
ax.set_xlabel("Years")
ax.set_ylabel("Overall Survival (%)")
ax.set_xlim(0, 40)
ax.set_ylim(0, 100)
ax.legend(fontsize=8, loc="lower left")
ax.grid(True, alpha=0.3)

fig.text(0.5, 0.01, "GP-constrained curves (solid) use a mortality hazard floor from month 60. "
         "The constrained combo OS stays below the GP curve at all time points, "
         "reflecting the excess mortality experienced in the first 5 years. "
         "Unconstrained log-normal (dotted) produces implausibly high 40-year survival (76.4%).",
         ha="center", fontsize=8, fontstyle="italic")

plt.tight_layout(rect=[0, 0.03, 1, 1])
plt.savefig("/tmp/figures/fig_s4_os_vs_gp.png", dpi=300)
plt.savefig("/tmp/figures/fig_s4_os_vs_gp.pdf")
print("Saved Figure S4 (hazard floor, CEAModelV2)")