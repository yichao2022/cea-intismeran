"""Generate Figure S2: Digitized KEYNOTE-942 survival data + fitted parametric curves."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm
import sys
sys.path.insert(0, "/Users/cary/cea-intismeran")
from model import ModelParams, lognorm_surv

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10,
    "axes.titlesize": 12, "axes.labelsize": 11,
    "figure.dpi": 300, "savefig.dpi": 300,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.1,
})
C = {"combo": "#2166ac", "pembro": "#b2182b"}
LAB = {"combo": "Intismeran + Pembro", "pembro": "Pembrolizumab"}

p = ModelParams()
t_fine = np.arange(0, 61, 0.5)

# Digitized data points
data = {
    "combo": {
        "RFS": (p.rfs_time_mo, p.rfs_combo_pct),
        "DMFS": (p.dmfs_time_mo, p.dmfs_combo_pct),
        "OS":  (p.os_time_mo,  p.os_combo_pct),
    },
    "pembro": {
        "RFS": (p.rfs_time_mo, p.rfs_pembro_pct),
        "DMFS": (p.dmfs_time_mo, p.dmfs_pembro_pct),
        "OS":  (p.os_time_mo,  p.os_pembro_pct),
    },
}

# Recompute fitted curves
def pct(v): return np.array(v, dtype=float) / 100.0
t_os = p.os_time_mo
# OS
os_c = lognorm_surv(t_fine, p.os_mu_combo, p.os_sigma_combo) * 100
os_p = lognorm_surv(t_fine, p.os_mu_pembro, p.os_sigma_pembro) * 100
# r_DM = DMFS/OS
rd_c_fit = lognorm_surv(t_fine, p.rd_mu_combo, p.rd_sigma_combo)
rd_p_fit = lognorm_surv(t_fine, p.rd_mu_pembro, p.rd_sigma_pembro)
dmfs_c = os_c / 100 * rd_c_fit * 100
dmfs_p = os_p / 100 * rd_p_fit * 100
# r_LR = RFS/DMFS
rl_c_fit = lognorm_surv(t_fine, p.rl_mu_combo, p.rl_sigma_combo)
rl_p_fit = lognorm_surv(t_fine, p.rl_mu_pembro, p.rl_sigma_pembro)
rfs_c = dmfs_c / 100 * rl_c_fit * 100
rfs_p = dmfs_p / 100 * rl_p_fit * 100

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)

for ax, arm in [(ax1, "combo"), (ax2, "pembro")]:
    assert arm in ("combo", "pembro")
    ax.set_title(LAB[arm], fontweight="bold")
    ax.set_xlabel("Months")
    if arm == "combo":
        ax.set_ylabel("Survival Probability (%)")

    # Fitted curves
    ax.plot(t_fine, rfs_c if arm == "combo" else rfs_p, "-", color=C[arm], lw=1.8, label="RFS (fitted)")
    ax.plot(t_fine, dmfs_c if arm == "combo" else dmfs_p, "--", color=C[arm], lw=1.8, label="DMFS (fitted)")
    ax.plot(t_fine, os_c if arm == "combo" else os_p, ":", color=C[arm], lw=1.8, label="OS (fitted)")

    # Digitized points
    for curve, marker, ms in [("RFS", "o", 7), ("DMFS", "s", 7), ("OS", "^", 7)]:
        tx, vy = data[arm][curve]
        ax.scatter(tx, vy, marker=marker, s=ms**2, color=C[arm], zorder=5, edgecolors="white", linewidth=0.5)

    ax.set_xlim(0, 62)
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")

# Footnote
fig.text(0.5, 0.01, "Digitized from Figure 1 of Khattak et al. (2026) JCO. Points: digitized KM values. Lines: fitted log-normal parametric curves.",
         ha="center", fontsize=8, fontstyle="italic")

plt.tight_layout(rect=[0, 0.03, 1, 1])
plt.savefig("/tmp/figures/fig_s2_survival_data.png", dpi=300)
plt.savefig("/tmp/figures/fig_s2_survival_data.pdf")
print("Saved fig_s2")