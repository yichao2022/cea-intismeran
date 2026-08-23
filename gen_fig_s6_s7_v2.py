"""Regenerate Figure S6 (PSA convergence) and Figure S7 (CE plane) from psa_v2 draws."""
import sys, os
sys.path.insert(0, "/Users/cary/cea-intismeran")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.stats import norm, multivariate_normal
from model import ModelParams
from rerun_primary import CEAModelV2

plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "figure.dpi": 300, "savefig.dpi": 300})

# Load full draws from psa_v2.py
data = np.load("/tmp/figures/psa_v2_draws.npz")
dq = data["dq"]  # all 2000 draws
dc = data["dc"]

# ── Figure S6: PSA convergence ──
fig, axes = plt.subplots(2, 2, figsize=(10, 7))
cum_q = np.cumsum(dq) / np.arange(1, len(dq) + 1)
cum_c = np.cumsum(dc) / np.arange(1, len(dc) + 1)
cum_nmb = np.cumsum(-dc + 100_000 * dq) / np.arange(1, len(dc) + 1)
cum_pce = np.array([np.mean(-dc[:i+1] + 100_000 * dq[:i+1] > 0) for i in range(len(dq))])

for ax, data, label, unit in [
    (axes[0, 0], cum_q, "Mean ΔQALY", "QALYs"),
    (axes[0, 1], cum_c, "Mean ΔCost", "USD"),
    (axes[1, 0], cum_nmb, "Mean NMB ($100K WTP)", "USD"),
    (axes[1, 1], cum_pce, "P(CE) @ $100K", "Probability"),
]:
    ax.plot(range(1, len(data) + 1), data, "-", color="#2166ac", lw=1)
    ax.axvline(500, color="gray", ls=":", alpha=0.3)
    ax.axvline(1000, color="gray", ls=":", alpha=0.3)
    ax.axvline(2000, color="gray", ls=":", alpha=0.3)
    ax.set_xlabel("PSA Iterations")
    ax.set_ylabel(label)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("/tmp/figures/fig_s6_psa_convergence.png", dpi=300)
plt.savefig("/tmp/figures/fig_s6_psa_convergence.pdf")
print("Saved Figure S6")

# ── Figure S7: CE plane ──
fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(dq, dc, s=5, alpha=0.3, color="#2166ac", label="PSA iterations (n=2000)")
ax.axhline(0, color="gray", ls=":", alpha=0.5)
ax.axvline(0, color="gray", ls=":", alpha=0.5)
for wtp, color, ls in [(50_000, "#e41a1c", "-"), (100_000, "#4daf4a", "-"), (150_000, "#984ea3", "-")]:
    x = np.linspace(-2, 10, 100)
    ax.plot(x, wtp * x, color=color, ls=ls, lw=1, label=f"WTP=${wtp:,}/QALY")
# Base case
p = ModelParams(constraint_general_pop=True)
m = CEAModelV2(p)
r = m.run()
base_dq = r["combo"]["qaly"] - r["pembro"]["qaly"]
base_dc = r["combo"]["cost"] - r["pembro"]["cost"]
ax.scatter([base_dq], [base_dc], s=80, color="black", marker="D", zorder=10, label="Base case")
ax.set_xlabel("Incremental QALYs")
ax.set_ylabel("Incremental Cost (USD)")
ax.legend(fontsize=8, loc="upper left")
ax.grid(True, alpha=0.3)
ax.set_xlim(-2, 10)
ax.set_ylim(-500_000, 800_000)
plt.tight_layout()
plt.savefig("/tmp/figures/fig_s7_ce_plane.png", dpi=300)
plt.savefig("/tmp/figures/fig_s7_ce_plane.pdf")
print("Saved Figure S7")