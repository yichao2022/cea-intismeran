"""Regenerate Figure S6 (convergence) and Figure S7 (CE plane) from 6,000-run PSA draws."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

with open("output/psa_v2_draws.json") as f:
    draws = json.load(f)
dq = np.array([d["dq"] for d in draws])
dc = np.array([d["dc"] for d in draws])
n = len(dq)
print(f"{n} draws loaded")

# ── Figure S6: convergence ──
cum_dq = np.cumsum(dq) / np.arange(1, n + 1)
cum_dc = np.cumsum(dc) / np.arange(1, n + 1)
nmbs = {100_000: -dc + 100_000 * dq}
cum_nmb = np.cumsum(nmbs[100_000]) / np.arange(1, n + 1)
cum_pce = np.array([np.mean(-dc[:i] + 100_000 * dq[:i] > 0) for i in range(1, n + 1)])

fig, axes = plt.subplots(2, 2, figsize=(10, 7))
x = np.arange(1, n + 1)
axes[0, 0].plot(x, cum_dq, "b-")
axes[0, 0].set_ylabel("Cumulative mean ΔQALY")
axes[0, 1].plot(x, cum_dc / 1000, "b-")
axes[0, 1].set_ylabel("Cumulative mean ΔCost ($K)")
axes[1, 0].plot(x, cum_nmb / 1000, "b-")
axes[1, 0].set_ylabel("Cumulative mean NMB@$100K ($K)")
axes[1, 1].plot(x, cum_pce * 100, "b-")
axes[1, 1].set_ylabel("Cumulative P(CE)@$100K (%)")
for ax in axes.flat:
    ax.set_xlabel("Iterations")
    ax.grid(alpha=0.3)
    for v in [500, 1000, 2000, 4000]:
        ax.axvline(v, color="gray", linestyle=":", linewidth=0.6)
fig.suptitle("PSA Convergence (6,000 Iterations, 0 Rejections)", fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig("fig_s6_psa_convergence.png", dpi=300)
plt.savefig("output/fig_s6_psa_convergence.png", dpi=300)
print("Saved fig_s6_psa_convergence.png")

# ── Figure S7: CE plane ──
fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(dq, dc / 1000, s=4, alpha=0.4, c="#2166ac", edgecolors="none")
# WTP lines: dc = wtp * dq  (in $K)
for wtp, color, ls in [(50_000, "gray", "--"), (100_000, "gray", "-."), (150_000, "gray", ":")]:
    xx = np.linspace(0, max(dq) * 1.05, 2)
    ax.plot(xx, wtp * xx / 1000, color=color, linestyle=ls, linewidth=1,
            label=f"WTP ${wtp/1000:.0f}K/QALY")
ax.axhline(0, color="black", linewidth=0.8)
ax.axvline(0, color="black", linewidth=0.8)
# Deterministic base case
ax.plot(3.76, 146.961, "kD", markersize=9, label="Deterministic base case")
ax.set_xlabel("Incremental QALYs", fontsize=12)
ax.set_ylabel("Incremental Cost ($1,000)", fontsize=12)
ax.set_title("Cost-Effectiveness Plane (6,000 PSA Iterations)", fontsize=13, fontweight="bold")
ax.legend(fontsize=9, loc="upper left")
ax.set_xlim(-0.5, max(dq) * 1.05)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("fig_s7_ce_plane.png", dpi=300)
plt.savefig("output/fig_s7_ce_plane.png", dpi=300)
print("Saved fig_s7_ce_plane.png")

# sanity: fraction SE quadrant
frac_se = np.mean((dc < 0) & (dq > 0))
print(f"Sanity: n={n}, frac cost-saving={frac_se:.3f}, P(CE)@100K={np.mean(-dc+100_000*dq>0):.3f}")