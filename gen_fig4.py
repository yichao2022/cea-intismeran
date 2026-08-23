"""Regenerate fig4_price_icer.png with $0-800K x-axis, clear threshold markers."""
import sys
sys.path.insert(0, "/Users/cary/cea-intismeran")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model import ModelParams
from rerun_primary import CEAModelV2

def icer_at_price(price):
    p = ModelParams(constraint_general_pop=True, cost_intismeran=price)
    r = CEAModelV2(p).run()
    dc = r["combo"]["cost"] - r["pembro"]["cost"]
    dq = r["combo"]["qaly"] - r["pembro"]["qaly"]
    return dc, dq, dc / dq

# Price ladder: 0-800K, dense
prices = np.linspace(0, 800_000, 161)
icers = []
for pr in prices:
    _, _, icer = icer_at_price(pr)
    icers.append(icer)

# Threshold values
thresholds = {
    "$100K/QALY": 100_000,
    "$150K/QALY": 150_000,
}
prices_v = {
    "Base price ($200K)": 200_000,
    "Threshold @ $100K WTP (~$429K)": 429_000,
    "Threshold @ $150K WTP (~$616K)": 616_000,
}

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(prices / 1000, icers, "b-", linewidth=2.5, label="ICER")

# Horizontal WTP lines
for label, val in thresholds.items():
    ax.axhline(val, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.text(10, val + 3000, label, fontsize=9, color="gray")

# Vertical price lines
colors = ["#d62728", "#2ca02c", "#2ca02c"]
ls = ["-", "--", ":"]
for (label, val), c, l in zip(prices_v.items(), colors, ls):
    ax.axvline(val / 1000, color=c, linestyle=l, linewidth=1.2, alpha=0.7)
    y_pos = 470_000 if val < 500_000 else 420_000
    ax.text(val / 1000 + 3, y_pos, label, fontsize=8, color=c, rotation=90)

ax.set_xlabel("Intismeran Course Price ($1,000)", fontsize=12)
ax.set_ylabel("ICER ($/QALY)", fontsize=12)
ax.set_title("Price–ICER Curve and Deterministic Threshold Pricing", fontsize=13, fontweight="bold")
ax.set_xlim(0, 800)
ax.set_ylim(0, 500_000)
ax.set_xticks(np.arange(0, 801, 100))
ax.grid(alpha=0.3)
ax.legend(fontsize=9, loc="upper left")
plt.tight_layout()
plt.savefig("fig4_price_icer.png", dpi=300)
plt.savefig("output/fig4_price_icer.png", dpi=300)
print("Saved fig4_price_icer.png (0-800K, 3 vertical markers)")