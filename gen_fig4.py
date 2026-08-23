"""Regenerate fig4_price_icer.png — journal style, no in-figure title."""
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

prices = np.linspace(0, 800_000, 161)
icers = np.array([icer_at_price(pr)[2] for pr in prices])

fig, ax = plt.subplots(figsize=(8, 5.4))

# Main ICER curve (thick, no legend entry needed — it's the only main curve)
ax.plot(prices / 1000, icers, "k-", linewidth=2.2)

# Horizontal WTP lines (thin gray dashed; small labels at left)
ax.axhline(100_000, color="gray", linestyle="--", linewidth=0.9, alpha=0.8)
ax.axhline(150_000, color="gray", linestyle="--", linewidth=0.9, alpha=0.8)
ax.text(14, 100_000 + 6000, "$100K/QALY", fontsize=8, color="dimgray")
ax.text(14, 150_000 + 6000, "$150K/QALY", fontsize=8, color="dimgray")

# Vertical reference lines: solid / dashed / dotted
# Labels placed in three horizontal rows at different heights to avoid overlap
vlines = [
    (200_000, "Base case price: $200K",         "solid",  "black",   204, 512_000),
    (429_000, "$100K/QALY threshold price: $429K", "dashed", "dimgray", 433, 489_000),
    (616_000, "$150K/QALY threshold price: $616K", "dotted", "dimgray", 620, 458_000),
]
for x, label, ls, c, lx, ly in vlines:
    ax.axvline(x / 1000, color=c, linestyle=ls, linewidth=1.3, alpha=0.85)
    ax.text(lx + 2, ly, label, fontsize=7.5, color=c, ha="left", va="top")

ax.set_xlabel("Intismeran acquisition price (US$000/course)", fontsize=11)
ax.set_ylabel("ICER (US$/QALY)", fontsize=11)
ax.set_xlim(0, 800)
ax.set_ylim(0, 540_000)
ax.set_xticks(np.arange(0, 801, 100))
ax.yaxis.set_major_formatter(lambda x, pos: f"{x/1000:.0f}K" if x else "0")
ax.grid(alpha=0.25, linewidth=0.5)
plt.tight_layout()
plt.savefig("fig4_price_icer.png", dpi=300)
plt.savefig("output/fig4_price_icer.png", dpi=300)
print("Saved fig4_price_icer.png (journal style)")