"""Regenerate fig4_price_icer.png with corrected title."""
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

prices = np.linspace(0, 1_500_000, 151)
icers = []
for pr in prices:
    _, _, icer = icer_at_price(pr)
    icers.append(icer)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(prices / 1000, icers, "b-", linewidth=2)
ax.axhline(100_000, color="gray", linestyle="--", linewidth=0.8, alpha=0.5, label="$100K/QALY")
ax.axhline(150_000, color="gray", linestyle=":", linewidth=0.8, alpha=0.5, label="$150K/QALY")
ax.axvline(200_000 / 1000, color="red", linestyle="--", linewidth=0.8, alpha=0.5, label="Base price ($200K)")
ax.set_xlabel("Intismeran Course Price ($1,000)", fontsize=12)
ax.set_ylabel("ICER ($/QALY)", fontsize=12)
ax.set_title("Price–ICER Curve and Deterministic Threshold Pricing", fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.set_xlim(0, 1500)
ax.set_ylim(0, 500_000)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("fig4_price_icer.png", dpi=300)
plt.savefig("output/fig4_price_icer.png", dpi=300)
print("Saved fig4_price_icer.png (root + output/)")